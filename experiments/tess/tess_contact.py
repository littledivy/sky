"""Validate an alternating TESS signal with a fixed Fourier light-curve template."""

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from scipy.optimize import minimize_scalar

from sky.tess import download, photometry
from sky.tess_pixels import localize
from tess_vet import aperture_photometry, write_phase


def design(time, period, epoch, cadence=0, *, even_only=False):
    harmonics = np.arange(2, 9, 2) if even_only else np.arange(1, 9)
    angles = 2 * np.pi * ((time - epoch) / period)[:, None] * harmonics
    smoothing = np.sinc(cadence * harmonics / period)
    return np.column_stack((np.ones(len(time)), np.cos(angles) * smoothing, np.sin(angles) * smoothing))


def template(time, flux, period, epoch, cadence=0, *, even_only=False):
    matrix = design(time, period, epoch, cadence, even_only=even_only)
    coefficients = np.linalg.lstsq(matrix, flux, rcond=None)[0]
    return coefficients, float(np.mean((flux - matrix @ coefficients)**2))


def fit_period(curve, seed, epoch):
    time, flux, cadence = curve
    return float(minimize_scalar(
        lambda period: template(time, flux, period, epoch, cadence)[1],
        bounds=(seed * 0.9997, seed * 1.0003), method="bounded",
        options={"xatol": 1e-12}).x)


def score(curve, period, epoch, coefficients, shift=0):
    time, flux, cadence = curve
    predicted = design(time, period, epoch + shift * period, cadence) @ coefficients
    predicted -= predicted.mean()
    centered = flux - flux.mean()
    scale = max(0, predicted @ centered / (predicted @ predicted))
    return float(1 - np.mean((centered - scale * predicted)**2) / np.var(centered))


def fit_shift(curve, period, epoch, coefficients):
    grid = np.linspace(-.5, .5, 101)
    phase = grid[np.argmax([score(curve, period, epoch, coefficients, shift) for shift in grid])]
    result = minimize_scalar(lambda shift: -score(curve, period, epoch, coefficients, shift),
                             bounds=(phase - .01, phase + .01), method="bounded")
    return float(result.x), float(-result.fun)


def archive_curve(cache, tic, sector, ra, dec):
    query = urlencode({"ra": ra, "dec": dec, "x": 9, "y": 9, "units": "px", "sector": sector})
    url = f"https://mast.stsci.edu/tesscut/api/v0.1/astrocut?{query}"
    path = download(url, cache / f"{tic}-cut-s{sector:04d}.zip")
    time, flux, cadence = aperture_photometry(path)
    selected = np.isfinite(time + flux) & (abs(flux) < .4)
    return url, (time[selected], flux[selected], cadence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("tic", type=int)
    parser.add_argument("--orbital-period", required=True, type=float)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ffi-sectors", nargs="*", type=int, default=[])
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    candidate = next(item for item in report["candidates"] if item["tic"] == args.tic
                     and item["status"] == "repeated")
    epoch = candidate["epoch_btjd"]
    ra, dec = candidate["RA_OBJ"], candidate["DEC_OBJ"]
    discovery, validation = report["discovery_sector"], report["validation_sector"]
    paths, curves, sectors = {}, {}, {}
    for sector, key in ((discovery, "discovery_source"), (validation, "validation_source")):
        url = candidate[key]
        paths[sector] = download(url, args.cache / "fits" / url.rsplit("/", 1)[-1], verify_fits=True)
        _, data = photometry(paths[sector], native=True)
        curves[sector] = (*data[:2], 2/1440)
        sectors[sector] = {"light_curve_url": url}
    p = fit_period(curves[discovery], args.orbital_period, epoch)
    coefficients, mse = template(*curves[discovery][:2], p, epoch, curves[discovery][2])
    _, mse_half = template(*curves[discovery][:2], p, epoch, curves[discovery][2], even_only=True)
    phases = np.linspace(-.2, .8, 6001)
    model = design(epoch + phases * p, p, epoch) @ coefficients
    groups = (phases < .25, phases >= .25)
    minima = [float(phases[g][np.argmin(model[g])]) for g in groups]
    for sector in args.ffi_sectors:
        url, curves[sector] = archive_curve(args.cache, args.tic, sector, ra, dec)
        sectors[sector] = {"full_frame_cutout_url": url}
    shifts = []
    for sector, curve in curves.items():
        shift, adjusted = fit_shift(curve, p, epoch, coefficients)
        mid = float(np.median(curve[0]))
        shifts.append((mid, shift))
        sectors[sector].update(points=len(curve[0]), cadence_minutes=curve[2] * 1440,
                               independent_period_days=fit_period(curve, p, epoch),
                               fixed_phase_r_squared=score(curve, p, epoch, coefficients),
                               adjusted_phase_r_squared=adjusted, phase_shift_cycles=shift)
    from scipy.stats import linregress

    times, offsets = np.array(shifts).T
    drift = linregress(times - epoch, offsets)
    refined = 1 / (1/p - drift.slope)
    refined_epoch = epoch + (drift.intercept + minima[0]) * refined
    for sector in paths:
        url = sectors[sector]["light_curve_url"].replace("_lc.fits", "_tp.fits")
        pixel = download(url, args.cache / "pixels" / url.rsplit("/", 1)[-1], verify_fits=True)
        sectors[sector]["pixel_url"] = url
        sectors[sector]["minimum_localizations"] = [
            localize(pixel, p, epoch + phase*p, p/5, target_ra=ra, target_dec=dec)
            for phase in minima
        ]
    results = {
        "aperture_tic": args.tic, "aperture_coordinate_icrs_deg": [ra, dec],
        "discovery_period_days": p, "discovery_mse_full_template": mse,
        "discovery_mse_half_period_template": mse_half,
        "primary_and_secondary_relative_flux": [float(model[g].min()) for g in groups],
        "orbital_period_days": refined, "mid_eclipse_btjd": refined_epoch,
        "long_baseline_sector_timing_residual_rms_minutes": float(
            np.sqrt(np.mean((offsets - drift.intercept - drift.slope*(times-epoch))**2)) * refined * 1440),
        "sectors": sectors,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "fit.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    recent = {sector: photometry(path, native=True)[1] for sector, path in paths.items()}
    write_phase(recent, {"period_days": refined, "mid_eclipse_btjd": refined_epoch}, args.output)


if __name__ == "__main__":
    main()
