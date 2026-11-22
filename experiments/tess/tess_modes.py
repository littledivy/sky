"""Find additional TESS frequencies and check them in other sectors."""

import argparse
import json
from pathlib import Path
import warnings

from astropy.io import fits
from astropy.timeseries import LombScargle
from astropy.wcs import FITSFixedWarning, WCS
import numpy as np
from scipy.optimize import minimize

from sky.tess import download, photometry
from sky.tess_pixels import positive_centroid
from tess_contact import archive_curve


def basis(time, frequencies, cadence):
    harmonics = (np.asarray(frequencies)[:, None] * [1, 2]).ravel()
    phase = 2 * np.pi * (time - np.median(time))[:, None] * harmonics
    attenuation = np.sinc(cadence * harmonics)
    return np.column_stack((np.cos(phase) * attenuation, np.sin(phase) * attenuation))


def fit(time, flux, cadence, frequencies):
    matrix = np.column_stack((basis(time, frequencies, cadence), np.ones(len(time))))
    coefficients = np.linalg.lstsq(matrix, flux, rcond=None)[0]
    return coefficients, flux - matrix @ coefficients


def refine(time, flux, cadence, seeds):
    variance = np.var(flux)
    return minimize(lambda f: np.mean(fit(time, flux, cadence, f)[1]**2) / variance,
                    seeds, method="Nelder-Mead", bounds=[(f - .035, f + .035) for f in seeds],
                    options={"xatol": 1e-9}).x


def discover(curve, seed, components):
    time, flux, cadence = curve
    trials = [refine(time, flux, cadence, [frequency]) for frequency in (seed, 2*seed)]
    frequencies = min(trials, key=lambda f: np.mean(fit(time, flux, cadence, f)[1]**2))
    for _ in range(components - 1):
        _, residual = fit(time, flux, cadence, frequencies)
        trial, power = LombScargle(time, residual).autopower(
            minimum_frequency=1, maximum_frequency=12, samples_per_peak=10, method="fast")
        harmonics = (frequencies[:, None] * [1, 2]).ravel()
        separate = np.all(abs(trial[:, None] - harmonics) > .1, axis=1)
        strongest = trial[separate][np.argmax(power[separate])]
        frequencies = refine(time, flux, cadence, np.append(frequencies, strongest))
    return frequencies


def summarize(curve, frequencies):
    time, flux, cadence = curve
    coefficients, residual = fit(time, flux, cadence, frequencies)
    components = len(frequencies) * 2
    amplitude = np.hypot(coefficients[:components], coefficients[components:2*components])
    rms = np.sqrt(np.mean(residual**2))
    residual_power = np.mean(residual**2)
    variance_reductions = []
    for i in range(len(frequencies)):
        single_residual = fit(time, flux, cadence, np.delete(frequencies, i))[1]
        variance_reductions.append(float(1 - residual_power / np.mean(single_residual**2)))
    return {
        "cadences": len(time),
        "exposure_minutes": cadence * 1440,
        "independently_fitted_frequencies_per_day": refine(*curve, frequencies).tolist(),
        "fixed_frequency_semi_amplitudes_fraction": amplitude.reshape(-1, 2).tolist(),
        "fixed_frequency_r_squared": float(1 - rms**2 / np.var(flux)),
        "residual_variance_reduction_per_frequency": variance_reductions,
    }


def localize_frequencies(path, frequencies):
    with fits.open(path) as hdu:
        table = hdu[1].data
        selected = np.isfinite(table["TIME"]) & (table["QUALITY"] == 0)
        time = table["TIME"][selected]
        pixels = np.asarray(table["FLUX"][selected], dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(hdu[2].header).celestial
    rows, columns = pixels.shape[1:]
    flux = pixels.reshape(len(time), -1)
    segments = np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > .35) + 1)
    trend = np.zeros((len(time), 3 * len(segments)))
    for i, indices in enumerate(segments):
        offset = time[indices] - np.median(time[indices])
        offset /= max(1, np.max(abs(offset)))
        trend[indices, 3*i:3*i+3] = np.column_stack((np.ones(len(indices)), offset, offset**2))
    matrix = np.column_stack((basis(time, frequencies, 2/1440), trend))
    coefficients = np.zeros((matrix.shape[1], flux.shape[1]))
    for i in np.flatnonzero(np.isfinite(flux).mean(axis=0) > .9):
        selected = np.isfinite(flux[:, i])
        coefficients[:, i] = np.linalg.lstsq(matrix[selected], flux[selected, i], rcond=None)[0]
    locations = []
    for i in range(len(frequencies)):
        amplitude = np.hypot(coefficients[2*i], coefficients[2*i+2*len(frequencies)])
        amplitude = amplitude.reshape(rows, columns)
        peak = np.unravel_index(amplitude.argmax(), amplitude.shape)
        position = wcs.pixel_to_world(*positive_centroid(amplitude, peak, radius=1.5))
        locations.append({"ra_deg": float(position.ra.deg), "dec_deg": float(position.dec.deg)})
    return locations


def plot_curves(curves, frequencies, output):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(curves), len(frequencies), figsize=(4 * len(frequencies), 2.1 * len(curves)),
                             layout="constrained", sharex=True, squeeze=False)
    for row, (sector, curve) in enumerate(curves.items()):
        time, flux, cadence = curve
        coefficients, _ = fit(*curve, frequencies)
        matrix = basis(time, frequencies, cadence)
        for col, frequency in enumerate(frequencies):
            ax = axes[row, col]
            phase = ((time - np.median(time)) * frequency) % 1
            cosine = np.arange(2*col, 2*col + 2)
            indices = np.concatenate((cosine, cosine + 2*len(frequencies)))
            signal = matrix[:, indices] @ coefficients[indices]
            separated = flux - (matrix @ coefficients[:-1] - signal) - coefficients[-1]
            bins = np.minimum((phase * 30).astype(int), 29)
            points = np.array([np.median(separated[bins == i]) if (bins == i).any() else np.nan
                               for i in range(30)])
            points -= np.nanmedian(points)
            ax.plot((np.arange(30) + .5) / 30, points * 100, ".", color="#215485", ms=4)
            order = np.argsort(phase)
            ax.plot(phase[order], signal[order] * 100, color="#b54a22", lw=1)
            if col == 0:
                ax.set_ylabel(f"S{sector} flux (%)")
            if row == 0:
                ax.set_title(f"{24 / frequency:.2f} hours")
            if row == len(curves) - 1:
                ax.set_xlabel("Phase (cycles)")
    fig.savefig(output, dpi=155)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("tic", type=int)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ffi-sectors", nargs="*", type=int, default=[])
    parser.add_argument("--components", type=int, choices=(2, 3), default=2)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    candidate = next(c for c in report["candidates"] if c["tic"] == args.tic and c["status"] == "repeated")
    ra, dec = candidate["RA_OBJ"], candidate["DEC_OBJ"]
    discovery, validation = report["discovery_sector"], report["validation_sector"]
    curves, details = {}, {}
    for sector, key in ((discovery, "discovery_source"), (validation, "validation_source")):
        url = candidate[key]
        lightcurve = download(url, args.cache / "fits" / url.rsplit("/", 1)[-1], verify_fits=True)
        _, values = photometry(lightcurve, native=True)
        curves[sector] = (*values[:2], 2/1440)
        details[sector] = {"light_curve_url": url, "pixel_url": url.replace("_lc.fits", "_tp.fits")}
    frequencies = discover(curves[discovery], 1 / candidate["spacing_days"], args.components)
    for sector in args.ffi_sectors:
        url, curves[sector] = archive_curve(args.cache, args.tic, sector, ra, dec)
        details[sector] = {"full_frame_cutout_url": url}
    for sector, curve in curves.items():
        details[sector].update(summarize(curve, frequencies))
        url = details[sector].get("pixel_url")
        if url:
            pixels = download(url, args.cache / "pixels" / url.rsplit("/", 1)[-1], verify_fits=True)
            details[sector]["frequency_amplitude_centroids_icrs"] = localize_frequencies(pixels, frequencies)
    result = {"tic": args.tic, "aperture_icrs_deg": [ra, dec], "discovery_sector": discovery,
              "discovery_frequencies_per_day": frequencies.tolist(),
              "discovery_periods_hours": (24 / frequencies).tolist(), "sectors": details}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "fit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    plot_curves(curves, frequencies, args.output / "phase.png")


if __name__ == "__main__":
    main()
