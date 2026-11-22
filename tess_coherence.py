"""Test whether one TESS frequency predicts phase in held-out sectors."""

import argparse
import json
from pathlib import Path

from astropy.io import fits
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

from sky.tess import download, photometry
from tess_modes import basis


def project_out(time, flux, cadence, frequency, other_frequencies):
    harmonic = 4 * np.pi * frequency * time
    nuisance = np.column_stack((basis(time, other_frequencies, cadence), np.cos(harmonic),
                                np.sin(harmonic), np.ones(len(time)), time - np.median(time)))
    q = np.linalg.qr(nuisance, mode="reduced")[0]
    residual = flux - q @ (q.T @ flux)
    return time, residual, q, cadence


def measure(frequency, curve):
    time, flux, q, cadence = curve
    phase = 2 * np.pi * frequency * time
    matrix = np.column_stack((np.cos(phase), np.sin(phase))) * np.sinc(cadence * frequency)
    matrix -= q @ (q.T @ matrix)
    return matrix.T @ matrix, matrix.T @ flux


def solve(frequency, curves):
    matrices, vectors = zip(*(measure(frequency, curve) for curve in curves))
    matrix, vector = sum(matrices), sum(vectors)
    coefficients = np.linalg.solve(matrix, vector)
    return coefficients, float(vector @ coefficients)


def search(curves, frequency, half_width):
    time_min = min(curve[0][0] for curve in curves)
    time_max = max(curve[0][-1] for curve in curves)
    samples = max(31, int(np.ceil(20 * half_width * (time_max - time_min))))
    grid = frequency + np.linspace(-half_width, half_width, samples)
    scores = np.array([solve(f, curves)[1] for f in grid])
    peaks, _ = find_peaks(scores)
    if not len(peaks) or scores.max() > scores[peaks].max():
        raise ValueError("Strongest frequency lies at the edge of the search interval")
    aliases = []
    for sample in peaks[np.argsort(scores[peaks])[-4:][::-1]]:
        fit = minimize_scalar(lambda f: -solve(f, curves)[1], bounds=grid[[sample - 1, sample + 1]],
                              method="bounded", options={"xatol": 1e-12})
        aliases.append({"frequency_per_day": float(fit.x), "fit_power": float(-fit.fun)})
    return sorted(aliases, key=lambda peak: peak["fit_power"], reverse=True)


def read_curves(fit, products, cache):
    curves = {}
    for product in products:
        if int(product["target_name"]) != fit["tic"]:
            continue
        url = product["url"]
        path = download(url, cache / product["dataURL"].rsplit("/", 1)[-1], verify_fits=True)
        with fits.open(path) as hdu:
            sector = int(hdu[0].header["SECTOR"])
            time, flux = (np.asarray(hdu[1].data[name], dtype=float) for name in ("TIME", "FLUX"))
        selected = np.isfinite(time + flux)
        time, flux = time[selected], flux[selected]
        curves[sector] = (time, flux, float(np.median(np.diff(time))))
    for sector, info in fit["sectors"].items():
        url = info.get("light_curve_url")
        if url:
            path = download(url, cache / "fits" / url.rsplit("/", 1)[-1], verify_fits=True)
            _, (time, flux, _) = photometry(path, native=True)
            curves[int(sector)] = (time, flux, 2/1440)
    return curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument("tars_products", type=Path)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--held-out", type=int, nargs="+", required=True)
    parser.add_argument("--component", type=int, default=2, help="One-based frequency index")
    parser.add_argument("--half-width", type=float, default=.0025, help="Frequency window in cycles per day")
    args = parser.parse_args()
    fit = json.loads(args.fit.read_text())
    curves = read_curves(fit, json.loads(args.tars_products.read_text()), args.cache)
    held_out = set(args.held_out)
    if not held_out < curves.keys():
        parser.error("Holdout sectors must be present and leave at least one training sector")
    frequencies = np.array(fit["discovery_frequencies_per_day"])
    index = args.component - 1
    if not 0 <= index < len(frequencies):
        parser.error("Component is outside the fitted frequency list")
    frequency, others = frequencies[index], np.delete(frequencies, index)
    reference = float(np.median(curves[fit["discovery_sector"]][0]))
    standardized = {sector: project_out(time - reference, flux, cadence, frequency, others)
                    for sector, (time, flux, cadence) in curves.items()}
    training = [curve for s, curve in standardized.items() if s not in held_out]
    aliases = search(training, frequency, args.half_width)
    best = aliases[0]["frequency_per_day"]
    train_coefficients, _ = solve(best, training)
    train_phase = np.arctan2(train_coefficients[1], train_coefficients[0]) / (2*np.pi)
    sectors = {}
    for sector, curve in sorted(standardized.items()):
        coefficients, _ = solve(best, [curve])
        phase = np.arctan2(coefficients[1], coefficients[0]) / (2*np.pi)
        offset = (phase - train_phase + .5) % 1 - .5
        sectors[sector] = {"held_out": sector in held_out, "phase_offset_minutes": float(offset * 1440 / best),
                           "semi_amplitude_ppm": float(np.linalg.norm(coefficients) * 1e6)}
    result = {
        "tic": fit["tic"], "reference_btjd": reference, "training_aliases": aliases,
        "training_period_hours": 24 / best, "sectors": sectors,
        "joint_fit_frequency_per_day": search(list(standardized.values()), frequency, args.half_width)[0][
            "frequency_per_day"],
    }
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
