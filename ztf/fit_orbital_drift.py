#!/usr/bin/env python3
"""Fit a constant or changing period to candidate photometry from ZTF DR24."""

import argparse
import json

from astropy import units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
import numpy as np
from pyarrow.parquet import read_table
from scipy.optimize import minimize, minimize_scalar

if __package__:
    from .periodogram import periodogram
    from .short_period_eclipses import BAD_FLAG, MINUTES_PER_DAY
else:
    from periodogram import periodogram
    from short_period_eclipses import BAD_FLAG, MINUTES_PER_DAY


def read_photometry(source):
    path, separator, objectid = source.rpartition(":")
    if not separator:
        raise ValueError(f"Expected PARQUET:OBJECT_ID, got {source}")
    columns = ["objra", "objdec", "hmjd", "mag", "magerr", "catflags"]
    rows = read_table(path, columns=columns, filters=[("objectid", "=", int(objectid))])
    if rows.num_rows != 1:
        raise ValueError(f"Expected one light curve in {source}, found {rows.num_rows}")
    row, = rows.to_pylist()
    hmjd, magnitude, error = (np.asarray(row[name], dtype=float)
                             for name in ("hmjd", "mag", "magerr"))
    good = np.isfinite(hmjd) & np.isfinite(magnitude) & np.isfinite(error)
    good &= (error > 0) & ((np.asarray(row["catflags"]) & BAD_FLAG) == 0)
    if np.count_nonzero(good) < 10:
        raise ValueError(f"Too few usable observations in {source}")
    heliocentric = Time(hmjd[good], format="mjd", scale="utc")
    # HMJD already has the observatory-to-Sun correction. Add the projected
    # Sun-to-barycentre term and express the times in TDB.
    direction = SkyCoord(row["objra"], row["objdec"], unit="deg").cartesian
    delay = (get_body_barycentric("sun", heliocentric.tdb).dot(direction) / c).to_value(u.day)
    return heliocentric.tdb.mjd + delay, magnitude[good], error[good]


def fit_drift(photometry, period_minutes):
    start = min(time.min() for time, _, _ in photometry)
    end = max(time.max() for time, _, _ in photometry)
    midpoint = (start + end) / 2
    half_baseline = (end - start) / 2
    if half_baseline <= 0:
        raise ValueError("Observations must span more than one timestamp")
    bands = [((time - midpoint) / half_baseline, magnitude, error)
             for time, magnitude, error in photometry]
    frequency_guess = MINUTES_PER_DAY / period_minutes

    # Search the seed frequency and its daily aliases. Work in cycles across
    # the half-baseline so frequency and curvature have comparable scales.
    frequencies = np.linspace(frequency_guess - 2, frequency_guess + 2,
                              max(3, int(40 * half_baseline) + 1))

    def chi_squared(frequency, curvature):
        total = 0.0
        for time, magnitude, error in bands:
            phase = 2 * np.pi * (frequency * half_baseline * time + curvature * time**2)
            design = np.column_stack((np.ones_like(time), np.cos(phase), np.sin(phase),
                                      np.cos(2 * phase), np.sin(2 * phase)))
            design /= error[:, None]
            scaled = magnitude / error
            coefficients = np.linalg.lstsq(design, scaled, rcond=None)[0]
            total += np.sum((scaled - design @ coefficients)**2)
        return total

    def frequency_seeds(curvature):
        power = np.zeros_like(frequencies)
        for time, magnitude, error in bands:
            mean = np.average(magnitude, weights=error**-2)
            flat_chi_squared = np.sum(((magnitude - mean) / error)**2)
            adjusted_time = half_baseline * time + curvature * time**2 / frequency_guess
            power += flat_chi_squared * periodogram(adjusted_time, magnitude, error, frequencies)
        peaks = np.flatnonzero(
            np.isfinite(power)
            & (power > np.r_[-np.inf, power[:-1]])
            & (power > np.r_[power[1:], -np.inf])
        )
        return frequencies[peaks[np.argsort(power[peaks])[-12:]]]

    constant_fits = []
    for frequency in frequency_seeds(0):
        result = minimize_scalar(
            lambda shift: chi_squared(frequency + shift / half_baseline, 0),
            bounds=(-0.10, 0.10), method="bounded",
        )
        constant_fits.append((result.fun, frequency + result.x / half_baseline))
    constant_chi_squared, constant_frequency = min(constant_fits)

    trials = [(chi_squared(frequency, curvature), frequency, curvature)
              for curvature in np.linspace(-3, 3, 31)
              for frequency in frequency_seeds(curvature)]
    drift_fits = []
    for _, frequency, curvature in sorted(trials)[:24]:
        result = minimize(
            lambda parameters: chi_squared(frequency + parameters[0] / half_baseline, parameters[1]),
            [0, curvature], method="Nelder-Mead",
            options={
                "initial_simplex": [[0, curvature], [0.04, curvature], [0, curvature + 0.04]],
                "xatol": 1e-7, "fatol": 1e-7,
            },
        )
        drift_fits.append((result.fun, frequency + result.x[0] / half_baseline, result.x[1]))
    drift_chi_squared, drift_frequency, drift_curvature = min(drift_fits)

    def local_chi_squared(shift, curve_offset):
        return chi_squared(drift_frequency + shift / half_baseline, drift_curvature + curve_offset)

    step = 0.001
    mixed = (local_chi_squared(step, step) + local_chi_squared(-step, -step)
             - local_chi_squared(step, -step) - local_chi_squared(-step, step)) / 4
    hessian = np.array([
        [local_chi_squared(step, 0) + local_chi_squared(-step, 0) - 2 * drift_chi_squared, mixed],
        [mixed, local_chi_squared(0, step) + local_chi_squared(0, -step) - 2 * drift_chi_squared],
    ]) / step**2
    dof = sum(len(time) for time, _, _ in bands) - 5 * len(bands) - 2
    covariance = 2 * np.linalg.inv(hessian) * max(1, drift_chi_squared / dof)
    derivative_scale = -2 / (half_baseline**2 * drift_frequency**2)
    return {
        "epochs": [len(time) for time, _, _ in bands],
        "baseline_days": end - start,
        "epoch_bmjd_tdb": midpoint,
        "input_period_minutes": period_minutes,
        "best_period_minutes": MINUTES_PER_DAY / drift_frequency,
        "pdot_s_per_s": derivative_scale * drift_curvature,
        "local_pdot_stderr": abs(derivative_scale) * np.sqrt(covariance[1, 1]),
        "constant_period_minutes": MINUTES_PER_DAY / constant_frequency,
        "constant_chi2": constant_chi_squared,
        "drifting_chi2": drift_chi_squared,
        "dof": dof,
        "delta_chi2": constant_chi_squared - drift_chi_squared,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="+", metavar="PARQUET:OBJECT_ID")
    parser.add_argument("--period-minutes", type=float, required=True, help="Period from the survey scan")
    args = parser.parse_args()
    if not 0 < args.period_minutes < MINUTES_PER_DAY / 2:
        parser.error("Pass a positive period shorter than 12 hours.")
    try:
        photometry = [read_photometry(source) for source in args.sources]
        result = fit_drift(photometry, args.period_minutes)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (ValueError, np.linalg.LinAlgError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
