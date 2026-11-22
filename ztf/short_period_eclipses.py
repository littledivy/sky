#!/usr/bin/env python3
"""Search ZTF DR24 bulk Parquet for short-period eclipses."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import heapq
import json
from pathlib import Path
import sys

from astropy.timeseries import BoxLeastSquares
import numpy as np
from pyarrow.parquet import ParquetFile

if __package__:
    from .periodogram import periodogram
else:
    from periodogram import periodogram


COLUMNS = (
    "objectid", "objra", "objdec", "filterid", "hmjd", "mag", "magerr", "catflags",
)
MINUTES_PER_DAY = 1440
BAD_FLAG = 0x8000


@dataclass(frozen=True)
class Search:
    min_period: float
    max_period: float
    window_days: float
    min_epochs: int
    top: int
    limit: int | None


def search_curve(row, settings):
    time = np.asarray(row["hmjd"], dtype=float)
    mag = np.asarray(row["mag"], dtype=float)
    error = np.asarray(row["magerr"], dtype=float)
    flags = np.asarray(row["catflags"])
    valid = np.isfinite(time) & np.isfinite(mag) & np.isfinite(error)
    valid &= (error > 0) & ((flags & BAD_FLAG) == 0)
    order = np.argsort(time[valid])
    time, mag, error = (values[valid][order] for values in (time, mag, error))
    _, unique = np.unique(time, return_index=True)
    time, mag, error = (values[unique] for values in (time, mag, error))
    if len(time) < settings.min_epochs:
        return "insufficient_good_epochs", None

    # Choose a window using observation times only. A short time baseline keeps
    # the blind frequency grid tractable and limits smearing from orbital drift.
    ends = np.searchsorted(time, time + settings.window_days, side="right")
    start = np.argmax(ends - np.arange(len(time)))
    window = slice(start, ends[start])
    epoch = time[start]
    time, mag, error = time[window] - epoch, mag[window], error[window]
    if len(time) < settings.min_epochs:
        return "insufficient_window_epochs", None

    # Split observations before choosing any frequency using photometry.
    held_out = np.arange(len(time)) % 3 == 2
    train_time = time[~held_out]
    # Negating magnitudes makes eclipses dips for BoxLeastSquares.
    train_y, train_error = -mag[~held_out], error[~held_out]
    baseline = np.ptp(train_time)
    if baseline < 3 * settings.max_period / MINUTES_PER_DAY:
        return "insufficient_baseline", None

    min_freq = MINUTES_PER_DAY / settings.max_period
    max_freq = MINUTES_PER_DAY / settings.min_period
    step = 1 / (5 * baseline)
    intervals = int(np.ceil((max_freq - min_freq) / step))
    frequencies = np.linspace(min_freq, max_freq, max(3, intervals + 1))
    step = frequencies[1] - frequencies[0]
    power = periodogram(train_time, train_y, train_error, frequencies)
    peaks = np.flatnonzero(
        np.isfinite(power)
        & (power >= np.r_[-np.inf, power[:-1]])
        & (power >= np.r_[power[1:], -np.inf])
    )
    if not len(peaks):
        return "no_frequency_peaks", None

    centers = frequencies[peaks[np.argsort(power[peaks])[-8:]]]
    # A sinusoidal peak can land at the first harmonic of the orbital frequency.
    centers = np.r_[centers, centers / 2, centers * 2]
    trial = (centers[:, None] + step * np.linspace(-2, 2, 17)).ravel()
    trial = np.unique(trial[(trial >= min_freq) & (trial <= max_freq)])
    bls = BoxLeastSquares(train_time, train_y, train_error)
    durations = np.array([15, 30, 60, 90]) / 86400
    durations = durations[durations < settings.min_period / MINUTES_PER_DAY]
    fits = bls.power(1 / trial, durations, objective="likelihood")
    eligible = np.flatnonzero(np.isfinite(fits.power) & (fits.depth > 0))
    if not len(eligible):
        return "no_eclipse_fit", None

    best = eligible[np.argmax(fits.power[eligible])]
    period, duration, mid = (array[best] for array in
                             (fits.period, fits.duration, fits.transit_time))
    validation_time = time[held_out]
    observed = -mag[held_out]
    predicted = bls.model(validation_time, period, duration, mid)
    flat = np.average(train_y, weights=train_error ** -2)
    sigma = error[held_out]
    gain = np.sum(((observed - flat) / sigma) ** 2
                  - ((observed - predicted) / sigma) ** 2)
    cycles = np.rint((validation_time - mid) / period)
    in_eclipse = np.abs(validation_time - mid - cycles * period) < duration / 2
    eclipse_count = len(np.unique(cycles[in_eclipse]))
    if eclipse_count < 2 or gain <= 0:
        return "failed_holdout", None

    return "candidate", {
        **({"source": row["_file"]} if "_file" in row else {}),
        "objectid": row["objectid"],
        "ra_deg": row["objra"],
        "dec_deg": row["objdec"],
        "filterid": row["filterid"],
        "period_minutes": float(period * MINUTES_PER_DAY),
        "eclipse_duration_seconds": float(duration * 86400),
        "eclipse_epoch_hmjd": float(epoch + mid),
        "window_start_hmjd": float(epoch),
        "training_epochs": len(train_time),
        "holdout_epochs": len(validation_time),
        "holdout_eclipses": eclipse_count,
        "holdout_delta_chi2": float(gain),
    }


def scan_rows(rows, total, partition, settings):
    counts = {"rows_total": total, "rows_examined": 0, "eligible_searched": 0}
    candidates = []
    try:
        for row in rows:
            status, candidate = search_curve(row, settings)
            counts["rows_examined"] += 1
            counts[status] = counts.get(status, 0) + 1
            if not status.startswith("insufficient_"):
                counts["eligible_searched"] += 1
            if candidate:
                key = (candidate["holdout_delta_chi2"], counts["rows_examined"])
                heapq.heappush(candidates, (key, candidate))
                if len(candidates) > settings.top:
                    heapq.heappop(candidates)
            if settings.limit is not None and counts["eligible_searched"] >= settings.limit:
                break
    except Exception as error:
        return {"partition": partition, "error": str(error), "coverage": counts, "candidates": []}
    return {
        "partition": partition,
        "coverage": counts,
        "candidates": [candidate for _, candidate in sorted(candidates, reverse=True)],
    }


def scan_partition(path, settings):
    partition = str(path)
    try:
        with ParquetFile(path) as parquet:
            rows = (row for batch in parquet.iter_batches(batch_size=128, columns=COLUMNS)
                    for row in batch.to_pylist())
            return scan_rows(rows, parquet.metadata.num_rows, partition, settings)
    except Exception as error:
        counts = {"rows_total": 0, "rows_examined": 0, "eligible_searched": 0}
        return {"partition": partition, "error": str(error), "coverage": counts, "candidates": []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("partitions", nargs="+", type=Path, help="Parquet files or directories")
    parser.add_argument("--min-period", type=float, default=5, help="Minimum trial period in minutes")
    parser.add_argument("--max-period", type=float, default=60, help="Maximum trial period in minutes")
    parser.add_argument("--window-days", type=float, default=60)
    parser.add_argument("--min-epochs", type=int, default=51)
    parser.add_argument("--top", type=int, default=20, help="Candidates to retain per partition")
    parser.add_argument("--limit-per-partition", type=int, help="Stop after this many eligible curves")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if not (1.5 < args.min_period < args.max_period and args.window_days > 0
            and args.min_epochs >= 10 and args.top >= 1 and args.workers >= 1
            and (args.limit_per_partition is None or args.limit_per_partition >= 1)):
        parser.error("Invalid search bounds, epoch minimum, candidate limit, or worker count.")

    files = set()
    for path in args.partitions:
        if path.is_dir():
            partitions = list(path.rglob("*.parquet"))
            if not partitions:
                parser.error(f"No Parquet partitions in {path}")
            files.update(file.resolve() for file in partitions)
        elif path.is_file():
            files.add(path.resolve())
        else:
            parser.error(f"Input path does not exist: {path}")

    settings = Search(args.min_period, args.max_period, args.window_days, args.min_epochs,
                      args.top, args.limit_per_partition)
    paths = sorted(files)
    failed = False

    def emit(results):
        nonlocal failed
        for result in results:
            result["search"] = vars(settings)
            print(json.dumps(result, allow_nan=False), flush=True)
            failed |= "error" in result

    if args.workers == 1:
        emit(scan_partition(path, settings) for path in paths)
    else:
        from itertools import repeat

        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            emit(pool.map(scan_partition, paths, repeat(settings),
                          buffersize=2 * args.workers))
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
