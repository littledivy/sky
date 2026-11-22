"""Search public TESS light curves with independent-sector eclipse validation."""

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import random
import re
from tempfile import NamedTemporaryFile
import threading
from time import monotonic

from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import numpy as np
import requests

from ztf.periodogram import periodogram
from . import tess_sparse


MANIFEST = "https://archive.stsci.edu/missions/tess/download_scripts/sector"
PUBLIC_TESS = "https://stpubdata.s3.us-east-1.amazonaws.com/tess/public/tid"
FILENAME = re.compile(r"tess\d+-s(?P<sector>\d{4})-(?P<tic>\d{16})-\d+-s_lc\.fits")
DOWNLOADS = threading.local()
DEFAULT_MIN = 0.12
DEFAULT_MAX = 1.5
BIN_DAYS = 6 / 1440
DURATIONS = np.array([0.010, 0.018, 0.03, 0.05, 0.08])


def download(url, destination, *, verify_fits=False):
    if destination.exists():
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not hasattr(DOWNLOADS, "session"):
        DOWNLOADS.session = requests.Session()
    temporary = None
    for attempt in range(3):
        try:
            with DOWNLOADS.session.get(url, timeout=(8, 35), stream=True) as response:
                response.raise_for_status()
                with NamedTemporaryFile(dir=destination.parent, delete=False) as output:
                    temporary = Path(output.name)
                    chunks = response.iter_content(chunk_size=256 * 1024)
                    initial = next(chunks, b"")
                    if not initial or (verify_fits and not initial.startswith(b"SIMPLE")):
                        raise ValueError(f"Unexpected archive response at {url}")
                    output.write(initial)
                    for chunk in chunks:
                        output.write(chunk)
            temporary.replace(destination)
            return destination
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
        finally:
            if temporary:
                temporary.unlink(missing_ok=True)
    raise RuntimeError("Download attempts exhausted.")


def sector_products(sector, cache):
    filename = f"tesscurl_sector_{sector}_lc.sh"
    manifest = download(f"{MANIFEST}/{filename}", cache / f"sector-{sector}.txt")
    products = {
        int(match["tic"]): match.group()
        for match in FILENAME.finditer(manifest.read_text())
        if int(match["sector"]) == sector
    }
    if not products:
        raise ValueError(f"No 2-minute light curves in the sector {sector} manifest.")
    return products


def product_url(filename):
    match = FILENAME.fullmatch(filename)
    if not match:
        raise ValueError(f"Invalid TESS light-curve filename: {filename}")
    tic = match["tic"]
    ids = "/".join(tic[i:i+4] for i in range(0, 16, 4))
    return f"{PUBLIC_TESS}/s{match['sector']}/{ids}/{filename}"


def photometry(path, *, native=False):
    with fits.open(path, memmap=True) as hdu:
        table = hdu[1].data
        values = {
            key: float(hdu[0].header[key])
            for key in ("RA_OBJ", "DEC_OBJ", "TESSMAG") if key in hdu[0].header
        }
        values["crowdsap"] = float(hdu[1].header.get("CROWDSAP", 0))
        time, flux, error = [np.asarray(table[name], dtype=np.float64)
                             for name in ("TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR")]
        good = np.isfinite(time + flux + error) & (error > 0) & (table["QUALITY"] == 0)
    if good.sum() < 600:
        return values, None
    order = np.argsort(time[good])
    time, flux, error = (column[good][order] for column in (time, flux, error))
    median = np.median(flux)
    if median <= 0:
        return values, None
    if native:
        return values, (time, flux / median - 1, error / median)
    slots = np.floor((time - time[0]) / BIN_DAYS + 1e-8).astype(int)
    weight = (error / median) ** -2
    binned_weight = np.bincount(slots, weights=weight)
    filled = binned_weight > 0
    slot_weight = binned_weight[filled]
    time = np.bincount(slots, (time - time[0]) * weight)[filled] / slot_weight + time[0]
    flux = np.bincount(slots, ((flux / median) - 1) * weight)[filled] / slot_weight
    error = 1 / np.sqrt(slot_weight)
    adjacent = np.diff(flux)[np.diff(time) < 1.6 * BIN_DAYS]
    if len(adjacent) > 20:
        floor = 1.4826 / np.sqrt(2) * np.median(np.abs(adjacent - np.median(adjacent)))
        error = np.maximum(error, floor)
    return values, (time, flux, error)


def cycle_medians(cycles, flux, selected, minimum):
    cycles, flux = cycles[selected], flux[selected]
    order = np.lexsort((flux, cycles))
    ids, start, count = np.unique(cycles[order], return_index=True, return_counts=True)
    good = count >= minimum
    start, count = start[good], count[good]
    lower = order[start + (count - 1) // 2]
    upper = order[start + count // 2]
    return ids[good], (flux[lower] + flux[upper]) / 2


def predicted_eclipses(photometry_data, period, epoch, duration):
    time, flux, _ = photometry_data
    times = time - epoch
    cycles = np.rint(times / period).astype(int)
    phases = abs(times - cycles * period)
    in_eclipse = phases < duration / 2
    nearby = (phases > 0.7 * duration) & (phases < min(1.8 * duration, 0.45 * period))
    in_ids, in_flux = cycle_medians(cycles, flux, in_eclipse, 2)
    out_ids, out_flux = cycle_medians(cycles, flux, nearby, 3)
    cycle_ids, in_index, out_index = np.intersect1d(in_ids, out_ids, assume_unique=True,
                                                   return_indices=True)
    depths = out_flux[out_index] - in_flux[in_index]
    if len(depths) < 4:
        return {"cycles": len(depths), "status": "insufficient_eclipse_coverage"}
    depth = float(np.median(depths))
    scatter = 1.4826 * np.median(abs(depths - depth))
    uncertainty = max(float(scatter) * 1.2533 / np.sqrt(len(depths)), 1e-7)
    return {"status": "measured", "cycles": len(depths), "depth_fraction": depth,
            "depth_sigma": depth / uncertainty, "dimming_fraction": float(np.mean(depths > 0)),
            "even_depth_fraction": median_or_none(depths[cycle_ids % 2 == 0]),
            "odd_depth_fraction": median_or_none(depths[cycle_ids % 2 == 1])}


def median_or_none(values):
    return float(np.median(values)) if len(values) else None


def passes(eclipses):
    return (eclipses["status"] == "measured" and eclipses["cycles"] >= 6
            and eclipses["depth_fraction"] >= 0.001
            and eclipses["depth_sigma"] >= 7 and eclipses["dimming_fraction"] >= 0.8)


def find_eclipse(data, minimum=DEFAULT_MIN, maximum=DEFAULT_MAX):
    time, flux, error = data
    t = time - time[0]
    baseline = np.ptp(t)
    if baseline < 8:
        return None
    low, high = 1 / maximum, 1 / minimum
    count = max(3, int(np.ceil((high - low) * 10 * baseline)) + 1)
    frequencies = np.linspace(low, high, count)
    step = frequencies[1] - frequencies[0]
    power = periodogram(t, flux, error, frequencies)
    peaks = np.flatnonzero(np.isfinite(power) & (power >= np.r_[-np.inf, power[:-1]])
                          & (power >= np.r_[power[1:], -np.inf]))
    if not len(peaks):
        return None
    best = peaks[np.argsort(power[peaks])[-12:]]
    centers = frequencies[best]
    centers = np.r_[centers, centers / 2, centers * 2]
    trials = (centers[:, None] + step * np.linspace(-1, 1, 9)).ravel()
    trials = np.unique(trials[(trials >= low) & (trials <= high)])
    model = BoxLeastSquares(t, flux, error)
    durations = DURATIONS[DURATIONS < minimum]

    def evaluate(frequencies):
        result = model.power(1 / frequencies, durations, objective="snr", oversample=8)
        eligible = np.flatnonzero(np.isfinite(result.power) & (result.duration <= 0.25 * result.period))
        return result, eligible

    coarse, eligible = evaluate(trials)
    if not len(eligible):
        return None
    # Refine the strongest separated periods before predicting another sector.
    seeds = []
    for index in eligible[np.argsort(coarse.power[eligible])[::-1]]:
        frequency = 1 / coarse.period[index]
        if all(abs(frequency - seen) > 3 * step for seen in seeds):
            seeds.append(frequency)
        if len(seeds) == 4:
            break
    refined = (np.asarray(seeds)[:, None] + step * np.linspace(-1.4, 1.4, 81)).ravel()
    refined = np.unique(refined[(refined >= low) & (refined <= high)])
    fit, eligible = evaluate(refined)
    for index in eligible[np.argsort(fit.power[eligible])[::-1]][:15]:
        period, duration = float(fit.period[index]), float(fit.duration[index])
        epoch = float(time[0] + fit.transit_time[index])
        eclipses = predicted_eclipses(data, period, epoch, duration)
        if passes(eclipses):
            return {"spacing_days": period, "duration_days": duration, "epoch_btjd": epoch,
                    "bls_depth_snr": float(fit.depth_snr[index]), "discovery_eclipses": eclipses}
    return None


def evaluate_target(tic, first, second, cache, search):
    filename = first[tic]
    source = product_url(filename)
    sparse = search == "sparse"
    try:
        meta, data = photometry(download(source, cache / "fits" / filename, verify_fits=True), native=sparse)
        if data is None:
            return "insufficient_discovery", None
        if meta["crowdsap"] < 0.60:
            return "crowded_discovery", None
        if sparse:
            candidates = tess_sparse.ephemerides(tess_sparse.dip_events(data[0], data[1]), data[0])
        else:
            eclipse = find_eclipse(data)
            candidates = [eclipse] if eclipse else []
        if not candidates:
            return "no_discovery_eclipses", None
        target = {"tic": tic, **meta, "discovery_source": source, **candidates[0]}
        other = second[tic]
        validation_source = product_url(other)
        validation_meta, validation = photometry(
            download(validation_source, cache / "fits" / other, verify_fits=True), native=sparse)
        if validation is None or validation_meta["crowdsap"] < 0.60:
            return "insufficient_validation", target
        if sparse:
            target.update(validation_source=validation_source,
                          validation_crowdsap=validation_meta["crowdsap"])
            events = tess_sparse.dip_events(validation[0], validation[1])
            for candidate in candidates:
                matches = tess_sparse.match_ephemeris(events, candidate)
                coverage = tess_sparse.covered_predictions(
                    validation[0], candidate["spacing_days"],
                    candidate["epoch_btjd"], candidate["duration_days"])
                if len(matches) >= 2 and coverage - len(matches) <= 1:
                    target.update(candidate, validation_events=matches, validation_predictions_covered=coverage)
                    return "repeated", target
            return "failed_validation", target
        check = predicted_eclipses(validation, target["spacing_days"],
                                   target["epoch_btjd"], target["duration_days"])
        target.update(validation_source=validation_source, validation_eclipses=check,
                      validation_crowdsap=validation_meta["crowdsap"])
        return ("repeated" if passes(check) else "failed_validation"), target
    except (OSError, ValueError, requests.RequestException) as error:
        return "error", {"tic": tic, "discovery_source": source, "message": str(error)}


def run_pair(discovery, validation, cache, *, offset=0, limit=384, workers=8, seed=196,
             prioritize_new_since=(), search="eclipses"):
    if not (0 < discovery < validation and offset >= 0 and limit > 0 and 1 <= workers <= 32):
        raise ValueError("Use increasing sector numbers and positive workers and limit.")
    started = monotonic()
    first, second = (sector_products(sector, cache) for sector in (discovery, validation))
    common = sorted(first.keys() & second.keys())
    random.Random(seed).shuffle(common)
    earlier = set().union(*(sector_products(sector, cache) for sector in prioritize_new_since))
    common.sort(key=earlier.__contains__)
    counts = Counter()
    candidates = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = pool.map(lambda tic: evaluate_target(tic, first, second, cache, search),
                        common[offset:offset + limit], buffersize=workers * 2)
        for status, candidate in jobs:
            counts[status] += 1
            if candidate:
                candidates.append({**candidate, "status": status})
    if search == "sparse":
        candidates.sort(key=lambda item: sum(event["robust_score"] for event in item.get("validation_events", [])),
                        reverse=True)
    else:
        candidates.sort(key=lambda item: item.get("validation_eclipses", {}).get("depth_sigma", -1),
                        reverse=True)
    return {"status": "partial" if counts["error"] else "ok", "discovery_sector": discovery,
            "validation_sector": validation, "common_targets": len(common),
            "seed": seed, "offset": offset, "limit": limit, "search": search,
            "prioritized_absence_in_sectors": list(prioritize_new_since),
            "coverage": dict(counts), "elapsed_seconds": monotonic() - started,
            "candidates": candidates}
