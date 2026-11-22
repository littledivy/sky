"""Search native-cadence TESS photometry for short, infrequent dips."""

import numpy as np
from scipy.signal import find_peaks


CADENCE_DAYS = 2 / 1440
WIDTHS = (2, 4, 8, 16, 32)
MIN_PERIOD = 2
MAX_PERIOD = 12


def dip_events(time, flux):
    """Find dips relative to a local baseline using five box widths."""
    sums = np.r_[0, np.cumsum(flux)]
    candidates = []
    for width in WIDTHS:
        shoulder = max(8, width)
        padding = shoulder + width
        center = np.arange(padding, len(time) - padding)
        start = center - width // 2
        end = start + width
        left = center - width
        right = center + width
        covered = time[right + shoulder - 1] - time[left - shoulder] < 2.25 * padding * CADENCE_DAYS
        baseline = (sums[left] - sums[left - shoulder] + sums[right + shoulder] - sums[right]) / (2 * shoulder)
        depth = baseline - (sums[end] - sums[start]) / width
        valid_depth = depth[covered]
        if len(valid_depth) < 100:
            continue
        middle = np.median(valid_depth)
        scatter = 1.4826 * np.median(abs(valid_depth - middle))
        if scatter <= 0:
            continue
        score = np.where(covered, (depth - middle) / scatter, 0)
        for i in find_peaks(score, height=4.2, distance=max(width, 4))[0]:
            candidates.append((time[center[i]], depth[i], score[i], width * CADENCE_DAYS))
    events = []
    for event in sorted(candidates, key=lambda event: event[2], reverse=True):
        if all(abs(event[0] - previous[0]) > max(event[3], previous[3]) for previous in events):
            events.append(event)
        if len(events) == 48:
            break
    return np.array(sorted(events)).reshape(-1, 4)


def covered_predictions(time, period, epoch, width):
    cycles = np.arange(np.ceil((time[0] - epoch) / period), np.floor((time[-1] - epoch) / period) + 1)
    predicted = epoch + period * cycles
    radius = max(0.01, width)
    before = np.searchsorted(time, predicted - radius)
    after = np.searchsorted(time, predicted + radius)
    return int(np.count_nonzero(after - before >= 1.5 * radius / CADENCE_DAYS))


def ephemerides(events, time):
    """Seed periods from event pairs and verify them against discovery events."""
    if len(events) < 3:
        return []
    times, depths, scores, widths = events.T
    first, last = np.triu_indices(len(events), 1)
    gaps = times[last] - times[first]
    cycles = np.arange(1, 1 + int(np.ptp(times) / MIN_PERIOD))
    periods = gaps[:, None] / cycles
    pair, alias = np.where((periods >= MIN_PERIOD) & (periods <= MAX_PERIOD))
    start = first[pair]
    periods = periods[pair, alias]
    if not len(periods):
        return []
    delta = times - times[start, None]
    labels = np.rint(delta / periods[:, None])
    residuals = abs(delta - labels * periods[:, None])
    matched = (residuals <= np.maximum(0.006, 0.5 * widths)) & (labels >= 0)
    ranks = (matched * np.minimum(scores, 12)).sum(axis=1)
    results = []
    for i in np.lexsort((-periods, -ranks)):
        selected = matched[i]
        labels_i = labels[i, selected]
        if len(np.unique(labels_i)) != len(labels_i) or len(labels_i) < 3:
            continue
        measured = times[selected]
        period, epoch = np.polynomial.polynomial.polyfit(labels_i, measured, 1)[::-1]
        width = float(np.median(widths[selected]))
        coverage = covered_predictions(time, period, epoch, width)
        if coverage - len(measured) > 1:
            continue
        if any(abs(period - result["spacing_days"]) < 0.005 and
               abs((epoch - result["epoch_btjd"] + period / 2) % period - period / 2) < 0.02
               for result in results):
            continue
        results.append({
            "spacing_days": float(period), "epoch_btjd": float(epoch), "duration_days": width,
            "discovery_events": summarize(events[selected]), "discovery_predictions_covered": coverage,
        })
        if len(results) == 8:
            break
    return results


def summarize(events):
    return [{"time_btjd": float(t), "depth_fraction": float(d), "robust_score": float(s),
             "duration_days": float(w)} for t, d, s, w in events]


def match_ephemeris(events, candidate):
    period = candidate["spacing_days"]
    delta = events[:, 0] - candidate["epoch_btjd"]
    cycles = np.rint(delta / period)
    residual = abs(delta - cycles * period)
    tolerance = max(0.013, 0.6 * candidate["duration_days"])
    depth = np.median([event["depth_fraction"] for event in candidate["discovery_events"]])
    selected = ((residual <= tolerance) & (events[:, 3] <= 3 * candidate["duration_days"])
                & (events[:, 1] >= depth / 3) & (events[:, 1] <= depth * 3))
    found = []
    for cycle in np.unique(cycles[selected]):
        group = events[selected & (cycles == cycle)]
        found.append(group[group[:, 2].argmax()])
    return summarize(found)
