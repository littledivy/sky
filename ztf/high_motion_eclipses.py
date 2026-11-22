"""Search ZTF light curves of Gaia sources with large proper motions."""

from collections import Counter, defaultdict
from itertools import chain
from time import monotonic

from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from pyarrow.parquet import ParquetFile, read_table
from scipy.spatial import cKDTree

from sky.archives import tap_csv
from ztf.short_period_eclipses import COLUMNS, Search, search_curve


def alternate_signal(candidate, row):
    t, y, error, flags = (np.asarray(row[column]) for column in ("hmjd", "mag", "magerr", "catflags"))
    start = candidate["window_start_hmjd"]
    valid = np.isfinite(t) & np.isfinite(y) & np.isfinite(error) & (error > 0)
    valid &= ((flags & 0x8000) == 0) & (t >= start) & (t <= start + 60)
    t, y, error = t[valid], y[valid], error[valid]
    period = candidate["period_minutes"] / 1440
    offset = t - candidate["eclipse_epoch_hmjd"]
    cycles = np.rint(offset / period)
    eclipse = abs(offset - period * cycles) < candidate["eclipse_duration_seconds"] / 172800
    count = len(np.unique(cycles[eclipse]))
    if count < 2 or (~eclipse).sum() < 10:
        return {"status": "insufficient_coverage", "epochs": len(t), "eclipses": count}
    baseline = y[~eclipse]
    noise = 1.4826 * np.median(abs(baseline - np.median(baseline)))
    weights = 1 / np.maximum(noise, error)**2
    inside, outside = weights[eclipse], weights[~eclipse]
    depth = np.average(y[eclipse], weights=inside) - np.average(baseline, weights=outside)
    scale = np.sqrt(1/inside.sum() + 1/outside.sum())
    return {"status": "tested", "epochs": len(t), "eclipses": count, "depth_mag": float(depth),
            "nominal_depth_sigma": float(depth / scale)}


def gaia_sample(ra, dec, radius):
    query = f"""SELECT TOP 501 source_id,ra,dec,ref_epoch,pmra,pmdec,phot_g_mean_mag,bp_rp
      FROM gaiadr3.gaia_source
      WHERE CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{ra},{dec},{radius}))=1
      AND (pmra*pmra+pmdec*pmdec)>=22500
      AND phot_g_mean_mag BETWEEN 12 AND 21.5 ORDER BY source_id"""
    result = tap_csv(query, 501, ("source_id", "pmra", "pmdec"))
    if len(result) > 500:
        raise ValueError("Gaia sample was truncated")
    return result


def run_patch(paths):
    catalogs = {}
    footprints = []
    for path in paths:
        table = read_table(path, columns=["objectid", "objra", "objdec"])
        coords = SkyCoord(table["objra"].to_numpy(), table["objdec"].to_numpy(), unit="deg")
        footprints.append(coords)
        catalogs[path] = (cKDTree(coords.cartesian.xyz.value.T), table["objectid"].to_numpy())
    xyz = np.concatenate([coords.cartesian.xyz.value for coords in footprints], axis=1)
    center = SkyCoord(x=xyz[0].mean(), y=xyz[1].mean(), z=xyz[2].mean(),
                      unit=(u.one, u.one, u.one), representation_type="cartesian")
    radius = max(center.separation(coords).deg.max() for coords in footprints) + 0.015
    center = center.icrs.spherical
    gaia = gaia_sample(center.lon.deg, center.lat.deg, radius)
    counts = Counter()
    results = []
    settings = Search(5, 60, 60, 51, 20, None)
    tolerance = 0.85 * u.arcsec
    chord = 2 * np.sin(tolerance.to_value(u.rad) / 2)
    for source in gaia:
        counts["gaia_sample"] += 1
        pmra, pmdec = (float(source[key]) for key in ("pmra", "pmdec"))
        motion = np.hypot(pmra, pmdec)
        reference = SkyCoord(float(source["ra"]), float(source["dec"]), unit="deg")
        dates = np.linspace(2018, 2025.9, max(81, int(motion / 12)))
        track = reference.directional_offset_by(np.arctan2(pmra, pmdec) * u.rad,
                                                 (dates - float(source["ref_epoch"])) * motion * u.mas)
        matches = []
        for path, (tree, ids) in catalogs.items():
            indices = set(chain.from_iterable(tree.query_ball_point(track.cartesian.xyz.value.T, chord)))
            if len(indices) > 1:
                counts["ambiguous_band_matches"] += 1
            if len(indices) == 1:
                matches.append((path, int(ids[indices.pop()])))
        if not matches:
            continue
        counts["gaia_sources_with_unique_band_match"] += 1
        curves = []
        for path, oid in matches:
            row, = read_table(path, columns=list(COLUMNS),
                              filters=[("objectid", "=", oid)]).to_pylist()
            outcome, candidate = search_curve(row,settings)
            counts["curves_examined"] += 1
            counts[outcome] += 1
            curves.append((row, candidate))
        candidates = [candidate for _, candidate in curves if candidate]
        record = {"gaia": source, "proper_motion_mas_yr": float(motion),
                  "ztf_objects": [oid for _, oid in matches],
                  "candidates": []}
        for candidate in candidates:
            other = next((row for row, _ in curves if row["filterid"] != candidate["filterid"]), None)
            confirmation = {"status": "no_other_band"} if other is None else alternate_signal(candidate,other)
            record["candidates"].append({**candidate, "other_band": confirmation})
        results.append(record)
    return {"files": [str(path) for path in paths], "radius_deg": radius,
            "center_ra_deg":float(center.lon.deg), "center_dec_deg":float(center.lat.deg),
            "counts": dict(counts), "results": results}



def run(context):
    if not {"ztf", "gaia"}.issubset(context.surveys):
        raise ValueError("Enable ztf and gaia to run the high-motion eclipse search.")
    started = monotonic()
    patches = defaultdict(list)
    for path in context.partitions:
        with ParquetFile(path) as parquet:
            batch = next(parquet.iter_batches(batch_size=1, columns=["fieldid", "rcid"]), None)
        if batch is None:
            continue
        row, = batch.to_pylist()
        patches[row["fieldid"], row["rcid"]].append(path)
    results = []
    for (field, channel), paths in sorted(patches.items()):
        results.append({"field": field, "channel": channel, **run_patch(paths)})
    return {
        "experiment": "high-proper-motion-eclipses",
        "selection": {
            "minimum_proper_motion_mas_yr": 150,
            "gaia_g_magnitude_range": [12, 21.5],
            "position_track_years": [2018, 2025.9],
            "position_tolerance_arcsec": 0.85,
            "period_range_minutes": [5, 60],
            "confirmation_threshold_nominal_sigma": 3,
        },
        "seconds": monotonic() - started,
        "patches": results,
    }
