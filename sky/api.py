"""Shared Python API for CLI commands and agent-authored experiments."""

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import csv
from dataclasses import dataclass
from itertools import repeat
import math
from pathlib import Path

import numpy as np

from ztf.fit_orbital_drift import fit_drift, read_photometry
from ztf.short_period_eclipses import Search, scan_partition

from .archives import (
    SURVEYS, gaia_epoch_photometry, gaia_sources, roman_observations,
    roman_simulations, validate_region,
)


def survey_names(value, allowed=SURVEYS):
    names = tuple(dict.fromkeys(name.strip() for name in value.split(",") if name.strip()))
    if not names or any(name not in allowed for name in names):
        raise ValueError(f"Choose surveys from {', '.join(allowed)}.")
    return names


def archive_call(name, operation, *args):
    try:
        return operation(*args)
    except Exception as error:
        return {"status": "error", "source": SURVEYS[name]["source"], "error": str(error), "rows": []}


def cone(ra, dec, *, radius_arcsec=60, limit=20, surveys=("gaia", "roman"), program="gbtds"):
    validate_region(ra, dec, radius_arcsec, limit)
    operations = {"gaia": gaia_sources, "roman": roman_observations,
                  "roman-sim": roman_simulations}
    if not surveys or any(name not in operations for name in surveys):
        raise ValueError("Cone searches support gaia, roman and roman-sim.")
    with ThreadPoolExecutor(max_workers=len(surveys)) as pool:
        pending = {
            name: pool.submit(archive_call, name, operations[name], ra, dec, radius_arcsec, limit,
                              *((program,) if name == "roman-sim" else ()))
            for name in surveys
        }
        return {name: future.result() for name, future in pending.items()}


@dataclass(frozen=True)
class RunContext:
    partitions: tuple[Path, ...]
    surveys: tuple[str, ...] = ("ztf", "gaia", "roman")
    top: int = 20
    workers: int = 1
    limit_per_partition: int | None = None
    roman_photometry: Path | None = None

    def scan_eclipses(self, *, min_period=5, max_period=60, window_days=60, min_epochs=51):
        settings = Search(min_period, max_period, window_days, min_epochs, self.top,
                          self.limit_per_partition)
        if self.workers == 1:
            return [scan_partition(path, settings) for path in self.partitions]
        with ProcessPoolExecutor(max_workers=self.workers) as pool:
            return list(pool.map(scan_partition, self.partitions, repeat(settings),
                                 buffersize=2 * self.workers))

    def cone(self, ra, dec, radius_arcsec=2, limit=5):
        surveys = tuple(name for name in self.surveys if name != "ztf")
        return cone(ra, dec, radius_arcsec=radius_arcsec, limit=limit, surveys=surveys) if surveys else {}


def local_partitions(paths):
    files = set()
    for path in map(Path, paths):
        if path.is_dir():
            matches = list(path.rglob("*.parquet"))
            if not matches:
                raise ValueError(f"No Parquet files in {path}")
            files.update(file.resolve() for file in matches)
        elif path.is_file():
            files.add(path.resolve())
        else:
            raise ValueError(f"Input path does not exist: {path}")
    if not files:
        raise ValueError("Provide at least one ZTF Parquet input.")
    return tuple(sorted(files))


def read_roman_photometry(path):
    if path is None:
        return {}
    sources = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"ztf_objectid", "band", "bmjd_tdb", "magnitude", "magnitude_error", "dataset_id"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"Roman photometry CSV requires {', '.join(sorted(required))}.")
        for row in reader:
            band, dataset = row["band"].strip(), row["dataset_id"].strip()
            time, mag, err = (float(row[key]) for key in ("bmjd_tdb", "magnitude", "magnitude_error"))
            if not band or not dataset or not all(map(math.isfinite, (time, mag, err))) or err <= 0:
                raise ValueError(f"Invalid Roman photometry on CSV line {reader.line_num}.")
            sources.setdefault(int(row["ztf_objectid"]), {}).setdefault(band, []).append(
                {"bmjd_tdb": time, "magnitude": mag, "magnitude_error": err, "dataset_id": dataset})
    return sources


def photometry_array(rows):
    return tuple(np.array([row[key] for row in rows])
                 for key in ("bmjd_tdb", "magnitude", "magnitude_error"))


def orbital_decay(context):
    if "ztf" not in context.surveys:
        raise ValueError("The orbital-decay experiment uses ZTF for its discovery scan.")
    if context.roman_photometry and "roman" not in context.surveys:
        raise ValueError("Enable roman when supplying Roman photometry.")
    if not (1 <= context.top <= 100 and context.workers >= 1
            and (context.limit_per_partition is None or context.limit_per_partition >= 1)):
        raise ValueError("Invalid top, worker count, or partition limit.")
    roman_curves = read_roman_photometry(context.roman_photometry)
    partitions = context.scan_eclipses()
    ranking = sorted(
        ((row, partition["partition"]) for partition in partitions for row in partition["candidates"]),
        key=lambda item: item[0]["holdout_delta_chi2"], reverse=True,
    )
    outcomes = []
    seen = set()
    incomplete = any("error" in entry for entry in partitions)
    for candidate, path in ranking:
        oid = candidate["objectid"]
        if oid in seen:
            continue
        seen.add(oid)
        ra, dec = candidate["ra_deg"], candidate["dec_deg"]
        response = context.cone(ra, dec)
        evidence = []
        data = []
        if gaia := response.get("gaia"):
            if len(gaia["rows"]) == 1 and not gaia["truncated"]:
                match = gaia["rows"][0]
                gaia["match"] = "single_positional_candidate"
                if match["has_epoch_photometry"]:
                    epoch = archive_call("gaia", gaia_epoch_photometry, match["source_id"])
                    gaia["epoch_photometry"] = {**epoch, "rows": len(epoch["rows"])}
                    for band in sorted({row["band"] for row in epoch["rows"]}):
                        rows = [row for row in epoch["rows"] if row["band"] == band]
                        used = len(rows) >= 10
                        evidence.append({"survey": "gaia", "band": band, "epochs": len(rows),
                                         "used": used, "source_id": match["source_id"]})
                        if used:
                            data.append(photometry_array(rows))
                    incomplete |= epoch["status"] == "error"
                else:
                    gaia["epoch_photometry"] = {"status": "empty", "reason": "no_epoch_photometry"}
            elif gaia["rows"]:
                gaia["match"] = "ambiguous"
        for band, rows in sorted(roman_curves.get(oid, {}).items()):
            used = len(rows) >= 10
            evidence.append({"survey": "roman", "band": band, "epochs": len(rows), "used": used,
                             "dataset_ids": sorted({row["dataset_id"] for row in rows})})
            if used:
                data.append(photometry_array(rows))
        incomplete |= any(archive["status"] == "error" for archive in response.values())
        try:
            ztf_photometry = read_photometry(f"{path}:{oid}")
            evidence.insert(0, {"survey": "ztf", "filterid": candidate["filterid"],
                                "epochs": len(ztf_photometry[0]), "used": True, "path": path})
            timing = {"status": "ok", **fit_drift([ztf_photometry, *data], candidate["period_minutes"])}
            if not all(math.isfinite(value) for value in timing.values()
                       if isinstance(value, (int, float))):
                raise ValueError("The timing fit produced nonfinite statistics.")
        except (ValueError, np.linalg.LinAlgError) as error:
            timing = {"status": "fit_unavailable", "reason": str(error)}
        outcomes.append({**candidate, "discovery_rank": len(outcomes) + 1,
                         "archives": response, "photometry": evidence, "timing": timing})
        if len(outcomes) == context.top:
            break
    outcomes.sort(key=lambda row: row["timing"].get("delta_chi2", -math.inf), reverse=True)
    return {
        "experiment": "orbital-decay", "status": "partial" if incomplete else "ok",
        "surveys": list(context.surveys),
        "coverage": [{key: value for key, value in partition.items() if key != "candidates"}
                     for partition in partitions],
        "follow_up": {"rule": "highest ZTF holdout scores", "limit": context.top,
                      "discovery_curve_limit_per_partition": context.limit_per_partition,
                      "ranking": "timing delta chi-squared"},
        "candidates": outcomes,
    }
