"""Plan, execute, and merge independent ZTF DR24 discovery shards."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import csv
from functools import lru_cache
import hashlib
import heapq
from io import StringIO
from itertools import repeat
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
from time import monotonic

from pyarrow.fs import FileSelector, S3FileSystem
from pyarrow.parquet import ParquetFile

from ztf.short_period_eclipses import Search, scan_rows

from .archives import request


CATALOG = "ipac-irsa-ztf/ztf/enhanced/dr24/lc/hats/ztf_dr24_lc-hats"
MANIFEST = "https://ipac-irsa-ztf.s3.us-east-1.amazonaws.com/" + CATALOG.split("/", 1)[1] + "/partition_info.csv"
REGION = "us-east-1"
LARGE_TILE = 1536 * (1 << 20)
ROW_GROUP_TARGET = 512 * (1 << 20)
COLUMNS = ("objectid", "objra", "objdec", "filterid", "lightcurve.hmjd",
           "lightcurve.mag", "lightcurve.magerr", "lightcurve.catflags")


@lru_cache(maxsize=1)
def public_s3():
    return S3FileSystem(anonymous=True, region=REGION, connect_timeout=10, request_timeout=60)


def digest(value):
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as temp:
        temporary = Path(temp.name)
        try:
            json.dump(value, temp, allow_nan=False, separators=(",", ":"))
            temp.flush()
        except BaseException:
            temporary.unlink()
            raise
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def split_large_tile(tile):
    if tile["bytes"] <= LARGE_TILE:
        return [tile]
    fs = public_s3()
    pieces = []
    for file in tile["files"]:
        with ParquetFile(file["path"], filesystem=fs) as parquet:
            groups = []
            byte_count = row_count = 0
            for index in range(parquet.metadata.num_row_groups):
                group = parquet.metadata.row_group(index)
                size = sum(group.column(column).total_compressed_size
                           for column in range(group.num_columns))
                if groups and byte_count + size > ROW_GROUP_TARGET:
                    pieces.append({"rows": row_count, "bytes": byte_count,
                                   "files": [{**file, "row_groups": groups}]})
                    groups, byte_count, row_count = [], 0, 0
                groups.append(index)
                byte_count += size
                row_count += group.num_rows
            if groups:
                pieces.append({"rows": row_count, "bytes": byte_count,
                               "files": [{**file, "row_groups": groups}]})
    if sum(piece["rows"] for piece in pieces) != tile["rows"]:
        raise ValueError(f"Parquet metadata disagrees with manifest tile {tile['id']}.")
    return [{**piece, "id": f"{tile['id']}.{index}", "tile_id": tile["id"]}
            for index, piece in enumerate(pieces)]


def create_plan(shards, top=20, sample=None):
    if shards < 1 or not 1 <= top <= 100 or (sample is not None and sample < 1):
        raise ValueError("Use positive shard and sample counts and a top count from 1 to 100.")
    with request(MANIFEST) as response:
        records = list(csv.DictReader(StringIO(response.read().decode("utf-8"))))
    tiles = {}
    for row in records:
        order, directory, pixel, count = (int(row[key]) for key in ("Norder", "Dir", "Npix", "num_rows"))
        path = f"{CATALOG}/dataset/Norder={order}/Dir={directory}/Npix={pixel}"
        if path in tiles or count < 1:
            raise ValueError(f"Invalid or duplicate tile in the archive manifest: {path}")
        tiles[path] = {"id": f"{order}-{pixel}", "rows": count, "files": [], "bytes": 0}
    if not tiles:
        raise ValueError("The archive returned an empty partition manifest.")
    selected = (dict(sorted(tiles.items(), key=lambda item: (item[1]["rows"], item[0]))[:sample])
                if sample is not None else tiles)
    fs = public_s3()
    selectors = (
        [FileSelector(path) for path in selected]
        if sample is not None else [FileSelector(f"{CATALOG}/dataset", recursive=True)]
    )
    for selector in selectors:
        for file in fs.get_file_info(selector):
            if not file.is_file or not file.path.endswith(".parquet"):
                continue
            tile_path, _, _ = file.path.rpartition("/")
            if tile_path not in selected:
                raise ValueError(f"Parquet file is absent from the partition manifest: {file.path}")
            tile = selected[tile_path]
            tile["files"].append({"path": file.path, "bytes": file.size})
            tile["bytes"] += file.size
    for path, tile in selected.items():
        if not tile["files"]:
            raise ValueError(f"No Parquet files for manifest tile {path}")
        tile["files"].sort(key=lambda file: file["path"])

    # Compressed bytes are a proxy for the number of observations to process.
    with ThreadPoolExecutor(max_workers=16) as pool:
        pieces = [piece for group in pool.map(split_large_tile, selected.values()) for piece in group]
    work = [[] for _ in range(min(shards, len(pieces)))]
    pending = [(0, 0, index) for index in range(len(work))]
    heapq.heapify(pending)
    for tile in sorted(pieces, key=lambda tile: (-tile["bytes"], tile["id"])):
        size, rows, index = heapq.heappop(pending)
        work[index].append(tile)
        heapq.heappush(pending, (size + tile["bytes"], rows + tile["rows"], index))

    payload = {
        "schema": 1, "catalog": CATALOG, "region": REGION, "manifest": MANIFEST,
        "sample": sample is not None, "tiles": len(selected), "rows": sum(tile["rows"] for tile in selected.values()),
        "bytes": sum(tile["bytes"] for tile in selected.values()),
        "settings": vars(Search(5, 60, 60, 51, top, None)), "shards": work,
    }
    return {"plan_id": digest(payload), **payload}


def load_plan(path):
    plan = json.loads(path.read_text())
    plan_id = plan.pop("plan_id", None)
    if plan.get("schema") != 1 or plan.get("catalog") != CATALOG or digest(plan) != plan_id:
        raise ValueError("Invalid ZTF plan or changed manifest.")
    return {"plan_id": plan_id, **plan}


def tile_rows(tile):
    fs = public_s3()
    for file in tile["files"]:
        with ParquetFile(file["path"], filesystem=fs, pre_buffer=True) as parquet:
            groups = file.get("row_groups", range(parquet.metadata.num_row_groups))
            for group in groups:
                for batch in parquet.iter_batches(batch_size=256, columns=COLUMNS, use_threads=False,
                                                  row_groups=[group]):
                    for row in batch.to_pylist():
                        row.update(row.pop("lightcurve"))
                        row["_file"] = f"s3://{file['path']}"
                        yield row


def scan_tile(tile, settings):
    started = monotonic()
    result = scan_rows(tile_rows(tile), tile["rows"], tile["id"], settings)
    examined = result["coverage"]["rows_examined"]
    if "error" not in result and examined != tile["rows"]:
        result["error"] = f"Manifest expected {tile['rows']} rows; scanned {examined}."
    result["seconds"] = monotonic() - started
    return result


def checkpoint_path(directory, plan, tile):
    return directory / plan["plan_id"] / f"{tile['id']}.json"


def read_checkpoint(directory, plan, tile):
    path = checkpoint_path(directory, plan, tile)
    try:
        data = json.loads(path.read_text())
    except (FileNotFoundError, ValueError):
        return None
    if (not isinstance(data, dict) or data.get("plan_id") != plan["plan_id"] or data.get("partition") != tile["id"]
            or "error" in data or data.get("coverage", {}).get("rows_examined") != tile["rows"]):
        return None
    return data


def scan_shard(plan, index, directory, workers=1):
    if not 0 <= index < len(plan["shards"]) or workers < 1:
        raise ValueError("Invalid shard index or worker count.")
    tiles = plan["shards"][index]
    pending = [tile for tile in tiles if read_checkpoint(directory, plan, tile) is None]
    failures = []
    settings = Search(**plan["settings"])

    def save(results):
        for tile, result in zip(pending, results, strict=True):
            write_json(checkpoint_path(directory, plan, tile), {"plan_id": plan["plan_id"], **result})
            if "error" in result:
                failures.append({"tile": tile["id"], "error": result["error"]})

    if workers == 1 or len(pending) < 2:
        save(scan_tile(tile, settings) for tile in pending)
    else:
        with ProcessPoolExecutor(max_workers=min(workers, len(pending))) as pool:
            save(pool.map(scan_tile, pending, repeat(settings), buffersize=workers * 2))
    return {"plan_id": plan["plan_id"], "shard": index, "status": "partial" if failures else "ok",
            "tiles": len(tiles), "reused": len(tiles) - len(pending), "failures": failures}


def merge_shards(plan, directory):
    missing = []
    coverage = Counter()
    candidates = []
    sequence = 0
    for shard in plan["shards"]:
        for tile in shard:
            result = read_checkpoint(directory, plan, tile)
            if result is None:
                missing.append(tile["id"])
                continue
            coverage.update(result["coverage"])
            for candidate in result["candidates"]:
                sequence += 1
                score = candidate["holdout_delta_chi2"]
                heapq.heappush(candidates, ((score, -sequence), {**candidate, "tile": tile.get("tile_id", tile["id"])}))
                if len(candidates) > plan["settings"]["top"]:
                    heapq.heappop(candidates)
    return {"plan_id": plan["plan_id"], "status": "partial" if missing else "ok",
            "sample": plan["sample"], "expected_rows": plan["rows"], "coverage": dict(coverage),
            "missing_tiles": missing, "candidates": [candidate for _, candidate in sorted(candidates, reverse=True)]}
