"""Cross-match the repeated signals in a Sky TESS scan with published TIC tables."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from urllib.parse import urlencode

from sky.tess import download


CATALOGS = {
    "ijspeert-2024": ("J/A+A/691/A242/obaf-eb1", 69057),
    "prsa-2022": ("J/ApJS/258/16/tess-ebs", 4584),
    "kostov-2025-new": ("J/ApJS/279/50/table3", 7936),
    "kostov-2025-known": ("J/ApJS/279/50/table4", 2065),
    "olmschenk-2024": ("J/AJ/168/83", 14156),
}


def catalog(cache, name, table, minimum_rows):
    query = urlencode({"-source": table, "-out": "TIC", "-out.max": "unlimited"})
    url = f"https://vizier.idia.ac.za/viz-bin/asu-tsv?{query}"
    path = download(url, cache / "catalogs" / f"{name}.tsv")
    ids = [int(line) for row in path.read_text().splitlines() if (line := row.strip()).isdigit()]
    if len(ids) < minimum_rows:
        raise ValueError(f"Incomplete TIC table {table}: {len(ids)} rows at {path}")
    return name, {"table": table, "url": url, "rows": len(ids)}, set(ids)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--cache", type=Path, default=Path("data/tess"))
    args = parser.parse_args()
    with ThreadPoolExecutor(max_workers=len(CATALOGS)) as pool:
        index = list(pool.map(lambda entry: catalog(args.cache, entry[0], *entry[1]), CATALOGS.items()))
    results = []
    for path in args.reports:
        report = json.loads(path.read_text())
        for candidate in report["candidates"]:
            if candidate["status"] != "repeated":
                continue
            matches = [name for name, _, ids in index if candidate["tic"] in ids]
            results.append({"report": str(path), "tic": candidate["tic"],
                            "spacing_days": candidate["spacing_days"], "catalog_matches": matches})
    print(json.dumps({"catalogs": {name: meta for name, meta, _ in index}, "results": results}, indent=2))


if __name__ == "__main__":
    main()
