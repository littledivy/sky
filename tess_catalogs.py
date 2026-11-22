"""Cross-match the repeated signals in a Sky TESS scan with published TIC tables."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import re
from urllib.parse import urlencode, urljoin

import requests

from sky.tess import download


CATALOGS = {
    "ijspeert-2024": ("J/A+A/691/A242/obaf-eb1", 69057),
    "prsa-2022": ("J/ApJS/258/16/tess-ebs", 4584),
    "kostov-2025-new": ("J/ApJS/279/50/table3", 7936),
    "kostov-2025-known": ("J/ApJS/279/50/table4", 2065),
    "olmschenk-2024": ("J/AJ/168/83", 14156),
}
TCE_INDEX = "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"


def catalog(cache, name, table, minimum_rows):
    query = urlencode({"-source": table, "-out": "TIC", "-out.max": "unlimited"})
    url = f"https://vizier.idia.ac.za/viz-bin/asu-tsv?{query}"
    path = download(url, cache / "catalogs" / f"{name}.tsv")
    ids = [int(line) for row in path.read_text().splitlines() if (line := row.strip()).isdigit()]
    if len(ids) < minimum_rows:
        raise ValueError(f"Incomplete TIC table {table}: {len(ids)} rows at {path}")
    return name, {"table": table, "url": url, "rows": len(ids)}, set(ids)


def talk_search(project, section, query):
    response = requests.get("https://talk.zooniverse.org/searches",
                            params={"section": section, "query": query}, timeout=15)
    response.raise_for_status()
    payload = response.json()
    listings = payload["searches"]
    matches = {}
    for row in listings:
        if row["type"] == "Discussion":
            discussion_id, title = row["id"], row["title"]
        elif row["type"] == "Comment":
            discussion_id, title = row["discussion_id"], row["discussion_title"]
        else:
            continue
        matches[discussion_id] = {
            "title": title,
            "url": f"https://www.zooniverse.org/projects/{project}/talk/{row['board_id']}/{discussion_id}",
        }
    return {"query": query, "search_url": response.url,
            "total_hits": payload["meta"]["searches"]["count"], "threads_on_first_page": list(matches.values())}


def tce_tables(cache, sectors, tics):
    response = requests.get(TCE_INDEX, timeout=20)
    response.raise_for_status()
    paths = set(re.findall(r'href="([^"]+_dvr-tcestats\.csv)"', response.text))
    urls = {}
    for sector in sorted(set(sectors)):
        suffix = f"-s{sector:04d}-s{sector:04d}_dvr-tcestats.csv"
        matches = [urljoin(TCE_INDEX, path) for path in paths if path.endswith(suffix)]
        urls[sector] = sorted(matches)

    def scan(url):
        path = download(url, cache / url.rsplit("/", 1)[-1])
        with path.open(newline="") as source:
            reader = csv.DictReader(line for line in source if not line.startswith("#"))
            required = {"ticid", "tceid", "tce_period", "tce_time0bt"}
            if not required.issubset(reader.fieldnames or ()):
                raise ValueError(f"Missing TCE columns in {url}")
            return [{"tic": int(row["ticid"]), "event_id": row["tceid"],
                     "period_days": float(row["tce_period"]) if row["tce_period"] else None,
                     "epoch_btjd": float(row["tce_time0bt"]) if row["tce_time0bt"] else None}
                    for row in reader if int(row["ticid"]) in tics]

    with ThreadPoolExecutor(max_workers=8) as pool:
        distinct_urls = list(dict.fromkeys(url for sector_urls in urls.values() for url in sector_urls))
        tables = dict(zip(distinct_urls, pool.map(scan, distinct_urls)))
    return {sector: [{"url": url, "matches": tables[url]} for url in sector_urls]
            for sector, sector_urls in urls.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--cache", type=Path, default=Path("data/tess"))
    parser.add_argument("--source-tic", type=int, action="append", default=[],
                        help="Also check a source identified by pixel localization")
    parser.add_argument("--talk-query", action="append", default=[],
                        help="Search Planet Hunters TESS and Eclipsing Binary Patrol Talk")
    parser.add_argument("--tce-sector", action="append", type=int, default=[],
                        help="Check the public SPOC threshold-crossing events for this single sector")
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
    sources = [{"tic": tic, "catalog_matches": [name for name, _, ids in index if tic in ids]}
               for tic in args.source_tic]
    projects = {"planet_hunters_tess": ("nora-dot-eisner/planet-hunters-tess", "project-7929"),
                "eclipsing_binary_patrol": ("vbkostov/eclipsing-binary-patrol", "project-22939")}
    queries = dict.fromkeys(args.talk_query)
    with ThreadPoolExecutor(max_workers=4) as pool:
        talk = {key: list(pool.map(lambda query: talk_search(*project, query), queries))
                for key, project in projects.items()}
    tics = {row["tic"] for row in results + sources}
    tces = tce_tables(args.cache, args.tce_sector, tics) if args.tce_sector else {}
    print(json.dumps({"catalogs": {name: meta for name, meta, _ in index}, "results": results,
                      "localized_sources": sources, "zooniverse_talk": talk, "spoc_single_sector_tces": tces},
                     indent=2))


if __name__ == "__main__":
    main()
