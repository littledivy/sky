"""Command-line entry point for sky."""

import argparse
import importlib.util
import json
from pathlib import Path
import sys

from .api import RunContext, archive_call, cone, local_partitions, orbital_decay, survey_names
from .archives import SURVEYS, fetch_ztf, gaia_epoch_photometry


def trusted_experiment(name):
    if name == "orbital-decay":
        return orbital_decay
    path = Path(name)
    if path.suffix != ".py" or not path.is_file():
        raise ValueError("Use orbital-decay or the path to a local Python experiment.")
    spec = importlib.util.spec_from_file_location("sky_user_experiment", path)
    if not spec or not spec.loader:
        raise ValueError(f"Cannot load {path}.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    run = getattr(module, "run", None)
    if not callable(run):
        raise ValueError(f"{path} must export a callable run(context).")
    return run


def main():
    parser = argparse.ArgumentParser(prog="sky", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("surveys", help="List survey adapters and products")

    search = commands.add_parser("cone", help="Search Gaia or Roman by position")
    search.add_argument("ra", type=float)
    search.add_argument("dec", type=float)
    search.add_argument("--radius-arcsec", type=float, default=60)
    search.add_argument("--limit", type=int, default=20)
    search.add_argument("--surveys", default="gaia,roman")
    search.add_argument("--roman-program", choices=("gbtds", "hlwas"), default="gbtds")

    gaia = commands.add_parser("gaia-photometry", help="Retrieve Gaia DR3 epoch photometry")
    gaia.add_argument("source_id")

    download = commands.add_parser("fetch-ztf", help="Download one official ZTF DR24 Parquet file")
    download.add_argument("field", type=int)
    download.add_argument("ccd", type=int)
    download.add_argument("quadrant", type=int)
    download.add_argument("band", choices=("g", "r", "i"))
    download.add_argument("--directory", type=Path, default=Path("data/ztf"))

    run = commands.add_parser("run", help="Run orbital-decay or a trusted local Python experiment")
    run.add_argument("experiment")
    run.add_argument("--ztf", nargs="+", type=Path, required=True, help="Parquet files or directories")
    run.add_argument("--surveys", default="ztf,gaia,roman")
    run.add_argument("--top", type=int, default=20)
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--limit-per-partition", type=int, help="Limit eligible curves for a pilot run")
    run.add_argument("--roman-photometry", type=Path, help="CSV of extracted, barycentric photometry")

    tess = commands.add_parser("tess", help="Search TESS eclipses and test them in another sector")
    tess.add_argument("discovery", type=int, help="Sector used to find the ephemeris")
    tess.add_argument("validation", type=int, help="Later sector used for validation")
    tess.add_argument("--cache", type=Path, default=Path("data/tess"))
    tess.add_argument("--offset", type=int, default=0)
    tess.add_argument("--limit", type=int, default=384)
    tess.add_argument("--workers", type=int, default=8)
    tess.add_argument("--prioritize-new-since", type=int, nargs="*", default=[])
    tess.add_argument("--search", choices=("eclipses", "sparse"), default="eclipses")

    survey = commands.add_parser("hats", help="Plan and scan the public ZTF DR24 S3 archive")
    actions = survey.add_subparsers(dest="hats_command", required=True)
    plan = actions.add_parser("plan", help="Build a manifest of independently runnable shards")
    plan.add_argument("file", type=Path)
    plan.add_argument("--shards", type=int, default=4096)
    plan.add_argument("--top", type=int, default=20)
    plan.add_argument("--sample", type=int, help="Use the N smallest tiles for a smoke test")
    scan = actions.add_parser("scan", help="Run or resume one shard")
    scan.add_argument("file", type=Path)
    scan.add_argument("index", type=int, help="Zero-based shard index")
    scan.add_argument("--output", type=Path, required=True, help="Checkpoint directory")
    scan.add_argument("--workers", type=int, default=1)
    merge = actions.add_parser("merge", help="Check coverage and merge ranked candidates")
    merge.add_argument("file", type=Path)
    merge.add_argument("--output", type=Path, required=True, help="Checkpoint directory")

    server = commands.add_parser("serve", help="Start the local HTTP API")
    server.add_argument("--port", type=int, default=8765)
    server.add_argument("--workers", type=int, default=1, help="Parallel ZTF scan workers per run")
    server.add_argument("--ztf", nargs="+", type=Path, help="Enable experiment runs with these inputs")

    args = parser.parse_args()
    try:
        match args.command:
            case "surveys":
                output = SURVEYS
            case "cone":
                output = cone(args.ra, args.dec, radius_arcsec=args.radius_arcsec, limit=args.limit,
                              surveys=survey_names(args.surveys, ("gaia", "roman", "roman-sim")),
                              program=args.roman_program)
            case "gaia-photometry":
                output = archive_call("gaia", gaia_epoch_photometry, args.source_id)
            case "fetch-ztf":
                output = fetch_ztf(args.field, args.ccd, args.quadrant, args.band, args.directory)
            case "run":
                context = RunContext(local_partitions(args.ztf), survey_names(args.surveys),
                                     args.top, args.workers, args.limit_per_partition,
                                     args.roman_photometry)
                output = trusted_experiment(args.experiment)(context)
            case "tess":
                from .tess import run_pair

                output = run_pair(
                    args.discovery, args.validation, args.cache, offset=args.offset,
                    limit=args.limit, workers=args.workers,
                    prioritize_new_since=tuple(args.prioritize_new_since),
                    search=args.search,
                )
            case "hats":
                from .hats import create_plan, load_plan, merge_shards, scan_shard, write_json

                match args.hats_command:
                    case "plan":
                        plan = create_plan(args.shards, args.top, args.sample)
                        write_json(args.file, plan)
                        output = {"path": str(args.file), "plan_id": plan["plan_id"],
                                  "region": plan["region"], "shards": len(plan["shards"]),
                                  "sample": plan["sample"], "tiles": plan["tiles"],
                                  "rows": plan["rows"], "bytes": plan["bytes"]}
                    case "scan":
                        output = scan_shard(load_plan(args.file), args.index, args.output, args.workers)
                    case "merge":
                        output = merge_shards(load_plan(args.file), args.output)
            case "serve":
                import uvicorn
                from .server import create_app

                if args.workers < 1:
                    parser.error("--workers must be positive")
                partitions = local_partitions(args.ztf) if args.ztf else ()
                uvicorn.run(create_app(partitions, workers=args.workers), host="127.0.0.1", port=args.port)
                return 0
        print(json.dumps(output, indent=2, allow_nan=False))
    except (ValueError, OSError) as error:
        print(f"sky: {error}", file=sys.stderr)
        return 1
    unavailable = ("error", "partial")
    if isinstance(output, dict):
        return int(output.get("status") in unavailable or
                   (args.command == "cone" and any(value["status"] in unavailable for value in output.values())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
