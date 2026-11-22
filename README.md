# Sky

Run astronomy experiments across ZTF, Gaia, Roman, and TESS. Sky provides a local CLI, a Python experiment interface, and a local HTTP API.

The [ZTF experiments](ztf/README.md) include the J1539 orbital-decay reproduction and a pilot search for eclipses around fast-moving stars. [Saved reports](output/ztf/) include ranked candidates, survey coverage, and timing fits.

From the repository root with Python 3.14:

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
sky surveys
```

Download ZTF DR24, query Gaia and Roman near J1539, and run the experiment:

```sh
sky fetch-ztf 760 14 3 g --directory data/ztf
sky cone 234.88396 50.4608 --radius-arcsec 2 --surveys gaia,roman
sky gaia-photometry 246756973652292992
sky run orbital-decay --ztf data/ztf --surveys ztf,gaia,roman --top 20 > run.json
```

For a smoke test, pass `--limit-per-partition 8 --top 1`. Omit the limit to
scan every eligible curve in the supplied files. After the ZTF scan, the
experiment takes the highest-ranked candidates under the same rule, queries
their positions, and attempts an orbital timing fit on each. The output
includes survey responses, partition coverage, discovery rank, timing
statistics and the photometry that went into each fit. Timing scores are
exploratory rankings.

Gaia sources come from ESA TAP. A unique positional candidate within two
arcseconds can contribute qualifying DR3 epoch photometry from the official
Gaia@AIP archive. Positional matches use the catalog coordinates at their
reference epoch and are not proper-motion propagated. Gaia epochs are converted
from barycentric TCB to barycentric TDB. ZTF timestamps also get converted to
barycentric TDB. The timing fit gives each instrument and band independent
amplitudes and a zero point.

Roman observation searches use MAST's public ROMAN collection and preserve
archive provenance, including test-product labels. An exposure is metadata.
To fit Roman observations, pass `--roman-photometry FILE` with already
extracted and vetted photometry. The CSV header is:

```csv
ztf_objectid,band,bmjd_tdb,magnitude,magnitude_error,dataset_id
```

Match the Roman photometry to the ZTF object ID before importing. Sky imports
these measurements by ID and band and records each Roman dataset ID. It does
not extract photometry from Roman images. Simulated exposure metadata is
available as an explicit separate survey. Results are matches to detector
centers and do not establish that a source falls on the detector.

```sh
sky cone 8.1 -43.7 --radius-arcsec 360 --surveys roman-sim --roman-program hlwas
```

## TESS repeatability search

```sh
sky tess 107 108 --cache data/tess --limit 1000 --workers 10 \
  --prioritize-new-since {1..106} --search sparse > tess-107-108.json
```

Sky intersects the two public SPOC two-minute sector manifests and picks a
reproducible target order. The optional prior-sector list puts targets absent
from those short-cadence manifests first. Use `--offset` and `--limit` to
divide a target list across runs with the same prior-sector list.

The default search uses an accelerated periodogram and a box search on
six-minute averages for 0.12–1.5 day recurrence. It downloads holdout
photometry only for targets with enough discoveries to make a prediction.
The sparse search instead uses native two-minute observations, tests dip
widths of 4–64 minutes, and looks for recurrence every 2–12 days. It needs
at least three events to estimate an ephemeris in the first sector and two
events at predicted times in the second. A period is penalized when it
predicts events in adequately sampled windows with no corresponding dip.
Both searches report the input product URLs, coverage counts, and timings.

The output is a shortlist for pixel localization and catalogue and literature
checks. Run `python experiments/tess/tess_catalogs.py tess-107-108.json` to match
repeaters against five published TIC tables, including the IJspeert et al.
eclipsing-binary table. For box-shaped repeaters, `sky.tess_pixels.localize` estimates the
position of the variable component from public target pixel files using the
ephemeris and paired in-event and out-of-event flux measurements.
Use `--source-tic` when the variable has a different TIC identifier from the
SPOC aperture target. The follow-up script `experiments/tess/tess_vet.py` can fit and localize
a second eclipse. `experiments/tess/tess_ztf_check.py` checks a fixed two-eclipse ephemeris in
resolved ZTF light curves.

The [Gaia DR3 5073603089619768704 analysis](output/tess/88481083/report.md)
includes the proposed VSX revision, phase and localization plots, TESS and ZTF
measurements, and the earlier Planet Hunters TESS report that the Sky analysis
independently reproduces. The [initial aperture fit](output/tess/88481081/fit.json)
records the signal in the light curve of its brighter neighbor.

The [HD 96135 recovery](output/tess/80845742/report.md) includes the scan,
independent-sector predictions, archival measurements, pixel checks, and
prior-publication match. [TIC 390471207](output/tess/390471207/catalog-checks.json)
also matched a previously published eclipsing signal.

The [HD 107126 analysis](output/tess/287893983/report.md) documents a 6.368-hour
photometric signal across seven years of TESS data and the checks for prior reports.

## Local HTTP API

```sh
sky serve --ztf data/ztf --workers 2
curl 'http://127.0.0.1:8765/v1/cone?ra=234.88396&dec=50.4608&surveys=gaia,roman'
curl -X POST 'http://127.0.0.1:8765/v1/experiments/orbital-decay' \
  -H 'Content-Type: application/json' \
  -d '{"surveys":"ztf,gaia,roman","top":1,"limit_per_partition":8}'
```

Use `--workers` on either `sky run` or `sky serve` to scan multiple input
files in parallel. Each light curve's Lomb–Scargle search uses one FINUFFT
thread so workers can use separate CPU cores. The box eclipse search and
holdout scoring follow the periodogram.

The server listens on loopback. OpenAPI docs are at
[127.0.0.1:8765/docs](http://127.0.0.1:8765/docs). The built-in experiment
endpoint uses only inputs configured at server start. Archive errors and empty
results have different status codes in each archive's result object. CLI
commands exit nonzero on archive errors while still returning available results.

## Local performance

On the M1 Pro, the first 96 eligible green-band curves in the J1539 quadrant
took 1.46 seconds with FINUFFT, down from 10.40 seconds with the previous
Astropy backend. Both commands produced identical output. Scanning the full
green and red files with eight workers took 24.0 seconds for 56,171 input
rows, of which 8,370 met the scan cuts. Both bands reproduced their previous
coverage counts and top 20 candidates. The joint timing fit also reproduced
the measured J1539 period derivative.

Extrapolating this patch's row throughput to the full DR24 release gives
about 26 days of local discovery scanning. This estimate depends on cadence
and source density and excludes input transfer and follow-up. The bulk
Parquet archive is 8.3 TB.

## Distributed ZTF discovery scan

The public DR24 HATS catalog on S3 in `us-east-1` contains 12,485 sky tiles,
5,401,134,101 light curves and 9.47 TB of Parquet. Sky can scan it directly.
A full plan records the input files, expected row counts, algorithm settings
and an assignment of work to shards. Large tiles are divided at existing
Parquet row groups. Shards apply the same 5–60 minute discovery search and
holdout to every eligible light curve.

```sh
sky hats plan dr24.json --shards 4096
sky hats scan dr24.json 0 --output /shared/dr24-results
sky hats merge dr24.json --output /shared/dr24-results > dr24-candidates.json
```

Run one worker per shard index, starting at zero. Give all workers the same
plan and a shared checkpoint directory. A repeated shard resumes from its
completed checkpoints. The merge reports missing work and exits nonzero
until it has every expected checkpoint. For a small public-archive smoke
test, add `--sample 2 --shards 2` when creating the plan and run both indices.

One hour requires about 1.50 million input rows per second and 2.63 GB/s of
aggregate Parquet throughput if the archive is read in full. The M1 Pro
benchmark consumed 115.9 CPU-seconds over 56,171 input rows, equivalent to
about 3,100 of those measured CPU cores at sustained throughput. A
4,096-shard plan leaves nominal headroom, subject to worker throughput, S3
reads, startup time and the slowest partitions. This is a sizing estimate.
A full timed cluster run and downstream timing fits remain unmeasured.

## Agent-written experiments

Save a trusted local Python file exposing `run(context)` and returning
JSON-compatible output. For example:

```python
def run(context):
    scans = context.scan_eclipses(min_period=5, max_period=40)
    tiles = []
    for scan in scans:
        for candidate in scan["candidates"][:1]:
            tiles.append({
                "candidate": candidate,
                "surveys": context.cone(candidate["ra_deg"], candidate["dec_deg"]),
            })
    return {"tiles": tiles}
```

```sh
sky run ./my_experiment.py --ztf data/ztf --surveys ztf,gaia,roman
```

[ZTF archive](https://irsa.ipac.caltech.edu/Missions/ztf.html) ·
[Gaia epoch photometry](https://gaia.aip.de/metadata/gaiadr3/epoch_photometry/) ·
[Roman archive](https://roman-docs.stsci.edu/data-handbook/accessing-wfi-data/the-roman-archive-in-mast) ·
[Roman simulations](https://spacetelescope.github.io/roman_notebooks/notebooks/ccs_simulations/ccs_simulations.html)
