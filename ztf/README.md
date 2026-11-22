# Blind ZTF short-period eclipse search

The scanner reads ZTF DR24 bulk light-curve Parquet files and outputs one JSON
record per partition. It uses only observation times, magnitudes, uncertainties,
and quality flags to select and score periods. The sky position is copied to the
candidate record for follow-up.

## Run a pilot

With Python 3.14, run these commands from the repository root:

```sh
python3 -m venv /tmp/ztf-experiment-env
/tmp/ztf-experiment-env/bin/pip install -r ztf/requirements.txt
mkdir -p /tmp/ztf-dr24-pilot
curl -fL 'https://irsa.ipac.caltech.edu/data/ZTF/lc/lc_dr24/0/field000697/ztf_000697_zg_c01_q1_dr24.parquet' \
  -o /tmp/ztf-dr24-pilot/ztf_000697_zg_c01_q1_dr24.parquet
/tmp/ztf-experiment-env/bin/python ztf/short_period_eclipses.py \
  /tmp/ztf-dr24-pilot --limit-per-partition 256 > /tmp/ztf-dr24-pilot/results.jsonl
```

The pilot uses the first green-band file in the previously selected field 697
and takes the first 256 eligible curves in archive order. For each curve, the
script selects its densest 60-day window by cadence, scans 5–60 minutes with
Lomb–Scargle, then fits box eclipses around the strongest frequency peaks and
their harmonics. Every third observation is kept out of the fit. Candidate scores
measure improvement over a flat model on those held-out points. Scores are
ranking statistics. They need injection/recovery tests and false-alarm
calibration before any candidates can be interpreted as discoveries.

Pass a directory of bulk-release partitions and omit the pilot limit to scan
every light curve in that directory. Use `--workers` to process several files
concurrently. Each record reports the file row count, the rows examined, the
number searched, and the outcomes. The default output retains the top 20 entries
per partition. A partition error is reported in the JSON and the command exits
with a nonzero status.

## Local result

The full green-band file contained 21,102 light curves. Of those, 2,748 met
the cadence and data-quality cuts within the search window, and 223 passed
the initial holdout filter. The report retains the 20 highest-scoring entries
in `output/ztf/field000697-green-q1.json`. The full file was divided into
eight worker inputs and the concatenated inputs were checked against the
original Parquet table before reporting coverage.

On one real ZTF cadence, a synthetic 17.39-minute signal was recovered and
a 27.43-minute signal was missed. Two of 32 Gaussian noise-only simulations
on that cadence passed the same filter. These checks illustrate why the
holdout score alone cannot establish a discovery.

## Reproduction scope

Burdge et al. used conditional entropy on a June 2018 high-cadence ZTF sample.
Their orbital-decay measurement used separate PTF/iPTF, CHIMERA, and KPED
observations. This script implements a blind discovery-stage search on DR24 with
a different period-finding method. Published timings, ephemerides and object
coordinates are unnecessary for the search. A reproduction of the decay
measurement would also need to extract eclipse epochs from the observational
data and fit the cycle counts and the changing period.

[ZTF DR24 release notes](https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/dr24/ztf_release_notes_dr24.pdf)

[Burdge et al. methods and data availability](https://arxiv.org/pdf/1907.11291)

## J1539 recovery in field 760

The positive-control run used the full green and red archive files for
field 760, CCD 14, quadrant 3. The sky position identified the archive files
and was used to identify J1539 after ranking. The scanner used the default
5–60 minute range and took no source ID or reference period.

```sh
mkdir -p /tmp/ztf-field760
for band in zg zr; do
  file="ztf_000760_${band}_c14_q3_dr24.parquet"
  curl -fL "https://irsa.ipac.caltech.edu/data/ZTF/lc/lc_dr24/0/field000760/$file" \
    -o "/tmp/ztf-field760/$file"
done
/tmp/ztf-experiment-env/bin/python ztf/short_period_eclipses.py \
  /tmp/ztf-field760 --workers 2 > /tmp/ztf-field760/ranked.jsonl
```

J1539 ranked 16th of 2,865 searched green-band light curves at
6.913144745659857 minutes. Its red-band light curve failed the blind holdout
filter. The red-band measurements were included after source identification
to fit the timing drift.

```sh
/tmp/ztf-experiment-env/bin/python ztf/fit_orbital_drift.py \
  --period-minutes 6.913144745659857 \
  /tmp/ztf-field760/ztf_000760_zg_c14_q3_dr24.parquet:760114300005872 \
  /tmp/ztf-field760/ztf_000760_zr_c14_q3_dr24.parquet:760214300009158
```

The timing code compares constant and quadratic phase models with two
harmonics and per-band amplitudes. It searched around the recovered period,
including daily frequency aliases. The joint fit used 1,537 observations
spanning 2,756 days and measured a period derivative of
`(-2.390 ± 0.028) × 10^-11 s/s`. The uncertainty is local and statistical.
The report in `output/ztf/j1539-field760.json` includes the fits from each
band, the sample coverage, the highest-ranked candidates, and the published
value for comparison.

## Fast-moving stars

Run the targeted search on original DR24 Parquet files with both filters
available for each detector quadrant.

```sh
sky run ./ztf/high_motion_eclipses.py --surveys ztf,gaia --ztf data/ztf
```

The experiment selects Gaia sources with total proper motion of at least
150 mas per year. It matches their projected 2018–2025 positions to ZTF
catalog centroids, searches each matched light curve, and checks predicted
eclipses in the other filter. These data are existing catalog light curves.
The position track is used for matching and the search uses the photometry
already extracted by ZTF.
