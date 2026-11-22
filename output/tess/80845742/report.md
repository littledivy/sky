# Recovered eclipses of HD 96135

HD 96135 (TIC 80845742, HIP 54129, Gaia DR3 5347442439131706624) exhibits a recurring primary eclipse in public TESS photometry. The eclipse is about **0.595% deep** and **11.28 hours long**. TESS recorded six eclipses in its two-minute photometry and 15 more in independently extracted full-frame photometry from five earlier sectors, spanning 2019–2026. Gaia DR3 previously published a single-lined spectroscopic orbit. IJspeert et al. (2024) had already identified the eclipses. The Sky scan independently recovered them and refined the ephemeris. [IJspeert et al. (2024)](https://doi.org/10.1051/0004-6361/202450507) [SIMBAD](https://simbad.idia.ac.za/simbad/sim-id?Ident=HD%2096135) · [Gaia DR3 spectroscopic orbits](https://cdsarc.cds.unistra.fr/viz-bin/cat/I/357#/browse)

## Measured parameters

| Field | Value |
|---|---|
| Primary name | HD 96135 |
| Other identifiers | TIC 80845742, HIP 54129, Gaia DR3 5347442439131706624, TYC 8615-1089-1, 2MASS J11042879-5306306 |
| Coordinates, ICRS at J2000 | 11h 04m 28.78s, −53° 06′ 30.65″ |
| Variability type | EA |
| Period | 7.73629 d |
| Primary minimum | BJD TDB 2461046.55638, equivalent to geocentric HJD UTC 2461046.55556 |
| Eclipse duration | 0.4698 d, about 6.1% of a cycle |
| Catalog magnitude and measured depth | TESS magnitude 7.726. The primary eclipse is about 0.0065 mag below the local continuum in the TESS passband. |

The estimated epoch is the center of the fitted primary eclipse. The HJD value uses the geocenter and the Astropy built-in Solar System ephemeris. The magnitude is from the TESS Input Catalog and the depth is differential. [AAVSO submission guidelines](https://vsx.aavso.org/index.php?view=about.notice) · [TESS Input Catalog](https://tess.mit.edu/science/tess-input-catalogue/)

## Evidence

1. A blind native-cadence search of 900 stars with two-minute light curves in Sectors 99 and 100 found three repeating flux changes in Sector 99. The sparse filter picked up the sharp ingress. A trapezoid fitted to Sector 99 gave a period of 7.73648 d and predicted all three usable eclipse windows in Sector 100. With that ephemeris and the eclipse shape fixed, Sector 100 depths were 0.597%, 0.596%, and 0.559%. The [discovery scan](discovery-scan.json) and [holdout measurements](discovery_ephemeris_checks.csv) preserve both steps.
2. The six-event fit across Sectors 99 and 100 gave 7.73609 ± 0.00016 d and a depth of 0.5946 ± 0.0062%. The uncertainties are leave-one-event-out estimates and do not include systematic model uncertainty. TESS difference images placed the signal 3.2″ and 4.7″ from HD 96135. The brightest other Gaia source within 90″ is 23.2″ away and 6.64 mag fainter in Gaia G. [Fit and pixel provenance](fit.json) · [Individual eclipses](eclipses.png)
3. Direct aperture photometry on archival TESS full-frame cutouts recovered the eclipse in all 15 usable windows across Sectors 10, 37, 63, 64, and 90. All 15 depths were positive with the recent two-sector ephemeris and shape held fixed. The median fixed-template depth was 0.564%. Allowing the old eclipse centers to shift and fitting their sector medians across the longer baseline gave a photometric period of 7.7362925 d and about 0.41 minutes RMS across the fitted sector medians. The latter is a fit diagnostic. [Archival measurements](ffi-events.csv) · [MAST TESSCut](https://mast.stsci.edu/tesscut/)
4. Gaia DR3 gives an SB1 period of 7.736533 ± 0.000166 d and a radial-velocity semi-amplitude of 16.802 ± 0.094 km/s for the same source. The period agreement and eclipse localization support identifying the photometric eclipse with the known spectroscopic orbit. [Gaia DR3 SB1 data](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I%2F357%2Ftbosb1) · [Gaia SB1 processing paper](https://doi.org/10.1051/0004-6361/202450600)

[Phase-folded light curve and time series](phase.png) · [Two-minute photometry](photometry.csv)

## Prior reports

Checked the target in the AAVSO VSX, live Villanova TESS Eclipsing Binary Catalog, TESS Ten Thousand, ASAS-SN variable catalog, Gaia DR3 photometric-variability tables, NASA TOI and ExoFOP community TOI lists, and the published TESS periodic-signal and short-period variable tables listed in [catalog checks](catalog-checks.json). The VSX query returned two other variables with different coordinates and variability. SIMBAD identifies HD 96135 as a spectroscopic binary. Its bibliography includes the Gaia spectroscopic solution, proper-motion and acceleration studies, radial-velocity measurements, and spectral and stellar-parameter catalogs. An additional machine-readable catalog check identified HD 96135 in the IJspeert et al. (2024) eclipse table, with a 7.73435 d orbital period from TESS-SPOC and 7.73432 d from QLP. This record establishes prior publication of the eclipses. [VizieR catalog](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J%2FA%2BA%2F691%2FA242%2Fobaf-eb1)

## Reproduce

From the repository root, with the Sky virtual environment active and Matplotlib installed, run:

```sh
sky tess 99 100 --cache data/tess --prioritize-new-since {1..98} \
  --limit 900 --workers 10 --search sparse > tess-99-100.json

MPLBACKEND=Agg python experiments/tess/tess_vet.py tess-99-100.json 80845742 \
  --cache data/tess --output output/tess/80845742 \
  --shift-days .125 --duration-days .42 --window-days .62 \
  --ffi-sectors 10 37 63 64 90
```

The follow-up fit uses a trapezoid and a local quadratic continuum per eclipse. The target-specific shift and duration initialize the fit around the eclipse seen in the discovery light curve. The scan covered the complete 4,453-star overlap of Sectors 107 and 108 in about 244 seconds on an Apple M1 Pro with a partly populated cache. The Sectors 99 and 100 candidate batch took 7.7 seconds with cached photometry.
