# Eclipses in Gaia DR3 5073603089619768704

Gaia DR3 5073603089619768704 (TIC 88481083, 2MASS J03090374−2611523) shows two eclipses per **3.46071-day orbit**. Sky found the signal in TESS photometry extracted around the brighter star CD−26 1162 (TIC 88481081). Difference images from four sectors place both eclipses on the fainter Gaia source, about 21 arcseconds away from the center of the original aperture. Gaia and VSX previously listed the fainter star as an RS-type variable. Planet Hunters TESS volunteers reported the signal, its location, and a 3.4608-day period on November 19, 2026. This analysis independently recovers their result. [Planet Hunters TESS discussion](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/2110/4185529) [Existing VSX record](https://vsx.aavso.org/index.php?view=detail.top&oid=2868729) · [Gaia DR3 variable classification](https://doi.org/10.1051/0004-6361/202245591)

## Proposed revision to VSX entry 2868729

| Field | Proposed value |
|---|---|
| Designation | Gaia DR3 5073603089619768704 |
| Aliases | TIC 88481083, 2MASS J03090374−2611523 |
| Position, ICRS at J2000 | 03h 09m 03.75s, −26° 11′ 52.51″ |
| Variability type | EA, with out-of-eclipse modulation |
| Orbital period | 3.46071 ± 0.00003 d |
| Primary minimum | BJD TDB 2461210.49074, equivalent to geocentric HJD UTC 2461210.48996 |
| Eclipse duration | 2.12 h primary, 2.10 h secondary, about 2.5% of an orbit each |
| Secondary phase | Approximately 0.500 |
| Photometric range | Existing VSX range 13.38–13.53 in Gaia G. The blended TESS aperture shows depths of 0.542% and 0.482% for the two eclipses. |

Suggested remarks. Two narrow eclipses, with their centers separated by half an orbital cycle, recur across TESS Sectors 105–108. Difference images localize both eclipse series to Gaia DR3 5073603089619768704. The TESS light curves were extracted for TIC 88481081 and include both stars. The reported TESS depths are relative to the blended aperture flux. Resolved ZTF photometry provides supporting detections near the predicted eclipse windows.

[Phase plot and time series](phase.png) · [TESS eclipse localization](localizations.png) · [Individual primary eclipses](eclipses.png)

## Measurement record

The blind sparse search found the signal in Sector 105. A fit to Sector 105 alone gave a period of 3.46078 d for the primary. Both sets of eclipse times and shapes fitted to Sector 105 were held fixed for checks in Sectors 106–108. All 16 usable primary windows and all 14 usable secondary windows in those three sectors showed dimming. The four-sector primary fit gives a period of 3.4607089 d with a leave-one-event-out standard error of 0.0000264 d. The secondary-only solution gives 3.4607828 d. The period uncertainty reflects the scatter in the TESS measurement and excludes model systematics. [Blind scan](discovery-scan.json) · [Primary holdout](discovery_ephemeris_checks.csv) · [Secondary holdout](secondary_ephemeris_checks.csv)

The calibrated target pixel files contain 24 primary and 21 secondary eclipses with sufficient before-event, in-event, and after-event sampling for paired difference images. Across the four camera/detector combinations, all eight difference centroids lie within 2.9″ of the fainter Gaia source. The primary and secondary each recur at the same position. The two Gaia catalog positions have similar parallaxes and proper motions, suggesting a wide comoving pair. Binding of the wider pair remains unconfirmed. [Pixel measurements and source FITS URLs](fit.json) · [Gaia neighborhood](catalog-checks.json)

The public IRSA ZTF light curve is position-matched to the fainter star. After selecting unflagged exposures and taking one median measurement per night and band, six of 68 groups fall within 0.045 d of the TESS-predicted primary or secondary centers. All six are fainter than their band medians, with a mean difference of 0.039 mag. The TESS period and epochs were held fixed in this check. ZTF's sparse sampling gives supporting evidence for the assignment and ephemeris. [Source ZTF table](ztf-irsa.csv) · [Resolved ZTF measurements](ztf-eclipse-check.csv) · [IRSA light-curve API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html)

## Prior-publication checks

The VSX entry identifies an RS-type variable and has empty period, epoch, and eclipse-duration fields. The Gaia DR3 general classifier gives RS and the Gaia eclipsing-binary and orbital-solution tables lack a matching entry. Exact TIC matching and source-position checks covered the live TESS eclipsing-binary catalog, the IJspeert OBAF eclipsing table, TESS Ten Thousand, TESS prime-mission variability, the TESS All-Sky Rotation Survey, ASAS, ASAS-SN, ATLAS, published ZTF variables, and the NASA TOI and ExoFOP CTOI tables. In the November 19, 2026 Planet Hunters TESS discussion, DAAN11 reported dips in TIC 88481081 spaced by about 1.73 d. Later that day, moosehunter21 reported alternating eclipse depths and used pixel data to identify TIC 88481083 as the source, with an orbital period of 3.4608 d. That public report predates the Sky analysis. The pending VSX revision was visible during the catalog check, and its scope remains unknown. [Planet Hunters TESS discussion](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/2110/4185529) · [Earlier eclipse analysis and plots](https://github.com/moosesteve21/PlanetCandidates/tree/main/checks/TIC88481081) [Catalog and identifier checks](catalog-checks.json) · [Published TIC crossmatch](tic-catalog-crossmatch.json)

## Reproduce

From the repository root, with the Sky virtual environment active and Matplotlib installed:

```sh
sky tess 105 106 --cache data/tess --prioritize-new-since {1..104} \
  --limit 900 --workers 10 --search sparse > tess-105-106.json

python experiments/tess/tess_catalogs.py tess-105-106.json --cache data/tess \
  --source-tic 88481083

MPLBACKEND=Agg python experiments/tess/tess_vet.py tess-105-106.json 88481081 \
  --cache data/tess --output output/tess/88481083 \
  --extra-sectors 107 108 --duration-days .10 --window-days .25 \
  --source-coordinate 47.26535016 -26.19797613 --secondary-phase .5

python experiments/tess/tess_ztf_check.py output/tess/88481083/fit.json \
  output/tess/88481083/ztf-irsa.csv --half-window-days .045 \
  --output output/tess/88481083/ztf-eclipse-check.csv
```

The `88481081` argument selects the source of the downloaded photometry. The `--source-coordinate` argument specifies the position of the variable identified in the difference images. The optional phase seeds the fit to the secondary seen in Sector 105. [TESS two-minute photometry](photometry.csv) · [TESS archive documentation](https://archive.stsci.edu/missions-and-data/tess/data-products)
