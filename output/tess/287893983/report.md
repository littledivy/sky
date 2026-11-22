# A 6.368-hour photometric signal in the HD 107126 aperture

The TESS light curve of TIC 287893983 contains a signal at **3.768821 cycles per day** (6.36804 hours). Its flux semi-amplitude in seven independently reduced sector light curves is 578–651 ppm. HD 107126 is the known visual double WDS J12193−5917 AB. TESS pixel fits placed the signal within 4.6 arcsec of the Gaia DR3 source 6059332375114497792 in both sectors with two-minute target-pixel data. This localizes it to the vicinity of the pair. The individual source and physical mechanism are unresolved. [SIMBAD bibliography](https://simbad.idia.ac.za/simbad/sim-id-refs?Ident=HD+107126) · [Mason et al.](https://arxiv.org/abs/0811.0492)

## Independent checks

The discovery search in sector 99 found the dominant 5.927-hour signal and recovered the 6.368-hour signal after subtracting it. Removing the 6.368-hour term increases the residual variance in sector 100. Adding it back reduces that variance by 40% with the discovery frequency held fixed.

For the long-baseline test I used sectors 10, 37, 64 and 99 to select a frequency and phase. Sectors 11, 38 and 100 had independent timing estimates within 5.6 minutes of that prediction. The seven sectors span March 2019 to February 2026. The five earlier curves use the TESS All-Sky Rotation Survey's independent photometry reduction. [Long-baseline timing](coherence.json) · [Sector fits](fit.json) · [Phase plots](phase.png)

## Prior reports

Boyle, Bouma and Mann's TESS All-Sky Rotation Survey reports the 5.927-hour period for this target. Its sector products also list a peak near 7.075 hours. The 6.368-hour frequency is absent from all five archived per-sector peak lists for this target. It is recoverable from those same public light curves. The AAVSO VSX cone search, Gaia DR3 variability tables, public single-sector SPOC event tables, and both Zooniverse Talk searches yielded no prior report of this signal. These are checks for prior identification, with no claim of exhaustive priority. [TARS catalog](https://zenodo.org/records/19917941) · [Per-sector FITS headers](tars-products.json) · [Catalog queries](catalog-checks.json)

## Reproduce

Run from the repository root in the Sky Python environment, with Matplotlib installed.

```sh
python tess_modes.py output/tess/287893983/scan-candidate.json 287893983 \
  --cache data/tess --output output/tess/287893983 \
  --ffi-sectors 10 11 37 38 64 101 --components 3
python tess_coherence.py output/tess/287893983/fit.json \
  output/tess/287893983/tars-products.json --cache data/tess --held-out 11 38 100
```
