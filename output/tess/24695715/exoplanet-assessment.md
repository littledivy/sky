# Is TIC 24695715 a planet host?

**Current read:** possible in principle, but the data now favor an eclipsing-binary interpretation over a single planet transiting every 5.7766 hours. Nothing here confirms either interpretation.

## The strongest clue: alternating eclipse depths

The discovery period is the interval between successive dips, 0.240691392 days. If that is half the binary's orbital period, the true period would be 0.481382784 days (11.5532 hours), with alternating primary and secondary eclipses.

Paired event-depth measurements support checking that alias:

| Data | Adjacent event pairs | Mean even-minus-odd depth | Simple standard error |
| --- | ---: | ---: | ---: |
| Sectors 95/96 | 38 | 0.090 percentage points | 0.023 percentage points |
| Older full-frame cutouts | 40 | 0.136 percentage points | 0.025 percentage points |

Even-numbered events are deeper in both independent time baselines. These are paired differences from the trapezoid event fits in `events.csv`; the quoted errors describe scatter across event pairs and are not a full transit-model posterior. The repeating parity is consistent with two similar, but not identical, eclipses in a binary orbit. A planet on a 5.7766-hour orbit would normally give the same primary transit on every cycle; it would not naturally produce this alternating transit depth pattern.

The absence of a dip at phase 0.5 when folded on 5.7766 hours does not rule out the double-period binary: the two eclipse types land on top of each other at phase 0 in that fold. A fold at 11.5532 hours should show the two types separately.

## Host-star and contamination checks

The TESS Input Catalog v8.2 lists TIC 24695715 as a 5313 K star with radius 0.957 solar radii. Its estimated contamination ratio is 0.000354; the TESS-SPOC CROWDSAP values are 0.982 and 0.986 in sectors 95 and 96. The nearest other TIC source in a two-arcminute cone is 43 arcseconds away and 7.1 TESS magnitudes fainter. These checks make a bright, widely separated contaminant less likely, but do not exclude a close unresolved companion.

The matched Gaia DR3 source is `4697634881027710464`. Gaia reports parallax 1.908 ± 0.102 mas, G=12.723, BP−RP=0.870, and Teff=5612 K. Its radius estimates differ: FLAME gives 1.451 (reported interval 1.376–1.531) solar radii, while GSP-Phot gives 1.826 (1.753–1.940). The solution also has RUWE=8.918, astrometric excess noise 1.215 mas with significance 1573, and no Gaia variability classification. No row was returned for this source in Gaia DR3's two-body-orbit or astrometric-acceleration tables. The poor astrometric fit and discrepant TIC/Gaia radii mean the host-star radius is not secure; unresolved multiplicity is one possible explanation, not a confirmed one.

If the 1.31% depth were an undiluted planetary transit, the simple radius-ratio estimate gives about 1.06 Jupiter radii using the TIC radius, 1.61 using Gaia FLAME, or 2.03 using Gaia GSP-Phot. Those estimates assume a single host star and a non-grazing transit; the alternating depths and radius disagreement invalidate using them as a planet measurement. For the Gaia FLAME mass (0.960 solar masses), Kepler's law gives a semimajor axis of about 1.61 solar radii at 5.7766 hours and 2.55 solar radii at 11.5532 hours. The Gaia single-star radius is therefore especially difficult to reconcile with a 5.7766-hour planet interpretation, but its high RUWE cautions against treating that radius as definitive.

## Catalog status and what would settle it

There is no NASA Exoplanet Archive TOI row and no TIC 24695715 match in the public SPOC TCE statistics for sectors 95 or 96. No match was found in the five checked published TESS eclipsing-binary tables or in the VSX search made before submission. Catalog absence does not establish that the signal is a planet.

The next useful work is to fit the two eclipse types separately at 11.5532 hours with a physical binary/transit model, compare odd/even timing and shapes across the older cutouts, and determine whether Gaia's astrometric anomaly is consistent with an unresolved companion. High-resolution imaging or multi-band ground-based photometry would test for a nearby eclipsing source; radial-velocity monitoring would distinguish a planetary companion from a stellar companion. Until then, call this an eclipsing-object candidate and do not claim an exoplanet.

## Data sources

- [MAST TIC cone-search service](https://mast.stsci.edu/vo-conesearch/api/v0.1/tic?RA=25.0617122455082&DEC=-67.8150662187333&SR=0.0004)
- [TIC v8.2 field definitions](https://outerspace.stsci.edu/spaces/TESS/pages/111023212/TIC%2Bv8.2%2Band%2BCTL%2Bv8.xx%2BData%2BRelease%2BNotes)
- [Gaia DR3 source query](https://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT%20source_id%2Cra%2Cdec%2Cparallax%2Cparallax_error%2Cphot_g_mean_mag%2Cbp_rp%2Cteff_gspphot%2Clogg_gspphot%2Cruwe%2Castrometric_excess_noise%20FROM%20gaiadr3.gaia_source%20WHERE%20source_id%3D4697634881027710464)
- [NASA Exoplanet Archive TOI table](https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select%20tid%2Ctoi%20from%20toi%20where%20tid%3D24695715&format=csv)
- [MAST public TESS SPOC TCE tables](https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html)
