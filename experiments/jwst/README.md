# JWST public-data discovery workflow

Sky can search public JWST science observations through MAST, list the calibrated
FITS science products for an observation, and download those products on request.
Start with observation metadata and a small product shortlist; NIRCam mosaics
and NIRSpec spectra can be large. The current adapter does not yet identify
sources in images or extract spectra/photometry.

## Finding overlooked compact red sources

Little red dots (LRDs) are an active topic, so reproducing an existing LRD
catalog is not a discovery. A useful search should target a clearly defined
gap, such as candidate transition objects, unusual variability, or sources
whose imaging colors and public spectra disagree. The field should have public
multi-band NIRCam coverage and, ideally, NIRSpec or independent HST/Chandra data.

For each lead:

1. Record the exact MAST program, observation IDs, filters, and calibrated
   product URIs; only use public data.
2. Detect and measure the source in the science mosaics with their weight and
   segmentation products. Check artifacts, blending, and astrometric offsets.
3. Fit photometric redshift and compare with any public spectroscopic redshift.
   Inspect spectra when available; a red color alone does not establish an
   LRD, AGN, black hole, or high redshift.
4. Cross-match catalogues and papers before describing a source as new. Check
   alternate explanations, including dusty star formation and foreground
   contamination, and seek an independent data product to validate it.
5. Save a reproducible record: query, source coordinates, data URIs, quality
   checks, measurements, plots, and literature/catalogue search date.

The first integration milestone is MAST discovery and data provenance. Full
image-level source finding and spectral fitting are follow-up work, to be built
only after selecting a tractable field and reproducible data products.
