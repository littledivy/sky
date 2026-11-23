# Gaia DR4 discovery readiness

Gaia DR4 is expected on 2 December 2026. Sky currently queries Gaia DR3. The
official archive's DR4 table names and column semantics must be checked against
the live TAP schema before adding production queries; DR3 names and time
conversions must not be assumed to carry over.

## First archive pass

After release, save the live schema and column definitions:

```sh
python experiments/gaia/dr4_schema_probe.py > output/gaia-dr4-schema.json
```

The probe discovers DR4 schemas and tables through ESA TAP's `TAP_SCHEMA` and
records their column metadata. An empty DR4 section means the archive has not
published discoverable DR4 tables yet; retry later rather than treating that
as a scientific null result.

## Search priorities

1. **Astrometric binaries and dormant black holes.** Inspect DR4's non-single-
   star and black-hole binary solutions, then independently vet high-value
   systems using Gaia astrometry, spectra, and available TESS/ZTF photometry.
2. **Unusual variables and eclipsing systems.** Use released epoch photometry
   and variability classifications to surface unclassified or poorly modeled
   sources; verify periods and eclipse shapes with TESS, ZTF, and other public
   time-domain surveys.
3. **Planet candidates.** Cross-match DR4's exoplanet products with known
   planet and candidate lists, then seek independent transit or radial-velocity
   evidence. A catalog candidate alone is not a new discovery.
4. **Rare outliers.** Rank objects with internally inconsistent astrometry,
   photometry, and source classifications, while checking crowding, blends,
   quality flags, and known literature before interpreting them physically.

## Discovery standard

Every shortlist item should carry reproducible query text, source IDs, data
release/table provenance, quality cuts, known-object checks, plots, and an
independent-data check. Prioritize a small number of high-information leads
over a large unvetted anomaly list. Do not describe an object as new until
catalog and literature checks are complete and the evidence supports the claim.

## Scope note: JWST little red dots

Gaia DR4 is not a practical direct survey for the high-redshift JWST “little red
dots”: they are extremely distant, compact galaxies, while Gaia's strength is
precise astrometry and time-domain measurements of sources it can resolve.
Sky can still search for nearby black-hole binaries and unusual stellar
systems, which are a distinct and accessible discovery space.
