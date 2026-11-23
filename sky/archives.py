"""Public survey archive adapters."""

import csv
from contextlib import contextmanager
from functools import lru_cache
from io import BytesIO, StringIO
import math
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse
from astropy.table import Table
from astropy.time import Time
import numpy as np
import requests
from pyarrow.parquet import ParquetFile


GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_DATA = "https://gaia.aip.de/tap/sync"
MAST_API = "https://mast.stsci.edu/api/v0/invoke"
ROMAN_SIMULATIONS = "https://raw.githubusercontent.com/spacetelescope/roman_notebooks/b940b033f0474b28b9fe63973d72ce01d903683d/notebooks/ccs_simulations"
ZTF_BULK = "https://irsa.ipac.caltech.edu/data/ZTF/lc/lc_dr24"
SURVEYS = {
    "ztf": {"release": "DR24", "source": ZTF_BULK, "products": ["local light curves"]},
    "gaia": {"release": "DR3", "source": GAIA_TAP, "products": ["source catalog", "epoch photometry"]},
    "jwst": {"source": MAST_API, "products": ["public observation metadata", "calibrated FITS products"]},
    "roman": {"source": MAST_API, "products": ["WFI dataset metadata", "local extracted photometry"]},
    "roman-sim": {"source": ROMAN_SIMULATIONS, "products": ["simulated WFI exposure metadata"]},
}


@contextmanager
def request(url, **params):
    with requests.get(url, params=params, headers={"User-Agent": "sky-survey/0.1"},
                      timeout=(10, 30), stream=True) as response:
        response.raise_for_status()
        response.raw.decode_content = True
        yield response.raw


def validate_region(ra, dec, radius_arcsec, limit):
    if not (0 <= ra < 360 and -90 <= dec <= 90 and 0 < radius_arcsec <= 3600):
        raise ValueError("Use decimal coordinates and a radius of at most 3600 arcseconds.")
    if not 1 <= limit <= 1000:
        raise ValueError("Limit must be between 1 and 1000.")


def result(rows, source, limit, **metadata):
    return {
        "status": "ok" if rows else "empty", "source": source,
        "rows": rows[:limit], "truncated": len(rows) > limit, **metadata,
    }


def tap_csv(query, limit, required):
    with request(GAIA_TAP, REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                 MAXREC=limit, QUERY=query) as response:
        table = csv.DictReader(StringIO(response.read().decode("utf-8-sig")))
        if not set(required).issubset(table.fieldnames or ()):
            raise ValueError("Gaia TAP did not return the expected CSV columns.")
        return list(table)


def gaia_sources(ra, dec, radius_arcsec=2, limit=20):
    validate_region(ra, dec, radius_arcsec, limit)
    point = "POINT('ICRS', ra, dec)"
    query = f"""
        SELECT TOP {limit + 1} source_id, ra, dec, ref_epoch, pmra, pmdec,
            parallax, phot_g_mean_mag, bp_rp, phot_variable_flag, has_epoch_photometry,
            DISTANCE({point}, POINT('ICRS', {ra:.12f}, {dec:.12f})) AS separation_deg
        FROM gaiadr3.gaia_source
        WHERE CONTAINS({point}, CIRCLE('ICRS', {ra:.12f}, {dec:.12f}, {radius_arcsec / 3600:.12f})) = 1
        ORDER BY separation_deg
    """
    rows = tap_csv(query, limit + 1, ("source_id", "separation_deg", "has_epoch_photometry"))
    for row in rows:
        # JSON clients must be able to retain Gaia's 64-bit identifiers exactly.
        for name in ("ra", "dec", "ref_epoch", "pmra", "pmdec", "parallax",
                     "phot_g_mean_mag", "bp_rp", "separation_deg"):
            row[name] = float(row[name]) if row[name] else None
        row["has_epoch_photometry"] = row["has_epoch_photometry"].lower() == "true"
    return result(rows, GAIA_TAP, limit, data_release="DR3", product="sources")


def gaia_epoch_photometry(source_id):
    source_id = str(source_id)
    if not (source_id.isascii() and source_id.isdecimal() and 0 < int(source_id) < 2**63):
        raise ValueError("Expected a Gaia DR3 source_id.")
    info = tap_csv(
        f"SELECT source_id, has_epoch_photometry FROM gaiadr3.gaia_source WHERE source_id={source_id}",
        1, ("source_id", "has_epoch_photometry"),
    )
    if not info or info[0]["has_epoch_photometry"].lower() != "true":
        return result([], GAIA_DATA, 1, source_id=source_id, time_system="BMJD_TDB",
                      reason="no_epoch_photometry" if info else "source_absent")
    query = f"SELECT * FROM gaiadr3.epoch_photometry WHERE source_id={source_id}"
    with request(GAIA_DATA, REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=query) as response:
        votable = parse(BytesIO(response.read()))
    for resource in votable.resources:
        for item in resource.infos:
            if item.name == "QUERY_STATUS" and item.value != "OK":
                raise ValueError(item.content or "Gaia AIP TAP epoch query was incomplete.")
    table = votable.get_first_table().to_table(use_names_over_ids=True)
    rows = []
    for entry in table:
        for band, prefix, time_column in (("G", "g_transit", "g_transit_time"),
                                           ("BP", "bp", "bp_obs_time"), ("RP", "rp", "rp_obs_time")):
            arrays = [np.ma.asarray(entry[key]).filled(np.nan)
                      for key in (time_column, f"{prefix}_mag", f"{prefix}_flux", f"{prefix}_flux_error")]
            time, mag, flux, flux_error = arrays
            rejected = np.ma.asarray(entry[f"variability_flag_{band.lower()}_reject"]).filled(True)
            valid = np.isfinite(time + mag + flux + flux_error) & (flux > 0) & (flux_error > 0) & ~rejected
            if not valid.any():
                continue
            bmjd = Time(2455197.5, time[valid], format="jd", scale="tcb").tdb.mjd
            errors = 2.5 / math.log(10) * flux_error[valid] / flux[valid]
            rows.extend({"band": band, "bmjd_tdb": float(t), "magnitude": float(m),
                         "magnitude_error": float(e)}
                        for t, m, e in zip(bmjd, mag[valid], errors, strict=True))
    return result(rows, GAIA_DATA, max(len(rows), 1), source_id=source_id,
                  time_system="BMJD_TDB", product="epoch_photometry")


def roman_observations(ra, dec, radius_arcsec=60, limit=20):
    from astroquery.mast.observations import ObservationsClass

    validate_region(ra, dec, radius_arcsec, limit)
    table = ObservationsClass().query_criteria(
        coordinates=SkyCoord(ra, dec, unit="deg"), radius=radius_arcsec * u.arcsec,
        obs_collection="ROMAN", pagesize=limit + 1, page=1,
    )
    rows = []
    for row in table:
        values = {}
        for name in table.colnames:
            value = row[name]
            if np.ma.is_masked(value):
                values[name] = None
            elif isinstance(value, (float, np.floating)) and not np.isfinite(value):
                values[name] = None
            elif isinstance(value, np.generic):
                values[name] = value.item()
            elif isinstance(value, (bool, float, int, str)) or value is None:
                values[name] = value
            else:
                values[name] = str(value)
        rows.append(values)
    return result(rows, MAST_API, limit, product="wfi_observation_metadata",
                  archive_collection="ROMAN", search_geometry="MAST footprint intersection")


def jwst_observations(ra, dec, radius_arcsec=60, limit=20):
    """Find public JWST science observations intersecting a sky position."""
    from astroquery.mast import Observations

    validate_region(ra, dec, radius_arcsec, limit)
    table = Observations.query_criteria(
        coordinates=SkyCoord(ra, dec, unit="deg"), radius=radius_arcsec * u.arcsec,
        obs_collection="JWST", intentType="science", dataRights="PUBLIC",
        pagesize=limit + 1, page=1,
    )
    rows = []
    for row in table:
        values = {}
        for name in table.colnames:
            value = row[name]
            if np.ma.is_masked(value):
                values[name] = None
            elif isinstance(value, (float, np.floating)) and not np.isfinite(value):
                values[name] = None
            elif isinstance(value, np.generic):
                values[name] = value.item()
            elif isinstance(value, (bool, float, int, str)) or value is None:
                values[name] = value
            else:
                values[name] = str(value)
        rows.append(values)
    rows.sort(key=lambda value: value.get("distance") if value.get("distance") is not None else float("inf"))
    return result(rows, MAST_API, limit, product="jwst_observation_metadata",
                  archive_collection="JWST", search_geometry="MAST footprint intersection")


def jwst_products(obsid, *, download=False, directory="data/jwst"):
    """List public calibrated JWST FITS science products for one MAST obsid."""
    from astroquery.mast import Observations

    try:
        obsid = int(obsid)
    except (TypeError, ValueError) as error:
        raise ValueError("Expected a numeric MAST obsid from sky cone --surveys jwst.") from error
    if obsid <= 0:
        raise ValueError("Expected a positive MAST obsid.")
    products = Observations.get_product_list(str(obsid))
    products = Observations.filter_products(
        products, extension="fits", calib_level=[2, 3], productType="SCIENCE",
    )
    rows = []
    for row in products:
        rows.append({
            "obsid": int(row["obsID"]),
            "dataURI": str(row["dataURI"]),
            "productFilename": str(row["productFilename"]),
            "productSubGroupDescription": str(row["productSubGroupDescription"]),
            "calib_level": int(row["calib_level"]),
            "size": int(row["size"]),
            "filters": str(row["filters"]),
            "instrument_name": str(row["instrument_name"]),
        })
    manifest = None
    if download and len(products):
        manifest_table = Observations.download_products(products, download_dir=str(directory))
        manifest = [
            {name: (None if np.ma.is_masked(row[name]) else str(row[name]))
             for name in manifest_table.colnames}
            for row in manifest_table
        ]
    return {"status": "ok" if rows else "empty", "source": MAST_API,
            "rows": rows, "truncated": False, "product": "calibrated_jwst_fits",
            "obsid": obsid, "downloaded": bool(download), "download_directory": str(directory),
            "manifest": manifest}


@lru_cache(maxsize=2)
def roman_simulation_manifest(program):
    url = f"{ROMAN_SIMULATIONS}/{program}_manifest.ecsv"
    with request(url) as response:
        table = Table.read(BytesIO(response.read()), format="ascii.ecsv")
    centers = SkyCoord(table["detector_ra"], table["detector_dec"], unit="deg")
    return url, table, centers


def roman_simulations(ra, dec, radius_arcsec=360, limit=20, program="gbtds"):
    validate_region(ra, dec, radius_arcsec, limit)
    if program not in ("gbtds", "hlwas"):
        raise ValueError("Roman simulation program must be gbtds or hlwas.")
    url, table, centers = roman_simulation_manifest(program)
    distances = SkyCoord(ra, dec, unit="deg").separation(centers).arcsec
    indices = np.flatnonzero(distances <= radius_arcsec)
    indices = indices[np.argsort(distances[indices], kind="stable")[:limit + 1]]
    rows = [
        {**{key: str(table[key][i]) for key in ("detector", "optical_element", "start_time", "l2_s3_uri")},
         "detector_ra": float(table["detector_ra"][i]), "detector_dec": float(table["detector_dec"][i]),
         "center_separation_arcsec": float(distances[i])} for i in indices
    ]
    return result(rows, url, limit, product="detector_centers", simulated=True, program=program)


def fetch_ztf(field, ccd, quadrant, band, directory):
    if not (1 <= field <= 1999 and 1 <= ccd <= 16 and 1 <= quadrant <= 4 and band in ("g", "r", "i")):
        raise ValueError("Use a ZTF field, CCD 1–16, quadrant 1–4, and band g, r, or i.")
    filename = f"ztf_{field:06d}_z{band}_c{ccd:02d}_q{quadrant}_dr24.parquet"
    url = f"{ZTF_BULK}/{field // 1000}/field{field:06d}/{filename}"
    path = Path(directory) / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with NamedTemporaryFile(dir=path.parent, prefix=".ztf-", delete=False) as temporary:
            temp_path = Path(temporary.name)
        try:
            with request(url) as response, temp_path.open("wb") as output:
                shutil.copyfileobj(response, output)
            with ParquetFile(temp_path):
                pass
            temp_path.replace(path)
        finally:
            temp_path.unlink(missing_ok=True)
    with ParquetFile(path) as parquet:
        return {"path": str(path.resolve()), "source": url, "rows": parquet.metadata.num_rows}
