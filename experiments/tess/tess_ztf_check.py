"""Check a fixed two-eclipse TESS ephemeris against resolved ZTF photometry."""

import argparse
import csv
import json
from pathlib import Path

from astropy.constants import c
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
import astropy.units as u
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument("ztf", type=Path, help="Public IRSA light-curve CSV for one source")
    parser.add_argument("--half-window-days", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fit = json.loads(args.fit.read_text())
    with args.ztf.open() as stream:
        observations = [row for row in csv.DictReader(stream) if int(row["catflags"]) == 0]
    if not observations:
        raise ValueError("The ZTF light curve has no unflagged observations.")
    time = Time([float(row["hjd"]) for row in observations], format="jd", scale="utc")
    source = SkyCoord(*fit["proposed_source_icrs_deg"], unit="deg").cartesian
    delay = (get_body_barycentric("sun", time.tdb).dot(source) / c).to_value(u.day)
    btjd = time.tdb.jd + delay - 2457000
    night = np.floor(btjd)
    band = np.array([row["filtercode"] for row in observations])
    magnitude = np.array([float(row["mag"]) for row in observations])
    first, second = fit["joint_fit"], fit["secondary_fit"]
    period = first["period_days"]
    rows = []
    for name in np.unique(band):
        same_band = band == name
        groups = [same_band & (night == day) for day in np.unique(night[same_band])]
        reference = float(np.median([np.median(magnitude[group]) for group in groups]))
        for group in groups:
            t = float(np.median(btjd[group]))
            mag = float(np.median(magnitude[group]))
            offsets = [(t - eclipse["mid_eclipse_btjd"] + period / 2) % period - period / 2
                       for eclipse in (first, second)]
            order = int(np.argmin(np.abs(offsets)))
            rows.append({"btjd": t, "band": name, "magnitude": mag, "median_band_magnitude": reference,
                         "depth_mag": mag - reference, "closest_eclipse": ("primary", "secondary")[order],
                         "offset_days": offsets[order],
                         "inside": abs(offsets[order]) < args.half_window_days})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    selected = [row["depth_mag"] for row in rows if row["inside"]]
    print(json.dumps({"good_observations": len(observations), "night_band_groups": len(rows),
                      "in_window_groups": len(selected), "groups_fainter_than_band_median": sum(x > 0 for x in selected),
                      "mean_dimming_mag": float(np.mean(selected)) if selected else None}, indent=2))


if __name__ == "__main__":
    main()
