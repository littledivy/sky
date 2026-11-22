"""Fit and localize one candidate from a Sky TESS scan."""

import argparse
import csv
import hashlib
from io import BytesIO
import json
from pathlib import Path
from urllib.parse import urlencode
from zipfile import ZipFile

from astropy.io import fits
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from sky.tess import download, photometry, product_url, sector_products
from sky.tess_pixels import localize


def events_for(data, sector, period, epoch, window):
    time, flux, _ = data
    start = int(np.ceil((time[0] - epoch) / period))
    end = int(np.floor((time[-1] - epoch) / period))
    events = []
    for cycle in range(start, end + 1):
        offset = time - (epoch + cycle * period)
        inside = abs(offset) < window
        x, y = offset[inside], flux[inside]
        if len(x) < 120 or min(sum(x < -window / 2), sum(x > window / 2)) < 20:
            continue
        adjacent = np.diff(y)[np.diff(x) < 3 / 1440]
        scatter = 1.4826 * np.median(abs(adjacent - np.median(adjacent))) / np.sqrt(2)
        if not np.isfinite(scatter) or scatter <= 0:
            continue
        baseline = np.column_stack((np.ones_like(x), x / window, (x / window)**2))
        basis, _ = np.linalg.qr(baseline)
        events.append((sector, cycle, x, y, scatter, basis))
    return events


def shape(offset, duration, ingress):
    return np.clip((duration / 2 - abs(offset)) / (duration * ingress), 0, 1)


def fit_eclipse(events, period, epoch, duration, window):
    # Epoch and period corrections are fitted relative to the discovery estimate.
    def residual(parameters):
        shift, drift, span, ingress, depth = parameters
        chunks = []
        for _, cycle, time, flux, noise, basis in events:
            offset = time - shift - drift * cycle
            model = -depth * shape(offset, span, ingress)
            delta = flux - model
            chunks.append((delta - basis @ (basis.T @ delta)) / noise)
        return np.concatenate(chunks)

    guess = [0, 0, duration, 0.25, 0.005]
    bounds = ([-window / 3, -0.02, duration / 3, 0.04, 0],
              [window / 3, 0.02, min(3 * duration, window * 1.45), 0.5, 0.25])
    fit = least_squares(residual, guess, bounds=bounds, x_scale="jac", loss="soft_l1")
    shift, drift, span, ingress, depth = fit.x
    return {"period_days": period + drift, "mid_eclipse_btjd": epoch + shift,
            "duration_days": span, "ingress_fraction": ingress, "depth_fraction": depth}


def event_measurements(events, fit, initial_period, initial_epoch):
    rows = []
    for sector, cycle, time, flux, noise, basis in events:
        mid = fit["mid_eclipse_btjd"] + cycle * fit["period_days"]
        offset = time + initial_epoch + cycle * initial_period - mid
        profile = shape(offset, fit["duration_days"], fit["ingress_fraction"])
        projected = profile - basis @ (basis.T @ profile)
        norm = np.linalg.norm(projected)
        if norm <= 0:
            continue
        depth = float(-projected @ flux / norm**2)
        baseline = basis @ (basis.T @ (flux + depth * profile))
        rows.append((sector, cycle, mid, depth, float(noise / norm), offset, flux - baseline))
    return rows


def aperture_photometry(path):
    with ZipFile(path) as archive:
        filename = next(name for name in archive.namelist() if name.endswith(".fits"))
        with fits.open(BytesIO(archive.read(filename))) as hdu:
            pixels = np.asarray(hdu[1].data["FLUX"], dtype=float)
            time = np.asarray(hdu[1].data["TIME"], dtype=float)
            quality = hdu[1].data["QUALITY"]
            cadence = float(hdu[1].header["TIMEDEL"])
    yy, xx = np.indices(pixels.shape[1:])
    radius = np.hypot(xx - (xx.shape[1] - 1) / 2, yy - (yy.shape[0] - 1) / 2)
    background = np.median(pixels[:, radius > 3.5], axis=1)
    flux = np.sum(pixels[:, radius <= 2.3] - background[:, None], axis=1)
    keep = (quality == 0) & np.isfinite(time + flux)
    return time[keep], flux[keep] / np.median(flux[keep]) - 1, cadence


def check_ffi(candidate, sectors, cache, output, fit):
    rows, products = [], {}
    for sector in sectors:
        query = urlencode({"ra": f'{candidate["RA_OBJ"]:.6f}', "dec": f'{candidate["DEC_OBJ"]:.6f}',
                           "x": 9, "y": 9, "units": "px", "sector": sector})
        url = f"https://mast.stsci.edu/tesscut/api/v0.1/astrocut?{query}"
        path = download(url, cache / f'{candidate["tic"]}-cut-s{sector:04d}.zip')
        time, flux, cadence = aperture_photometry(path)
        products[sector] = {"url": url, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                            "cadence_minutes": cadence * 1440}
        period, epoch = fit["period_days"], fit["mid_eclipse_btjd"]
        cycles = np.rint((time - epoch) / period).astype(int)
        for cycle in np.unique(cycles):
            middle = epoch + period * cycle
            window = (cycles == cycle) & (abs(time - middle) < 0.75)
            x, y = time[window] - middle, flux[window]
            if min(sum(x < -.4), sum(x > .4)) < 6 or sum(abs(x) < .15) < 6:
                continue
            basis, _ = np.linalg.qr(np.column_stack((np.ones_like(x), x, x*x)))
            residual = y - basis @ (basis.T @ y)
            exposure = np.linspace(-0.4, 0.4, 5) * cadence

            def measure(shift):
                profile = shape(x[:, None] - shift + exposure, fit["duration_days"],
                                fit["ingress_fraction"]).mean(axis=1)
                profile -= basis @ (basis.T @ profile)
                depth = -profile @ residual / (profile @ profile)
                return depth, np.sum((residual + depth * profile)**2)

            best = minimize_scalar(lambda shift: measure(shift)[1], bounds=(-.15, .15), method="bounded")
            rows.append((sector, cycle, middle, measure(0)[0], best.x, measure(best.x)[0]))
    with (output / "ffi-events.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("sector", "cycle", "predicted_mid_btjd", "fixed_depth_fraction",
                         "timing_offset_days", "recentered_depth_fraction"))
        writer.writerows(rows)
    if len({row[0] for row in rows}) < 2:
        raise ValueError("Multiple populated archival sectors are required for the long-baseline fit.")
    sector_centers = [
        np.median([(row[1], row[4]) for row in rows if row[0] == sector], axis=0)
        for sector in sorted({row[0] for row in rows})
    ]
    cycle, offset = np.asarray([*sector_centers, (0, 0)]).T
    drift, shift = np.polynomial.polynomial.polyfit(cycle, offset, 1)[::-1]
    rms = np.sqrt(np.mean((offset - shift - drift * cycle)**2))
    return products, {
        "period_days": fit["period_days"] + float(drift),
        "mid_eclipse_btjd": fit["mid_eclipse_btjd"] + float(shift),
        "sector_median_timing_residual_rms_minutes": float(rms * 1440),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("tic", type=int)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--extra-sectors", type=int, nargs="*", default=[])
    parser.add_argument("--ffi-sectors", type=int, nargs="*", default=[])
    parser.add_argument("--shift-days", type=float, default=0)
    parser.add_argument("--duration-days", type=float, required=True)
    parser.add_argument("--window-days", type=float, required=True)
    parser.add_argument("--source-coordinate", nargs=2, type=float, metavar=("RA", "DEC"),
                        help="ICRS degrees of the proposed source of the eclipse")
    parser.add_argument("--secondary-phase", type=float, help="Initial phase of a second eclipse")
    args = parser.parse_args()
    if not 0 < args.duration_days < args.window_days:
        parser.error("Use a positive eclipse duration shorter than the fit window.")
    report = json.loads(args.report.read_text())
    candidate = next(c for c in report["candidates"] if c["tic"] == args.tic and c["status"] == "repeated")
    first, second = report["discovery_sector"], report["validation_sector"]
    sources = {first: candidate["discovery_source"], second: candidate["validation_source"]}
    for sector in args.extra_sectors:
        if sector not in sources:
            sources[sector] = product_url(sector_products(sector, args.cache)[args.tic])
    events, files, curves = [], {}, {}
    period = candidate["spacing_days"]
    epoch = candidate["epoch_btjd"] + args.shift_days
    for sector, source in sources.items():
        path = download(source, args.cache / "fits" / source.rsplit("/", 1)[-1], verify_fits=True)
        _, data = photometry(path, native=True)
        if data is None:
            raise ValueError(f"No usable photometry in sector {sector}.")
        curves[sector] = data
        events.extend(events_for(data, sector, period, epoch, args.window_days))
        files[sector] = {"light_curve_url": source, "light_curve_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    discovery = [event for event in events if event[0] == first]
    if len(discovery) < 3:
        raise ValueError("At least three usable events are needed in the discovery sector.")
    first_fit = fit_eclipse(discovery, period, epoch, args.duration_days, args.window_days)
    holdout = event_measurements(events, first_fit, period, epoch)
    fit = fit_eclipse(events, period, epoch, args.duration_days, args.window_days)
    measurement = event_measurements(events, fit, period, epoch)
    secondary_holdout = None
    secondary_fit = None
    if args.secondary_phase is not None:
        second_epoch = epoch + args.secondary_phase * period
        secondary_events = [
            event for sector, curve in curves.items()
            for event in events_for(curve, sector, period, second_epoch, args.window_days)
        ]
        secondary_discovery = [event for event in secondary_events if event[0] == first]
        if len(secondary_discovery) < 3:
            raise ValueError("At least three usable secondary events are needed in the discovery sector.")
        secondary_discovery_fit = fit_eclipse(
            secondary_discovery, period, second_epoch, args.duration_days, args.window_days)
        secondary_holdout = event_measurements(
            secondary_events, secondary_discovery_fit, period, second_epoch)
        secondary_fit = fit_eclipse(
            secondary_events, period, second_epoch, args.duration_days, args.window_days)
    jackknife = [
        fit_eclipse(events[:i] + events[i+1:], period, epoch, args.duration_days, args.window_days)
        for i in range(len(events))
    ]
    result = {"tic": args.tic, "discovery_fit": first_fit, "joint_fit": fit, "files": files,
              "jackknife_se": {
                  key: float(np.std([row[key] for row in jackknife]) * np.sqrt(len(jackknife) - 1))
                  for key in fit
              }}
    if secondary_fit:
        result.update(secondary_fit=secondary_fit, secondary_discovery_fit=secondary_discovery_fit)
    catalog_coordinate = (candidate["RA_OBJ"], candidate["DEC_OBJ"])
    source_coordinate = args.source_coordinate or catalog_coordinate
    result.update(aperture_target_icrs_deg=catalog_coordinate, proposed_source_icrs_deg=source_coordinate)
    for sector, source in sources.items():
        pixel_url = source.replace("_lc.fits", "_tp.fits")
        path = download(pixel_url, args.cache / "pixels" / pixel_url.rsplit("/", 1)[-1], verify_fits=True)
        files[sector].update(pixel_url=pixel_url, pixel_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        eclipse_fits = {"localization": fit}
        if secondary_fit:
            eclipse_fits["secondary_localization"] = secondary_fit
        for key, eclipse in eclipse_fits.items():
            try:
                files[sector][key] = localize(
                    path, eclipse["period_days"], eclipse["mid_eclipse_btjd"], eclipse["duration_days"],
                    target_ra=source_coordinate[0], target_dec=source_coordinate[1], minimum_events=3)
            except ValueError as error:
                files[sector][f"{key}_error"] = str(error)
    args.output.mkdir(parents=True, exist_ok=True)
    if args.ffi_sectors:
        result["ffi_products"], result["long_baseline_ephemeris"] = check_ffi(
            candidate, args.ffi_sectors, args.cache, args.output, fit)
    (args.output / "fit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    evidence = {"events": measurement, "discovery_ephemeris_checks": holdout}
    if secondary_holdout:
        evidence["secondary_ephemeris_checks"] = secondary_holdout
    for name, rows in evidence.items():
        with (args.output / f"{name}.csv").open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(("sector", "cycle", "predicted_mid_btjd", "depth_fraction", "white_noise_depth_se"))
            writer.writerows(row[:5] for row in rows)
    plot(measurement, fit, args.output / "eclipses.png")
    write_phase(curves, result.get("long_baseline_ephemeris", fit), args.output)
    if args.source_coordinate:
        plot_localizations(files, source_coordinate, catalog_coordinate, args.output / "localizations.png")


def plot(rows, fit, output):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(rows), figsize=(2.7 * len(rows), 3.2),
                            sharex=True, sharey=True, layout="constrained")
    for axis, (sector, cycle, mid, _, _, offset, flux) in zip(np.atleast_1d(axes), rows):
        x = offset * 24
        axis.scatter(x, flux * 1000, s=3, alpha=.45, color="#266b9d", rasterized=True)
        y = -fit["depth_fraction"] * shape(offset, fit["duration_days"], fit["ingress_fraction"])
        axis.plot(x, y * 1000, color="#aa4f27", lw=1.2)
        axis.set_title(f"Sector {sector}\nBTJD {mid:.3f}", fontsize=9)
        axis.set_xlabel("Hours from mid-eclipse")
    np.atleast_1d(axes)[0].set_ylabel("Relative flux (ppt)")
    fig.savefig(output, dpi=160)
    plt.close(fig)


def write_phase(curves, fit, output):
    import matplotlib.pyplot as plt

    fig, (axis, timeline) = plt.subplots(2, 1, figsize=(8, 6.5), layout="constrained")
    phase_edges = np.linspace(-0.5, 0.5, 181)
    centers = (phase_edges[1:] + phase_edges[:-1]) / 2
    with (output / "photometry.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("sector", "btjd", "relative_flux", "relative_flux_error"))
        for sector, (time, flux, error) in curves.items():
            writer.writerows(zip(np.repeat(sector, len(time)), time, 1 + flux, error))
            phase = ((time - fit["mid_eclipse_btjd"]) / fit["period_days"] + .5) % 1 - .5
            index = np.searchsorted(phase_edges, phase, side="right") - 1
            mag = -2.5 * np.log10(np.maximum(1 + flux, np.finfo(float).tiny))
            medians = np.array([np.median(mag[index == i]) if np.any(index == i) else np.nan
                                for i in range(len(centers))])
            phases = np.r_[centers, centers + 1]
            differential_mmag = 1000 * np.tile(medians, 2)
            show = (phases > -.25) & (phases < 1.25)
            line, = axis.plot(phases[show], differential_mmag[show], lw=1, label=f"Sector {sector}")
            timeline.scatter(time, 1000 * mag, s=2, alpha=.24, color=line.get_color(), rasterized=True)
    axis.set(xlabel="Orbital phase", ylabel="Differential TESS magnitude (mmag)")
    axis.invert_yaxis()
    axis.legend()
    timeline.set(xlabel="BTJD (BJD TDB − 2,457,000)", ylabel="Differential TESS magnitude (mmag)")
    timeline.invert_yaxis()
    fig.savefig(output / "phase.png", dpi=180)
    plt.close(fig)


def plot_localizations(files, source, aperture, output):
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(7.4, 3.7), layout="constrained")
    projection = np.array([np.cos(np.deg2rad(source[1])), 1]) * 3600
    aperture_offset = (np.asarray(aperture) - source) * projection
    axis.scatter(*aperture_offset, marker="x", s=100, color="#794837", label="SPOC aperture target")
    axis.scatter(0, 0, marker="*", s=150, color="black", label="Proposed Gaia source")
    for key, label, marker, color in (
        ("localization", "Primary eclipse", "o", "#286a9d"),
        ("secondary_localization", "Secondary eclipse", "s", "#dc6b37"),
    ):
        localized = [file[key] for file in files.values() if key in file]
        if not localized:
            continue
        positions = np.asarray([(row["ra_deg"], row["dec_deg"]) for row in localized])
        x, y = ((positions - source) * projection).T
        axis.scatter(x, y, marker=marker, color=color, s=36, alpha=.8, label=label)
    axis.set(xlabel="Eastward offset (arcsec)", ylabel="Northward offset (arcsec)", aspect="equal")
    axis.set_ylim(-6, 6)
    axis.legend(loc="upper center", bbox_to_anchor=(.5, 1.27), ncols=2, fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
