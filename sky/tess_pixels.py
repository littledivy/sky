"""Check whether an eclipse falls on the target's TESS pixels."""

import warnings

from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from astropy.wcs import FITSFixedWarning, WCS
import numpy as np


def eclipse_image(path, period, epoch, duration, *, minimum_events=6):
    """Build a difference image using paired flux measurements around each eclipse."""
    with fits.open(path, memmap=True) as hdu:
        table = hdu[1].data
        time = np.asarray(table["TIME"], dtype=float)
        quality = table["QUALITY"] == 0
        pixels = table["FLUX"]
        header = hdu[2].header.copy()
        valid = np.flatnonzero(np.isfinite(time) & quality)
        offset = time[valid] - epoch
        cycle = np.rint(offset / period).astype(int)
        phase_time = offset - cycle * period
        inner = abs(phase_time) < duration * 0.35
        outer = (abs(phase_time) > duration * 0.85) & (
            abs(phase_time) < min(duration * 2, period * 0.45))
        events = []
        for index in np.unique(cycle[inner]):
            is_cycle = cycle == index
            before = is_cycle & outer & (phase_time < 0)
            after = is_cycle & outer & (phase_time > 0)
            eclipse = is_cycle & inner
            if eclipse.sum() < 3 or before.sum() < 3 or after.sum() < 3:
                continue
            left = np.nanmedian(pixels[valid[before]], axis=0)
            right = np.nanmedian(pixels[valid[after]], axis=0)
            inside = np.nanmedian(pixels[valid[eclipse]], axis=0)
            events.append((left + right) / 2 - inside)
        if len(events) < minimum_events:
            raise ValueError("Too few eclipses with pixel coverage on both sides.")
    events = np.stack(events)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        image = np.nanmedian(events, axis=0)
        noise = 1.4826 * np.nanmedian(abs(events - image), axis=0) / np.sqrt(len(events))
    return image, noise, events, header


def positive_centroid(image, peak, radius=2.0):
    y, x = np.indices(image.shape)
    selected = (x - peak[1])**2 + (y - peak[0])**2 <= radius**2
    weight = np.where(selected & np.isfinite(image), np.maximum(image, 0), 0)
    total = weight.sum()
    if total <= 0:
        raise ValueError("No positive signal near the difference-image peak.")
    return np.array([(x * weight).sum(), (y * weight).sum()]) / total


def localize(path, period, epoch, duration, *, target_ra, target_dec, minimum_events=6):
    image, noise, events, header = eclipse_image(
        path, period, epoch, duration, minimum_events=minimum_events)
    peak = np.unravel_index(np.nanargmax(image), image.shape)
    if image[peak] <= 0 or not np.isfinite(noise[peak]) or noise[peak] <= 0:
        raise ValueError("The eclipse difference image has no measurable positive peak.")
    center = positive_centroid(image, peak)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header).celestial
    sky = wcs.pixel_to_world(*center)
    target = SkyCoord(target_ra, target_dec, unit="deg")
    rng = np.random.default_rng(107)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        samples = np.array([
            positive_centroid(np.nanmedian(events[rng.integers(len(events), size=len(events))], axis=0), peak)
            for _ in range(80)
        ])
    scatter_pixels = np.sqrt(np.sum(np.var(samples, axis=0, ddof=1)))
    return {
        "eclipses": len(events),
        "ra_deg": float(sky.ra.deg),
        "dec_deg": float(sky.dec.deg),
        "distance_from_target_arcsec": float(sky.separation(target).arcsec),
        "centroid_bootstrap_scatter_pixels": float(scatter_pixels),
        "peak_signal_over_scatter": float(image[peak] / noise[peak]),
        "centroid_xy": center.tolist(),
        "difference_image": [[float(value) if np.isfinite(value) else None for value in row]
                             for row in image],
    }
