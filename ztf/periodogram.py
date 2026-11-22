"""Evaluate a floating-mean Lomb–Scargle periodogram on a regular grid."""

from nifty_ls import lombscargle


def periodogram(time, magnitude, error, frequencies):
    # Each scan worker owns a light curve, so let the scheduler parallelize
    # across curves instead of spawning FFT threads for each one.
    return lombscargle(
        time, magnitude, error,
        fmin=frequencies[0], fmax=frequencies[-1], Nf=len(frequencies),
        backend="finufft", nthreads=1, finufft_kwargs={"eps": 1e-12, "upsampfac": 2.0},
    ).power
