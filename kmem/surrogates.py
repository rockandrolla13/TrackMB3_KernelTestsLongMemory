"""Surrogate series for a known independence null on real data (EXP8, DESIGN.md).

iaaft keeps a series' exact values and, approximately, its periodogram, so its memory. Its phases come from rng,
so it carries no time alignment with any other series. iaaft_drift also keeps the series' slow scale drift,
fitted as a linear trend in log segment variance.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.legacy import mmdcal as M

IAAFT_ITERS = 100
SEGMENTS = 8
WHITTLE_POWER = 0.65


def iaaft(y: np.ndarray, rng: np.random.Generator, iters: int = IAAFT_ITERS) -> np.ndarray:
    """Iterated amplitude-adjusted Fourier transform surrogate. Ends on the rank remap, so the output is a
    permutation of y (exact values); the amplitude spectrum matches approximately."""
    y = np.asarray(y, dtype=float)
    n = y.size
    values = np.sort(y)
    amp = np.abs(np.fft.rfft(y))
    s = rng.permutation(y)
    for _ in range(iters):
        spec = np.fft.rfft(s)
        s = np.fft.irfft(amp * np.exp(1j * np.angle(spec)), n)
        s = values[np.argsort(np.argsort(s, kind="stable"), kind="stable")]
    return s


def segment_positions(n: int, segments: int = SEGMENTS) -> np.ndarray:
    """Position of each time point on the 0..1 segment scale used by variance_trend (segment j sits at j/(S-1))."""
    return (((np.arange(n) + 0.5) * segments / n) - 0.5) / (segments - 1)


def variance_trend(w: np.ndarray, segments: int = SEGMENTS) -> tuple[float, float]:
    """OLS slope and intercept of log segment variance on segment position 0..1.

    Same arithmetic as experiments/STOCK-DIURNAL-20261025/diagnostics.py segment_stats (the slope). When n is not
    a multiple of segments (n = 500), the last n mod segments points are dropped from the fit."""
    w = np.asarray(w, dtype=float)
    v = w[: w.size - w.size % segments].reshape(segments, -1).var(axis=1)
    slope, intercept = np.polyfit(np.linspace(0, 1, segments), np.log(v), 1)
    return float(slope), float(intercept)


def iaaft_drift(y: np.ndarray, rng: np.random.Generator, iters: int = IAAFT_ITERS,
                segments: int = SEGMENTS) -> np.ndarray:
    """IAAFT of y divided by its fitted local scale, then the scale multiplied back (DESIGN.md S2)."""
    y = np.asarray(y, dtype=float)
    slope, intercept = variance_trend(y, segments)
    scale = np.exp(0.5 * (intercept + slope * segment_positions(y.size, segments)))
    m = y.mean()
    return m + scale * iaaft((y - m) / scale, rng, iters)


def hurst_lw(x: np.ndarray) -> float:
    """Local Whittle estimate of H with bandwidth floor(n^0.65), as in the EXP7 layer check."""
    return M.local_whittle_d(np.asarray(x, dtype=float), math.floor(x.size ** WHITTLE_POWER)) + 0.5
