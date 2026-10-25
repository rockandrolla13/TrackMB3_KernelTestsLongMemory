"""EXP4 seasonal adjustments, each estimated from one observed series only (see DESIGN.md).

A1: subtract period-bin means (bin = t mod P), P = 50 known.
A2: A1, then divide by period-bin standard deviations (ddof = 1), P = 50 known.
A3: estimate P from the series (periodogram peak of x and of x^2), then A2 with that P.
Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import numpy as np

P_TRUE = 50


def bin_mean(x: np.ndarray, p: int) -> np.ndarray:
    b = np.arange(x.size) % p
    m = np.bincount(b, weights=x, minlength=p) / np.bincount(b, minlength=p)
    return x - m[b]


def bin_scale(r: np.ndarray, p: int) -> np.ndarray:
    """Divide residuals r by their bin sd (ddof = 1). r already has bin mean zero."""
    b = np.arange(r.size) % p
    cnt = np.bincount(b, minlength=p)
    ss = np.bincount(b, weights=r * r, minlength=p)
    sd = np.sqrt(ss / np.maximum(cnt - 1, 1))
    sd[sd <= 0] = 1.0
    return r / sd[b]


def estimate_period(x: np.ndarray) -> tuple[int, str]:
    """Period from the largest normalised periodogram ordinate of x or of x^2.

    Candidate Fourier indices j = ceil(10)...floor(n/2), i.e. periods n/j from n/10 down to 2.
    Each periodogram is divided by its sum over the candidate range, so the mean and the
    square compete on the same scale. Returns (round(n / j*), 'mean' or 'square').
    """
    n = x.size
    xc = x - x.mean()
    sq = xc * xc
    sq = sq - sq.mean()
    j = np.arange(10, n // 2 + 1)
    i1 = np.abs(np.fft.rfft(xc))[j] ** 2
    i2 = np.abs(np.fft.rfft(sq))[j] ** 2
    s1, s2 = i1 / i1.sum(), i2 / i2.sum()
    k1, k2 = int(np.argmax(s1)), int(np.argmax(s2))
    if s1[k1] >= s2[k2]:
        return int(round(n / j[k1])), "mean"
    return int(round(n / j[k2])), "square"


def adjust(x: np.ndarray, adj: str) -> tuple[np.ndarray, int]:
    """Return (adjusted series, period used). A0 returns x unchanged with period 0."""
    if adj == "A0":
        return x, 0
    if adj == "A1":
        return bin_mean(x, P_TRUE), P_TRUE
    if adj == "A2":
        return bin_scale(bin_mean(x, P_TRUE), P_TRUE), P_TRUE
    if adj == "A3":
        p, _ = estimate_period(x)
        return bin_scale(bin_mean(x, p), p), p
    raise ValueError(adj)
