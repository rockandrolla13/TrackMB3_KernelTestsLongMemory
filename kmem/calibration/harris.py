"""Harris (2020) non-circular shift test, conservative and approximate forms, and the literal N = 19 form."""
from __future__ import annotations

import math

import numpy as np

from kmem.bandwidth import median_bw_fast
from kmem.calibration.common import ALPHA
from kmem.features import gauss_features, select_pairs

HARRIS_N_MIN = 19


def harris_big_n(n: int) -> int:
    """N near 0.05 n with N + 1 a multiple of 20, so alpha (N + 1) is a whole number at alpha = 5%. At least 19."""
    return max(HARRIS_N_MIN, 20 * round(0.05 * n / 20) - 1)


def harris_shifts(fx, fy, big_n: int) -> np.ndarray:
    """V_s for s = -N..N (index s + N). V_s = (1/D) sum_ij Kc_ij Lc^{(s)}_ij on the segments X[N:m-N], Y[s+N:s+m-N].

    fx, fy = compress(...) on all m points. Centring X within its segment makes centring Y unnecessary. Feature
    pairs are pruned as in rules_fast. Linear correlation by FFT, without wrap-around.
    """
    px, lx = fx
    py, ly = fy
    m = px.shape[0]
    d = m - 2 * big_n
    i, j = select_pairs(lx, ly)
    ux, iu = np.unique(i, return_inverse=True)
    uy, ju = np.unique(j, return_inverse=True)
    a = px[big_n:m - big_n][:, ux]
    a = a - a.mean(axis=0, keepdims=True)
    size = m + d
    c = np.fft.irfft(np.conj(np.fft.rfft(a, size, axis=0)[:, iu]) * np.fft.rfft(py[:, uy], size, axis=0)[:, ju], n=size, axis=0)
    c = c[:2 * big_n + 1]
    return (c * c).sum(axis=1) / d


def harris_ratio(v: np.ndarray, big_n: int, m_limit: float) -> float:
    """R = V_0 / (the M-th largest V_s, s != 0), M = ceil(m_limit). R > 1 iff m < m_limit, m = #{s != 0 : V_s >= V_0}."""
    k = math.ceil(m_limit - 1e-9)
    rest = np.delete(v, big_n)
    return float(v[big_n] / np.partition(rest, rest.size - k)[rest.size - k])


def harris_rules(fx, fy, n: int) -> dict:
    """rhc: conservative at N near 0.05 n (primary). rh19: conservative at N = 19. rha: approximate at N near 0.05 n."""
    big = harris_big_n(n)
    vb = harris_shifts(fx, fy, big)
    v19 = vb if big == HARRIS_N_MIN else harris_shifts(fx, fy, HARRIS_N_MIN)
    return {"rhc": harris_ratio(vb, big, ALPHA * (big + 1)), "rh19": harris_ratio(v19, HARRIS_N_MIN, ALPHA * (HARRIS_N_MIN + 1)),
            "rha": harris_ratio(vb, big, ALPHA * (2 * big + 1))}


def harris_strict19(x: np.ndarray, y: np.ndarray) -> dict:
    """The literal form at N = 19 for raw and inc1: every bandwidth is a function of the segment, so V_s is a stationary sequence.

    The main path takes the bandwidths from the whole series, which is much cheaper but not literally a function of the segment.
    """
    n, big = x.size, HARRIS_N_MIN
    xs = x[big:n - big]
    fxs = gauss_features(xs, median_bw_fast(xs))
    out = {}
    vals = {"raw": [], "inc1": []}
    ax = {"raw": fxs - fxs.mean(axis=0, keepdims=True), "inc1": (fxs[1:] - fxs[:-1]) - (fxs[1:] - fxs[:-1]).mean(axis=0, keepdims=True)}
    for s in range(-big, big + 1):
        ys = y[s + big:s + n - big]
        f = gauss_features(ys, median_bw_fast(ys))
        for k, by in (("raw", f), ("inc1", f[1:] - f[:-1])):
            c = ax[k].T @ by
            vals[k].append(float((c * c).sum() / by.shape[0]))
    for k in vals:
        out[f"xh19_{k}"] = harris_ratio(np.array(vals[k]), big, ALPHA * (big + 1))
    return out
