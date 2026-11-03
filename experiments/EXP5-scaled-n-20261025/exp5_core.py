"""EXP5 rule bundles. Imports the shared EXP2 core from kmem.legacy (unmodified).

Part 1 (split MMD): mmd_ext is the EXP3 power_core.mmd_ext code (same formulas, same rng
consumption order), plus the half-mean difference and a naive two-sample t statistic for the
linear benchmark. The shift is scaled in null-SD units: delta_H = c * sd_half(H, n).

Part 2 (independence): ind_rep is the EXP1 E3 code path (mbcore.gen_series with proc G4,
mbcore.kind_tests, mbcore.lind_tests), unchanged.
Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.calibration.common import ALPHA, pboot, q95  # noqa: F401  (re-exported)
from kmem.generators.split import split_series
from kmem.legacy import mbcore as mb
from kmem.legacy import mmdcal as M
from kmem.mmd import split_mmd_ext


def sd_half(H: float, n: int) -> float:
    """Null SD of mean(second half) - mean(first half) for unit-variance fGn(H), n = 2m."""
    m = n // 2
    return math.sqrt(4.0 * m ** (2 * H) - (2.0 * m) ** (2 * H)) / m


def sd_half_exact(H: float, n: int) -> float:
    """Same quantity from the fGn autocovariance (check of the closed form)."""
    m = n // 2
    g = mb.fgn_acov(H, n)
    s = np.r_[np.full(m, -1.0 / m), np.full(n - m, 1.0 / (n - m))]
    i = np.arange(n)
    S = g[np.abs(i[:, None] - i[None, :])]
    return float(math.sqrt(s @ S @ s))


# ---------------------------------------------------------------- Part 1: split MMD

def mmd_series(H: float, n: int, delta: float, rng: np.random.Generator) -> np.ndarray:
    """G0 fGn (EXP2/EXP3 generator); delta added to the second half (delta = 0 is the null)."""
    return split_series("G0", "mean" if delta != 0.0 else "null", H, n, delta, rng)


def mmd_ext(x: np.ndarray, rng: np.random.Generator, B: int = M.B_WILD) -> dict:
    """EXP3 mmd_ext (uncentred wild l = floor(n^0.5), sqrt subsampling, self-normalised block)
    plus the half-mean difference hm_diff and the naive Welch t statistic t_naive."""
    n = x.size
    h1 = n // 2
    a1, a2 = x[:h1], x[h1:]
    hm = float(a2.mean() - a1.mean())
    t_naive = hm / math.sqrt(a1.var(ddof=1) / a1.size + a2.var(ddof=1) / a2.size)
    out: dict = {"hm_diff": hm, "t_naive": t_naive}
    out.update(split_mmd_ext(x, rng, B))  # no random draw above, so the order does not matter
    return out


# ---------------------------------------------------------------- Part 2: independence (E3 path)

def ind_rep(H: float, n: int, rng: np.random.Generator, B: int = 500) -> dict:
    x, y = mb.gen_series({"proc": "G4", "HX": H, "HY": H, "lam": None}, n, rng)
    out = mb.kind_tests(x, y, rng, B=B)
    out.update(mb.lind_tests(x, y))
    return out
