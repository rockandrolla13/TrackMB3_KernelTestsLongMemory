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
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from kmem.generators.split import split_series  # noqa: E402
from kmem.legacy import mbcore as mb  # noqa: E402
from kmem.legacy import mmdcal as M  # noqa: E402

ALPHA = 0.05


def q95(v: np.ndarray) -> float:
    return float(np.quantile(v, 1.0 - ALPHA))


def pboot(t: float, b: np.ndarray) -> float:
    return float((1.0 + np.sum(b >= t)) / (1.0 + b.size))


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

    d = np.abs(x[:, None] - x[None, :])
    sigma = M.median_heuristic(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    s = np.empty(n)
    s[:h1] = 1.0 / h1
    s[h1:] = -1.0 / (n - h1)
    stat = float(s @ k @ s)
    t_n = n * stat
    out: dict = {"mmd2": stat, "sigma": sigma, "t_n": t_n, "hm_diff": hm, "t_naive": t_naive}

    kc = k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()
    l = max(1, math.floor(n ** 0.5))
    w = M.multipliers("unc", n, h1, l, rng, B)
    a = w * s[:, None]
    boot = np.einsum("ij,ij->j", a, kc @ a)
    del kc, a, w
    out["p_unc_05"] = pboot(stat, boot)
    out["r_unc_05"] = stat / q95(boot)

    p = M.prefix2(k)
    del k
    v = np.sort(M.window_mmd(p, n, l))
    q_sub = float(v[int(math.ceil((1 - ALPHA) * v.size)) - 1])  # rep_e2 convention
    out["l_sqrt"] = l
    out["q95_sqrt"] = q_sub
    out["r_sub_sqrt"] = t_n / q_sub
    out["sn_T"] = float(M.sn_stat(p, np.array([0]), n)[0])
    v = np.sort(M.sn_stat(p, np.arange(n - l + 1), l))
    out["sn_q95_05"] = float(v[int(math.ceil((1 - ALPHA) * v.size)) - 1])
    out["r_sn_05"] = out["sn_T"] / out["sn_q95_05"]
    return out


# ---------------------------------------------------------------- Part 2: independence (E3 path)

def ind_rep(H: float, n: int, rng: np.random.Generator, B: int = 500) -> dict:
    x, y = mb.gen_series({"proc": "G4", "HX": H, "HY": H, "lam": None}, n, rng)
    out = mb.kind_tests(x, y, rng, B=B)
    out.update(mb.lind_tests(x, y))
    return out
