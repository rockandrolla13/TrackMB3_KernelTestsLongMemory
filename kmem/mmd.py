"""Split MMD at the midpoint: the contrast, the uncentred wild rule, and EXP3's three-rule bundle (CR-001).

split_mmd_ext is EXP3 power_core.mmd_ext, which EXP5 re-typed. split_mmd_wild is the statistic and wild draws
shared with EXP7's exact path. The centring is centre_cols_first, as in all three experiments.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.bandwidth import median_mmdcal
from kmem.calibration.common import ALPHA, pboot, q95
from kmem.calibration.wild import multipliers_mmdcal
from kmem.gram import centre_cols_first
from kmem.legacy import mmdcal as M


def contrast(m: int, h: int) -> np.ndarray:
    s = np.empty(m)
    s[:h] = 1.0 / h
    s[h:] = -1.0 / (m - h)
    return s


def split_mmd_wild(k: np.ndarray, s: np.ndarray, a: np.ndarray) -> tuple[float, np.ndarray]:
    """MMD^2_V = s' K s and its uncentred wild draws a_b' Kc a_b, where a = w * s and Kc is K centred cols-first."""
    stat = float(s @ k @ s)
    kc = centre_cols_first(k)
    return stat, np.einsum("ij,ij->j", a, kc @ a)


def split_mmd_ext(x: np.ndarray, rng: np.random.Generator, B: int = M.B_WILD) -> dict:
    """Uncentred wild (l = floor(n^0.5)), sqrt subsampling, self-normalised block (l = floor(n^0.5)).

    Same formulas as mmdcal.rep_e2 for these three rules (rep_e2 also draws other multipliers first, so per-seed
    values differ).
    """
    n = x.size
    h1 = n // 2
    d = np.abs(x[:, None] - x[None, :])
    sigma = median_mmdcal(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    s = contrast(n, h1)
    l = max(1, math.floor(n ** 0.5))
    w = multipliers_mmdcal("unc", n, h1, l, rng, B)
    stat, boot = split_mmd_wild(k, s, w * s[:, None])
    del w
    t_n = n * stat
    out: dict = {"mmd2": stat, "sigma": sigma, "t_n": t_n, "p_unc_05": pboot(stat, boot), "r_unc_05": stat / q95(boot)}

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
