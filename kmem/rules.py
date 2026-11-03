"""Independence rules on one pair of kernels: bare statistic, wild bootstrap, circular shift, subsampling.

Three paths compute the same rules: rules_fast on compressed low-rank features, rules_full on transformed n x n
Gram matrices with the mbcore functions, and rules_distance on n x n matrices with two exact shortcuts.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.calibration.common import pboot, q95
from kmem.calibration.shift import shift_stats_sym
from kmem.calibration.subsample import sub_fast
from kmem.calibration.wild import wild_prefix
from kmem.features import select_pairs
from kmem.legacy import mbcore as mb


def rules_fast(fx, fy, w: np.ndarray, sub: bool = False) -> dict:
    """Bare statistic, wild bootstrap, circular shift, optional subsampling, from compressed features.

    fx, fy = compress(...). T = (1/m) sum_ij Kc_ij Lc_ij; w is m x B and centred.
    """
    px, lx = fx
    py, ly = fy
    m = px.shape[0]
    A = math.floor(0.1 * m)
    i, j = select_pairs(lx, ly)
    p = px[:, i] * py[:, j]
    t = float((p.sum(axis=0) ** 2).sum() / m)
    a = p.T @ w
    bw = (a * a).sum(axis=0) / m
    del a, p
    ux, iu = np.unique(i, return_inverse=True)
    uy, ju = np.unique(j, return_inverse=True)
    xc = np.fft.irfft(np.conj(np.fft.rfft(px[:, ux], axis=0)[:, iu]) * np.fft.rfft(py[:, uy], axis=0)[:, ju], n=m, axis=0)
    s = (xc * xc).sum(axis=1) / m
    del xc
    shifts = s[A:m - A + 1]
    q_w, q_s = q95(bw), q95(shifts)
    out = {"t": t, "m": m, "q": i.size, "p_wb": pboot(t, bw), "r_wb": t / q_w,
           "p_shift": pboot(s[0], shifts), "r_shift": s[0] / q_s, "check": float(s[0] - t),
           "r975_wb": t / float(np.quantile(bw, 0.975)), "r975_shift": s[0] / float(np.quantile(shifts, 0.975))}
    if sub:
        sb = sub_fast(px, py, i, j, math.floor(m ** 0.5), math.floor(0.1 * m))
        q_b = q95(sb)
        out.update({"p_sub": float(np.mean(sb >= t)), "rej_sub": int(t > q_b), "r_sub": t / q_b})
    return out


def rules_full(kc: np.ndarray, lc: np.ndarray, w: np.ndarray) -> dict:
    """Full path on transformed Gram matrices. Same statistic and rules as rules_fast (no subsampling)."""
    m = kc.shape[0]
    A = math.floor(0.1 * m)
    mh = kc * lc
    t = float(mh.sum() / m)
    bw = mb.wild_boot_stats(mh, w)
    del mh
    s = mb.shift_stats(kc, lc)
    shifts = s[A:m - A + 1]
    q_w, q_s = q95(bw), q95(shifts)
    return {"t": t, "m": m, "q": -1, "p_wb": pboot(t, bw), "r_wb": t / q_w,
            "p_shift": pboot(s[0], shifts), "r_shift": s[0] / q_s, "check": float(s[0] - t)}


def rules_distance(kc: np.ndarray, lc: np.ndarray, w: np.ndarray) -> dict:
    """rules_full with the two shortcuts. Same values up to rounding (unit_checks.py)."""
    m = kc.shape[0]
    A = math.floor(0.1 * m)
    mh = kc * lc
    t = float(mh.sum() / m)
    bw = wild_prefix(mh, w)
    del mh
    s = shift_stats_sym(kc, lc)
    shifts = s[A:m - A + 1]
    q_w, q_s = q95(bw), q95(shifts)
    return {"t": t, "m": m, "q": -1, "p_wb": pboot(t, bw), "r_wb": t / q_w,
            "p_shift": pboot(s[0], shifts), "r_shift": s[0] / q_s, "check": float(s[0] - t)}
