"""Independence rules on one pair of kernels: the bare statistic plus calibration rules from kmem.calibration.registry.

Three paths compute the same rules: rules_fast on compressed low-rank features, rules_full on transformed n x n
Gram matrices with the mbcore functions, and rules_distance on n x n matrices with two exact shortcuts.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.calibration import registry as cal
from kmem.calibration.shift import shift_stats_sym
from kmem.calibration.subsample import sub_fast
from kmem.calibration.wild import wild_prefix
from kmem.features import select_pairs
from kmem.legacy import mbcore as mb


def _decide(out: dict, names, draws: dict) -> dict:
    """Add each named calibration rule's columns, computed from the draws it names."""
    for name in names:
        rule = cal.get(name)
        out.update(rule.decide(*draws[rule.draws]))
    return out


def rules_fast(fx, fy, w: np.ndarray, sub: bool = False, rules=None) -> dict:
    """Bare statistic and calibration rules from compressed features.

    fx, fy = compress(...). T = (1/m) sum_ij Kc_ij Lc_ij; w is m x B and centred. rules: names in
    kmem.calibration.registry; the default is wild, shift and their 97.5% ratios, plus subsampling if sub.
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
    names = tuple(rules) if rules is not None else cal.FAST_RULES + (("subsample",) if sub else ())
    draws = {"wild": (t, bw), "shift": (s[0], shifts)}
    if any(cal.get(nm).draws == "subsample" for nm in names):
        draws["subsample"] = (t, sub_fast(px, py, i, j, math.floor(m ** 0.5), math.floor(0.1 * m)))
    return _decide({"t": t, "m": m, "q": i.size, "check": float(s[0] - t)}, names, draws)


def rules_full(kc: np.ndarray, lc: np.ndarray, w: np.ndarray, rules=cal.FULL_RULES) -> dict:
    """Full path on transformed Gram matrices. Same statistic and rules as rules_fast (no subsampling)."""
    m = kc.shape[0]
    A = math.floor(0.1 * m)
    mh = kc * lc
    t = float(mh.sum() / m)
    bw = mb.wild_boot_stats(mh, w)
    del mh
    s = mb.shift_stats(kc, lc)
    shifts = s[A:m - A + 1]
    return _decide({"t": t, "m": m, "q": -1, "check": float(s[0] - t)}, rules, {"wild": (t, bw), "shift": (s[0], shifts)})


def rules_distance(kc: np.ndarray, lc: np.ndarray, w: np.ndarray, rules=cal.FULL_RULES) -> dict:
    """rules_full with the two shortcuts. Same values up to rounding (unit_checks.py)."""
    m = kc.shape[0]
    A = math.floor(0.1 * m)
    mh = kc * lc
    t = float(mh.sum() / m)
    bw = wild_prefix(mh, w)
    del mh
    s = shift_stats_sym(kc, lc)
    shifts = s[A:m - A + 1]
    return _decide({"t": t, "m": m, "q": -1, "check": float(s[0] - t)}, rules, {"wild": (t, bw), "shift": (s[0], shifts)})
