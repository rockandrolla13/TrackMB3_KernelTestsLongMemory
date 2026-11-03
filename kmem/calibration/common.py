"""Shared pieces of the calibration rules: the level, bootstrap p-values and 95% reference quantiles."""
from __future__ import annotations

import numpy as np

from kmem.legacy import mbcore as mb

ALPHA = 0.05
NORMAL_975 = 1.959963984540054  # the two-sided 5% normal cut-off, as written in mbcore.lind_tests


def t975(nu: int) -> float:
    """Two-sided 5% Student-t cut-off with nu degrees of freedom: the value mbcore.lind_tests caches in its
    private _TQ dict (same function, same argument), without depending on lind_tests having run first."""
    return mb.t_quantile(0.975, nu)


def q95(v: np.ndarray) -> float:
    return float(np.quantile(v, 1.0 - ALPHA))


def pboot(t: float, b: np.ndarray) -> float:
    return float((1.0 + np.sum(b >= t)) / (1.0 + b.size))
