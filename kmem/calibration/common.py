"""Shared pieces of the calibration rules: the level, bootstrap p-values and 95% reference quantiles."""
from __future__ import annotations

import numpy as np

ALPHA = 0.05


def q95(v: np.ndarray) -> float:
    return float(np.quantile(v, 1.0 - ALPHA))


def pboot(t: float, b: np.ndarray) -> float:
    return float((1.0 + np.sum(b >= t)) / (1.0 + b.size))
