"""Hermite layers of a series: centred He_1(z), ..., He_r(z) of its normal scores, and each layer's memory."""
from __future__ import annotations

import math

import numpy as np

from kmem import memory_kernels as mk
from kmem.calibration.selfnorm import SN_GRID
from kmem.legacy import mmdcal as M


def hermite_layers(x: np.ndarray, r: int = 2) -> np.ndarray:
    """Columns He_1(z), ..., He_r(z) of the normal scores of x, each centred. n x r."""
    h = mk.hermite_design(mk.normal_scores(x), r)[:, 1:]
    return h - h.mean(axis=0, keepdims=True)


def layer_memory(v: np.ndarray) -> float:
    """Local Whittle H of one layer series, m = floor(n^0.65) (mmdcal), clipped to the table range."""
    d = M.local_whittle_d(v, math.floor(v.size ** 0.65))
    return float(min(max(d + 0.5, SN_GRID[0]), SN_GRID[-1]))
