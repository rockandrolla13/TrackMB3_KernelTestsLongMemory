"""Gram-matrix helpers, with the two double-centring orders named apart (code review CR-005).

Floating-point subtraction depends on order, so the two forms differ in the last bits. The saved HSIC and dCov
rows of EXP2 to EXP7 used centre_rows_first (legacy mbcore.double_center, also kmem.memory_kernels.double_center).
The saved split-MMD rows of EXP3, EXP5 and EXP7 used centre_cols_first. Keeping both reproduces both.
"""
from __future__ import annotations

import numpy as np

from kmem.legacy import mbcore as mb

gauss_gram = mb.gauss_gram
centre_rows_first = mb.double_center  # a - row means - column means + grand mean


def centre_cols_first(k: np.ndarray) -> np.ndarray:
    """H K H with the column means subtracted first: the inline form of the EXP3, EXP5 and EXP7 split-MMD code."""
    return k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()


def abs_diff(x: np.ndarray) -> np.ndarray:
    """|x_i - x_j| for all i, j: the distance matrix of one series."""
    return np.abs(x[:, None] - x[None, :])
