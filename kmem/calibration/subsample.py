"""Subsampling (Betken-Dehling sampling windows) on low-rank features."""
from __future__ import annotations

import numpy as np

from kmem.features import window_sum


def sub_fast(px: np.ndarray, py: np.ndarray, i: np.ndarray, j: np.ndarray, l: int, dn: int) -> np.ndarray:
    """Sampling-window statistics T_{l,k}, k = 0..m-l-dn-1, X block at k, Y block at k + dn (mbcore rule)."""
    m = px.shape[0]
    nw = m - l - dn
    a, b = px[:m - dn], py[dn:]
    ux, iu = np.unique(i, return_inverse=True)
    uy, ju = np.unique(j, return_inverse=True)
    sa, sb = window_sum(a[:, ux], l), window_sum(b[:, uy], l)
    cab = window_sum(a[:, i] * b[:, j], l) - sa[:, iu] * sb[:, ju] / l
    return ((cab * cab).sum(axis=1) / l)[:nw]
