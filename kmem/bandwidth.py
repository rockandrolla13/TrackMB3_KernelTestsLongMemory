"""Bandwidths of the Gaussian kernel, named by algorithm.

Three algorithms share the name "median bandwidth" in the code base (code review CR-020):
- median_full_matrix: mbcore.median_bandwidth, the off-diagonal median of a full n x n distance matrix,
  with a positive-distance fallback. Used by the HSIC and dCov rules of EXP2 to EXP7.
- median_mmdcal: mmdcal.median_heuristic, used by the split-MMD rules of EXP2, EXP3, EXP5 and EXP7.
- median_bw_fast: median pairwise distance by bisection on the sorted series (the EXP7 fast path). It equals
  median_full_matrix to about 1e-13 relative when the median is positive, and falls back to it otherwise.
"""
from __future__ import annotations

import numpy as np

from kmem.legacy import mbcore as mb
from kmem.legacy import mmdcal as M

median_full_matrix = mb.median_bandwidth
median_mmdcal = M.median_heuristic


def median_bw_fast(x: np.ndarray) -> float:
    """Median pairwise distance by bisection on the sorted series, about 1e-13 relative to the range.

    Equals mbcore.median_bandwidth up to that tolerance when the median is positive. O(n log n).
    """
    s = np.sort(x)
    n = s.size
    pairs = n * (n - 1) // 2
    base = n * (n + 1) // 2

    def kth(k: int) -> float:
        lo, hi = 0.0, float(s[-1] - s[0])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if int(np.searchsorted(s, s + mid, side="right").sum()) - base >= k:
                hi = mid
            else:
                lo = mid
        return hi

    if pairs % 2:
        med = kth(pairs // 2 + 1)
    else:
        med = 0.5 * (kth(pairs // 2) + kth(pairs // 2 + 1))
    if med <= 0.0:
        return mb.median_bandwidth(np.abs(x[:, None] - x[None, :]))
    return med
