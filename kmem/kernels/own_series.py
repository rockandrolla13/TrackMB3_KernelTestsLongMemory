"""Own-series kernels: a Gaussian kernel on a transformed series with that series' own median bandwidth.

sym uses |x - median(x)|. datadiff uses the first difference x_{t+1} - x_t, so it loses one point.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from kmem import memory_kernels as mk
from kmem.legacy import mbcore as mb


@dataclass(frozen=True)
class OwnSeries:
    name: str
    series: Callable[[np.ndarray], np.ndarray]
    lag: int = 0
    kind: str = "own"

    def gram(self, x: np.ndarray, sigma: float | None = None, H: float | None = None) -> np.ndarray:
        """sigma is ignored: the bandwidth is the median heuristic of the own series."""
        u = self.series(x)
        d = np.abs(u[:, None] - u[None, :])
        return mk.double_center(mb.gauss_gram(d, mb.median_bandwidth(d)))


SYM = OwnSeries("sym", mk.even_map)
DATADIFF = OwnSeries("datadiff", mk.first_difference, lag=1)
