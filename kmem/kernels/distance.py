"""Distance-kernel maps: the increment kernel built on K = -|x_i - x_j| (dcov_inc1 in EXP7)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kmem import memory_kernels as mk


@dataclass(frozen=True)
class DistanceIncrement:
    lag: int
    name: str = ""
    kind: str = "distance"

    def __post_init__(self):
        if not self.name:
            object.__setattr__(self, "name", f"dcov_inc{self.lag}")

    def from_distances(self, d: np.ndarray) -> np.ndarray:
        """H D (-d) D' H from a precomputed distance matrix d = |x_i - x_j|."""
        return mk.increment(-d, self.lag)

    def gram(self, x: np.ndarray, sigma: float | None = None, H: float | None = None) -> np.ndarray:
        return mk.increment(mk.distance_kernel(x), self.lag)
