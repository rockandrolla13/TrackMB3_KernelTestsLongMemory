"""The kernel protocol and the registry (plan step 3.2, findings AR-EXT-001 and CR-006).

A kernel turns one series x into a centred Gram matrix on n - lag points. It has one of three kinds:

- "feature": a linear map A on the Gaussian features F of x (bandwidth sigma from x), so A K A' = (A F)(A F)'.
  It supplies features(f, x, H) for the fast path and gram(x, sigma, H) for the exact path.
- "own": a Gaussian kernel on its own series u = series(x), with u's own median bandwidth.
- "distance": a map on the distance kernel -|x_i - x_j|, supplied as from_distances(d) with d = |x_i - x_j|.

Which calibration rules run for which kernel is an experiment's choice, not a kernel property.

To add a kernel: write a small class with these attributes in a new module and call register(MyKernel(...)).
"""
from __future__ import annotations

from typing import Protocol

import numpy as np


class Kernel(Protocol):
    name: str
    kind: str  # "feature", "own" or "distance"
    lag: int  # number of points lost: the Gram matrix is (n - lag) x (n - lag)

    def gram(self, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
        """Exact-path transformed Gram matrix of one series."""
        ...


_REGISTRY: dict[str, Kernel] = {}


def register(kernel: Kernel) -> Kernel:
    """Add a kernel under its name. Re-registering a name with a different object is an error."""
    old = _REGISTRY.get(kernel.name)
    if old is not None and old != kernel:
        raise ValueError(f"kernel {kernel.name!r} is already registered as {old!r}")
    if kernel.kind not in ("feature", "own", "distance"):
        raise ValueError(f"kernel {kernel.name!r}: unknown kind {kernel.kind!r}")
    _REGISTRY[kernel.name] = kernel
    return kernel


def get(name: str) -> Kernel:
    try:
        return _REGISTRY[name]
    except KeyError:
        raise ValueError(f"unknown kernel {name!r}; registered: {sorted(_REGISTRY)}") from None


def names() -> tuple[str, ...]:
    return tuple(_REGISTRY)
