"""Split-MMD series: a null series, or one whose second half changes (EXP3, EXP5, EXP7)."""
from __future__ import annotations

import math

import numpy as np

from kmem.legacy import gens

ALIASES = {"shift": "mean"}  # EXP3 and EXP5 called the mean change "shift"


def split_series(gen: str, alt: str, H: float, n: int, par: float, rng: np.random.Generator) -> np.ndarray:
    """Split-MMD series. mean: add delta to the second half (EXP3 'shift'). scale: multiply the second
    half by par. shape: second half becomes sqrt(1 - a^2) xi + a (xi^2 - 1)/sqrt(2), same mean and variance."""
    alt = ALIASES.get(alt, alt)
    x = gens.generate(gen, H, n, rng)
    if alt == "null":
        return x
    x = x.copy()
    h = n // 2
    if alt == "mean":
        x[h:] += par
    elif alt == "scale":
        x[h:] *= par
    elif alt == "shape":
        xi = x[h:]
        x[h:] = math.sqrt(1.0 - par * par) * xi + par * (xi * xi - 1.0) / math.sqrt(2.0)
    else:
        raise ValueError(alt)
    return x
