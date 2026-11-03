"""Split-MMD series: a null series, or one whose second half changes (EXP3, EXP5, EXP7)."""
from __future__ import annotations

import numpy as np

from kmem.generators.registry import get_change, get_process

ALIASES = {"shift": "mean"}  # EXP3 and EXP5 called the mean change "shift"


def split_series(gen: str, alt: str, H: float, n: int, par: float, rng: np.random.Generator) -> np.ndarray:
    """Series of process gen. alt "null" leaves it; otherwise the registered change alt (mean: add par; scale: multiply
    by par; shape: sqrt(1 - a^2) xi + a (xi^2 - 1)/sqrt(2), same mean and variance) acts on the second half."""
    alt = ALIASES.get(alt, alt)
    x = get_process(gen)(H, n, rng)
    if alt == "null":
        return x
    change = get_change(alt)
    x = x.copy()
    h = n // 2
    x[h:] = change(x[h:], par)
    return x
