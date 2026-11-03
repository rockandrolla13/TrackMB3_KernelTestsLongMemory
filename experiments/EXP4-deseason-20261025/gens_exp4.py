"""EXP4 extra generator. SEAS-V-slow: seasonal scale whose period is 4n, so a window sees a quarter cycle.

X_t = (1 + 0.5 sin(2 pi (t + U) / (4n))) xi_t, U ~ Uniform{0, ..., 4n - 1}, independent per series.
Diagnostic for the stock arm, where a 2,000-duration window covers a small part of the trading day.
All other generators come unchanged from gens.py.
"""
from __future__ import annotations

import numpy as np

from kmem.generators.registry import get_process


def generate(gen: str, H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """One series; SEAS-V-slow and the EXP2 generators come from kmem.generators.registry (same draws)."""
    return get_process(gen)(H, n, rng)


def pair(gen: str, H: float, n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    return generate(gen, H, n, rng), generate(gen, H, n, rng)
