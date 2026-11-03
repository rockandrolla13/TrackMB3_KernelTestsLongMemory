"""Seasonal scales: the EXP2 seasonal scale with a random phase, and EXP4's slow scale with period 4n."""
from __future__ import annotations

import numpy as np

from kmem.legacy import gens
from kmem.legacy import mbcore as mb


def seas_scale(n: int, rng: np.random.Generator) -> np.ndarray:
    u = int(rng.integers(gens.PERIOD))
    return 1.0 + gens.SEAS_B * gens.season(n, u)


def seas_v_slow(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """EXP4 SEAS-V-slow, same draw order as gens_exp4.generate: period 4n, phase independent per series."""
    p = 4 * n
    u = int(rng.integers(p))
    return (1.0 + gens.SEAS_B * np.sin(2.0 * np.pi * (np.arange(n) + u) / p)) * mb.fgn(H, n, rng)
