"""Analysis helpers shared by the experiment analyses (plan step 3.3, AR-DRY-003, CR-007).

The numerics are those every EXP3 to EXP7 analysis used, so their saved outputs regenerate byte for byte.
"""
from __future__ import annotations

import math
import os

import numpy as np

from kmem.calibration.common import ALPHA, NORMAL_975

Z = NORMAL_975
NBOOT = 2000  # bootstrap resamples of replications
SEED = 42  # the bootstrap seed every saved analysis used


def require_input(path: str, what: str) -> str:
    """A declared input from another experiment (experiments/EXPn/spec.json, inputs). Missing is an error,
    never a silent skip (code review CR-009)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"declared input missing: {os.path.normpath(path)} ({what})")
    return path


def wilson(k: int, m: int, z: float = Z) -> tuple[float, float]:
    """Wilson 95% interval for k rejections in m replications; (nan, nan) when m = 0."""
    if m == 0:
        return float("nan"), float("nan")
    p = k / m
    d = 1 + z * z / m
    c = (p + z * z / (2 * m)) / d
    h = z * math.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / d
    return c - h, c + h


def size_class(k: int, m: int) -> str:
    """valid (Wilson interval holds 5%, or upper bound <= 7.5%), mild (rate <= 10%) or over."""
    p = k / m
    lo, hi = wilson(k, m)
    if lo <= ALPHA <= hi or hi <= 0.075:
        return "valid"
    return "mild" if p <= 0.10 else "over"


def cstar(v: np.ndarray) -> float:
    """Size-adjusted cut-off: the order statistic with floor(0.05 m) null values above it (no ties)."""
    v = np.sort(v)
    return float(v[v.size - int(math.floor(ALPHA * v.size)) - 1])


def cell_rng(seed_index: int, seed: int = SEED) -> np.random.Generator:
    """Bootstrap stream of one cell, for new analyses (decision D12): SeedSequence(seed, spawn_key=(seed_index,)).

    The saved EXP3, EXP6 and EXP7 analyses drew every cell's bootstrap from one default_rng(SEED) in cell order,
    so adding or removing a cell shifted every later interval (AR-PAR-001, CR-010). They keep that stream, so
    their saved intervals still reproduce. With cell_rng, a cell's intervals depend on that cell alone.
    """
    return np.random.default_rng(np.random.SeedSequence(entropy=seed, spawn_key=(seed_index,)))


def boot_idx(m: int, rng: np.random.Generator, nboot: int = NBOOT) -> np.ndarray:
    return rng.integers(0, m, size=(nboot, m))


def boot_power(nv: np.ndarray, av: np.ndarray, bi_n: np.ndarray, bi_a: np.ndarray) -> np.ndarray:
    """Size-adjusted power on each bootstrap resample: c* recomputed from the resampled null values."""
    vn = np.sort(nv[bi_n], axis=1)
    kk = vn.shape[1] - int(math.floor(ALPHA * vn.shape[1])) - 1
    return np.mean(av[bi_a] > vn[:, kk][:, None], axis=1)


def mcnemar_p(a: np.ndarray, b: np.ndarray) -> tuple[float, int, int]:
    """Exact two-sided McNemar p-value on paired 0/1 vectors, with the discordant counts."""
    n10 = int(np.sum(a & ~b))
    n01 = int(np.sum(~a & b))
    m = n10 + n01
    if m == 0:
        return 1.0, n10, n01
    k = min(n10, n01)
    return min(1.0, 2 * sum(math.comb(m, i) for i in range(k + 1)) / 2 ** m), n10, n01
