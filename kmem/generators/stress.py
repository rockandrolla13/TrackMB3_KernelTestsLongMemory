"""Building blocks of the EXP7 stress alternatives: a slow common factor and the scale-coupling normaliser."""
from __future__ import annotations

import math

import numpy as np

from kmem.legacy import mbcore as mb

FACTOR_PHI = 0.99
FACTOR_VAR = 3.0 / 7.0  # corr(X, Y) = v / (v + 1) = 0.3 with unit-variance noise


def factor(n: int, rng: np.random.Generator) -> np.ndarray:
    """Stationary AR(1), coefficient 0.99, variance FACTOR_VAR, Gaussian iid innovations."""
    e = rng.standard_normal(n)
    out = np.empty(n)
    prev = e[0]
    for t in range(n):
        if t:
            prev = FACTOR_PHI * prev + math.sqrt(1.0 - FACTOR_PHI ** 2) * e[t]
        out[t] = prev
    return math.sqrt(FACTOR_VAR) * out


def scale_norm(gamma: float) -> float:
    """sd of exp(gamma |xi|) Z with xi, Z independent standard normal: E exp(2 gamma |xi|) = 2 e^{2 gamma^2} Phi(2 gamma)."""
    return math.sqrt(2.0 * math.exp(2.0 * gamma * gamma) * float(mb.norm_cdf(np.array(2.0 * gamma))))


def flat_ma_norm(H: float) -> float:
    """sd of xi_{t-1} + xi_t + xi_{t+1} for unit fGn(H)."""
    g = mb.fgn_acov(H, 3)
    return math.sqrt(3.0 + 4.0 * g[1] + 2.0 * g[2])
