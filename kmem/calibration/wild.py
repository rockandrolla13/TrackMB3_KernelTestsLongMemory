"""Wild bootstrap: block-Markov Rademacher multipliers shared across kernels, and w' M w for piecewise-constant w."""
from __future__ import annotations

import math

import numpy as np

from kmem.legacy import gens
from kmem.legacy import mbcore as mb
from kmem.legacy import mmdcal as M

INC_LONG = gens.PERIOD  # lag of the long increment kernel inc50; its chain may need its own ln

# Two multiplier laws share the name markov_rademacher in the legacy core (code review CR-020).
markov_rademacher_mbcore = mb.markov_rademacher  # HSIC and dCov wild bootstrap of EXP2 to EXP7
multipliers_mmdcal = M.multipliers  # split-MMD wild bootstrap forms of mmdcal (EXP2, EXP3, EXP5, EXP7)


def wild_prefix(mh: np.ndarray, w: np.ndarray) -> np.ndarray:
    """(1/m) w_b' M w_b per column, for piecewise-constant multipliers. Equals mbcore.wild_boot_stats.

    With jumps J_a at positions f_a (the first row counts as a jump), w' M w = sum_ab J_a J_b Q(f_a, f_b),
    where Q(i, j) is the sum of M over s >= i, t >= j. The chains flip about m / (2 ln) times.
    """
    m, b = w.shape
    q = np.cumsum(np.cumsum(mh[::-1, ::-1], axis=0), axis=1)[::-1, ::-1]
    jumps = np.empty_like(w)
    jumps[0] = w[0]
    jumps[1:] = w[1:] - w[:-1]
    chain, pos = np.nonzero(jumps.T)
    cuts = np.searchsorted(chain, np.arange(b + 1))
    out = np.empty(b)
    for c in range(b):
        f = pos[cuts[c]:cuts[c + 1]]
        j = jumps[f, c]
        out[c] = j @ q[np.ix_(f, f)] @ j
    return out / m


class Multipliers:
    """One block-Markov chain per replication. Shorter kernels use its first m rows, re-centred.

    The main chain is drawn first, with ln = floor(n^0.5), as in EXP6, so raw reproduces EXP6.
    ln = floor(m^0.5) can differ only for inc50, and then one extra chain is drawn. Draw order does not
    depend on the kernel set, so any kernel subset sees the same multipliers.
    """

    def __init__(self, n: int, B: int, rng: np.random.Generator):
        self.n = n
        self.ln = math.floor(n ** 0.5)
        self.w = mb.markov_rademacher(n, B, self.ln, rng)
        self.ln_long = math.floor((n - INC_LONG) ** 0.5)
        self.w_long = self.w if self.ln_long == self.ln else mb.markov_rademacher(n, B, self.ln_long, rng)
        if math.floor((n - 1) ** 0.5) != self.ln:
            raise ValueError(f"n = {n}: floor(sqrt(n - 1)) differs from floor(sqrt(n))")

    def get(self, m: int) -> np.ndarray:
        if m == self.n:
            return self.w
        base = self.w_long if m == self.n - INC_LONG else self.w
        v = base[:m]
        return v - v.mean(axis=0, keepdims=True)
