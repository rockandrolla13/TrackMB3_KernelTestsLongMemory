"""Independence-arm pairs (X, Y): nulls and alternatives of EXP2, EXP3, EXP6 and EXP7, in their original draw order."""
from __future__ import annotations

import math

import numpy as np

from kmem.generators.seasonal import seas_scale, seas_v_slow
from kmem.generators.stress import factor, flat_ma_norm, scale_norm
from kmem.legacy import gens
from kmem.legacy import mbcore as mb


def pair(gen: str, alt: str, HX: float, HY: float, n: int, par: float, rng: np.random.Generator):
    """Independence-arm pair. Null generators and lin/sq follow EXP2 and EXP3 draw orders.

    gen: G0, an EXP2 process arm (null only), SEAS-V-slow, FACTOR or SCALE (matched nulls).
    alt: null, lin, sq (gen G0 or SEAS-V, as in EXP3 and EXP6), lag, flat, factor, scale (gen G0, EXP7).
    """
    if alt == "null":
        if gen == "G0":
            x = mb.fgn(HX, n, rng)
            return x, mb.fgn(HY, n, rng)
        if gen == "SEAS-V-slow":
            x = seas_v_slow(HX, n, rng)
            return x, seas_v_slow(HY, n, rng)
        if gen == "FACTOR":  # matched null of the common slow factor: two independent factors
            x = factor(n, rng) + mb.fgn(HX, n, rng)
            return x, factor(n, rng) + mb.fgn(HY, n, rng)
        if gen == "SCALE":  # matched null of the scale coupling: Y's scale comes from an independent copy
            x = mb.fgn(HX, n, rng)
            xp = mb.fgn(HX, n, rng)
            return x, np.exp(par * np.abs(xp)) * mb.fgn(HY, n, rng) / scale_norm(par)
        if HX != HY:
            raise ValueError("process arms use equal H")
        return gens.pair(gen, HX, n, rng)
    if HX != HY:
        raise ValueError("alternatives use equal H")
    H = HX
    if alt in ("lin", "sq"):  # EXP3 power_core.ind_pair draw order; for G0 it is also EXP7's
        if gen not in ("G0", "SEAS-V"):
            raise ValueError(f"{alt} is defined for G0 and SEAS-V, got {gen}")
        if gen == "SEAS-V":
            sx = seas_scale(n, rng)
        xi = mb.fgn(H, n, rng)
        if gen == "SEAS-V":
            sy = seas_scale(n, rng)
        w = mb.fgn(H, n, rng)
        if alt == "lin":
            eta = par * xi + math.sqrt(1.0 - par * par) * w
        else:
            eta = par * (xi * xi - 1.0) / math.sqrt(2.0) + math.sqrt(1.0 - par * par) * w
        return (xi, eta) if gen == "G0" else (sx * xi, sy * eta)
    if alt == "lag":  # Y_t = rho X_{t-1} + noise
        xi = mb.fgn(H, n + 1, rng)
        w = mb.fgn(H, n, rng)
        return xi[1:], par * xi[:-1] + math.sqrt(1.0 - par * par) * w
    if alt == "flat":  # Y_t = rho (X_{t-1} + X_t + X_{t+1}) / sd + noise
        xi = mb.fgn(H, n + 2, rng)
        w = mb.fgn(H, n, rng)
        ma = (xi[:-2] + xi[1:-1] + xi[2:]) / flat_ma_norm(H)
        return xi[1:-1], par * ma + math.sqrt(1.0 - par * par) * w
    if alt == "factor":  # X = m + e, Y = m + u, corr 0.3
        m = factor(n, rng)
        return m + mb.fgn(H, n, rng), m + mb.fgn(H, n, rng)
    if alt == "scale":  # Y = exp(gamma |X|) Z / sd
        xi = mb.fgn(H, n, rng)
        return xi, np.exp(par * np.abs(xi)) * mb.fgn(H, n, rng) / scale_norm(par)
    raise ValueError(alt)
