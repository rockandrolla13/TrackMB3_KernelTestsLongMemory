"""Independence-arm pairs (X, Y): nulls and alternatives of EXP2, EXP3, EXP6 and EXP7, in their original draw order."""
from __future__ import annotations

import numpy as np

from kmem.generators.registry import get_alternative, get_pair_null, get_process


def pair(gen: str, alt: str, HX: float, HY: float, n: int, par: float, rng: np.random.Generator):
    """Independence-arm pair, looked up in kmem.generators.registry.

    alt "null": a registered pair null (FACTOR, SCALE), else two independent draws of process gen at HX then HY.
    Any other alt is a registered alternative at equal memory: lin, sq (gen G0 or SEAS-V, as in EXP3 and EXP6),
    lag, flat, factor, scale (gen G0, EXP7).
    """
    if alt == "null":
        null = get_pair_null(gen)
        if null is not None:
            return null(HX, HY, n, par, rng)
        proc = get_process(gen)
        x = proc(HX, n, rng)
        return x, proc(HY, n, rng)
    if HX != HY:
        raise ValueError("alternatives use equal H")
    return get_alternative(alt)(gen, HX, n, par, rng)
