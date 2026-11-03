"""Generator registries (plan step 3.4, findings AR-EXT-001 and CR-004).

Four small tables, each name -> function:

- processes:     (H, n, rng) -> one series at memory H. The EXP2 generators (wrapped, legacy untouched),
                 and SEAS-V-slow from EXP4.
- pair nulls:    (HX, HY, n, par, rng) -> (X, Y) independent but not two copies of one process (FACTOR, SCALE).
- alternatives:  (gen, H, n, par, rng) -> dependent (X, Y): lin, sq, lag, flat, factor, scale.
- changes:       (second_half, par) -> new second half, for split-MMD series: mean, scale, shape.

A null pair of a process is two independent draws, X first. To add a generator, write its function in a new
module and call the matching register_* once. Each built-in keeps the draw order of the code it replaced.
"""
from __future__ import annotations

import math
from typing import Callable

import numpy as np

from kmem.generators.seasonal import seas_scale, seas_v_slow
from kmem.generators.stress import factor, flat_ma_norm, scale_norm
from kmem.legacy import gens
from kmem.legacy import mbcore as mb

Process = Callable[[float, int, np.random.Generator], np.ndarray]
PairNull = Callable[[float, float, int, float, np.random.Generator], tuple]
Alternative = Callable[[str, float, int, float, np.random.Generator], tuple]
Change = Callable[[np.ndarray, float], np.ndarray]

_TABLES: dict[str, dict[str, Callable]] = {"process": {}, "pair null": {}, "alternative": {}, "change": {}}


def _register(kind: str, name: str, fn: Callable) -> Callable:
    table = _TABLES[kind]
    if name in table and table[name] is not fn:
        raise ValueError(f"{kind} {name!r} is already registered")
    table[name] = fn
    return fn


def _get(kind: str, name: str) -> Callable:
    try:
        return _TABLES[kind][name]
    except KeyError:
        raise ValueError(f"unknown {kind} {name!r}; registered: {sorted(_TABLES[kind])}") from None


def register_process(name: str, fn: Process) -> Process:
    return _register("process", name, fn)


def register_pair_null(name: str, fn: PairNull) -> PairNull:
    return _register("pair null", name, fn)


def register_alternative(name: str, fn: Alternative) -> Alternative:
    return _register("alternative", name, fn)


def register_change(name: str, fn: Change) -> Change:
    return _register("change", name, fn)


def get_process(name: str) -> Process:
    return _get("process", name)


def get_pair_null(name: str) -> PairNull | None:
    return _TABLES["pair null"].get(name)


def get_alternative(name: str) -> Alternative:
    return _get("alternative", name)


def get_change(name: str) -> Change:
    return _get("change", name)


def names(kind: str) -> tuple[str, ...]:
    return tuple(_TABLES[kind])


# ---------------------------------------------------------------- built-in processes

def _legacy(gen: str) -> Process:
    def draw(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
        return gens.generate(gen, H, n, rng)
    draw.__name__ = f"legacy_{gen}"
    return draw


for _g in gens.MAIN + ("SEAS-M-shared",):
    register_process(_g, _legacy(_g))
register_process("SEAS-V-slow", seas_v_slow)


# ---------------------------------------------------------------- built-in pair nulls (EXP7 matched nulls)

def factor_null(HX: float, HY: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    """Matched null of the common slow factor: two independent factors."""
    x = factor(n, rng) + mb.fgn(HX, n, rng)
    return x, factor(n, rng) + mb.fgn(HY, n, rng)


def scale_null(HX: float, HY: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    """Matched null of the scale coupling: Y's scale comes from an independent copy."""
    x = mb.fgn(HX, n, rng)
    xp = mb.fgn(HX, n, rng)
    return x, np.exp(par * np.abs(xp)) * mb.fgn(HY, n, rng) / scale_norm(par)


register_pair_null("FACTOR", factor_null)
register_pair_null("SCALE", scale_null)


# ---------------------------------------------------------------- built-in alternatives

def _latent(gen: str, H: float, n: int, par: float, rng: np.random.Generator, eta: Callable) -> tuple:
    """EXP3 power_core.ind_pair draw order (seasonal scale of X, xi, scale of Y, w); for G0 it is also EXP7's."""
    if gen not in ("G0", "SEAS-V"):
        raise ValueError(f"lin and sq are defined for G0 and SEAS-V, got {gen}")
    if gen == "SEAS-V":
        sx = seas_scale(n, rng)
    xi = mb.fgn(H, n, rng)
    if gen == "SEAS-V":
        sy = seas_scale(n, rng)
    w = mb.fgn(H, n, rng)
    e = eta(xi, w, par)
    return (xi, e) if gen == "G0" else (sx * xi, sy * e)


def lin(gen: str, H: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    return _latent(gen, H, n, par, rng, lambda xi, w, p: p * xi + math.sqrt(1.0 - p * p) * w)


def sq(gen: str, H: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    return _latent(gen, H, n, par, rng, lambda xi, w, p: p * (xi * xi - 1.0) / math.sqrt(2.0) + math.sqrt(1.0 - p * p) * w)


def _g0_only(gen: str, alt: str) -> None:
    if gen != "G0":
        raise ValueError(f"{alt} is defined for G0, got {gen}")


def lag(gen: str, H: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    """Y_t = rho X_{t-1} + noise."""
    _g0_only(gen, "lag")
    xi = mb.fgn(H, n + 1, rng)
    w = mb.fgn(H, n, rng)
    return xi[1:], par * xi[:-1] + math.sqrt(1.0 - par * par) * w


def flat(gen: str, H: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    """Y_t = rho (X_{t-1} + X_t + X_{t+1}) / sd + noise."""
    _g0_only(gen, "flat")
    xi = mb.fgn(H, n + 2, rng)
    w = mb.fgn(H, n, rng)
    ma = (xi[:-2] + xi[1:-1] + xi[2:]) / flat_ma_norm(H)
    return xi[1:-1], par * ma + math.sqrt(1.0 - par * par) * w


def common_factor(gen: str, H: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    """X = m + e, Y = m + u, corr 0.3."""
    _g0_only(gen, "factor")
    m = factor(n, rng)
    return m + mb.fgn(H, n, rng), m + mb.fgn(H, n, rng)


def scale_coupling(gen: str, H: float, n: int, par: float, rng: np.random.Generator) -> tuple:
    """Y = exp(gamma |X|) Z / sd."""
    _g0_only(gen, "scale")
    xi = mb.fgn(H, n, rng)
    return xi, np.exp(par * np.abs(xi)) * mb.fgn(H, n, rng) / scale_norm(par)


for _name, _fn in (("lin", lin), ("sq", sq), ("lag", lag), ("flat", flat), ("factor", common_factor),
                   ("scale", scale_coupling)):
    register_alternative(_name, _fn)


# ---------------------------------------------------------------- built-in split-MMD changes

def mean_change(t: np.ndarray, par: float) -> np.ndarray:
    return t + par


def scale_change(t: np.ndarray, par: float) -> np.ndarray:
    return t * par


def shape_change(t: np.ndarray, par: float) -> np.ndarray:
    """Same mean and variance, different shape: sqrt(1 - a^2) xi + a (xi^2 - 1) / sqrt(2)."""
    return math.sqrt(1.0 - par * par) * t + par * (t * t - 1.0) / math.sqrt(2.0)


for _name, _fn in (("mean", mean_change), ("scale", scale_change), ("shape", shape_change)):
    register_change(_name, _fn)
