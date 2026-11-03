"""Calibration rules by name, shared by the runner and the analysis (plan step 3.3, AR-EXT-002, AR-DRY-003, CR-007).

A rule has two halves:
- decide(stat, draws) -> {column stem: value}: p-value, ratio or rejection from the statistic and the null draws
  it names ("wild" multiplier draws, circular "shift" statistics, or "subsample" window statistics).
  Rules decided elsewhere (Harris, Pearson, self-normalised) have draws = None and only the second half.
- reject(data, column): the rejection read back from a saved column (a value or a numpy column).

The runner writes what decide returns; the analysis calls reject on the same object, so the two cannot disagree.
To add a rule: build a CalibrationRule in a new module and call register(rule).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from kmem.calibration.common import ALPHA, pboot, q95


@dataclass(frozen=True)
class CalibrationRule:
    name: str
    draws: str | None  # "wild", "shift", "subsample", or None
    decide: Callable[[float, np.ndarray], dict] | None
    decision: str  # stem of the column reject() reads, e.g. "p_wb"
    reject_value: Callable[[Any], Any]  # decision column value(s) -> rejection(s)
    ratio: str | None = None  # stem of the continuous ratio column (size adjustment), e.g. "r_wb"

    def reject(self, data, column: str | None = None, suffix: str = ""):
        """Rejection(s) from data[column]; column defaults to the decision stem plus suffix (e.g. "_raw")."""
        return self.reject_value(data[column if column is not None else self.decision + suffix])


_REGISTRY: dict[str, CalibrationRule] = {}


def register(rule: CalibrationRule) -> CalibrationRule:
    old = _REGISTRY.get(rule.name)
    if old is not None and old != rule:
        raise ValueError(f"calibration rule {rule.name!r} is already registered")
    if rule.draws not in (None, "wild", "shift", "subsample"):
        raise ValueError(f"rule {rule.name!r}: unknown draws {rule.draws!r}")
    _REGISTRY[rule.name] = rule
    return rule


def get(name: str) -> CalibrationRule:
    try:
        return _REGISTRY[name]
    except KeyError:
        raise ValueError(f"unknown calibration rule {name!r}; registered: {sorted(_REGISTRY)}") from None


def names() -> tuple[str, ...]:
    return tuple(_REGISTRY)


# ---------------------------------------------------------------- built-in rules

def _p_at_level(v):
    return v <= ALPHA


def _above_one(v):
    return v > 1


def _flag_set(v):
    return v == 1


def _subsample(t: float, sb: np.ndarray) -> dict:
    q_b = q95(sb)
    return {"p_sub": float(np.mean(sb >= t)), "rej_sub": int(t > q_b), "r_sub": t / q_b}


WILD = register(CalibrationRule("wild", "wild", lambda t, b: {"p_wb": pboot(t, b), "r_wb": t / q95(b)},
                                "p_wb", _p_at_level, "r_wb"))
SHIFT = register(CalibrationRule("shift", "shift", lambda s0, sh: {"p_shift": pboot(s0, sh), "r_shift": s0 / q95(sh)},
                                 "p_shift", _p_at_level, "r_shift"))
SUBSAMPLE = register(CalibrationRule("subsample", "subsample", _subsample, "rej_sub", _flag_set, "r_sub"))
WILD975 = register(CalibrationRule("wild975", "wild", lambda t, b: {"r975_wb": t / float(np.quantile(b, 0.975))},
                                   "r975_wb", _above_one, "r975_wb"))
SHIFT975 = register(CalibrationRule("shift975", "shift",
                                    lambda s0, sh: {"r975_shift": s0 / float(np.quantile(sh, 0.975))},
                                    "r975_shift", _above_one, "r975_shift"))
HARRIS = register(CalibrationRule("harris", None, None, "rhc", _above_one, "rhc"))
HARRIS19 = register(CalibrationRule("harris19", None, None, "rh19", _above_one, "rh19"))
HARRIS_A = register(CalibrationRule("harrisA", None, None, "rha", _above_one, "rha"))
NEWEY_WEST = register(CalibrationRule("nw", None, None, "rej_nw", _flag_set, "r_nw"))
EWC = register(CalibrationRule("ewc", None, None, "rej_ewc", _flag_set, "r_ewc"))
SN_LOW = register(CalibrationRule("low", None, None, "rl", _above_one, "rl"))
SN_LOW_BLOCK = register(CalibrationRule("low_blk", None, None, "rb", _above_one, "rb"))

FAST_RULES = ("wild", "shift", "wild975", "shift975")  # what rules_fast reports by default
FULL_RULES = ("wild", "shift")  # what rules_full and rules_distance report by default
