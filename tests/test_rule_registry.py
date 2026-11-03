"""A new calibration rule is one object plus one register() call; the runner writes its columns and the
analysis rejects with the same object (plan step 3.3)."""
from __future__ import annotations

import numpy as np
import pytest

from kmem import features, rules
from kmem.calibration import registry as cal
from kmem.calibration.common import pboot
from kmem.calibration.wild import Multipliers
from kmem.legacy import mbcore as mb

# Toy rule: the wild bootstrap at the 10% level, with its own column names.
WILD10 = cal.register(cal.CalibrationRule(
    "wild10", "wild", lambda t, b: {"p_wb10": pboot(t, b), "r_wb10": t / float(np.quantile(b, 0.90))},
    "p_wb10", lambda v: v <= 0.10, "r_wb10"))


def _pair(n: int = 300):
    rng = np.random.default_rng(20261103)
    x, y = mb.fgn(0.7, n, rng), mb.fgn(0.7, n, rng)
    fx = features.compress(features.gauss_features(x, mb.median_bandwidth(np.abs(x[:, None] - x[None, :]))))
    fy = features.compress(features.gauss_features(y, mb.median_bandwidth(np.abs(y[:, None] - y[None, :]))))
    return fx, fy, Multipliers(n, 200, rng).w


def test_default_rules_give_the_documented_columns() -> None:
    fx, fy, w = _pair()
    out = rules.rules_fast(fx, fy, w, sub=True)
    assert set(out) == {"t", "m", "q", "check", "p_wb", "r_wb", "p_shift", "r_shift", "r975_wb", "r975_shift",
                        "p_sub", "rej_sub", "r_sub"}


def test_toy_rule_columns_and_rejection() -> None:
    fx, fy, w = _pair()
    base = rules.rules_fast(fx, fy, w)
    out = rules.rules_fast(fx, fy, w, rules=cal.FAST_RULES + ("wild10",))
    assert out["p_wb10"] == out["p_wb"]  # same draws, same statistic, same p-value
    assert all(out[k] == base[k] for k in base)  # adding a rule changes no other column
    assert cal.get("wild10").reject(out, suffix="") == (out["p_wb"] <= 0.10)
    col = {"p_wb10_raw": np.array([0.05, 0.2])}
    assert list(WILD10.reject(col, suffix="_raw")) == [True, False]


def test_unknown_rule_is_a_clear_error() -> None:
    with pytest.raises(ValueError, match="unknown calibration rule"):
        cal.get("bogus")
