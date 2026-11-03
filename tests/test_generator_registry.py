"""A new generator is one function plus one register_* call (plan step 3.4)."""
from __future__ import annotations

import numpy as np
import pytest

from kmem.generators import registry
from kmem.generators.independence import pair
from kmem.generators.split import split_series


def toy_white(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Short-memory toy process: iid standard normal, H ignored."""
    return rng.standard_normal(n)


registry.register_process("toy_white", toy_white)


def test_toy_process_gives_null_pairs_and_split_series() -> None:
    x, y = pair("toy_white", "null", 0.5, 0.5, 64, 0.0, np.random.default_rng(3))
    r = np.random.default_rng(3)
    assert np.array_equal(x, r.standard_normal(64)) and np.array_equal(y, r.standard_normal(64))
    s = split_series("toy_white", "mean", 0.5, 64, 0.25, np.random.default_rng(4))
    base = np.random.default_rng(4).standard_normal(64)
    assert np.array_equal(s[:32], base[:32]) and np.array_equal(s[32:], base[32:] + 0.25)


def test_every_exp7_name_is_registered() -> None:
    for g in ("G0", "HT-S", "HT-L", "NL-SV", "NL-TAR", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow"):
        assert registry.get_process(g)
    for a in ("lin", "sq", "lag", "flat", "factor", "scale"):
        assert registry.get_alternative(a)
    assert registry.get_pair_null("FACTOR") and registry.get_pair_null("SCALE")


def test_unknown_names_are_clear_errors() -> None:
    with pytest.raises(ValueError, match="unknown process"):
        pair("G9", "null", 0.5, 0.5, 64, 0.0, np.random.default_rng(1))
    with pytest.raises(ValueError, match="unknown alternative"):
        pair("G0", "cubic", 0.5, 0.5, 64, 0.1, np.random.default_rng(1))
    with pytest.raises(ValueError, match="already registered"):
        registry.register_process("G0", toy_white)
