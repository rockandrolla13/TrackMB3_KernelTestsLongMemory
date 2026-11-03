"""A new kernel is one class plus one register() call; the EXP7 bundles run it without edits (plan step 3.2).

The toy kernel negates the Gaussian features: (-F)(-F)' = F F', so every statistic must equal raw's up to rounding.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np
import pytest

from kmem import kernels
from kmem import memory_kernels as mk
from kmem.legacy import mbcore as mb

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "experiments", "EXP7-memory-kernel-20261101"))
import exp7_core as K  # noqa: E402  (the EXP7 bundle under test)


@dataclass(frozen=True)
class Negated:
    name: str = "toy_neg"
    kind: str = "feature"
    lag: int = 0

    def features(self, f, x, H=None):
        return -f

    def gram(self, x, sigma, H=None):
        return mk.double_center(mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sigma))


kernels.register(Negated())


def _xy(n: int = 200):
    rng = np.random.default_rng(20261103)
    return mb.fgn(0.8, n, rng), mb.fgn(0.8, n, rng)


def test_registry_holds_the_ten_exp7_kernels() -> None:
    assert set(K.IND_KERNELS) <= set(kernels.names())
    assert kernels.get("herm2").r == 2 and kernels.get("inc50").lag == 50


def test_unknown_kernel_is_a_clear_error() -> None:
    with pytest.raises(ValueError, match="unknown kernel"):
        kernels.get("herm10")


def test_toy_kernel_runs_in_the_independence_bundle() -> None:
    x, y = _xy()
    out = K.ind_bundle(x, y, np.random.default_rng(1), kernels=("raw", "toy_neg"), full=False, hsplit=False)
    cols = [c for c in K.kernel_columns(("toy_neg",))]
    assert all(c in out for c in cols), [c for c in cols if c not in out]
    for key in ("t", "r_wb", "r_shift"):
        assert out[f"{key}_toy_neg"] == pytest.approx(out[f"{key}_raw"], rel=1e-9)


def test_toy_kernel_runs_in_the_split_mmd_bundle_and_the_exact_path() -> None:
    x, _ = _xy()
    out = K.mmd_bundle(x, np.random.default_rng(2), kernels=("raw", "toy_neg"), hsplit=False)
    assert out["mmd2_toy_neg"] == pytest.approx(out["mmd2_raw"], rel=1e-9)
    sigma = K.median_bw_fast(x)
    assert np.array_equal(K.gram_full("toy_neg", x, sigma), K.gram_full("raw", x, sigma))
