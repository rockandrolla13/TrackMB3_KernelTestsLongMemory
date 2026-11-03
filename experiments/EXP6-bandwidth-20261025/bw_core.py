"""EXP6 bundle: EXP3 independence rules with the Gaussian bandwidth scaled by kappa.

sigma_X(kappa) = kappa * median heuristic, same for Y. The scaled bandwidths enter the full-sample
Gram matrices (wild, shift, oracle) and the subsampling windows. dCov and Pearson do not depend on
kappa and are computed once. All kappa share the data and the wild multipliers, so they are paired.
At kappa = 1 every field equals EXP3 power_core.kind_ext on the same seed (unit_checks.py).
Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.calibration.common import ALPHA, pboot, q95  # noqa: F401  (re-exported)
from kmem.generators.independence import pair
from kmem.generators.seasonal import seas_scale  # noqa: F401  (re-exported)
from kmem.hsic import dcov_full_rules, pearson_rules
from kmem.hsic import hsic_full_rules as hsic_fields  # noqa: F401  (EXP6 name, re-exported)
from kmem.legacy import mbcore as mb

KAPPAS = (0.5, 1.0, 2.0, 4.0)


def ktag(k: float) -> str:
    return f"k{k:g}".replace(".", "p")


# ---------------------------------------------------------------- generators (moved to kmem.generators)

def ind_pair(gen: str, alt: str, H: float, n: int, rho: float, rng: np.random.Generator):
    """EXP3 power_core.ind_pair: kmem.generators.independence.pair at equal memory. alt null, lin, sq; gen G0, SEAS-V."""
    return pair(gen, alt, H, H, n, rho, rng)


# ---------------------------------------------------------------- bundle


def kind_bw(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, kappas=KAPPAS, B: int = 500,
            oracle_kappas=()) -> dict:
    """HSIC rules per kappa (suffix _k0p5 etc.), dCov and Pearson once, mechanism diagnostics.

    oracle_kappas: extra kappas for which only the bare statistic n HSIC_V is computed (cheap).

    rng is used only for the wild multipliers, drawn once, as in kind_ext. The rule pieces live in kmem.hsic.
    """
    n = x.size
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    A = math.floor(0.1 * n)

    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    sx0 = mb.median_bandwidth(dx)
    sy0 = mb.median_bandwidth(dy)
    w = mb.markov_rademacher(n, B, ln, rng)
    out = dcov_full_rules(x, y, dx, dy, w, n, ln, dn)
    out.update({
        "sx_med": sx0, "sy_med": sy0,
        "corr_sq": float(np.corrcoef(x * x, y * y)[0, 1]),
        "corr_abs": float(np.corrcoef(np.abs(x), np.abs(y))[0, 1]),
    })
    for k in kappas:
        f = hsic_fields(x, y, dx, dy, k * sx0, k * sy0, w, n, ln, dn, A)
        tag = ktag(k)
        out.update({f"{key}_{tag}": v for key, v in f.items()})
    for k in oracle_kappas:
        kc = mb.double_center(mb.gauss_gram(dx, k * sx0))
        lc = mb.gauss_gram(dy, k * sy0)
        out[f"hsic_stat_{ktag(k)}"] = float(np.sum(kc * lc)) / n  # sum(Kc Lc) = sum(Kc L)
        del kc, lc
    del dx, dy, w
    out.update(pearson_rules(x, y))
    return out


HSIC_KEYS = ["hsic_stat", "sx", "sy", "p_hsic_wb", "r_hsic_wb", "p_hsic_shift", "r_hsic_shift",
             "p_hsic_sub", "rej_hsic_sub", "r_hsic_sub", "shift_check"]
SHARED_KEYS = ["dcov_stat", "p_dcov_wb", "r_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "r_dcov_sub",
               "sx_med", "sy_med", "corr_sq", "corr_abs", "pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw",
               "r_nw", "ewc_nu", "t_ewc", "rej_ewc", "r_ewc", "nw_lrv_nonpos"]


def fields(kappas=KAPPAS, oracle_kappas=()) -> list[str]:
    return (SHARED_KEYS + [f"{key}_{ktag(k)}" for k in kappas for key in HSIC_KEYS]
            + [f"hsic_stat_{ktag(k)}" for k in oracle_kappas])
