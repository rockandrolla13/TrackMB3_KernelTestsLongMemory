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

from kmem.legacy import gens
from kmem.legacy import mbcore as mb

ALPHA = 0.05
KAPPAS = (0.5, 1.0, 2.0, 4.0)


def q95(v: np.ndarray) -> float:
    return float(np.quantile(v, 1.0 - ALPHA))


def pboot(t: float, b: np.ndarray) -> float:
    return float((1.0 + np.sum(b >= t)) / (1.0 + b.size))


def ktag(k: float) -> str:
    return f"k{k:g}".replace(".", "p")


# ---------------------------------------------------------------- generators (EXP3 copy)

def seas_scale(n: int, rng: np.random.Generator) -> np.ndarray:
    u = int(rng.integers(gens.PERIOD))
    return 1.0 + gens.SEAS_B * gens.season(n, u)


def ind_pair(gen: str, alt: str, H: float, n: int, rho: float, rng: np.random.Generator):
    """EXP3 power_core.ind_pair, verbatim logic. alt in {null, lin, sq}; gen in {G0, SEAS-V}."""
    if alt == "null":
        return gens.pair(gen, H, n, rng)
    if gen == "SEAS-V":
        sx = seas_scale(n, rng)
    xi = mb.fgn(H, n, rng)
    if gen == "SEAS-V":
        sy = seas_scale(n, rng)
    w = mb.fgn(H, n, rng)
    if alt == "lin":
        eta = rho * xi + math.sqrt(1.0 - rho * rho) * w
    elif alt == "sq":
        eta = rho * (xi * xi - 1.0) / math.sqrt(2.0) + math.sqrt(1.0 - rho * rho) * w
    else:
        raise ValueError(alt)
    if gen == "G0":
        return xi, eta
    if gen == "SEAS-V":
        return sx * xi, sy * eta
    raise ValueError(gen)


# ---------------------------------------------------------------- bundle

def hsic_fields(x, y, dx, dy, sx, sy, w, n, ln, dn, A) -> dict:
    kc = mb.double_center(mb.gauss_gram(dx, sx))
    lc = mb.double_center(mb.gauss_gram(dy, sy))
    mh = kc * lc
    t_h = mh.sum() / n
    bh = mb.wild_boot_stats(mh, w)
    del mh
    s = mb.shift_stats(kc, lc)
    del kc, lc
    shifts = s[A:n - A + 1]
    sub_h = mb.subsample_stats(x, y, ln, dn, "hsic", sx, sy)
    q_hw, q_hs, q_hb = q95(bh), q95(shifts), q95(sub_h)
    return {
        "hsic_stat": t_h, "sx": sx, "sy": sy,
        "p_hsic_wb": pboot(t_h, bh), "r_hsic_wb": t_h / q_hw,
        "p_hsic_shift": pboot(s[0], shifts), "r_hsic_shift": s[0] / q_hs,
        "p_hsic_sub": float(np.mean(sub_h >= t_h)), "rej_hsic_sub": int(t_h > q_hb), "r_hsic_sub": t_h / q_hb,
        "shift_check": float(s[0] - t_h),
    }


def kind_bw(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, kappas=KAPPAS, B: int = 500,
            oracle_kappas=()) -> dict:
    """HSIC rules per kappa (suffix _k0p5 etc.), dCov and Pearson once, mechanism diagnostics.

    oracle_kappas: extra kappas for which only the bare statistic n HSIC_V is computed (cheap).

    rng is used only for the wild multipliers, drawn once, as in kind_ext.
    """
    n = x.size
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    A = math.floor(0.1 * n)

    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    sx0 = mb.median_bandwidth(dx)
    sy0 = mb.median_bandwidth(dy)

    ac = mb.double_center(dx)
    bc = mb.double_center(dy)
    md = ac * bc
    del ac, bc
    t_d = md.sum() / n
    w = mb.markov_rademacher(n, B, ln, rng)
    bd = mb.wild_boot_stats(md, w)
    del md
    sub_d = mb.subsample_stats(x, y, ln, dn, "dcov")
    q_dw, q_db = q95(bd), q95(sub_d)
    out = {
        "dcov_stat": t_d,
        "p_dcov_wb": pboot(t_d, bd), "r_dcov_wb": t_d / q_dw,
        "p_dcov_sub": float(np.mean(sub_d >= t_d)), "rej_dcov_sub": int(t_d > q_db), "r_dcov_sub": t_d / q_db,
        "sx_med": sx0, "sy_med": sy0,
        "corr_sq": float(np.corrcoef(x * x, y * y)[0, 1]),
        "corr_abs": float(np.corrcoef(np.abs(x), np.abs(y))[0, 1]),
    }
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
    lt = mb.lind_tests(x, y)
    nu = lt["ewc_nu"]
    out.update(lt)
    out["r_nw"] = abs(lt["t_nw"]) / 1.959963984540054 if np.isfinite(lt["t_nw"]) else float("nan")
    out["r_ewc"] = abs(lt["t_ewc"]) / mb._TQ[nu]
    return out


HSIC_KEYS = ["hsic_stat", "sx", "sy", "p_hsic_wb", "r_hsic_wb", "p_hsic_shift", "r_hsic_shift",
             "p_hsic_sub", "rej_hsic_sub", "r_hsic_sub", "shift_check"]
SHARED_KEYS = ["dcov_stat", "p_dcov_wb", "r_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "r_dcov_sub",
               "sx_med", "sy_med", "corr_sq", "corr_abs", "pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw",
               "r_nw", "ewc_nu", "t_ewc", "rej_ewc", "r_ewc", "nw_lrv_nonpos"]


def fields(kappas=KAPPAS, oracle_kappas=()) -> list[str]:
    return (SHARED_KEYS + [f"{key}_{ktag(k)}" for k in kappas for key in HSIC_KEYS]
            + [f"hsic_stat_{ktag(k)}" for k in oracle_kappas])
