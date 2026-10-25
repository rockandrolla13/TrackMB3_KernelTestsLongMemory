"""EXP3 rule bundles that return, for every rule, its raw decision input AND a continuous ratio.

The ratio of a rule is R = T / c_ref, where T is the rule's statistic and c_ref is the rule's own
95% reference quantile (wild draws, circular shifts, or sampling windows). The raw rule rejects
iff p <= 0.05 (wild, shift) or T > q95 (subsampling), exactly as in EXP2. Size adjustment
replaces the cut-off 1 on R by the 95% quantile of R over null replications (DESIGN.md).

The statistics and references are the EXP2 code (mbcore.kind_tests, mmdcal.rep_e2), with the
extra quantiles read off the same reference draws. unit_checks.py verifies the p-values match.
Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import math

import numpy as np

import gens
import mbcore as mb
import mmdcal as M

ALPHA = 0.05
SEAS_PERIOD = gens.PERIOD


def q95(v: np.ndarray) -> float:
    return float(np.quantile(v, 1.0 - ALPHA))


def pboot(t: float, b: np.ndarray) -> float:
    return float((1.0 + np.sum(b >= t)) / (1.0 + b.size))


# ---------------------------------------------------------------- generators for EXP3

def seas_scale(n: int, rng: np.random.Generator) -> np.ndarray:
    u = int(rng.integers(SEAS_PERIOD))
    return 1.0 + gens.SEAS_B * gens.season(n, u)


def ind_pair(gen: str, alt: str, H: float, n: int, rho: float, rng: np.random.Generator):
    """Independence-arm pair. alt in {null, lin, sq}; gen in {G0, SEAS-V}.

    null: X, Y independent (gens.pair, the EXP2 code path).
    lin:  latent eta = rho xi + sqrt(1 - rho^2) w, w independent fGn(H).
    sq:   eta = rho (xi^2 - 1)/sqrt(2) + sqrt(1 - rho^2) w  (uncorrelated, rank-two dependence).
    SEAS-V multiplies each latent series by its own seasonal scale with an independent phase.
    """
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


def mmd_series(gen: str, alt: str, H: float, n: int, delta: float, rng: np.random.Generator):
    """Split-MMD series. alt 'null' is the EXP2 generator; 'shift' adds delta to the second half."""
    x = gens.generate(gen, H, n, rng)
    if alt == "shift":
        x = x.copy()
        x[n // 2:] += delta
    elif alt != "null":
        raise ValueError(alt)
    return x


# ---------------------------------------------------------------- independence bundle

def kind_ext(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, B: int = 500) -> dict:
    """mbcore.kind_tests plus 95% reference quantiles. Same rng consumption order."""
    n = x.size
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    A = math.floor(0.1 * n)

    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    sx = mb.median_bandwidth(dx)
    sy = mb.median_bandwidth(dy)
    kc = mb.double_center(mb.gauss_gram(dx, sx))
    lc = mb.double_center(mb.gauss_gram(dy, sy))
    mh = kc * lc
    t_h = mh.sum() / n
    ac = mb.double_center(dx)
    bc = mb.double_center(dy)
    del dx, dy
    md = ac * bc
    del ac, bc
    t_d = md.sum() / n

    w = mb.markov_rademacher(n, B, ln, rng)
    bh = mb.wild_boot_stats(mh, w)
    del mh
    bd = mb.wild_boot_stats(md, w)
    del md, w
    s = mb.shift_stats(kc, lc)
    del kc, lc
    shifts = s[A:n - A + 1]
    sub_h = mb.subsample_stats(x, y, ln, dn, "hsic", sx, sy)
    sub_d = mb.subsample_stats(x, y, ln, dn, "dcov")

    q_hw, q_hs, q_hb = q95(bh), q95(shifts), q95(sub_h)
    q_dw, q_db = q95(bd), q95(sub_d)
    out = {
        "hsic_stat": t_h, "dcov_stat": t_d,
        "p_hsic_wb": pboot(t_h, bh), "r_hsic_wb": t_h / q_hw,
        "p_hsic_shift": pboot(s[0], shifts), "r_hsic_shift": s[0] / q_hs,
        "p_hsic_sub": float(np.mean(sub_h >= t_h)), "rej_hsic_sub": int(t_h > q_hb), "r_hsic_sub": t_h / q_hb,
        "p_dcov_wb": pboot(t_d, bd), "r_dcov_wb": t_d / q_dw,
        "p_dcov_sub": float(np.mean(sub_d >= t_d)), "rej_dcov_sub": int(t_d > q_db), "r_dcov_sub": t_d / q_db,
        "shift_check": float(s[0] - t_h),
    }
    lt = mb.lind_tests(x, y)
    nu = lt["ewc_nu"]
    out.update(lt)
    out["r_nw"] = abs(lt["t_nw"]) / 1.959963984540054 if np.isfinite(lt["t_nw"]) else float("nan")
    out["r_ewc"] = abs(lt["t_ewc"]) / mb._TQ[nu]
    return out


# ---------------------------------------------------------------- split-MMD bundle

def mmd_ext(x: np.ndarray, rng: np.random.Generator, B: int = M.B_WILD) -> dict:
    """Uncentred wild (l = floor(n^0.5)), sqrt subsampling, self-normalised block (l = floor(n^0.5)).

    Same formulas as mmdcal.rep_e2 for these three rules (rep_e2 also draws other multipliers
    first, so per-seed values differ; the null is re-run in this code path).
    """
    n = x.size
    h1 = n // 2
    d = np.abs(x[:, None] - x[None, :])
    sigma = M.median_heuristic(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    s = np.empty(n)
    s[:h1] = 1.0 / h1
    s[h1:] = -1.0 / (n - h1)
    stat = float(s @ k @ s)
    t_n = n * stat
    out: dict = {"mmd2": stat, "sigma": sigma, "t_n": t_n}

    kc = k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()
    l = max(1, math.floor(n ** 0.5))
    w = M.multipliers("unc", n, h1, l, rng, B)
    a = w * s[:, None]
    boot = np.einsum("ij,ij->j", a, kc @ a)
    del kc, a, w
    out["p_unc_05"] = pboot(stat, boot)
    out["r_unc_05"] = stat / q95(boot)

    p = M.prefix2(k)
    del k
    v = np.sort(M.window_mmd(p, n, l))
    q_sub = float(v[int(math.ceil((1 - ALPHA) * v.size)) - 1])  # rep_e2 convention
    out["l_sqrt"] = l
    out["q95_sqrt"] = q_sub
    out["r_sub_sqrt"] = t_n / q_sub
    out["sn_T"] = float(M.sn_stat(p, np.array([0]), n)[0])
    v = np.sort(M.sn_stat(p, np.arange(n - l + 1), l))
    out["sn_q95_05"] = float(v[int(math.ceil((1 - ALPHA) * v.size)) - 1])
    out["r_sn_05"] = out["sn_T"] / out["sn_q95_05"]
    return out
