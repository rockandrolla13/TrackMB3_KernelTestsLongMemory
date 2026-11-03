"""EXP7 core: generators, memory-kernel rule bundles, and a low-rank fast path (see DESIGN.md).

Three paths compute the same rules. unit_checks.py compares them on every kernel.

Exact path: mbcore functions on n x n Gram matrices (raw_exact, mmd_raw_exact). It reproduces EXP6
hsic_fields and EXP3 mmd_ext bit for bit on the same seed. It runs on the first replications of each
cell (columns xraw_*), so the fast path is also compared with it on the real cells.

Fast path: a Gaussian Gram matrix on one series has fast-decaying eigenvalues. A pivoted Cholesky
factor K = F F' gets the trace error below 1e-12 relative. Each memory kernel is a linear map A on
features, so A K A' = (A F)(A F)'. All three rules then reduce to sums over a few products of feature
columns.

Distance path: the distance kernel has slowly decaying eigenvalues, so dCov and dcov_inc1 keep their
n x n matrices. Two exact shortcuts remove most of their cost. The block-Markov multipliers are
piecewise constant, so w' M w needs only 2-d tail sums of M at the sign flips. Symmetric Gram matrices
repeat each circulant diagonal, so the shift statistic needs only half of them.

The general pieces (bandwidth, features, multipliers, rules, Harris, self-normalised tests, Hermite layers)
moved to the kmem package on 2026-11-03 with their code unchanged. They are re-exported here, so EXP7
scripts and checks keep their K.<name> references.

Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import math
import os

import numpy as np

from kmem import memory_kernels as mk
from kmem.bandwidth import median_bw_fast  # noqa: F401  (re-exported)
from kmem.calibration.common import ALPHA, pboot, q95  # noqa: F401
from kmem.calibration.harris import (HARRIS_N_MIN, harris_big_n, harris_ratio, harris_rules,  # noqa: F401
                                     harris_shifts, harris_strict19)
from kmem.calibration.selfnorm import (SN_GRID, SN_LEVELS, SNTable, betken_sup, bridge_sum,  # noqa: F401
                                       bridge_sums_ranges, bridge_sums_windows, sn_mean_stat, sn_mean_windows,
                                       sn_split_stat, sn_split_windows)
from kmem.calibration.shift import shift_stats_sym  # noqa: F401
from kmem.calibration.subsample import sub_fast  # noqa: F401
from kmem.calibration.wild import INC_LONG, Multipliers, wild_prefix  # noqa: F401
from kmem.generators.independence import pair as ind_pair  # noqa: F401
from kmem.generators.seasonal import seas_scale, seas_v_slow  # noqa: F401
from kmem.generators.split import split_series as mmd_series  # noqa: F401
from kmem.generators.stress import FACTOR_PHI, FACTOR_VAR, factor, flat_ma_norm, scale_norm  # noqa: F401
from kmem.hsic import dcov_full_rules, pearson_rules
from kmem.mmd import contrast, split_mmd_wild  # noqa: F401  (contrast re-exported)
from kmem.features import TOL_CHOL, TOL_PAIR, compress, gauss_features, select_pairs, window_sum  # noqa: F401
from kmem.layers import hermite_layers, layer_memory  # noqa: F401
from kmem.legacy import mbcore as mb
from kmem.legacy import mmdcal as M
from kmem.rules import rules_distance, rules_fast, rules_full  # noqa: F401

B_WILD = 500
HP_FRAC = 0.1

IND_KERNELS = ("raw", "inc1", "inc50", "herm1", "herm2", "sym", "datadiff", "hp10", "dcov_inc1", "frac_oracle")
SUB_KERNELS = ("raw", "inc1")
MMD_KERNELS = ("raw", "herm1", "herm2", "sym", "datadiff")

# hsplit: herm_r on the remainder plus a self-normalised test of the removed layers, Bonferroni (DESIGN.md A2)
HSPLIT_R = (1, 2)
IND_PAIRS = {1: ((1, 1),), 2: ((1, 1), (1, 2), (2, 1), (2, 2))}
MMD_LAYERS = {1: (1,), 2: (1, 2)}
SN_TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sn_table.json")
_SN = SNTable(SN_TABLE_PATH)

# Harris (2020) non-circular shift test (DESIGN.md A3)
HARRIS_KERNELS = ("raw", "inc1")
HARRIS_STRICT_REPS = {500: 40, 2000: 40, 4000: 10}


# ---------------------------------------------------------------- kernel feature maps (fast path)

def transform_features(kernel: str, f: np.ndarray, x: np.ndarray, H: float | None) -> np.ndarray:
    """Apply the kernel's linear feature map to the Gaussian features f of raw x."""
    if kernel == "raw":
        return f
    if kernel == "inc1":
        return f[1:] - f[:-1]
    if kernel == "inc50":
        return f[INC_LONG:] - f[:-INC_LONG]
    if kernel in ("herm1", "herm2"):
        u = mk.hermite_basis(x, int(kernel[-1]))
        return f - u @ (u.T @ f)
    if kernel == "hp10":
        u = mk.dct_basis(x.size, mk.n_low(x.size, HP_FRAC))
        return f - u @ (u.T @ f)
    if kernel == "frac_oracle":
        if H is None:
            raise ValueError("frac_oracle needs the true H")
        d = H - 0.5
        return f - f.mean(axis=0, keepdims=True) if d == 0.0 else mk.fracdiff_features(f - f.mean(axis=0, keepdims=True), d)
    raise ValueError(kernel)


def own_series(kernel: str, x: np.ndarray) -> np.ndarray:
    """Series on which sym and datadiff build their own Gaussian kernel."""
    return mk.even_map(x) if kernel == "sym" else mk.first_difference(x)


# ---------------------------------------------------------------- full Gram matrices of each kernel (checks)

def gram_full(kernel: str, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
    """Full-path transformed Gram matrix for one series. Used by unit checks and tests."""
    if kernel in ("sym", "datadiff"):
        u = own_series(kernel, x)
        d = np.abs(u[:, None] - u[None, :])
        return mk.double_center(mb.gauss_gram(d, mb.median_bandwidth(d)))
    if kernel == "dcov_inc1":
        return mk.increment(mk.distance_kernel(x), 1)
    k = mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sigma)
    if kernel == "raw":
        return mk.double_center(k)
    if kernel == "inc1":
        return mk.increment(k, 1)
    if kernel == "inc50":
        return mk.increment(k, INC_LONG)
    if kernel in ("herm1", "herm2"):
        return mk.hermite(k, x, int(kernel[-1]))
    if kernel == "hp10":
        return mk.highpass(k, HP_FRAC)
    if kernel == "frac_oracle":
        return mk.fracdiff(k, H - 0.5)
    raise ValueError(kernel)


# ---------------------------------------------------------------- hsplit: tests of the removed Hermite layers

def sn_table() -> dict:
    return _SN.data


def sn_cut(kind: str, n: int, level: float, h: float) -> float:
    """Critical value at memory h: interpolated in h on the log scale from the EXP7 simulated table."""
    return _SN.cut(kind, n, level, h)


def hsplit_ind(x: np.ndarray, y: np.ndarray) -> dict:
    """SN tests of the products He_a(z_x) He_b(z_y), a, b <= 2. rl_r and rb_r exceed 1 iff the low-layer part rejects.

    hsplit_r spends alpha/2 on the low layers, split equally over r^2 products. rl uses the plug-in table, keyed on
    H_u = max(1/2, H_a + H_b - 1) with H_a, H_b the local Whittle memories of the two layers. rb is the H-free
    variant: the cut-off is a quantile of the same statistic on windows of length floor(n^0.5).
    """
    n = x.size
    l = math.floor(n ** 0.5)
    lx, ly = hermite_layers(x), hermite_layers(y)
    hx = [layer_memory(lx[:, 0]), layer_memory(lx[:, 1])]
    hy = [layer_memory(ly[:, 0]), layer_memory(ly[:, 1])]
    out = {"hl_x1": hx[0], "hl_x2": hx[1], "hl_y1": hy[0], "hl_y2": hy[1]}
    rl2, rb2 = 0.0, 0.0
    for a, b in IND_PAIRS[2]:
        u = lx[:, a - 1] * ly[:, b - 1]
        t = sn_mean_stat(u)
        win = sn_mean_windows(u, l)
        hu = max(0.5, hx[a - 1] + hy[b - 1] - 1.0)
        out[f"T_{a}{b}"] = t
        rl2 = max(rl2, t / sn_cut("prod", n, 0.99375, hu))
        rb2 = max(rb2, t / float(np.quantile(win, 0.99375)))
        if (a, b) == (1, 1):
            out["rl1"] = t / sn_cut("prod", n, 0.975, hu)
            out["rb1"] = t / float(np.quantile(win, 0.975))
    out["rl2"], out["rb2"] = rl2, rb2
    return out


def hsplit_mmd(x: np.ndarray) -> dict:
    """SN tests of a change at n // 2 in He_1(z) and He_2(z). hsplit_r spends alpha/2 on the low layers, split over r."""
    n = x.size
    l = math.floor(n ** 0.5)
    v = hermite_layers(x)
    h1, h2 = layer_memory(v[:, 0]), layer_memory(v[:, 1])
    t1, t2 = sn_split_stat(v[:, 0]), sn_split_stat(v[:, 1])
    w1, w2 = sn_split_windows(v[:, 0], l), sn_split_windows(v[:, 1], l)
    return {"hl_1": h1, "hl_2": h2, "T_1": t1, "T_2": t2,
            "rl1": t1 / sn_cut("split1", n, 0.975, h1),
            "rl2": max(t1 / sn_cut("split1", n, 0.9875, h1), t2 / sn_cut("split2", n, 0.9875, h2)),
            "rb1": t1 / float(np.quantile(w1, 0.975)),
            "rb2": max(t1 / float(np.quantile(w1, 0.9875)), t2 / float(np.quantile(w2, 0.9875)))}


# ---------------------------------------------------------------- independence bundle

def tagged(kernel: str, rules: dict) -> dict:
    return {f"{key}_{kernel}": v for key, v in rules.items()}


def kernel_columns(kernels) -> list[str]:
    out = []
    for k in kernels:
        out += [f"{key}_{k}" for key in ("t", "m", "q", "p_wb", "r_wb", "p_shift", "r_shift", "check")]
        if k in SUB_KERNELS:
            out += [f"{key}_{k}" for key in ("p_sub", "rej_sub", "r_sub")]
        if k in ("herm1", "herm2"):
            out += [f"{key}_{k}" for key in ("r975_wb", "r975_shift")]
        if k in HARRIS_KERNELS:
            out += [f"{key}_{k}" for key in ("rhc", "rh19", "rha", "xh19")]
    return out


HSPLIT_IND_COLS = ["hl_x1", "hl_x2", "hl_y1", "hl_y2", "T_11", "T_12", "T_21", "T_22", "rl1", "rl2", "rb1", "rb2"]
HSPLIT_MMD_COLS = ["hl_1", "hl_2", "T_1", "T_2", "rl1", "rl2", "rb1", "rb2"]

BASE_COLS = ["sx", "sy", "rx_rank", "ry_rank", "dcov_stat", "p_dcov_wb", "r_dcov_wb", "p_dcov_sub", "rej_dcov_sub",
             "r_dcov_sub", "pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw", "r_nw", "ewc_nu", "t_ewc", "rej_ewc",
             "r_ewc", "nw_lrv_nonpos"]
EXACT_COLS = [f"xraw_{k}" for k in ("sx", "sy", "t", "p_wb", "r_wb", "p_shift", "r_shift", "check", "p_sub", "rej_sub",
                                    "r_sub")]


def ind_columns(kernels) -> list[str]:
    return BASE_COLS + EXACT_COLS + HSPLIT_IND_COLS + kernel_columns(kernels)


def raw_exact(x, y, dx, dy, w, n: int) -> dict:
    """EXP6 hsic_fields at kappa = 1: same functions, same order, same arithmetic."""
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    sx, sy = mb.median_bandwidth(dx), mb.median_bandwidth(dy)
    r = rules_full(mb.double_center(mb.gauss_gram(dx, sx)), mb.double_center(mb.gauss_gram(dy, sy)), w)
    sub_h = mb.subsample_stats(x, y, ln, dn, "hsic", sx, sy)
    q_hb = q95(sub_h)
    r.update({"sx": sx, "sy": sy, "p_sub": float(np.mean(sub_h >= r["t"])), "rej_sub": int(r["t"] > q_hb),
              "r_sub": r["t"] / q_hb})
    return r


def ind_bundle(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, kernels=IND_KERNELS, B: int = B_WILD,
               HX: float | None = None, HY: float | None = None, exact: bool = False, full: bool = True,
               hsplit: bool = True, dist: bool = True, strict: bool = False) -> dict:
    """All independence rules for one pair on one set of multipliers.

    exact: also compute the exact raw path (xraw_* columns). full = False skips the n x n matrices: no
    dCov, dcov_inc1, Pearson or exact columns. dist = False skips dCov and dcov_inc1 only. hsplit: also test
    the removed Hermite layers. strict: also compute the literal Harris form at N = 19. No random draw follows
    the multipliers, so these switches and kernel subsets never change any other column.
    """
    n = x.size
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    out: dict = {}
    mult = Multipliers(n, B, rng)
    w = mult.w
    sx, sy = median_bw_fast(x), median_bw_fast(y)
    if full:
        if dist or exact:
            dx = np.abs(x[:, None] - x[None, :])
            dy = np.abs(y[:, None] - y[None, :])
        if exact:
            out.update({f"xraw_{k}": v for k, v in raw_exact(x, y, dx, dy, w, n).items() if k not in ("m", "q", "r975_wb", "r975_shift")})
        if dist:
            out.update(dcov_full_rules(x, y, dx, dy, w, n, ln, dn, wild=wild_prefix))
            if "dcov_inc1" in kernels:
                ka = mk.increment(-dx, 1)
                la = mk.increment(-dy, 1)
                out.update(tagged("dcov_inc1", rules_distance(ka, la, mult.get(n - 1))))
                del ka, la
        if dist or exact:
            del dx, dy
        out.update(pearson_rules(x, y))
    out["sx"], out["sy"] = sx, sy
    fx = gauss_features(x, sx)
    fy = gauss_features(y, sy)
    out["rx_rank"], out["ry_rank"] = fx.shape[1], fy.shape[1]
    for k in kernels:
        if k == "dcov_inc1":
            continue
        if k in ("sym", "datadiff"):
            ux, uy = own_series(k, x), own_series(k, y)
            gx = gauss_features(ux, median_bw_fast(ux))
            gy = gauss_features(uy, median_bw_fast(uy))
        else:
            gx = transform_features(k, fx, x, HX)
            gy = transform_features(k, fy, y, HY)
        cx, cy = compress(gx), compress(gy)
        out.update(tagged(k, rules_fast(cx, cy, mult.get(gx.shape[0]), sub=k in SUB_KERNELS)))
        if k in HARRIS_KERNELS:
            out.update(tagged(k, harris_rules(cx, cy, n)))
    if hsplit:
        out.update(hsplit_ind(x, y))
    if strict:
        out.update(harris_strict19(x, y))
    return out


# ---------------------------------------------------------------- split-MMD bundle

MMD_EXACT_COLS = ["xraw_sigma", "xraw_mmd2", "xraw_p_wb", "xraw_r_wb"]


def mmd_columns(kernels) -> list[str]:
    out = ["sigma", "rank"] + MMD_EXACT_COLS + HSPLIT_MMD_COLS
    for k in kernels:
        out += [f"{key}_{k}" for key in ("mmd2", "p_wb", "r_wb")]
        if k in ("herm1", "herm2"):
            out.append(f"r975_wb_{k}")
    return out


def mmd_raw_exact(x: np.ndarray, s: np.ndarray, a: np.ndarray) -> dict:
    """EXP3 mmd_ext for the uncentred wild rule: same arithmetic (kmem.mmd.split_mmd_wild). a = w * s."""
    d = np.abs(x[:, None] - x[None, :])
    sigma = M.median_heuristic(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    stat, boot = split_mmd_wild(k, s, a)
    return {"xraw_sigma": sigma, "xraw_mmd2": stat, "xraw_p_wb": pboot(stat, boot), "xraw_r_wb": stat / q95(boot)}


def mmd_fast_rule(g: np.ndarray, s: np.ndarray, a: np.ndarray) -> tuple[float, float, float, float]:
    """MMD^2_V = s' K s and its uncentred wild draws a' Kc a, from features g of K. Centring does not change s' K s."""
    psi = g - g.mean(axis=0, keepdims=True)
    v = psi.T @ s
    stat = float(v @ v)
    b = psi.T @ a
    boot = (b * b).sum(axis=0)
    return stat, pboot(stat, boot), stat / q95(boot), stat / float(np.quantile(boot, 0.975))


def mmd_bundle(x: np.ndarray, rng: np.random.Generator, kernels=MMD_KERNELS, B: int = M.B_WILD,
               exact: bool = False, hsplit: bool = True) -> dict:
    """Split MMD^2_V with the tested uncentred wild rule (mmdcal form 'unc', l = floor(n^0.5)).

    All kernels use the same two chains. datadiff has one point fewer: it drops the last row of the first
    chain and splits at h - 1. herm_r uses w' Q K Q w with the contrast w.
    """
    n = x.size
    h1 = n // 2
    s = contrast(n, h1)
    l = max(1, math.floor(n ** 0.5))
    if math.floor((n - 1) ** 0.5) != l:
        raise ValueError(f"n = {n}: floor(sqrt(n - 1)) differs from floor(sqrt(n))")
    w = M.multipliers("unc", n, h1, l, rng, B)
    a = w * s[:, None]
    out: dict = mmd_raw_exact(x, s, a) if exact else {}
    sigma = median_bw_fast(x)
    f = gauss_features(x, sigma)
    out["sigma"], out["rank"] = sigma, f.shape[1]
    for kn in kernels:
        if kn == "raw":
            st, pv, rv, r975 = mmd_fast_rule(f, s, a)
        elif kn in ("herm1", "herm2"):
            u = mk.hermite_basis(x, int(kn[-1]))
            st, pv, rv, r975 = mmd_fast_rule(f - u @ (u.T @ f), s, a)
            out[f"r975_wb_{kn}"] = r975
        elif kn == "sym":
            v = mk.even_map(x)
            st, pv, rv, r975 = mmd_fast_rule(gauss_features(v, median_bw_fast(v)), s, a)
        elif kn == "datadiff":
            v = mk.first_difference(x)
            sd = contrast(n - 1, h1 - 1)
            st, pv, rv, r975 = mmd_fast_rule(gauss_features(v, median_bw_fast(v)), sd, np.delete(w, h1 - 1, axis=0) * sd[:, None])
        else:
            raise ValueError(kn)
        out.update({f"mmd2_{kn}": st, f"p_wb_{kn}": pv, f"r_wb_{kn}": rv})
    if hsplit:
        out.update(hsplit_mmd(x))
    return out


# ---------------------------------------------------------------- bare-statistic growth

GROWTH_COLS = ["t_raw", "t_inc1", "sx", "sy", "rx_rank", "ry_rank"]


def growth_bundle(x: np.ndarray, y: np.ndarray) -> dict:
    """T = (1/m) sum_ij Kc_ij Lc_ij for raw and inc1 only. Fast path, no pruning of feature pairs."""
    sx, sy = median_bw_fast(x), median_bw_fast(y)
    fx, fy = gauss_features(x, sx), gauss_features(y, sy)
    out = {"sx": sx, "sy": sy, "rx_rank": fx.shape[1], "ry_rank": fy.shape[1]}
    for name, gx, gy in (("raw", fx, fy), ("inc1", fx[1:] - fx[:-1], fy[1:] - fy[:-1])):
        px = gx - gx.mean(axis=0, keepdims=True)
        py = gy - gy.mean(axis=0, keepdims=True)
        c = px.T @ py
        out[f"t_{name}"] = float((c * c).sum() / px.shape[0])
    return out
