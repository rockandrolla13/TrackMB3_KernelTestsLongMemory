"""HSIC, dCov and Pearson rules on full n x n matrices: the pieces EXP3, EXP6 and EXP7 each re-typed (CR-001).

The arithmetic and call order are those of the EXP2 core (mbcore.kind_tests, mbcore.lind_tests), with each rule's
95% reference quantile read off the same draws, as EXP3 power_core.kind_ext added. The only random draw is the
wild multiplier matrix w, which the caller draws once and passes in, so these pieces can run in any order.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.bandwidth import median_full_matrix
from kmem.calibration.common import NORMAL_975, pboot, q95, t975
from kmem.calibration.wild import markov_rademacher_mbcore
from kmem.gram import abs_diff
from kmem.legacy import mbcore as mb


def hsic_full_rules(x, y, dx, dy, sx, sy, w, n, ln, dn, A) -> dict:
    """HSIC with Gaussian bandwidths sx, sy: statistic, wild bootstrap, circular shift, subsampling (EXP6 hsic_fields)."""
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


def dcov_full_rules(x, y, dx, dy, w, n, ln, dn, wild=mb.wild_boot_stats) -> dict:
    """Distance covariance: statistic, wild bootstrap, subsampling. EXP7 passes wild=kmem.calibration.wild.wild_prefix."""
    ac = mb.double_center(dx)
    bc = mb.double_center(dy)
    md = ac * bc
    del ac, bc
    t_d = md.sum() / n
    bd = wild(md, w)
    del md
    sub_d = mb.subsample_stats(x, y, ln, dn, "dcov")
    q_dw, q_db = q95(bd), q95(sub_d)
    return {"dcov_stat": t_d, "p_dcov_wb": pboot(t_d, bd), "r_dcov_wb": t_d / q_dw,
            "p_dcov_sub": float(np.mean(sub_d >= t_d)), "rej_dcov_sub": int(t_d > q_db), "r_dcov_sub": t_d / q_db}


def pearson_rules(x: np.ndarray, y: np.ndarray) -> dict:
    """mbcore.lind_tests plus the ratios r_nw and r_ewc of each |t| to its 5% two-sided cut-off."""
    lt = mb.lind_tests(x, y)
    out = dict(lt)
    out["r_nw"] = abs(lt["t_nw"]) / NORMAL_975 if np.isfinite(lt["t_nw"]) else float("nan")
    out["r_ewc"] = abs(lt["t_ewc"]) / t975(lt["ewc_nu"])
    return out


def hsic_dcov_ext(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, B: int = 500) -> dict:
    """EXP3 power_core.kind_ext: mbcore.kind_tests plus 95% reference quantiles, and the Pearson rules."""
    n = x.size
    ln = math.floor(n ** 0.5)  # wild-bootstrap ln and subsampling block length
    dn = math.floor(0.1 * n)
    A = math.floor(0.1 * n)  # shift range [A, n - A]
    dx, dy = abs_diff(x), abs_diff(y)
    sx, sy = median_full_matrix(dx), median_full_matrix(dy)
    w = markov_rademacher_mbcore(n, B, ln, rng)
    out = hsic_full_rules(x, y, dx, dy, sx, sy, w, n, ln, dn, A)
    out.update(dcov_full_rules(x, y, dx, dy, w, n, ln, dn))
    out.update(pearson_rules(x, y))
    return out
