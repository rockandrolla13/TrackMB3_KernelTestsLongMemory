"""Build the single data file for the examples report.

Reads EXP7 outputs (summary_main.csv, results.json, main/*.csv) and runs one seeded
simulation for example 1 with kmem code. Writes data/report_numbers.json.
Every number in the report text and every figure comes from that file.

Run from anywhere: python -I build_data.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EXP = REPO / "experiments" / "EXP7-memory-kernel-20261101"
OUT = HERE / "data" / "report_numbers.json"
sys.path.insert(0, str(REPO))

from kmem import memory_kernels as mk  # noqa: E402
from kmem.legacy import mbcore  # noqa: E402

SEED_EX1 = 20261103
H_EX1 = 0.9
N_EX1 = 2000
MAXLAG = 1000


def load_summary() -> pd.DataFrame:
    return pd.read_csv(EXP / "summary_main.csv", dtype={"H": str, "n": str})


def rows(d: pd.DataFrame, **kw) -> pd.DataFrame:
    q = d
    for k, v in kw.items():
        q = q[q[k] == v]
    return q


def one(d: pd.DataFrame, **kw) -> dict:
    q = rows(d, **kw)
    if len(q) != 1:
        raise KeyError(f"expected 1 row for {kw}, got {len(q)}")
    r = q.iloc[0]
    return {"v": float(r.value), "lo": float(r.lo), "hi": float(r.hi), "m": int(r.m),
            "note": None if pd.isna(r.note) else str(r.note)}


# ---------------------------------------------------------------- example 1

def sample_acf(x: np.ndarray, maxlag: int) -> np.ndarray:
    x = x - x.mean()
    n = x.size
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[: maxlag + 1] / n
    return ac / ac[0]


def increment_acf_theory(H: float, maxlag: int) -> np.ndarray:
    rho = mbcore.fgn_acov(H, maxlag + 2)
    rho_full = np.concatenate([[rho[1]], rho])  # rho(-1) = rho(1) at index 0
    j = np.arange(maxlag + 1)
    g = 2 * rho_full[j + 1] - rho_full[j + 2] - rho_full[j]  # 2 rho(j) - rho(j+1) - rho(j-1)
    return g / g[0]


def example1() -> dict:
    rng = np.random.default_rng(SEED_EX1)
    x = mbcore.fgn(H_EX1, N_EX1, rng)
    y = mbcore.fgn(H_EX1, N_EX1, rng)
    dx, dy = mk.first_difference(x), mk.first_difference(y)
    ax, ay = sample_acf(x, MAXLAG), sample_acf(y, MAXLAG)
    adx, ady = sample_acf(dx, MAXLAG), sample_acf(dy, MAXLAG)
    th_raw = mbcore.fgn_acov(H_EX1, MAXLAG + 1)
    th_inc = increment_acf_theory(H_EX1, MAXLAG)
    # Running pair sum over lags |j| <= L of acf_X(j) * acf_Y(j). This drives the null mean of HSIC.
    def pair_sum(a, b):
        p = a * b
        return (p[0] + 2 * np.cumsum(p[1:])).tolist()
    lags_keep = [1, 2, 5, 10, 50, 100, 500]
    return {
        "H": H_EX1, "n": N_EX1, "seed": SEED_EX1, "maxlag": MAXLAG,
        "generator": "kmem.legacy.mbcore.fgn (Davies-Harte), increments by kmem.memory_kernels.first_difference",
        "x": x.tolist(), "dx": dx.tolist(),
        "acf_x": ax.tolist(), "acf_dx": adx.tolist(),
        "acf_theory_raw": th_raw.tolist(), "acf_theory_inc": th_inc.tolist(),
        "pairsum_sample_raw": pair_sum(ax, ay), "pairsum_sample_inc": pair_sum(adx, ady),
        "pairsum_theory_raw": pair_sum(th_raw, th_raw), "pairsum_theory_inc": pair_sum(th_inc, th_inc),
        "acf_at": {str(L): {"raw": float(ax[L]), "inc": float(adx[L]),
                            "raw_theory": float(th_raw[L]), "inc_theory": float(th_inc[L])}
                   for L in lags_keep},
    }


# ---------------------------------------------------------------- examples 2 to 6

def growth(d: pd.DataFrame) -> dict:
    out = {}
    for H in ["0.5", "0.6", "0.7", "0.8", "0.9", "0.95"]:
        for k in ["raw", "inc1"]:
            q = {n: one(d, family="growth", H=H, n=n, kernel=k, kind="growth_q95")
                 for n in ["500", "1000", "2000", "4000", "8000"]}
            s = one(d, family="growth", H=H, n="all", kernel=k, kind="growth_slope")
            out[f"{k}|{H}"] = {"q95": q, "slope": s, "ratio_8000_500": q["8000"]["v"] / q["500"]["v"]}
    return out


SIZE_KERNELS = {
    "wild": ["raw", "inc1", "herm1", "datadiff", "dcov_inc1", "hsplit1", "hsplit2"],
    "shift": ["raw", "inc1", "herm1", "herm2", "sym", "datadiff", "dcov_inc1", "hsplit1", "hsplit2"],
    "sub": ["raw", "inc1"],
    "harris": ["raw", "inc1"],
}
HS = ["0.5", "0.6", "0.7", "0.8", "0.9", "0.95"]


def size_vs_h(d: pd.DataFrame) -> dict:
    out = {}
    for n in ["500", "2000", "4000"]:
        for rule, ks in SIZE_KERNELS.items():
            for k in ks:
                for H in HS:
                    out[f"{k}|{rule}|{H}|{n}"] = one(d, family="ind", gen="G0", H=H, n=n,
                                                     kernel=k, rule=rule, kind="size")
    return out


POWER_SERIES = [("raw", "shift"), ("inc1", "shift"), ("herm1", "shift"), ("dcov_inc1", "shift"),
                ("raw", "harris"), ("datadiff", "shift"), ("hsplit1", "shift"), ("hsplit2", "shift")]
ALTS = ["lin", "sq", "lag", "flat", "factor", "scale"]
PH = ["0.5", "0.7", "0.8", "0.9"]


def power(d: pd.DataFrame) -> dict:
    out = {}
    for a in ALTS:
        for k, r in POWER_SERIES:
            for H in PH:
                out[f"{a}|{k}|{r}|{H}"] = one(d, family="ind", alt=a, H=H, kernel=k, rule=r, kind="adj_power")
    return out


def mmd_power(d: pd.DataFrame) -> dict:
    out = {}
    for a, prm in [("mean", 0.1), ("mean", 0.5), ("scale", 1.1), ("shape", 0.1)]:
        for k in ["raw", "herm1", "hsplit1", "hsplit2", "datadiff"]:
            for H in ["0.5", "0.7", "0.8"]:
                q = rows(d, family="mmd", alt=a, H=H, kernel=k, rule="wild", kind="adj_power")
                q = q[np.isclose(q.param, prm)]
                r = q.iloc[0]
                out[f"{a}{prm}|{k}|{H}"] = {"v": float(r.value), "lo": float(r.lo), "hi": float(r.hi)}
    return out


ARMS = ["SEAS-M-shared", "SEAS-V-slow", "SEAS-M", "SEAS-V", "NL-SV", "FACTOR"]
ARM_KERNELS = [("raw", r) for r in ["wild", "shift", "sub", "harris"]] + \
              [("inc1", r) for r in ["wild", "shift", "sub", "harris"]] + \
              [(k, r) for k in ["herm1", "herm2", "sym", "datadiff", "hsplit1", "hsplit2"] for r in ["wild", "shift"]] + \
              [("hsplit1", "low"), ("hsplit2", "low")]


def arms(d: pd.DataFrame) -> dict:
    out = {}
    for g in ARMS:
        for H in (["0.5", "0.8"] if g != "FACTOR" else ["0.5", "0.7", "0.8", "0.9"]):
            for k, r in ARM_KERNELS:
                q = rows(d, family="ind", gen=g, H=H, kernel=k, rule=r, kind="size")
                if len(q) == 1:
                    out[f"{g}|{k}|{r}|{H}"] = one(d, family="ind", gen=g, H=H, kernel=k, rule=r, kind="size")
    return out


def harris_stationary_range(d: pd.DataFrame) -> dict:
    """Range of (raw, harris) size over stationary null cells (all except SEAS-M-shared and SEAS-V-slow)."""
    q = rows(d, family="ind", kernel="raw", rule="harris", kind="size")
    q = q[~q.gen.isin(["SEAS-M-shared", "SEAS-V-slow"])]
    return {"min": float(q.value.min()), "max": float(q.value.max()), "cells": int(len(q))}


HARRIS_GROUP = {"G0": "Gaussian memory", "HT-L": "process", "HT-S": "process", "NL-SV": "process",
                "NL-TAR": "process", "SEAS-M": "process", "SEAS-V": "process", "FACTOR": "stress null",
                "SCALE": "stress null", "SEAS-M-shared": "non-stationary", "SEAS-V-slow": "non-stationary"}


def harris_cells(d: pd.DataFrame) -> list:
    """Every independence null cell: raw and inc1 under the Harris rule, with its cell group."""
    q = rows(d, family="ind", kind="size", rule="harris")
    out = []
    for (g, H, n), grp in q.groupby(["gen", "H", "n"]):
        r = {k: float(grp[grp.kernel == k].value.iloc[0]) for k in ["raw", "inc1"]}
        out.append({"gen": g, "H": H, "n": n, "group": HARRIS_GROUP[g], **r})
    return out


def sv_layers() -> dict:
    m = pd.read_csv(EXP / "main" / "ind__NL-SV__null0__H0.8__n2000.csv")
    m = m[m.status == "ok"]
    return {"layer1": float(np.mean([m.hl_x1.mean(), m.hl_y1.mean()])),
            "layer2": float(np.mean([m.hl_x2.mean(), m.hl_y2.mean()])), "reps": int(len(m))}


# ---------------------------------------------------------------- review additions (REVIEW.md, 2026-11-03)
# These follow the reviewer's scripts (scratchpad/review/trend.py, binom_boot.py, theory_lambda.py,
# mmd_size.py), reimplemented here so every number lives in this one data file.

ZW = 1.959963984540054
MMD_H = ["0.5", "0.6", "0.7", "0.75", "0.8", "0.85", "0.9"]
MMD_K = ["raw", "herm1", "herm2", "sym", "datadiff", "hsplit1", "hsplit2"]
GAUSS_HP = [(0.5, 0.5), (0.6, 0.6), (0.7, 0.7), (0.8, 0.8), (0.9, 0.9), (0.95, 0.95), (0.6, 0.9), (0.7, 0.9)]
GAUSS_N = [500, 2000, 4000]
TREND_RULES = [(k, r) for k in ["inc1", "herm1", "herm2", "sym", "datadiff", "dcov_inc1", "hp10", "hsplit1", "hsplit2"]
               for r in ["wild", "shift"]] + [("raw", "wild"), ("raw", "shift"), ("frac_oracle", "wild")]


def wilson(k: int, m: int) -> tuple[float, float]:
    p = k / m
    dd = 1 + ZW * ZW / m
    c = (p + ZW * ZW / (2 * m)) / dd
    h = ZW * np.sqrt(p * (1 - p) / m + ZW * ZW / (4 * m * m)) / dd
    return c - h, c + h


def size_class(k: int, m: int) -> str:
    lo, hi = wilson(k, m)
    if (lo <= 0.05 <= hi) or hi <= 0.075:
        return "valid"
    return "mild" if k / m <= 0.10 else "over"


def mmd_size(d: pd.DataFrame) -> dict:
    """Split-MMD null rates, wild rule, from summary_main.csv, plus the reviewer's two memory checks."""
    out = {}
    for g in ["G0", "HT-S"]:
        for k in MMD_K:
            cells = {}
            for n in ["500", "2000"]:
                for H in MMD_H:
                    r = one(d, family="mmd", gen=g, H=H, n=n, kernel=k, rule="wild", kind="size")
                    cells[f"{H}|{n}"] = r
                    out[f"{g}|{k}|{H}|{n}"] = r
            above = [f"{H}|{n}" for n in ["500", "2000"] for H in MMD_H
                     if cells[f"{H}|{n}"]["lo"] > cells[f"0.5|{n}"]["hi"]]
            rises = [H for H in MMD_H if cells[f"{H}|2000"]["lo"] > cells[f"{H}|500"]["hi"]]
            out[f"{g}|{k}|flat"] = {"above_baseline": above, "rises_500_to_2000": rises,
                                    "memory_flat": not above and not rises}
        for n_ in ["500", "2000"]:  # block variant of hsplit2 (lead's request, 2026-11-03)
            for H in MMD_H:
                r = one(d, family="mmd", gen=g, H=H, n=n_, kernel="hsplit2", rule="wild_blk", kind="size")
                out[f"{g}|hsplit2_blk|{H}|{n_}"] = r
    return out


def _gname(h, n):
    hx, hy = h
    hs = f"H{hx}" if hx == hy else f"HX{hx}_HY{hy}"
    return f"ind__G0__null0__{hs}__n{n}"


def _rej(df: pd.DataFrame, k: str, r: str) -> np.ndarray:
    if k.startswith("hsplit"):
        q = k[-1]
        pc = "p_shift" if r == "shift" else "p_wb"
        return ((df[f"{pc}_herm{q}"] <= 0.025) | (df[f"rl{q}"] > 1)).values
    return (df[f"p_{'wb' if r == 'wild' else 'shift'}_{k}"] <= 0.05).values


def trend_and_pooled() -> dict:
    """Pooled rate over the 24 Gaussian cells and a logistic slope in H (n fixed effects, 18 equal-H cells)."""
    import statsmodels.api as sm

    frames = {(h, n): pd.read_csv(EXP / "main" / f"{_gname(h, n)}.csv") for h in GAUSS_HP for n in GAUSS_N}
    out = {}
    for k, r in TREND_RULES:
        kk = mm = 0
        ys, xs = [], []
        for (h, n), df in frames.items():
            v = _rej(df, k, r)
            kk += int(v.sum())
            mm += v.size
            if h[0] == h[1]:
                ys.append(v.astype(float))
                xs.append(np.column_stack([np.full(v.size, h[0]), np.full(v.size, n == 2000), np.full(v.size, n == 4000)]))
        fit = sm.Logit(np.concatenate(ys), sm.add_constant(np.vstack(xs).astype(float))).fit(disp=0)
        lo, hi = wilson(kk, mm)
        out[f"{k}|{r}"] = {"pooled": kk / mm, "lo": lo, "hi": hi, "m": mm,
                           "slope": float(fit.params[1]), "p": float(fit.pvalues[1])}
    return out


def binom_all24() -> dict:
    """Exact chance that a rule at exactly 5% is 'valid' in one cell, and in all 24 Gaussian cells."""
    from scipy import stats

    def fail(m):
        return float(sum(stats.binom.pmf(k, m, 0.05) for k in range(m + 1) if size_class(k, m) != "valid"))

    f1, f3 = fail(1000), fail(300)
    return {"fail_cell_1000": f1, "fail_cell_300": f3, "all24": (1 - f1) ** 16 * (1 - f3) ** 8}


def _et_inc(H: float, n: int, sig: float) -> float:
    """Fixed-bandwidth exact E T for the increment kernel (THEORY.md identity, as in theory_lambda.py)."""
    k = np.arange(n + 3.0)
    rho = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))
    kap = sig / np.sqrt(sig ** 2 + 2 * (1 - rho)) - sig / np.sqrt(sig ** 2 + 2)
    kd = np.empty(n)
    kd[0] = 2 * kap[0] - 2 * kap[1]
    j = np.arange(1, n)
    kd[1:] = 2 * kap[j] - kap[j + 1] - kap[j - 1]
    w = 1 - np.arange(n) / n
    return float(kd[0] ** 2 + 2 * np.sum(w[1:] * kd[1:] ** 2))


def bandwidth() -> dict:
    s0 = 0.6744897501960817 * np.sqrt(2)  # median heuristic for N(0, 1)
    out = {"sigma_limit": s0}
    for H in [0.5, 0.8, 0.9, 0.95]:
        out[f"lambda|{H}"] = _et_inc(H, 200000, s0)
        for n in [500, 1000, 2000, 4000, 8000]:
            g = pd.read_csv(EXP / "main" / f"growth__G0__null0__H{H}__n{n}.csv")
            sx = float(g.sx.mean())
            out[f"{H}|{n}"] = {"mean_T_inc1": float(g.t_inc1.mean()), "mean_sigma": sx,
                               "formula_at_mean_sigma": _et_inc(H, n - 1, sx)}
    bc = pd.read_csv(EXP / "bandwidth_check.csv")  # sim-kernel's fixed-bandwidth rerun, 2026-11-03
    for _, r in bc.iterrows():
        out[f"check|{r.H}|{int(r.n)}"] = {"inc1_ratio_fixed": float(r.ratio_to_n500_inc1_fixed),
                                         "inc1_ratio_sample": float(r.ratio_to_n500_inc1_sample)}
    return out


def main() -> None:
    d = load_summary()
    res = json.loads((EXP / "results.json").read_text())
    data = {
        "_about": "Generated by build_data.py from EXP7 outputs. Do not edit by hand. Rates are fractions.",
        "source_files": ["experiments/EXP7-memory-kernel-20261101/summary_main.csv",
                         "experiments/EXP7-memory-kernel-20261101/results.json",
                         "experiments/EXP7-memory-kernel-20261101/main/ind__NL-SV__null0__H0.8__n2000.csv"],
        "run": {k: res[k] for k in ["seed_master_main", "replications", "cells", "rows", "error_rows", "cpu_hours"]},
        "verdicts": res["verdicts_gaussian_memory_cells"],
        "verdicts_mmd_G0": res["verdicts_mmd_G0"],
        "audit": res["audit"],
        "ex1": example1(),
        "growth": growth(d),
        "size": size_vs_h(d),
        "power": power(d),
        "mmd_power": mmd_power(d),
        "arms": arms(d),
        "harris_stationary": harris_stationary_range(d),
        "sv_layers": sv_layers(),
        "layer_true_vs_est": res["layer_memory_true_vs_estimated"],
        "mmd_size": mmd_size(d),
        "trend": trend_and_pooled(),
        "binom": binom_all24(),
        "bandwidth": bandwidth(),
        "harris_cells": harris_cells(d),
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(data, indent=1))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
