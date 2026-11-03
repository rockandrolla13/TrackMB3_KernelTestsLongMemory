"""EXP7 post-review checks (2026-11-03), answering investigations/memory-kernel-2026-11-01/REVIEW.md.

Usage: python -I review_checks.py  -> trend_main.csv, binomial_valid.csv, bandwidth_check.csv

1. Trend test (REVIEW M1, M2). For each rule, a logistic regression of the 0/1 rejection on H with n fixed effects. Independence uses the 18
   equal-H Gaussian cells (6 H x n = 500, 2000, 4000). Split MMD uses the 14 cells of each generator (7 H x n = 500, 2000). Also the pooled
   rate over all cells with a Wilson interval, and the rates pooled over H <= 0.6 and H >= 0.8.
2. Exact binomial chance that a rule at a true rate p is "valid" in every cell (REVIEW M1), at p = 0.05 and at the rule's pooled rate.
3. Bandwidth check (REVIEW M4). The growth rows store the median bandwidths. A fresh run recomputes T_inc1 and T_raw at a fixed bandwidth
   sigma = sqrt(2) * Phi^{-1}(0.75) = 0.9539, the median heuristic of N(0, 1), on the same draws as the sample-bandwidth T.
   H in {0.9, 0.95}, n in {500, 2000, 8000}, 1,000 replications, master 2026110177.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import csv  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import numpy as np  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy import special, stats  # noqa: E402

import analyze_exp7 as A  # noqa: E402
import exp7_core as K  # noqa: E402
import run_exp7 as R  # noqa: E402
from kmem.legacy import mbcore as mb  # noqa: E402

SIGMA_FIXED = math.sqrt(2.0) * float(special.ndtri(0.75))


def fail_prob(p: float, m: int) -> float:
    """P(a cell is not 'valid') when the true rate is p, exact binomial."""
    k = np.arange(m + 1)
    pmf = stats.binom.pmf(k, m, p)
    return float(sum(pmf[i] for i in k if A.size_class(int(i), m) != "valid"))


def p_all_valid(p: float, ms: list[int]) -> float:
    return float(np.prod([1.0 - fail_prob(p, m) for m in ms]))


def trend_rows(data: dict, family: str, gen: str, rules: list[tuple]) -> list[dict]:
    cells = [(c, d) for c, d in data.values() if c["family"] == family and c["gen"] == gen and c["alt"] == "null"]
    ns = sorted({c["n"] for c, _ in cells})
    out = []
    for kern, rule, rc, raw in rules:
        if raw is None or any(rc not in d for _, d in cells):
            continue
        ys, xs, ms = [], [], []
        k_all = m_all = 0
        lo = [0, 0]
        hi = [0, 0]
        for c, d in cells:
            ok = np.isfinite(d[rc])
            v = raw(d)[ok].astype(float)
            k_all += int(v.sum())
            m_all += v.size
            ms.append(v.size)
            if c["HX"] == c["HY"]:
                ys.append(v)
                xs.append(np.column_stack([np.full(v.size, c["HX"])] + [np.full(v.size, float(c["n"] == n)) for n in ns[1:]]))
                if c["HX"] <= 0.6:
                    lo[0] += int(v.sum())
                    lo[1] += v.size
                if c["HX"] >= 0.8:
                    hi[0] += int(v.sum())
                    hi[1] += v.size
        fit = sm.Logit(np.concatenate(ys), sm.add_constant(np.vstack(xs))).fit(disp=0)
        p = k_all / m_all
        wlo, whi = A.wilson(k_all, m_all)
        out.append({"family": family, "gen": gen, "kernel": kern, "rule": rule, "cells": len(cells), "pooled_rate": p, "pooled_lo": wlo, "pooled_hi": whi,
                    "rate_H_le_0.6": lo[0] / lo[1], "rate_H_ge_0.8": hi[0] / hi[1], "logit_slope_per_unit_H": float(fit.params[1]),
                    "slope_se": float(fit.bse[1]), "slope_p": float(fit.pvalues[1]), "p_all_cells_valid_at_pooled_rate": p_all_valid(p, ms),
                    "p_all_cells_valid_at_5pc": p_all_valid(0.05, ms)})
    return out


def bw_chunk(args):
    ih, h, n, c, reps = args
    rng = np.random.default_rng(np.random.SeedSequence(entropy=2026110177, spawn_key=(ih, n, c)))
    out = np.empty((reps, 6))
    for i in range(reps):
        x, y = mb.fgn(h, n, rng), mb.fgn(h, n, rng)
        sx, sy = K.median_bw_fast(x), K.median_bw_fast(y)
        row = [sx, sy]
        for sig_x, sig_y in ((sx, sy), (SIGMA_FIXED, SIGMA_FIXED)):
            fx, fy = K.gauss_features(x, sig_x), K.gauss_features(y, sig_y)
            for gx, gy in ((fx, fy), (fx[1:] - fx[:-1], fy[1:] - fy[:-1])):
                px, py = gx - gx.mean(axis=0), gy - gy.mean(axis=0)
                cc = px.T @ py
                row.append(float((cc * cc).sum() / px.shape[0]))
        out[i] = row[:2] + row[2:6]
    return ih, h, n, out


def bandwidth_check(workers: int = 4, reps: int = 1000) -> list[dict]:
    hs, ns = (0.9, 0.95), (500, 2000, 8000)
    tasks = [(ih, h, n, c, 100) for ih, h in enumerate(hs) for n in ns for c in range(reps // 100)]
    store: dict = {}
    with mp.get_context("fork").Pool(workers) as pool:
        for ih, h, n, out in pool.imap_unordered(bw_chunk, tasks):
            store.setdefault((h, n), []).append(out)
    rows = []
    for h in hs:
        base = {}
        for n in ns:
            v = np.vstack(store[(h, n)])
            q = {name: float(np.quantile(v[:, j], 0.95)) for j, name in ((2, "raw_sample"), (3, "inc1_sample"), (4, "raw_fixed"), (5, "inc1_fixed"))}
            if n == ns[0]:
                base = q
            rows.append({"H": h, "n": n, "reps": v.shape[0], "mean_bandwidth": float(v[:, :2].mean()), "sigma_fixed": SIGMA_FIXED,
                         **{f"q95_{k}": val for k, val in q.items()}, **{f"ratio_to_n500_{k}": val / base[k] for k, val in q.items()}})
    return rows


def write(name: str, rows: list[dict]) -> None:
    with open(os.path.join(HERE, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    write("binomial_valid.csv", [{"true_rate": 0.05, "m": m, "p_cell_not_valid": fail_prob(0.05, m)} for m in (1000, 300)]
          + [{"true_rate": 0.05, "m": "16 x 1000 + 8 x 300", "p_cell_not_valid": 1.0 - p_all_valid(0.05, [1000] * 16 + [300] * 8)}])
    data = A.load_all("main", R.main_cells())
    rows = trend_rows(data, "ind", "G0", A.ind_rules())
    rows += trend_rows(data, "mmd", "G0", A.mmd_rules()) + trend_rows(data, "mmd", "HT-S", A.mmd_rules())
    write("trend_main.csv", rows)
    # the growth rows: mean median bandwidth by H and n
    gb = []
    for c, d in data.values():
        if c["family"] == "growth":
            gb.append({"H": c["HX"], "n": c["n"], "mean_bandwidth": float(np.mean(np.r_[d["sx"], d["sy"]])), "mean_T_inc1": float(np.mean(d["t_inc1"])),
                       "mean_T_raw": float(np.mean(d["t_raw"]))})
    write("growth_bandwidth.csv", sorted(gb, key=lambda r: (r["H"], r["n"])))
    write("bandwidth_check.csv", bandwidth_check())
    print("wrote binomial_valid.csv, trend_main.csv, growth_bandwidth.csv, bandwidth_check.csv", flush=True)


if __name__ == "__main__":
    main()
