"""Analysis for the kernel-linear arm, exactly as fixed in the parameters_fixed event.

Usage: python -I analyze.py --cells cells --out results_auto.md
Writes cell_summary.csv, knot_analysis.json and a markdown report of every cell.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
Z = 1.959963984540054
KIND_TESTS = {"hsic_wb": ("p_hsic_wb", "p"), "hsic_shift": ("p_hsic_shift", "p"), "hsic_sub": ("rej_hsic_sub", "rej"),
              "dcov_wb": ("p_dcov_wb", "p"), "dcov_sub": ("rej_dcov_sub", "rej")}
LIND_TESTS = {"pearson_nw": ("rej_nw", "rej"), "pearson_ewc": ("rej_ewc", "rej")}
KNOT_GRID = np.round(np.arange(1.10, 1.70 + 1e-9, 0.005), 3)
BOOT = 2000
BOOT_SEED = 20261015003


def wilson(k, R):
    if R == 0:
        return float("nan"), float("nan")
    p = k / R
    den = 1 + Z * Z / R
    c = (p + Z * Z / (2 * R)) / den
    h = Z * math.sqrt(p * (1 - p) / R + Z * Z / (4 * R * R)) / den
    return c - h, c + h


def load(cells_dir):
    man = json.load(open(os.path.join(cells_dir, "manifest.json")))
    out = []
    for c in man["cells"]:
        if c.get("dropped"):
            continue
        path = os.path.join(cells_dir, c["name"] + ".csv")
        rows = list(csv.DictReader(open(path))) if os.path.exists(path) else []
        tests = KIND_TESTS if c["family"] == "kind" else LIND_TESTS
        ok = [r for r in rows if r["status"] == "ok"]
        errs = [r for r in rows if r["status"] != "ok"]
        cpu = sum(float(r["cpu_s"]) for r in rows)
        for tname, (col, kind) in tests.items():
            vals = []
            for r in ok:
                v = r[col]
                if v in ("", "nan", "-1"):
                    continue
                vals.append(float(v) <= ALPHA if kind == "p" else int(float(v)) == 1)
            R = len(vals)
            k = int(sum(vals))
            lo, hi = wilson(k, R)
            out.append({**{kk: c[kk] for kk in ("name", "family", "proc", "sid", "HX", "HY", "lam", "n", "r")},
                        "sumH": round(c["HX"] + c["HY"], 4), "test": tname, "reps_target": c["reps"],
                        "reps_done": len(rows), "reps_ok": len(ok), "errors": len(errs), "R_used": R,
                        "rej": k, "size": k / R if R else float("nan"), "ci_lo": lo, "ci_hi": hi,
                        "cpu_s_cell": cpu,
                        "error_examples": "; ".join(sorted({e["error"] for e in errs})[:2])})
    return man, out


def fit_knot(x, y):
    best = (np.inf, np.nan, np.nan)
    for k in KNOT_GRID:
        z = np.clip(x - k, 0, None)
        if np.count_nonzero(z) < 1 or np.count_nonzero(x < k) < 2:
            continue
        X = np.column_stack([np.ones_like(x), x, z])
        if np.linalg.matrix_rank(X) < 3:
            continue
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        sse = float(np.sum((y - X @ beta) ** 2))
        if sse < best[0] - 1e-15:
            best = (sse, float(k), float(beta[2]))
    return best[1], best[2]


def knot_analysis(summ, test, nmax):
    pts = [s for s in summ if s["test"] == test and s["n"] == nmax and s["r"] == 0 and s["proc"] in ("G1", "G2", "G3", "G4")]
    x = np.array([p["sumH"] for p in pts])
    y = np.array([p["size"] - ALPHA for p in pts])
    k_hat, c_hat = fit_knot(x, y)
    rng = np.random.default_rng(BOOT_SEED)
    ks, cs = [], []
    for _ in range(BOOT):
        idx = rng.integers(0, len(x), len(x))
        kb, cb = fit_knot(x[idx], y[idx])
        if np.isfinite(kb):
            ks.append(kb)
            cs.append(cb)
    ks, cs = np.array(ks), np.array(cs)
    return {"test": test, "n_max": nmax, "n_points": len(x), "knot_hat": k_hat, "slope_change_hat": c_hat,
            "knot_ci": [float(np.quantile(ks, 0.025)), float(np.quantile(ks, 0.975))] if ks.size else [None, None],
            "slope_change_ci": [float(np.quantile(cs, 0.025)), float(np.quantile(cs, 0.975))] if cs.size else [None, None],
            "boot_valid": int(ks.size)}


def grows(summ, test, sid, nmin, nmax):
    a = next(s for s in summ if s["test"] == test and s["sid"] == sid and s["n"] == nmin and s["r"] == 0)
    b = next(s for s in summ if s["test"] == test and s["sid"] == sid and s["n"] == nmax and s["r"] == 0)
    se = lambda s: math.sqrt(max(s["size"] * (1 - s["size"]), 1e-12) / s["R_used"])  # noqa: E731
    diff = b["size"] - a["size"]
    return diff, diff > Z * math.sqrt(se(a) ** 2 + se(b) ** 2)


def fmt(s):
    return f"{100 * s['size']:.1f} [{100 * s['ci_lo']:.1f}, {100 * s['ci_hi']:.1f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cells_dir = os.path.join(HERE, a.cells)
    man, summ = load(cells_dir)
    with open(os.path.join(cells_dir, "cell_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0].keys()))
        w.writeheader()
        w.writerows(summ)

    ns = {"kind": sorted({s["n"] for s in summ if s["family"] == "kind"}),
          "lind": sorted({s["n"] for s in summ if s["family"] == "lind"})}
    tests = {"kind": list(KIND_TESTS), "lind": list(LIND_TESTS)}
    sids = []
    for s in summ:
        if s["sid"] not in sids:
            sids.append(s["sid"])
    sumh = {s["sid"]: s["sumH"] for s in summ}

    knots, p1, growth = {}, {}, {}
    for fam in ("kind", "lind"):
        nmin, nmax = ns[fam][0], ns[fam][-1]
        for t in tests[fam]:
            knots[t] = knot_analysis(summ, t, nmax)
            below_ok, above_ok, detail = True, True, []
            for sid in sids:
                if sid.startswith("R1"):
                    continue
                cell = next(s for s in summ if s["test"] == t and s["sid"] == sid and s["n"] == nmax and s["r"] == 0)
                ex = cell["size"] - ALPHA
                d, g = grows(summ, t, sid, nmin, nmax)
                growth[(t, sid)] = (d, g)
                sh = sumh[sid]
                if sh < 1.5 - 1e-9:
                    below_ok &= ex <= 0.02 + 1e-12
                elif sh > 1.5 + 1e-9:
                    above_ok &= g
                detail.append({"sid": sid, "sumH": sh, "excess_nmax": ex, "growth": d, "grows": g})
            p1[t] = {"below_ok": below_ok, "above_grows": above_ok, "supported": below_ok and above_ok, "cells": detail}

    # KC1
    nk = max(s["n"] for s in summ if s["family"] == "kind")
    kc1_cells = [s for s in summ if s["test"] in ("hsic_wb", "dcov_wb") and s["n"] == nk and s["r"] == 0
                 and s["proc"] in ("G1", "G2", "G3")]
    kc1 = bool(kc1_cells) and all(0.035 <= s["size"] <= 0.065 for s in kc1_cells)
    # KC2 per K-ind test
    kc2 = {}
    for t in tests["kind"]:
        kn = knots[t]
        excl = kn["knot_ci"][0] is not None and not (kn["knot_ci"][0] <= 1.5 <= kn["knot_ci"][1])
        h06 = [s["size"] - ALPHA for s in summ if s["test"] == t and s["n"] == ns["kind"][-1] and s["r"] == 0
               and s["proc"] in ("G1", "G2", "G3") and s["HX"] == 0.6]
        m06 = float(np.mean(h06))
        nochange = kn["slope_change_ci"][0] is not None and kn["slope_change_ci"][0] <= 0 <= kn["slope_change_ci"][1]
        kc2[t] = {"knot_ci_excludes_1.5": excl, "mean_excess_H0.6": m06, "slope_change_ci_contains_0": nochange,
                  "triggered": excl or (m06 > 0.05 and nochange)}

    json.dump({"knots": knots, "P1": p1, "KC1": {"triggered": kc1, "n_cells": len(kc1_cells),
               "sizes": {f"{s['test']}|{s['sid']}": s["size"] for s in kc1_cells}}, "KC2": kc2},
              open(os.path.join(cells_dir, "knot_analysis.json"), "w"), indent=1, default=float)

    L = []
    L.append("## Cell tables: empirical size at 5%, % [Wilson 95% MC interval]\n")
    L.append("Null cells (r = 0). Every cell that was run is listed. 'err' = replications that raised an error. "
             "MC s.e. at 5% nominal: 0.49 pp at 2,000 reps, 0.69 pp at 1,000.\n")
    dropped = [c["name"] for c in man["cells"] if c.get("dropped")]
    if dropped:
        L.append(f"Dropped by author decision (not run): {len(dropped)} cells: " + ", ".join(dropped) + "\n")
    for fam in ("kind", "lind"):
        for t in tests[fam]:
            L.append(f"\n### {t}\n")
            L.append("| setting | H_X+H_Y | " + " | ".join(f"n={n}" for n in ns[fam]) + " | reps (n_max) | err |")
            L.append("|---|---|" + "---|" * len(ns[fam]) + "---|---|")
            for sid in sids:
                row = [s for s in summ if s["test"] == t and s["sid"] == sid and s["r"] == 0]
                row.sort(key=lambda s: s["n"])
                last = row[-1]
                sh = "—" if sid.startswith("R1") else f"{sumh[sid]:.2f}"
                L.append(f"| {sid} | {sh} | " + " | ".join(fmt(s) for s in row) +
                         f" | {last['R_used']} | {sum(s['errors'] for s in row)} |")
    L.append("\n## Power check: G2, H = 0.7, Betken–Dehling 'linear' r = 0.25 (rejection %, [95% MC])\n")
    L.append("| test | " + " | ".join(f"n={n}" for n in sorted({s['n'] for s in summ if s['r']})) + " |")
    L.append("|---|" + "---|" * len({s['n'] for s in summ if s['r']}))
    for fam in ("kind", "lind"):
        for t in tests[fam]:
            row = sorted([s for s in summ if s["test"] == t and s["r"]], key=lambda s: s["n"])
            cols = {s["n"]: fmt(s) for s in row}
            L.append(f"| {t} | " + " | ".join(cols.get(n, "—") for n in sorted({s['n'] for s in summ if s['r']})) + " |")
    L.append("\n## One-knot spline on excess at the largest n (x = H_X + H_Y, 22 G1–G4 cells)\n")
    L.append("| test | n_max | knot | knot 95% CI | slope change c | c 95% CI | CI contains 1.5 |")
    L.append("|---|---|---|---|---|---|---|")
    for t, kn in knots.items():
        ci = kn["knot_ci"]
        cc = kn["slope_change_ci"]
        L.append(f"| {t} | {kn['n_max']} | {kn['knot_hat']:.3f} | [{ci[0]:.3f}, {ci[1]:.3f}] | "
                 f"{100 * kn['slope_change_hat']:.1f} pp/unit | [{100 * cc[0]:.1f}, {100 * cc[1]:.1f}] | "
                 f"{'yes' if ci[0] <= 1.5 <= ci[1] else 'no'} |")
    L.append("\n## P1 rule per test\n")
    L.append("| test | all sum<1.5 cells: excess at n_max ≤ 2 pp | all sum>1.5 cells: excess grows with n | P1 supported |")
    L.append("|---|---|---|---|")
    for t, v in p1.items():
        L.append(f"| {t} | {v['below_ok']} | {v['above_grows']} | **{v['supported']}** |")
    L.append("\nGrowth with n (size at n_max minus size at n_min, pp; * = significant by the fixed rule):\n")
    L.append("| setting | sum | " + " | ".join(p1) + " |")
    L.append("|---|---|" + "---|" * len(p1))
    for sid in sids:
        if sid.startswith("R1"):
            continue
        L.append(f"| {sid} | {sumh[sid]:.2f} | " + " | ".join(
            f"{100 * growth[(t, sid)][0]:+.1f}{'*' if growth[(t, sid)][1] else ''}" for t in p1) + " |")
    L.append("\n## Kill criteria\n")
    L.append(f"- KC1 (no bite): triggered = **{kc1}** ({len(kc1_cells)} wild-bootstrap G1–G3 cells at n = {nk}; "
             f"in band [3.5, 6.5]%: {sum(0.035 <= s['size'] <= 0.065 for s in kc1_cells)}).")
    for t, v in kc2.items():
        L.append(f"- KC2 ({t}): knot CI excludes 1.5 = {v['knot_ci_excludes_1.5']}; mean G1–G3 excess at H = 0.6 = "
                 f"{100 * v['mean_excess_H0.6']:.1f} pp; slope-change CI contains 0 = {v['slope_change_ci_contains_0']}; "
                 f"triggered = **{v['triggered']}**")
    tot_cpu = sum({s["name"]: s["cpu_s_cell"] for s in summ}.values())
    tot_err = sum({s["name"]: s["errors"] for s in summ}.values())
    L.append(f"\nTotal CPU in these cells: {tot_cpu / 3600:.2f} CPU-hours. Replications with errors: {tot_err}.")
    open(os.path.join(HERE, a.out), "w").write("\n".join(L) + "\n")
    print(f"wrote {a.out}; CPU {tot_cpu / 3600:.2f} h; errors {tot_err}")


if __name__ == "__main__":
    main()
