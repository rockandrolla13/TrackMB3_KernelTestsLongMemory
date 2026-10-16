"""Analysis fixed before the main run (parameters_fixed event): per-cell size / coverage minus nominal
with Wilson 95% Monte Carlo intervals; hockey-stick breakpoint at the largest n / T; P2, P3, P4, KC3,
KC4 verdicts.

usage: python -I analyze.py {pilot,main}
Writes <tag>/summary_khom.csv, <tag>/summary_dml.csv, <tag>/analysis.json, <tag>/tables.md.
"""
import csv
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

import common as C  # noqa: E402
import dml  # noqa: E402
import khom  # noqa: E402

BREAK_PP = 2.0  # "breaks" = excess beyond 2 pp (size > 7% or coverage < 93%), as in P1/P4/KC wording
KAPPA_GRID = np.round(np.arange(0.5, 0.8995, 0.001), 3)
N_BOOT = 2000


def read_cell(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def rate(rows, key):
    ok = [r for r in rows if r["status"] == "ok"]
    k = sum(int(float(r[key])) for r in ok)
    return k, len(ok)


def hockey_fit(h, y):
    """Profile least squares for y = a + c (h - kappa)_+ over KAPPA_GRID (first minimiser)."""
    z = np.maximum(h[None, :] - KAPPA_GRID[:, None], 0.0)
    zc = z - z.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    szz = np.sum(zc * zc, axis=1)
    c = np.sum(zc * yc[None, :], axis=1) / szz
    a = y.mean() - c * z.mean(axis=1)
    sse = np.sum((y[None, :] - a[:, None] - c[:, None] * z) ** 2, axis=1)
    i = int(np.argmin(sse))
    return float(sse[i]), float(KAPPA_GRID[i]), float(a[i]), float(c[i])


def breakpoint(h, k, n, sign, rng):
    """Excess in pp = sign * (k/n - nominal) handled by caller via y; here k, n per cell."""
    h = np.asarray(h, float)
    p = np.asarray(k, float) / np.asarray(n, float)
    y = sign(p)
    _, kap, a, c = hockey_fit(h, y)
    boots = []
    for _ in range(N_BOOT):
        pb = rng.binomial(np.asarray(n), p) / np.asarray(n, float)
        boots.append(hockey_fit(h, sign(pb))[1])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"kappa": kap, "kappa_ci95": [float(lo), float(hi)], "intercept_pp": a, "slope_pp_per_H": c}


def growth(k1, n1, k2, n2):
    """Change of rate from smallest to largest n with a normal-approximation z."""
    p1, p2 = k1 / n1, k2 / n2
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    return {"delta_pp": 100 * (p2 - p1), "z": (p2 - p1) / se if se > 0 else float("nan")}


def fmt_ci(k, n, nominal, scale=100.0):
    lo, hi = C.wilson(k, n)
    return f"{scale * (k / n - nominal):+.1f} [{scale * (lo - nominal):+.1f}, {scale * (hi - nominal):+.1f}]"


# ------------------------------------------------------------------ K-hom

def analyze_khom(tag, rng, md):
    d = os.path.join(HERE, tag, "khom")
    cells = khom.cell_list()
    out_rows, res = [], {}
    for c in cells:
        path = os.path.join(d, c["cell_id"] + ".csv")
        if not os.path.exists(path):
            continue
        rows = read_cell(path)
        n_err = sum(r["status"] != "ok" for r in rows)
        rec = dict(c, reps=len(rows), errors=n_err)
        for cal in ("wild", "sub"):
            k, n = rate(rows, f"rej_{cal}")
            lo, hi = C.wilson(k, n)
            rec.update({f"rej_{cal}": k / n, f"rej_{cal}_lo": lo, f"rej_{cal}_hi": hi, f"k_{cal}": k, f"n_{cal}": n})
        rec["cpu_s"] = sum(float(r["cpu_s"]) for r in rows)
        res[c["cell_id"]] = rec
        out_rows.append(rec)
    if not out_rows:
        return {}
    keys = list(out_rows[0].keys())
    with open(os.path.join(HERE, tag, "summary_khom.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out_rows)

    ns = khom.N_GRID
    md.append("## K-hom: MMD^2 first half vs second half — size minus 5% (pp) [Wilson 95% MC interval]\n")
    for cal, name in (("wild", "wild bootstrap"), ("sub", "subsampling")):
        md.append(f"### {name}\n")
        md.append("| process | param | " + " | ".join(f"n={n}" for n in ns) + " |")
        md.append("|---|---|" + "---|" * len(ns))
        for proc, params in (("G1", C.HURST_GRID), ("G2", C.HURST_GRID), ("G3", C.HURST_GRID), ("R1", (0.01, 0.1))):
            for prm in params:
                cellsrow = []
                for n in ns:
                    r = res.get(f"khom_{proc}_{prm}_n{n}")
                    cellsrow.append(fmt_ci(r[f"k_{cal}"], r[f"n_{cal}"], 0.05) + (f" E{r['errors']}" if r and r["errors"] else "") if r else "–")
                md.append(f"| {proc} | {'H=' if proc != 'R1' else 'λ='}{prm} | " + " | ".join(cellsrow) + " |")
        md.append("")
    md.append("### Power check (G1, H = 0.7, second half shifted +0.3): rejection rate\n")
    md.append("| calibration | " + " | ".join(f"n={n}" for n in ns) + " |")
    md.append("|---|" + "---|" * len(ns))
    for cal in ("wild", "sub"):
        row = []
        for n in ns:
            r = res.get(f"khom_G1_0.7_n{n}_power")
            row.append(f"{r[f'rej_{cal}']:.3f}" if r else "–")
        md.append(f"| {cal} | " + " | ".join(row) + " |")
    md.append("")

    # breakpoints and P2 for MMD
    nmax, nmin = max(ns), min(ns)
    analysis = {"breakpoints": {}, "P2_mmd": {}}
    for proc in ("G1", "G2", "G3"):
        for cal in ("wild", "sub"):
            rs = [res.get(f"khom_{proc}_{h}_n{nmax}") for h in C.HURST_GRID]
            if any(r is None for r in rs):
                continue
            bp = breakpoint(C.HURST_GRID, [r[f"k_{cal}"] for r in rs], [r[f"n_{cal}"] for r in rs],
                            lambda p: 100 * (p - 0.05), rng)
            per_h = {}
            for h, r in zip(C.HURST_GRID, rs):
                r0 = res[f"khom_{proc}_{h}_n{nmin}"]
                per_h[str(h)] = {
                    "excess_pp_nmax": 100 * (r[f"rej_{cal}"] - 0.05),
                    "ci_pp_nmax": [100 * (r[f"rej_{cal}_lo"] - 0.05), 100 * (r[f"rej_{cal}_hi"] - 0.05)],
                    "breaks": 100 * (r[f"rej_{cal}"] - 0.05) > BREAK_PP,
                    "growth_nmin_to_nmax": growth(r0[f"k_{cal}"], r0[f"n_{cal}"], r[f"k_{cal}"], r[f"n_{cal}"]),
                }
            bp["contains_0.5"] = bp["kappa_ci95"][0] <= 0.5 + 1e-9
            analysis["breakpoints"][f"{proc}_{cal}"] = bp
            breaks_all = all(per_h[str(h)]["breaks"] for h in C.HURST_GRID if h > 0.5)
            grows_all = all(per_h[str(h)]["growth_nmin_to_nmax"]["z"] > 2 for h in C.HURST_GRID if h > 0.5)
            analysis["P2_mmd"][f"{proc}_{cal}"] = {"per_H": per_h, "breaks_every_H_gt_0.5": breaks_all,
                                                   "excess_grows_with_n_every_H_gt_0.5": grows_all,
                                                   "H0.5_breaks": per_h["0.5"]["breaks"]}
    return analysis


# ------------------------------------------------------------------ DML

def analyze_dml(tag, rng, md):
    d = os.path.join(HERE, tag, "dml")
    cells = dml.cell_list()
    out_rows, res = [], {}
    for c in cells:
        path = os.path.join(d, c["cell_id"] + ".csv")
        if not os.path.exists(path):
            continue
        rows = [r for r in read_cell(path) if r["status"] == "ok"]
        n_err = sum(r["status"] != "ok" for r in read_cell(path))
        t = c["T"]
        crit = {"nw": C.Z975, "ewc": C.t_quantile_975(C.ewc_nu(t))}
        path_th = np.array([float(r["theta_path"]) for r in rows])
        base = dict(c, reps=len(rows), errors=n_err,
                    irls_nonconv=sum(int(r["irls_nonconv_iid"]) + int(r["irls_nonconv_blk"]) for r in rows),
                    n_clip=sum(int(r["n_clip_iid"]) + int(r["n_clip_blk"]) for r in rows))
        for arm in dml.ARMS:
            est = np.array([float(r[f"{arm}_est"]) for r in rows])
            for v in ("nw", "ewc"):
                se = np.array([float(r[f"{arm}_se_{v}"]) for r in rows])
                rec = dict(base, arm=arm, var=v)
                for estimand, target in (("pop", np.ones_like(est)), ("path", path_th)):
                    cov = np.abs(est - target) <= crit[v] * se
                    k, n = int(cov.sum()), len(cov)
                    lo, hi = C.wilson(k, n)
                    rec.update({f"cov_{estimand}": k / n, f"cov_{estimand}_lo": lo, f"cov_{estimand}_hi": hi,
                                f"k_{estimand}": k, f"n_{estimand}": n})
                rec["rej_theta1.1"] = float(np.mean(np.abs(est - 1.1) > crit[v] * se))
                rec["bias_pop"] = float(est.mean() - 1.0)
                rec["sd_est"] = float(est.std(ddof=1))
                rec["mean_se"] = float(se.mean())
                rec["se_over_sd_pop"] = float(se.mean() / est.std(ddof=1))
                rec["sd_est_minus_path"] = float((est - path_th).std(ddof=1))
                res[(c["cell_id"], arm, v)] = rec
                out_rows.append(rec)
    if not out_rows:
        return {}
    keys = list(out_rows[0].keys())
    with open(os.path.join(HERE, tag, "summary_dml.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out_rows)

    md.append("## DML (AIPW): coverage minus 95% (pp) [Wilson 95% MC interval]\n")
    for estimand, ename in (("pop", "population ATE = 1"), ("path", "path-conditional ATE")):
        for v in ("nw", "ewc"):
            md.append(f"### {ename}, {'HAC Newey–West (Andrews)' if v == 'nw' else 'EWC fixed-b'}\n")
            md.append("| τ | arm | H | " + " | ".join(f"T={t}" for t in dml.T_GRID) + " |")
            md.append("|---|---|---|" + "---|" * len(dml.T_GRID))
            for tau in dml.TAUS:
                for arm in dml.ARMS:
                    for h in C.HURST_GRID:
                        row = []
                        for t in dml.T_GRID:
                            r = res.get((f"dml_{tau}_H{h}_T{t}", arm, v))
                            row.append(fmt_ci(r[f"k_{estimand}"], r[f"n_{estimand}"], 0.95) if r else "–")
                        md.append(f"| {tau} | {arm} | {h} | " + " | ".join(row) + " |")
            md.append("")
    tmax, tmin = max(dml.T_GRID), min(dml.T_GRID)
    md.append(f"### Power check: rejection of θ0 = 1.1 (truth 1), τ const, H = 0.7\n")
    md.append("| arm | var | " + " | ".join(f"T={t}" for t in dml.T_GRID) + " |")
    md.append("|---|---|" + "---|" * len(dml.T_GRID))
    for arm in dml.ARMS:
        for v in ("nw", "ewc"):
            row = [f"{res[(f'dml_const_H0.7_T{t}', arm, v)]['rej_theta1.1']:.3f}"
                   if (f"dml_const_H0.7_T{t}", arm, v) in res else "–" for t in dml.T_GRID]
            md.append(f"| {arm} | {v} | " + " | ".join(row) + " |")
    md.append("")
    md.append(f"### Diagnostics at T = {tmax} (NW): mean SE / Monte Carlo SD of estimate around 1\n")
    md.append("| τ | arm | " + " | ".join(f"H={h}" for h in C.HURST_GRID) + " |")
    md.append("|---|---|" + "---|" * len(C.HURST_GRID))
    for tau in dml.TAUS:
        for arm in dml.ARMS:
            row = [f"{res[(f'dml_{tau}_H{h}_T{tmax}', arm, 'nw')]['se_over_sd_pop']:.2f}"
                   if (f"dml_{tau}_H{h}_T{tmax}", arm, "nw") in res else "–" for h in C.HURST_GRID]
            md.append(f"| {tau} | {arm} | " + " | ".join(row) + " |")
    md.append("")

    analysis = {"breakpoints": {}, "P2_dml": {}, "P3": {}, "P4": {}, "KC3": {}, "KC4": {}}
    for tau in ("rank1", "rank2", "const"):
        for arm in dml.ARMS:
            for v in ("nw", "ewc"):
                rs = [res.get((f"dml_{tau}_H{h}_T{tmax}", arm, v)) for h in C.HURST_GRID]
                if any(r is None for r in rs):
                    continue
                bp = breakpoint(C.HURST_GRID, [r["k_pop"] for r in rs], [r["n_pop"] for r in rs],
                                lambda p: 100 * (0.95 - p), rng)
                per_h = {}
                for h, r in zip(C.HURST_GRID, rs):
                    r0 = res[(f"dml_{tau}_H{h}_T{tmin}", arm, v)]
                    gr = growth(r0["k_pop"], r0["n_pop"], r["k_pop"], r["n_pop"])
                    per_h[str(h)] = {"deficit_pp_Tmax": 100 * (0.95 - r["cov_pop"]),
                                     "breaks": 100 * (0.95 - r["cov_pop"]) > BREAK_PP,
                                     "deficit_growth_Tmin_to_Tmax_pp": -gr["delta_pp"], "z": -gr["z"]}
                bp["contains_0.5"] = bp["kappa_ci95"][0] <= 0.5 + 1e-9
                bp["contains_0.75"] = bp["kappa_ci95"][0] <= 0.75 <= bp["kappa_ci95"][1]
                analysis["breakpoints"][f"{tau}_{arm}_{v}"] = bp
                if tau == "rank1":
                    analysis["P2_dml"][f"{arm}_{v}"] = {
                        "per_H": per_h,
                        "breaks_every_H_gt_0.5": all(per_h[str(h)]["breaks"] for h in C.HURST_GRID if h > 0.5),
                        "grows_with_T_every_H_gt_0.5": all(per_h[str(h)]["z"] > 2 for h in C.HURST_GRID if h > 0.5),
                        "H0.5_breaks": per_h["0.5"]["breaks"]}
                if tau == "rank2":
                    analysis["P3"][f"{arm}_{v}"] = {
                        "per_H": per_h,
                        "no_break_H_le_0.75": not any(per_h[str(h)]["breaks"] for h in C.HURST_GRID if h <= 0.75),
                        "breaks_H_gt_0.75": all(per_h[str(h)]["breaks"] for h in C.HURST_GRID if h > 0.75),
                        "grows_with_T_H_gt_0.75": all(per_h[str(h)]["z"] > 2 for h in C.HURST_GRID if h > 0.75)}
    for arm in dml.ARMS:
        for v in ("nw", "ewc"):
            worst = None
            for tau in dml.TAUS:
                for h in C.HURST_GRID:
                    for t in dml.T_GRID:
                        r = res.get((f"dml_{tau}_H{h}_T{t}", arm, v))
                        if r and (worst is None or r["cov_path"] < worst[0]):
                            worst = (r["cov_path"], tau, h, t, r["cov_path_lo"], r["cov_path_hi"])
            if worst:
                analysis["P4"][f"{arm}_{v}"] = {"min_path_coverage": worst[0], "at": worst[1:4],
                                                "wilson": worst[4:6], "holds_ge_0.93": worst[0] >= 0.93}
            r = res.get((f"dml_rank1_H0.9_T{tmax}", arm, v))
            if r:
                analysis["KC3"][f"{arm}_{v}"] = {"pop_coverage": r["cov_pop"], "wilson": [r["cov_pop_lo"], r["cov_pop_hi"]],
                                                 "triggered_C1_killed": r["cov_pop"] >= 0.93}
            if arm.startswith(("oracle", "kp")):
                kc4 = {}
                for tau in dml.TAUS:
                    for t in dml.T_GRID:
                        r = res.get((f"dml_{tau}_H0.9_T{t}", arm, v))
                        if r:
                            kc4[f"{tau}_T{t}"] = {"path_coverage": r["cov_path"],
                                                  "wilson": [r["cov_path_lo"], r["cov_path_hi"]],
                                                  "below_0.93": r["cov_path"] < 0.93}
                analysis["KC4"][f"{arm}_{v}"] = {"cells": kc4, "triggered": any(x["below_0.93"] for x in kc4.values())}
    return analysis


def main():
    tag = sys.argv[1]
    rng = np.random.default_rng(np.random.SeedSequence([C.MASTER_SEED, 7777]))
    md = [f"# Tables ({tag})\n", "Generated by analyze.py. Exploratory simulation; not evidence for any accepted claim.\n"]
    out = {"khom": analyze_khom(tag, rng, md), "dml": analyze_dml(tag, rng, md)}
    with open(os.path.join(HERE, tag, "analysis.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    with open(os.path.join(HERE, tag, "tables.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print("written", tag)


if __name__ == "__main__":
    main()
