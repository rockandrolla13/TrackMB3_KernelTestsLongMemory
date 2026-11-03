"""EXP5 analysis (see DESIGN.md).

Usage: python -I -B analyze_exp5.py --in main --tag main
  -> summary_<tag>.csv (Part 1 power), summary_ind_<tag>.csv (Part 2 rates),
     tables_<tag>.md, results_<tag>.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))  # kmem (D7: not installed)
import numpy as np  # noqa: E402

import run_exp5 as R  # noqa: E402
from kmem.legacy import common as C  # noqa: E402
from kmem.analysis import boot_power, cstar, require_input  # noqa: E402,F401

EXPS = os.path.join(HERE, "..")
EXP3_MAIN = os.path.join(EXPS, "EXP3-power-20261025", "main")
E3_SUMMARY = os.path.join(EXPS, "EXP1-20261016", "e3_summary.csv")
KL_SUMMARY = os.path.join(EXPS, "JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear", "cells", "cell_summary.csv")
ALPHA = 0.05
NBOOT = 2000
SEED = 42
Z975 = 1.959963984540054

MMD_RULES = {
    "mmd_unc_wild": ("r_unc_05", lambda r: float(r["p_unc_05"]) <= ALPHA),
    "mmd_sub_sqrt": ("r_sub_sqrt", lambda r: float(r["t_n"]) > float(r["q95_sqrt"])),
    "mmd_sn_block": ("r_sn_05", lambda r: float(r["sn_T"]) > float(r["sn_q95_05"])),
}
MMD_ORACLE = {"oracle_mmd": lambda r: float(r["mmd2"]),
              "lin_oracle": lambda r: abs(float(r["hm_diff"])),
              "lin_t": lambda r: abs(float(r["t_naive"]))}
IND_RULES = {"hsic_wb": ("p_hsic_wb", "p"), "hsic_shift": ("p_hsic_shift", "p"), "hsic_sub": ("rej_hsic_sub", "rej"),
             "dcov_wb": ("p_dcov_wb", "p"), "dcov_sub": ("rej_dcov_sub", "rej"),
             "pearson_nw": ("rej_nw", "rej"), "pearson_ewc": ("rej_ewc", "rej")}
IND_MAIN = ["hsic_wb", "hsic_shift", "hsic_sub", "dcov_sub", "pearson_nw", "pearson_ewc"]


def load(path):
    rows = list(csv.DictReader(open(path)))
    ok = [r for r in rows if r["status"] == "ok"]
    return ok, len(rows) - len(ok), sum(float(r["cpu_s"]) for r in rows)


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def lin_exact(c):
    return norm_cdf(c - Z975) + norm_cdf(-c - Z975)


def fmt(p, lo, hi):
    return f"{100 * p:.1f} [{100 * lo:.1f}, {100 * hi:.1f}]"


# ------------------------------------------------------------------ Part 1

def part1(indir, cells, rng):
    data = {}
    for c in cells:
        if c["family"] != "mmd":
            continue
        p = os.path.join(HERE, indir, c["name"] + ".csv")
        if os.path.exists(p):
            data[(c["alt"], c["c"], c["H"], c["n"])] = (c, *load(p))
    rows, boots, cpu, errs = [], {}, 0.0, 0
    for k, (c, ok, err, cp) in data.items():
        cpu += cp
        errs += err
    # null regression against EXP3 (same rules, same generator)
    reg = []
    for (alt, cc, h, n), (c, ok, err, cp) in data.items():
        if alt != "null":
            continue
        p3 = os.path.join(EXP3_MAIN, f"mmd__G0__null0__H{h}__n{n}.csv")
        old = [r for r in csv.DictReader(open(require_input(p3, "EXP3 main run, null regression"))) if r["status"] == "ok"]
        for rule, (rc, raw) in MMD_RULES.items():
            k1, m1 = sum(raw(r) for r in ok), len(ok)
            rec = {"H": h, "n": n, "rule": rule, "exp5_size": k1 / m1, "exp5_m": m1,
                   "exp5_cstar": cstar(np.array([float(r[rc]) for r in ok]))}
            if old:
                k2, m2 = sum(raw(r) for r in old), len(old)
                pp = (k1 + k2) / (m1 + m2)
                se = math.sqrt(pp * (1 - pp) * (1 / m1 + 1 / m2))
                rec.update({"exp3_size": k2 / m2, "exp3_m": m2, "z": (k1 / m1 - k2 / m2) / se,
                            "exp3_cstar": cstar(np.array([float(r[rc]) for r in old]))})
            reg.append(rec)
    for (alt, cc, h, n), (c, ok, err, cp) in sorted(data.items(), key=lambda t: str(t[0])):
        if alt == "null" or ("null", 0.0, h, n) not in data:
            continue
        nrows = data[("null", 0.0, h, n)][1]
        bn = rng.integers(0, len(nrows), size=(NBOOT, len(nrows)))
        ba = rng.integers(0, len(ok), size=(NBOOT, len(ok)))
        for rule, spec in {**MMD_RULES, **MMD_ORACLE}.items():
            if rule in MMD_RULES:
                rc, raw = spec
                nv = np.array([float(r[rc]) for r in nrows])
                av = np.array([float(r[rc]) for r in ok])
            else:
                raw = None
                nv = np.array([spec(r) for r in nrows])
                av = np.array([spec(r) for r in ok])
            cs = cstar(nv)
            k = int(np.sum(av > cs))
            lo, hi = C.wilson(k, av.size)
            bp = boot_power(nv, av, bn, ba)
            boots[(cc, h, n, rule)] = bp
            rec = {"c": cc, "H": h, "n": n, "delta": c["delta"], "rule": rule, "adj_k": k, "adj_m": av.size,
                   "adj_power": k / av.size, "adj_wilson_lo": lo, "adj_wilson_hi": hi,
                   "adj_boot_lo": float(np.quantile(bp, 0.025)), "adj_boot_hi": float(np.quantile(bp, 0.975)),
                   "cstar": cs, "null_rate_at_cstar": float(np.mean(nv > cs)), "null_reps": nv.size}
            if raw is not None:
                kr, ks = sum(raw(r) for r in ok), sum(raw(r) for r in nrows)
                rec.update({"raw_power": kr / len(ok), "null_raw_size": ks / len(nrows)})
            if rule == "lin_oracle":
                rec["exact_power"] = lin_exact(cc)
            rows.append(rec)
    # reading rule: H = 0.5 vs H = 0.8 at constant c
    verdicts = {}
    for cc in sorted({r["c"] for r in rows}):
        for n in R.N_MMD:
            for rule in list(MMD_RULES) + list(MMD_ORACLE):
                a = [r for r in rows if (r["c"], r["n"], r["rule"], r["H"]) == (cc, n, rule, 0.5)]
                b = [r for r in rows if (r["c"], r["n"], r["rule"], r["H"]) == (cc, n, rule, 0.8)]
                if not a or not b:
                    continue
                d = a[0]["adj_power"] - b[0]["adj_power"]
                db = boots[(cc, 0.5, n, rule)] - boots[(cc, 0.8, n, rule)]
                lo, hi = float(np.quantile(db, 0.025)), float(np.quantile(db, 0.975))
                excl = lo > 0 or hi < 0
                if abs(d) <= 0.10:
                    v = "flat"
                elif d > 0.10 and excl:
                    v = "collapse" if b[0]["adj_power"] < 0.5 * a[0]["adj_power"] else "partial loss"
                elif d < -0.10 and excl:
                    v = "gain"
                else:
                    v = "unclear"
                verdicts[f"c{cc:g}|n{n}|{rule}"] = {"p_H0.5": round(100 * a[0]["adj_power"], 1),
                                                     "p_H0.8": round(100 * b[0]["adj_power"], 1),
                                                     "drop_pp": round(100 * d, 1),
                                                     "drop_boot_pp": [round(100 * lo, 1), round(100 * hi, 1)],
                                                     "verdict": v}
    return rows, reg, verdicts, cpu, errs, len(data)


# ------------------------------------------------------------------ Part 2

def part2(indir, cells):
    out, cpu, errs = [], 0.0, 0
    for c in cells:
        if c["family"] != "ind":
            continue
        p = os.path.join(HERE, indir, c["name"] + ".csv")
        if not os.path.exists(p):
            continue
        ok, err, cp = load(p)
        cpu += cp
        errs += err
        for t, (col, kind) in IND_RULES.items():
            fl = [(float(r[col]) <= ALPHA) if kind == "p" else (int(float(r[col])) == 1)
                  for r in ok if r[col] not in ("", "nan", "-1")]
            k, m = sum(fl), len(fl)
            lo, hi = C.wilson(k, m)
            out.append({"source": "EXP5", "sumH": c["sumH"], "H": c["H"], "n": c["n"], "test": t, "rej": k,
                        "R": m, "size": k / m, "ci_lo": lo, "ci_hi": hi,
                        "cpu_s_rep": cp / max(1, len(ok) + err)})
    prior = []
    for r in csv.DictReader(open(require_input(E3_SUMMARY, "EXP1 E3 summary"))):
        if r["HX"] == r["HY"] and r["HX"] in ("0.7", "0.8"):
            prior.append({"source": "E3", "sumH": float(r["sumH"]), "H": float(r["HX"]), "n": int(r["n"]),
                          "test": r["test"], "rej": int(r["rej"]), "R": int(r["R"]), "size": float(r["size"]),
                          "ci_lo": float(r["ci_lo"]), "ci_hi": float(r["ci_hi"])})
    seen = set()
    for r in csv.DictReader(open(require_input(KL_SUMMARY, "kernel-linear job, KL baseline"))):
        if r["proc"] == "G1" and r["HX"] == "0.5" and r["HY"] == "0.5" and float(r["r"] or 0) == 0:
            key = (int(r["n"]), r["test"])
            if key in seen or r["test"] not in IND_RULES:
                continue
            seen.add(key)
            prior.append({"source": "KL", "sumH": 1.0, "H": 0.5, "n": int(r["n"]), "test": r["test"],
                          "rej": int(r["rej"]), "R": int(r["R_used"]), "size": float(r["size"]),
                          "ci_lo": float(r["ci_lo"]), "ci_hi": float(r["ci_hi"])})
    allr = prior + out
    # reads: z vs n = 2000, WLS slope on log2 n
    reads = {}
    for s in (1.0, 1.4, 1.6):
        for t in IND_RULES:
            ser = sorted([r for r in allr if r["sumH"] == s and r["test"] == t], key=lambda r: r["n"])
            ref = [r for r in ser if r["n"] == 2000 and r["source"] != "EXP5"]
            rec = {"rates": {str(r["n"]): round(100 * r["size"], 2) for r in ser}}
            for big in (4000, 8000):
                new = [r for r in ser if r["n"] == big and r["source"] == "EXP5"]
                if ref and new:
                    a, b = ref[0], new[0]
                    se = math.sqrt(a["size"] * (1 - a["size"]) / a["R"] + b["size"] * (1 - b["size"]) / b["R"])
                    se_b = math.sqrt(b["size"] * (1 - b["size"]) / b["R"])
                    rec[f"z_{big}_vs_2000"] = round((b["size"] - a["size"]) / se, 2) if se > 0 else float("nan")
                    rec[f"diff_{big}_vs_2000_pp"] = round(100 * (b["size"] - a["size"]), 2)
                    rec[f"excess_{big}_over5_z"] = round((b["size"] - ALPHA) / se_b, 2) if se_b > 0 else float("nan")
            pts = [r for r in ser if r["n"] >= 500 and 0 < r["size"] < 1]
            if len(pts) >= 3:
                x = np.log2([r["n"] for r in pts])
                y = np.array([r["size"] for r in pts])
                w = np.array([r["R"] / (r["size"] * (1 - r["size"])) for r in pts])
                X = np.c_[np.ones_like(x), x]
                A = X.T @ (w[:, None] * X)
                beta = np.linalg.solve(A, X.T @ (w * y))
                se = math.sqrt(np.linalg.inv(A)[1, 1])
                rec["slope_pp_per_doubling"] = round(100 * beta[1], 2)
                rec["slope_se_pp"] = round(100 * se, 2)
                rec["slope_ns"] = [r["n"] for r in pts]
            reads[f"sum{s}|{t}"] = rec
    return allr, reads, cpu, errs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    cells = R.cell_list()
    rng = np.random.default_rng(SEED)
    rows, reg, verdicts, cpu1, err1, ncell1 = part1(a.indir, cells, rng)
    allr, reads, cpu2, err2 = part2(a.indir, cells)

    def dump(path, recs):
        if not recs:
            return
        fields = []
        for r in recs:
            for k in r:
                if k not in fields:
                    fields.append(k)
        with open(os.path.join(HERE, path), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(recs)

    dump(f"summary_{a.tag}.csv", rows)
    dump(f"null_regression_{a.tag}.csv", reg)
    dump(f"summary_ind_{a.tag}.csv", allr)
    with open(os.path.join(HERE, f"results_{a.tag}.json"), "w") as f:
        json.dump({"part1_cells_loaded": ncell1, "errors": err1 + err2,
                   "cpu_hours_part1": round(cpu1 / 3600, 2), "cpu_hours_part2": round(cpu2 / 3600, 2),
                   "nboot": NBOOT, "seed": SEED, "lin_exact_power": {f"{c:g}": lin_exact(c) for c in R.C_GRID},
                   "part1_verdicts_H05_vs_H08": verdicts, "part2_reads": reads}, f, indent=1)
    # tables
    L = [f"# EXP5 tables ({a.tag})", "",
         "Percent. Exploratory simulation.", "",
         "## Part 1: split MMD, shift δ_H = c · sd_H (null-SD units), size-adjusted power", "",
         "`adj` = size-adjusted power [bootstrap 95% interval incl. threshold error]. `raw` = nominal-5% power.",
         "`size` = raw null rejection. lin_oracle = |half-mean difference|; lin_t = naive Welch t; both size-adjusted.", ""]
    for cc in sorted({r["c"] for r in rows}, reverse=True):
        for n in R.N_MMD:
            sub = [r for r in rows if r["c"] == cc and r["n"] == n]
            if not sub:
                continue
            ds = ", ".join(f"H={h}: {next(r['delta'] for r in sub if r['H'] == h):.3f}"
                           for h in R.H_GRID if any(r["H"] == h for r in sub))
            L += [f"### c = {cc:g}, n = {n} (δ: {ds}; exact linear power {100 * lin_exact(cc):.1f})", "",
                  "| rule | " + " | ".join(f"H={h} adj | raw | size" for h in R.H_GRID) + " |",
                  "|---|" + "---|---|---|" * len(R.H_GRID)]
            for rule in list(MMD_RULES) + list(MMD_ORACLE):
                cl = []
                for h in R.H_GRID:
                    rr = [r for r in sub if r["rule"] == rule and r["H"] == h]
                    if not rr:
                        cl.append("– | – | –")
                        continue
                    r = rr[0]
                    raw = f"{100 * r['raw_power']:.1f}" if "raw_power" in r else "–"
                    sz = f"{100 * r['null_raw_size']:.1f}" if "null_raw_size" in r else "–"
                    cl.append(f"{fmt(r['adj_power'], r['adj_boot_lo'], r['adj_boot_hi'])} | {raw} | {sz}")
                L.append(f"| {rule} | " + " | ".join(cl) + " |")
            L.append("")
    L += ["### Reading rule: adjusted power at H = 0.5 minus H = 0.8, same c", "",
          "| c, n, rule | H=0.5 | H=0.8 | drop pp [boot 95%] | verdict |", "|---|---|---|---|---|"]
    for k, v in verdicts.items():
        L.append(f"| {k} | {v['p_H0.5']} | {v['p_H0.8']} | {v['drop_pp']} [{v['drop_boot_pp'][0]}, "
                 f"{v['drop_boot_pp'][1]}] | {v['verdict']} |")
    L += ["", "### Null regression: EXP5 fresh MMD nulls vs EXP3 main nulls (raw size, c*)", "",
          "| H, n | rule | EXP5 % | EXP3 % | z | c* EXP5 | c* EXP3 |", "|---|---|---|---|---|---|---|"]
    for r in reg:
        L.append(f"| H{r['H']} n{r['n']} | {r['rule']} | {100 * r['exp5_size']:.1f} | "
                 f"{100 * r.get('exp3_size', float('nan')):.1f} | {r.get('z', float('nan')):.2f} | "
                 f"{r['exp5_cstar']:.3f} | {r.get('exp3_cstar', float('nan')):.3f} |")
    L += ["", "## Part 2: independent fGn pairs, rejection % at nominal 5% [Wilson 95%]", "",
          "Sources: KL = kernel-linear arm (H = 0.5 baseline, 2,000 reps), E3 = EXP1 E3 (2,000 reps), "
          "EXP5 = this run.", ""]
    ns = sorted({r["n"] for r in allr if r["n"] >= 500})
    for s in (1.0, 1.4, 1.6):
        L += [f"### memory sum {s} (H_X = H_Y = {s / 2:g})", "",
              "| rule | " + " | ".join(f"n={n}" for n in ns) + " | Δ(4000−2000) pp, z | slope pp/doubling (SE) |",
              "|---|" + "---|" * len(ns) + "---|---|"]
        for t in IND_RULES:
            cl = []
            for n in ns:
                rr = [r for r in allr if r["sumH"] == s and r["test"] == t and r["n"] == n]
                cl.append(f"{100 * rr[0]['size']:.1f} [{100 * rr[0]['ci_lo']:.1f}, {100 * rr[0]['ci_hi']:.1f}]"
                          if rr else "–")
            rd = reads.get(f"sum{s}|{t}", {})
            dz = (f"{rd['diff_4000_vs_2000_pp']:+.1f}, z={rd['z_4000_vs_2000']:+.2f}"
                  if "z_4000_vs_2000" in rd else "–")
            sl = f"{rd['slope_pp_per_doubling']:+.2f} ({rd['slope_se_pp']:.2f})" if "slope_pp_per_doubling" in rd else "–"
            L.append(f"| {t} | " + " | ".join(cl) + f" | {dz} | {sl} |")
        L.append("")
    L += [f"CPU: Part 1 {cpu1 / 3600:.2f} h, Part 2 {cpu2 / 3600:.2f} h. Errors: {err1 + err2}."]
    with open(os.path.join(HERE, f"tables_{a.tag}.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"part1 cells {ncell1}, cpu {cpu1 / 3600:.2f} h; part2 cpu {cpu2 / 3600:.2f} h; errors {err1 + err2}")


if __name__ == "__main__":
    main()
