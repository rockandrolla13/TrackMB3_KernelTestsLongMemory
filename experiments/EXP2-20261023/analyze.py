"""EXP2 analysis: rejection rates, Wilson intervals, and the pre-declared reading rule (DESIGN.md).

Usage: python -I analyze.py --in main --tag main   ->  summary_<tag>.csv, tables_<tag>.md, verdicts_<tag>.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))  # kmem (D7: not installed)
import numpy as np  # noqa: E402

from kmem.legacy import common as C  # noqa: E402
from kmem.legacy import gens  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
H_GRID = (0.5, 0.7, 0.8)
N_GRID = (500, 2000)
RULES = {
    "ind": {"hsic_wild": lambda r: float(r["p_hsic_wb"]) <= ALPHA,
            "hsic_shift": lambda r: float(r["p_hsic_shift"]) <= ALPHA,
            "hsic_sub": lambda r: int(float(r["rej_hsic_sub"])) == 1,
            "dcov_wild": lambda r: float(r["p_dcov_wb"]) <= ALPHA,
            "dcov_sub": lambda r: int(float(r["rej_dcov_sub"])) == 1,
            "pearson_nw": lambda r: int(float(r["rej_nw"])) == 1,
            "pearson_ewc": lambda r: int(float(r["rej_ewc"])) == 1},
    "mmd": {"mmd_unc_wild": lambda r: float(r["p_unc_05"]) <= ALPHA,
            "mmd_sub_sqrt": lambda r: float(r["t_n"]) > float(r["q95_sqrt"]),
            "mmd_sn_block": lambda r: float(r["sn_T"]) > float(r["sn_q95_05"])},
}
VERDICT_RULES = {"ind": ["hsic_wild", "hsic_shift", "hsic_sub", "dcov_wild", "dcov_sub"],
                 "mmd": ["mmd_unc_wild", "mmd_sub_sqrt", "mmd_sn_block"]}
FEATURES = {"heavy tails": ["HT-S", "HT-L"], "nonlinear dynamics": ["NL-SV", "NL-TAR"],
            "seasonality": ["SEAS-M", "SEAS-V"]}


def load(indir: str) -> dict:
    out = {}
    for fam in ("mmd", "ind"):
        gl = list(gens.MAIN) + (["SEAS-M-shared"] if fam == "ind" else [])
        for g in gl:
            for h in H_GRID:
                for n in N_GRID:
                    p = os.path.join(HERE, indir, f"{fam}__{g}__H{h}__n{n}.csv")
                    if not os.path.exists(p):
                        continue
                    rows = list(csv.DictReader(open(p)))
                    ok = [r for r in rows if r["status"] == "ok"]
                    rec = {"reps": len(rows), "errors": len(rows) - len(ok),
                           "cpu_s": sum(float(r["cpu_s"]) for r in rows)}
                    for rule, f in RULES[fam].items():
                        v = [f(r) for r in ok if not (rule == "pearson_nw" and r["rej_nw"] == "-1")]
                        k = int(sum(v))
                        lo, hi = C.wilson(k, len(v))
                        rec[rule] = {"k": k, "m": len(v), "p": k / len(v) if v else float("nan"), "lo": lo, "hi": hi}
                    out[(fam, g, h, n)] = rec
    return out


def interaction(d, fam, g, rule, n, h):
    cells = [d.get((fam, g, h, n)), d.get((fam, g, 0.5, n)), d.get((fam, "G0", h, n)), d.get((fam, "G0", 0.5, n))]
    if any(c is None for c in cells):
        return None
    p = [c[rule]["p"] for c in cells]
    m = [c[rule]["m"] for c in cells]
    i = (p[0] - p[1]) - (p[2] - p[3])
    se = math.sqrt(sum(pp * (1 - pp) / mm for pp, mm in zip(p, m)))
    return i, (i / se if se > 0 else 0.0)


def verdict(d, fam, g, rule, n):
    f0, b0 = d.get((fam, g, 0.5, n)), d.get((fam, "G0", 0.5, n))
    if f0 is None or b0 is None:
        return None
    a, b = f0[rule], b0[rule]
    alone = (a["p"] - b["p"] > 0.02) and (a["lo"] > b["hi"])
    below = (b["p"] - a["p"] > 0.02) and (a["hi"] < b["lo"])
    ints = {h: interaction(d, fam, g, rule, n, h) for h in (0.7, 0.8)}
    ip = any(v and v[0] > 0.02 and v[1] > 1.96 for v in ints.values())
    im = any(v and v[0] < -0.02 and v[1] < -1.96 for v in ints.values())
    if alone and ip:
        v = "breaks alone and amplifies memory"
    elif alone:
        v = "breaks alone"
    elif ip:
        v = "only with memory"
    elif im:
        v = "dampens memory"
    else:
        v = "no effect"
    if im and not ip and alone:
        v = "breaks alone; dampens memory"
    return {"verdict": v, "alone_pp": 100 * (a["p"] - b["p"]), "below_baseline_alone": below,
            "I_pp": {str(h): (None if ints[h] is None else round(100 * ints[h][0], 2)) for h in ints},
            "I_z": {str(h): (None if ints[h] is None else round(ints[h][1], 2)) for h in ints}}


def fmt(c):
    return f"{100 * c['p']:.1f} [{100 * c['lo']:.1f}, {100 * c['hi']:.1f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    d = load(a.indir)
    with open(os.path.join(HERE, f"summary_{a.tag}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["family", "gen", "subordinated", "H", "n", "rule", "rejections", "valid_reps", "errors",
                     "rate", "wilson_lo", "wilson_hi", "cell_cpu_s"])
        for (fam, g, h, n), rec in sorted(d.items()):
            for rule in RULES[fam]:
                c = rec[rule]
                w.writerow([fam, g, gens.SUBORDINATED[g], h, n, rule, c["k"], c["m"], rec["errors"],
                            round(c["p"], 5), round(c["lo"], 5), round(c["hi"], 5), round(rec["cpu_s"], 1)])
    verd = {}
    L = [f"# EXP2 tables ({a.tag})", "",
         "Rejection % at nominal 5% [Wilson 95%]. Rows = generator, columns = H (H_X = H_Y for independence).", ""]
    for fam in ("mmd", "ind"):
        for rule in RULES[fam]:
            for n in N_GRID:
                L += [f"## {rule}, n = {n}", "", "| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |",
                      "|---|---|---|---|---|---|"]
                gl = list(gens.MAIN) + (["SEAS-M-shared"] if fam == "ind" else [])
                for g in gl:
                    cs = [d.get((fam, g, h, n)) for h in H_GRID]
                    if all(c is None for c in cs):
                        continue
                    reps = "/".join(str(c[rule]["m"]) if c else "-" for c in cs)
                    L.append(f"| {g} | {gens.SUBORDINATED[g]} | " + " | ".join(fmt(c[rule]) if c else "-" for c in cs)
                             + f" | {reps} |")
                L.append("")
            if rule in VERDICT_RULES[fam]:
                for n in N_GRID:
                    for g in gens.MAIN[1:]:
                        v = verdict(d, fam, g, rule, n)
                        if v:
                            verd[f"{rule}|n{n}|{g}"] = v
    L += ["# Verdicts (pre-declared rule, DESIGN.md)", "",
          "alone = F(0.5) − G0(0.5) in pp; I(h) = (F(h) − F(0.5)) − (G0(h) − G0(0.5)) in pp, z in brackets.", ""]
    for n in N_GRID:
        L += [f"## n = {n}", "", "| feature | generator | rule | alone pp | I(0.7) | I(0.8) | verdict |",
              "|---|---|---|---|---|---|---|"]
        for feat, gl in FEATURES.items():
            for g in gl:
                for fam in ("mmd", "ind"):
                    for rule in VERDICT_RULES[fam]:
                        v = verd.get(f"{rule}|n{n}|{g}")
                        if not v:
                            continue
                        i7 = f"{v['I_pp']['0.7']:+.1f} ({v['I_z']['0.7']:+.1f})"
                        i8 = f"{v['I_pp']['0.8']:+.1f} ({v['I_z']['0.8']:+.1f})"
                        flag = " (below G0)" if v["below_baseline_alone"] else ""
                        L.append(f"| {feat} | {g} | {rule} | {v['alone_pp']:+.1f}{flag} | {i7} | {i8} | {v['verdict']} |")
        L.append("")
    tot_cpu = sum(r["cpu_s"] for r in d.values())
    errs = sum(r["errors"] for r in d.values())
    L += [f"Cells: {len(d)}. Errors: {errs}. CPU-hours (row sum): {tot_cpu / 3600:.2f}."]
    open(os.path.join(HERE, f"tables_{a.tag}.md"), "w").write("\n".join(L) + "\n")
    json.dump(verd, open(os.path.join(HERE, f"verdicts_{a.tag}.json"), "w"), indent=1)
    print(f"cells {len(d)}, errors {errs}, cpu_h {tot_cpu / 3600:.2f}")


if __name__ == "__main__":
    main()
