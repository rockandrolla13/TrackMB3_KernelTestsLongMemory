"""EXP4 analysis: rejection rates by generator x adjustment, Wilson intervals, reading rule (DESIGN.md).

Usage: python -I analyze_exp4.py --in main --tag main  ->  summary_<tag>.csv, tables_<tag>.md, verdicts_<tag>.json
Unadjusted n = 2000 rates come from ../EXP2-20261023/summary_main.csv (row adj = 'A0-EXP2').
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))  # kmem (D7: not installed)
from kmem.legacy import common as C  # noqa: E402
from kmem.analysis import require_input  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXP2 = os.path.join(HERE, "..", "EXP2-20261023", "summary_main.csv")
ALPHA = 0.05
GENS = ("G0", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow")
H_GRID = (0.5, 0.7, 0.8)
N_GRID = (500, 2000)
ADJ = ("A0", "A0-EXP2", "A1", "A2", "A3")
DIAGNOSTIC = ("SEAS-M-shared", "SEAS-V-slow")  # DESIGN.md rule 5 and amendment 3: not used in rules 1-3
RULES = {"hsic_wild": lambda r: float(r["p_hsic_wb"]) <= ALPHA,
         "hsic_shift": lambda r: float(r["p_hsic_shift"]) <= ALPHA,
         "hsic_sub": lambda r: int(float(r["rej_hsic_sub"])) == 1,
         "dcov_wild": lambda r: float(r["p_dcov_wb"]) <= ALPHA,
         "dcov_sub": lambda r: int(float(r["rej_dcov_sub"])) == 1,
         "pearson_nw": lambda r: int(float(r["rej_nw"])) == 1,
         "pearson_ewc": lambda r: int(float(r["rej_ewc"])) == 1}


def rate(k, m):
    lo, hi = C.wilson(k, m)
    return {"k": k, "m": m, "p": k / m if m else float("nan"), "lo": lo, "hi": hi}


def load(indir):
    d, flags, period, cpu = {}, {}, {}, defaultdict(float)
    for g in GENS:
        for h in H_GRID:
            for n in N_GRID:
                p = os.path.join(HERE, indir, f"ind__{g}__H{h}__n{n}.csv")
                if not os.path.exists(p):
                    continue
                ok, failed = {}, defaultdict(set)
                for r in csv.DictReader(open(p)):
                    if None in r.values():  # truncated row from a crash mid-write
                        continue
                    cpu[(g, h, n)] += float(r["cpu_s"] or 0)
                    if r["status"] != "ok":
                        failed[r["adj"]].add(r["rep"])
                        continue
                    ok.setdefault((r["adj"], r["rep"]), r)  # resume never writes a second ok row
                by = defaultdict(list)
                for (adj, _), r in ok.items():
                    by[adj].append(r)
                for a, rows in by.items():
                    # an error that was retried successfully is not an error
                    rec = {"errors": len(failed[a] - {r["rep"] for r in rows})}
                    for rule, f in RULES.items():
                        v = [f(r) for r in rows if not (rule == "pearson_nw" and r["rej_nw"] == "-1")]
                        rec[rule] = rate(int(sum(v)), len(v))
                    d[(g, h, n, a)] = rec
                    flags[(g, h, n, a)] = {r["rep"]: {k: RULES[k](r) for k in RULES} for r in rows}
                    if a == "A3":
                        px = [int(r["px"]) for r in rows] + [int(r["py"]) for r in rows]
                        src = [r["srcx"] for r in rows] + [r["srcy"] for r in rows]
                        hit = sum(1 for q in px if q == 50)
                        period[(g, h, n)] = {"series": len(px), "recovered": hit,
                                             "both_recovered": sum(1 for r in rows if r["px"] == r["py"] == "50"),
                                             "pairs": len(rows),
                                             "src_square": sum(1 for s in src if s == "square"),
                                             "median_P": sorted(px)[len(px) // 2]}
    # a declared input (spec inputs): missing is an error, not a silent skip
    for r in csv.DictReader(open(require_input(EXP2, "EXP2 main summary, unadjusted A0-EXP2 rows"))):
        if r["family"] != "ind" or r["gen"] not in GENS:
            continue
        key = (r["gen"], float(r["H"]), int(r["n"]), "A0-EXP2")
        d.setdefault(key, {"errors": int(r["errors"])})[r["rule"]] = rate(int(r["rejections"]), int(r["valid_reps"]))
    return d, flags, period, cpu


def fmt(c):
    return f"{100 * c['p']:.1f} [{100 * c['lo']:.1f}, {100 * c['hi']:.1f}]"


def z2(a, b):
    pp = (a["k"] + b["k"]) / (a["m"] + b["m"])
    se = math.sqrt(pp * (1 - pp) * (1 / a["m"] + 1 / b["m"]))
    return (a["p"] - b["p"]) / se if se > 0 else 0.0


def write_summary(d: dict, path: str) -> None:
    """One row per (generator, H, n, adjustment, rule): rejections, valid replications, rate, Wilson interval."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gen", "H", "n", "adj", "rule", "rejections", "valid_reps", "errors", "rate", "wilson_lo", "wilson_hi"])
        for (g, h, n, adj), rec in sorted(d.items()):
            for rule in RULES:
                if rule in rec:
                    c = rec[rule]
                    w.writerow([g, h, n, adj, rule, c["k"], c["m"], rec["errors"], round(c["p"], 5),
                                round(c["lo"], 5), round(c["hi"], 5)])


def report(tag: str, d: dict, flags: dict, period: dict, cpu: dict) -> tuple[str, dict, int, float]:
    """The tables, with the reading-rule verdicts computed alongside: (markdown, verdicts, errors, CPU seconds)."""
    L = [f"# EXP4 tables ({tag})", "",
         "Rejection % at nominal 5% [Wilson 95%]. A0-EXP2 = unadjusted EXP2 main run (independent seeds).", ""]
    for n in N_GRID:
        for h in H_GRID:
            rows = [(g, adj) for g in GENS for adj in ADJ if (g, h, n, adj) in d]
            if not rows:
                continue
            L += [f"## n = {n}, H = {h}", "", "| generator | adj | " + " | ".join(RULES) + " | reps |",
                  "|---|---|" + "---|" * len(RULES) + "---|"]
            for g, adj in rows:
                rec = d[(g, h, n, adj)]
                L.append(f"| {g} | {adj} | " + " | ".join(fmt(rec[r]) if r in rec else "-" for r in RULES)
                         + f" | {rec['hsic_wild']['m']} |")
            L.append("")
    L += ["# A3 period recovery", "", "| generator | H | n | series with P̂ = 50 | pairs with both = 50 | won by x² | median P̂ |",
          "|---|---|---|---|---|---|---|"]
    for (g, h, n), p in sorted(period.items()):
        L.append(f"| {g} | {h} | {n} | {p['recovered']}/{p['series']} | {p['both_recovered']}/{p['pairs']} "
                 f"| {p['src_square']}/{p['series']} | {p['median_P']} |")
    L.append("")
    verd = {}
    L += ["# Reading rule 1-2: does adjustment k restore the adjusted-G0 rate?", "",
          "diff = rate(g, A_k) − rate(G0, A_k) in pp. 'fixed' if diff ≤ 2 pp or intervals overlap.", "",
          "| n | H | generator | adj | " + " | ".join(RULES) + " |", "|---|---|---|---|" + "---|" * len(RULES)]
    for n in N_GRID:
        for h in H_GRID:
            for g in GENS[1:]:
                for adj in ("A0", "A0-EXP2", "A1", "A2", "A3"):
                    f_, b_ = d.get((g, h, n, adj)), d.get(("G0", h, n, adj))
                    base = adj
                    if f_ and not b_ and adj == "A0":  # e.g. SEAS-V-slow at n = 2000: G0 A0 exists only in EXP2
                        b_, base = d.get(("G0", h, n, "A0-EXP2")), "A0-EXP2"
                    if not f_ or not b_:
                        continue
                    cells = []
                    for r in RULES:
                        x, y = f_[r], b_[r]
                        diff = 100 * (x["p"] - y["p"])
                        ok = diff <= 2 or x["lo"] <= y["hi"]
                        verd[f"{g}|{h}|{n}|{adj}|{r}"] = {"diff_pp": round(diff, 2), "fixed": ok,
                                                          "diagnostic": g in DIAGNOSTIC, "baseline": base}
                        cells.append(f"{diff:+.1f} {'fixed' if ok else 'NOT'}")
                    label = g + (" (diagnostic)" if g in DIAGNOSTIC else "")
                    vs = "" if base == adj else f" vs G0 {base}"
                    L.append(f"| {n} | {h} | {label} | {adj}{vs} | " + " | ".join(cells) + " |")
    L.append("")
    L += ["# Reading rule 3: does adjustment change G0 size?", "",
          "n = 500: paired against A0 on the same draws (discordant counts b = adj-only, c = A0-only).",
          "n = 2000: against EXP2 G0 (independent samples), two-sample z.", "",
          "| n | H | adj | rule | G0 A_k | G0 A0 | diff pp | z or (b, c) | flag |", "|---|---|---|---|---|---|---|---|---|"]
    for n in N_GRID:
        base = "A0" if n == 500 else "A0-EXP2"
        for h in H_GRID:
            for adj in ("A1", "A2", "A3"):
                f_, b_ = d.get(("G0", h, n, adj)), d.get(("G0", h, n, base))
                if not f_ or not b_:
                    continue
                for r in RULES:
                    x, y = f_[r], b_[r]
                    diff = 100 * (x["p"] - y["p"])
                    if n == 500 and (("G0", h, n, adj) in flags):
                        fa, fb = flags[("G0", h, n, adj)], flags[("G0", h, n, "A0")]
                        common = set(fa) & set(fb)
                        bb = sum(1 for k in common if fa[k][r] and not fb[k][r])
                        cc = sum(1 for k in common if fb[k][r] and not fa[k][r])
                        zz = (bb - cc) / math.sqrt(bb + cc) if bb + cc else 0.0
                        stat = f"({bb}, {cc}) z={zz:+.1f}"
                    else:
                        zz = z2(x, y)
                        stat = f"z={zz:+.1f}"
                    flag = "FLAG" if abs(diff) > 2 and abs(zz) > 1.96 else ""
                    verd[f"size|G0|{h}|{n}|{adj}|{r}"] = {"diff_pp": round(diff, 2), "z": round(zz, 2), "flag": bool(flag)}
                    L.append(f"| {n} | {h} | {adj} | {r} | {100 * x['p']:.1f} | {100 * y['p']:.1f} | {diff:+.1f} | {stat} | {flag} |")
    L.append("")
    L += ["# Reading rule 6: suspicion check", "",
          "Flag: adjusted G0 more than 2 pp below its unadjusted rate, or an adjusted seasonal cell",
          "more than 2 pp below the adjusted G0 rate.", "",
          "| n | H | generator | adj | rule | rate | reference | diff pp |", "|---|---|---|---|---|---|---|---|"]
    nflag = 0
    for n in N_GRID:
        base0 = "A0" if n == 500 else "A0-EXP2"
        for h in H_GRID:
            for g in GENS:
                for adj in ("A1", "A2", "A3"):
                    f_ = d.get((g, h, n, adj))
                    ref = d.get(("G0", h, n, base0 if g == "G0" else adj))
                    if not f_ or not ref:
                        continue
                    for r in RULES:
                        diff = 100 * (f_[r]["p"] - ref[r]["p"])
                        if diff < -2:
                            nflag += 1
                            refname = f"G0 {base0 if g == 'G0' else adj}"
                            verd[f"suspicion|{g}|{h}|{n}|{adj}|{r}"] = {"diff_pp": round(diff, 2), "reference": refname}
                            L.append(f"| {n} | {h} | {g} | {adj} | {r} | {100 * f_[r]['p']:.1f} | {refname} | {diff:+.1f} |")
    L += ["", f"Suspicion flags: {nflag}.", ""]
    tot = sum(cpu.values())
    errs = sum(rec["errors"] for k, rec in d.items() if k[3] != "A0-EXP2")
    L += [f"Cells: {len(cpu)}. Unrecovered (replication, adjustment) errors: {errs}. CPU-hours (row sum): {tot / 3600:.2f}."]
    return "\n".join(L) + "\n", verd, errs, tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    d, flags, period, cpu = load(a.indir)
    write_summary(d, os.path.join(HERE, f"summary_{a.tag}.csv"))
    text, verd, errs, tot = report(a.tag, d, flags, period, cpu)
    open(os.path.join(HERE, f"tables_{a.tag}.md"), "w").write(text)
    json.dump({"verdicts": verd, "period": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in period.items()},
               "cpu_h": tot / 3600}, open(os.path.join(HERE, f"verdicts_{a.tag}.json"), "w"), indent=1)
    print(f"cells {len(cpu)}, error rows {errs}, cpu_h {tot / 3600:.2f}")


if __name__ == "__main__":
    main()
