"""EXP6 analysis: raw size, size-adjusted power per kappa, paired kappa gaps, mechanism diagnostics.

Usage: python -I analyze_exp6.py --in main --tag main
  -> summary_<tag>.csv, tables_<tag>.md, results_<tag>.json
Size adjustment, intervals and c* follow EXP3 analyze.py exactly (DESIGN.md).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import numpy as np  # noqa: E402

import bw_core as K  # noqa: E402
import run_exp6 as R  # noqa: E402
from kmem.analysis import boot_power, cstar, mcnemar_p, size_class, wilson  # noqa: E402,F401

ALPHA = 0.05
NBOOT = 2000
SEED = 42

# rule -> (ratio column, raw-rejection column test); HSIC rules get the kappa suffix
HSIC_RULES = {
    "hsic_wild": ("r_hsic_wb", lambda r, s: float(r["p_hsic_wb" + s]) <= ALPHA),
    "hsic_shift": ("r_hsic_shift", lambda r, s: float(r["p_hsic_shift" + s]) <= ALPHA),
    "hsic_sub": ("r_hsic_sub", lambda r, s: int(float(r["rej_hsic_sub" + s])) == 1),
    "oracle_hsic": ("hsic_stat", None),
}
REF_RULES = {
    "dcov_wild": ("r_dcov_wb", lambda r, s: float(r["p_dcov_wb"]) <= ALPHA),
    "dcov_sub": ("r_dcov_sub", lambda r, s: int(float(r["rej_dcov_sub"])) == 1),
    "pearson_nw": ("r_nw", lambda r, s: int(float(r["rej_nw"])) == 1),
    "pearson_ewc": ("r_ewc", lambda r, s: int(float(r["rej_ewc"])) == 1),
    "oracle_dcov": ("dcov_stat", None),
    "oracle_pearson": ("abs_r", None),
}
NULL_OF = {("G0", "lin"): ("G0", "null"), ("G0", "sq"): ("G0", "null"), ("SEAS-V", "lin"): ("SEAS-V", "null")}


def load(path):
    rows = list(csv.DictReader(open(path)))
    ok = [r for r in rows if r["status"] == "ok"]
    for r in ok:
        r["abs_r"] = abs(float(r["pearson_r"]))
    return ok, len(rows) - len(ok), sum(float(r["cpu_s"]) for r in rows)


def col(rows, c):
    return np.array([float(r[c]) for r in rows])


def rule_list(kappas, okap=()):
    out = [("oracle_hsic", k, "hsic_stat_" + K.ktag(k), None) for k in okap]
    for rule, (rc, raw) in HSIC_RULES.items():
        for k in kappas:
            s = "_" + K.ktag(k)
            out.append((rule, k, rc + s, (lambda r, raw=raw, s=s: raw(r, s)) if raw else None))
    for rule, (rc, raw) in REF_RULES.items():
        out.append((rule, None, rc, (lambda r, raw=raw: raw(r, "")) if raw else None))
    return out


def fmt(p, lo=None, hi=None):
    if lo is None:
        return f"{100 * p:.1f}"
    return f"{100 * p:.1f} [{100 * lo:.1f}, {100 * hi:.1f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--n", default="2000,500")
    a = ap.parse_args()
    n_grid = tuple(int(v) for v in a.n.split(","))
    cells = R.cells_for("main", n_grid)
    data, kappas, okap = {}, None, None
    for c in cells:
        p = os.path.join(HERE, a.indir, c["name"] + ".csv")
        if os.path.exists(p):
            data[(c["gen"], c["alt"], c["H"], c["n"])] = (c, *load(p))
            hdr = open(p).readline().strip().split(",")
            kv = lambda h: float(h.rsplit("_k", 1)[1].replace("p", "."))  # noqa: E731
            ks = sorted({kv(h) for h in hdr if h.startswith("r_hsic_wb_k")})
            ok = sorted({kv(h) for h in hdr if h.startswith("hsic_stat_k")} - set(ks))
            kappas = ks if kappas is None else [k for k in kappas if k in ks]
            okap = ok if okap is None else [k for k in okap if k in ok]
    # stop rule (DESIGN.md A1): keep replication indices 0..m-1 complete in every cell
    def prefix(rows):
        have = {int(r["rep"]) for r in rows}
        m = 0
        while m in have:
            m += 1
        return m
    m_common = min(prefix(v[1]) for v in data.values())
    for key, (c, rows, err, cpu_) in list(data.items()):
        data[key] = (c, [r for r in rows if int(r["rep"]) < m_common], err, cpu_)
    print(f"common replication prefix m = {m_common}")
    rules = rule_list(kappas, okap)
    allk = sorted(kappas + okap)
    rng = np.random.default_rng(SEED)
    summary, res = [], {"kappas": kappas, "size": {}, "power": {}, "gaps": {}, "mechanism": {}}
    cpu = sum(v[3] for v in data.values())
    errs = sum(v[2] for v in data.values())
    # ------------------------------------------------ size
    for (g, alt, h, n), (c, rows, err, _) in data.items():
        if alt != "null":
            continue
        m = len(rows)
        key = f"{g}__H{h}__n{n}"
        res["size"][key] = {}
        rej = {}
        for rule, k, rc, raw in rules:
            if raw is not None:
                rr = rows if rule != "pearson_nw" else [r for r in rows if r["rej_nw"] != "-1"]
                rej[(rule, k)] = (rr, np.array([raw(r) for r in rr], dtype=bool))
        for rule, k, rc, raw in rules:
            if raw is None:
                continue
            rr, v = rej[(rule, k)]
            kk = int(v.sum())
            lo, hi = wilson(kk, len(rr))
            ent = {"rate": kk / len(rr), "lo": lo, "hi": hi, "m": len(rr), "class": size_class(kk, len(rr))}
            if k is not None and k != 1.0:
                ent["mcnemar_p_vs_k1"], ent["n_only_this"], ent["n_only_k1"] = mcnemar_p(v, rej[(rule, 1.0)][1])
            res["size"][key][f"{rule}|{k}"] = ent
            summary.append({"cell": f"{g}__null__H{h}__n{n}", "rule": rule, "kappa": k, "kind": "size",
                            "value": ent["rate"], "lo": lo, "hi": hi, "m": len(rr),
                            "note": ent["class"] + (f"; McNemar p vs k1 = {ent['mcnemar_p_vs_k1']:.3g}" if "mcnemar_p_vs_k1" in ent else "")})
    # ------------------------------------------------ size-adjusted power, paired kappa gaps
    for (g, alt, h, n), (c, rows, err, _) in data.items():
        if alt == "null":
            continue
        ng, nalt = NULL_OF[(g, alt)]
        if (ng, nalt, h, n) not in data:
            continue
        nrows = data[(ng, nalt, h, n)][1]
        key = f"{g}__{alt}{c['param']:g}__H{h}__n{n}"
        bi_n = rng.integers(0, len(nrows), size=(NBOOT, len(nrows)))
        bi_a = rng.integers(0, len(rows), size=(NBOOT, len(rows)))
        res["power"][key], boots = {}, {}
        for rule, k, rc, raw in rules:
            nv, av = col(nrows, rc), col(rows, rc)
            okn, oka = np.isfinite(nv), np.isfinite(av)
            nv, av = nv[okn], av[oka]
            cs = cstar(nv)
            kk = int(np.sum(av > cs))
            lo, hi = wilson(kk, av.size)
            b = boot_power(nv, av, bi_n[:, :nv.size] % nv.size, bi_a[:, :av.size] % av.size)
            boots[(rule, k)] = b
            raw_p = float(np.mean([raw(r) for r in rows])) if raw else float("nan")
            ent = {"adj": kk / av.size, "lo": lo, "hi": hi, "boot_lo": float(np.quantile(b, 0.025)),
                   "boot_hi": float(np.quantile(b, 0.975)), "raw": raw_p, "cstar": cs, "m_alt": int(av.size),
                   "m_null": int(nv.size)}
            res["power"][key][f"{rule}|{k}"] = ent
            summary.append({"cell": key, "rule": rule, "kappa": k, "kind": "adj_power", "value": ent["adj"],
                            "lo": ent["boot_lo"], "hi": ent["boot_hi"], "m": av.size,
                            "note": f"raw {raw_p:.3f}; Wilson [{lo:.3f}, {hi:.3f}]; c* {cs:.4g}"})
        res["gaps"][key] = {}
        best_ref = max(res["power"][key]["dcov_wild|None"]["adj"], 0.0)
        for rule in HSIC_RULES:
            for k in (allk if rule == "oracle_hsic" else kappas):
                if k == 1.0:
                    continue
                d = boots[(rule, k)] - boots[(rule, 1.0)]
                gain = res["power"][key][f"{rule}|{k}"]["adj"] - res["power"][key][f"{rule}|1.0"]["adj"]
                lo, hi = float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))
                verdict = ("gain" if lo > 0 and gain > 0.10 else "loss" if hi < 0 and gain < -0.10 else
                           "small/unclear")
                if verdict == "gain" and res["power"][key][f"{rule}|{k}"]["adj"] >= best_ref - 0.05:
                    verdict = "full recovery (within 5 pp of dCov wild)"
                res["gaps"][key][f"{rule}|{k}"] = {"gain": gain, "lo": lo, "hi": hi, "verdict": verdict}
                summary.append({"cell": key, "rule": rule, "kappa": k, "kind": "gain_vs_k1", "value": gain,
                                "lo": lo, "hi": hi, "m": len(rows), "note": verdict})
    # ------------------------------------------------ mechanism diagnostics (null cells)
    for (g, alt, h, n), (c, rows, err, _) in data.items():
        key = f"{g}__{alt}__H{h}__n{n}"
        r2 = col(rows, "pearson_r") ** 2
        cs2 = col(rows, "corr_sq") ** 2
        ent = {"median_sx": float(np.median(col(rows, "sx_med"))), "mean_corr_sq2": float(np.mean(cs2)),
               "mean_r2": float(np.mean(r2))}
        for k in allk:
            t = col(rows, f"hsic_stat_{K.ktag(k)}")
            rk = lambda v: np.argsort(np.argsort(v))  # noqa: E731
            ent[f"spearman_hsic_corrsq2|{k}"] = float(np.corrcoef(rk(t), rk(cs2))[0, 1])
            ent[f"spearman_hsic_r2|{k}"] = float(np.corrcoef(rk(t), rk(r2))[0, 1])
            # least-squares share of the HSIC statistic explained by n r^2 and n corr_sq^2
            X = np.column_stack([np.ones_like(t), n * r2, n * cs2])
            beta, *_ = np.linalg.lstsq(X, t, rcond=None)
            fit = X @ beta
            ent[f"r2_fit|{k}"] = float(1 - np.sum((t - fit) ** 2) / np.sum((t - t.mean()) ** 2))
            ent[f"mean_hsic|{k}"] = float(np.mean(t))
        res["mechanism"][key] = ent
    res["cpu_hours"] = cpu / 3600
    res["errors"] = errs
    res["m_common"] = m_common
    res["reps"] = {f"{k[0]}__{k[1]}__H{k[2]}__n{k[3]}": len(v[1]) for k, v in data.items()}
    with open(os.path.join(HERE, f"summary_{a.tag}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    with open(os.path.join(HERE, f"results_{a.tag}.json"), "w") as f:
        json.dump(res, f, indent=1)
    res["oracle_only_kappas"] = okap
    write_tables(res, kappas, allk, a.tag)
    print(f"cells {len(data)}, errors {errs}, CPU-h {cpu / 3600:.2f}")


def write_tables(res, kappas, allk, tag):
    L = [f"# EXP6 tables ({tag})", "", f"κ grid: {kappas}. Rates in %. CPU-hours: {res['cpu_hours']:.2f}. "
         f"Errors: {res['errors']}. Reps: {res['reps']}", ""]
    L += ["## Raw null rejection rate (size), Wilson 95%", ""]
    for cell, d in res["size"].items():
        L += [f"### {cell}", "", "| rule | " + " | ".join(f"κ={k:g}" for k in kappas) + " |",
              "|---|" + "---|" * len(kappas)]
        for rule in ("hsic_wild", "hsic_shift", "hsic_sub"):
            cells = []
            for k in kappas:
                e = d[f"{rule}|{k}"]
                s = f"{fmt(e['rate'], e['lo'], e['hi'])} {e['class']}"
                if "mcnemar_p_vs_k1" in e:
                    s += f" (McN p={e['mcnemar_p_vs_k1']:.2g})"
                cells.append(s)
            L.append(f"| {rule} | " + " | ".join(cells) + " |")
        refs = ", ".join(f"{r} {fmt(d[f'{r}|None']['rate'])}" for r in ("dcov_wild", "dcov_sub", "pearson_nw", "pearson_ewc"))
        L += ["", f"κ-free references: {refs}.", ""]
    L += ["## Size-adjusted power, bootstrap 95% (null and alternative resampled)", "",
          f"Oracle-only κ (bare statistic, no rule): {res['oracle_only_kappas']}.", ""]
    for cell, d in res["power"].items():
        L += [f"### {cell}", "", "| rule | " + " | ".join(f"κ={k:g}" for k in allk) + " |",
              "|---|" + "---|" * len(allk)]
        for rule in HSIC_RULES:
            L.append(f"| {rule} | " + " | ".join(
                fmt(d[f'{rule}|{k}']['adj'], d[f'{rule}|{k}']['boot_lo'], d[f'{rule}|{k}']['boot_hi'])
                if f'{rule}|{k}' in d else "–" for k in allk) + " |")
        refs = ", ".join(f"{r} {fmt(d[f'{r}|None']['adj'])}" for r in REF_RULES)
        L += ["", f"κ-free references: {refs}.", "", "Paired gain vs κ = 1 (pp, bootstrap 95%):", ""]
        for rk, g in res["gaps"][cell].items():
            L.append(f"- {rk}: {100 * g['gain']:+.1f} [{100 * g['lo']:+.1f}, {100 * g['hi']:+.1f}] {g['verdict']}")
        L.append("")
    L += ["## Mechanism diagnostics", "",
          "Spearman correlation of the HSIC statistic with corr(X²,Y²)² and with r², and R² of a linear fit "
          "HSIC ~ n r² + n corr(X²,Y²)², per cell.", "",
          "| cell | median σ_X | " + " | ".join(f"ρS(·,c²) κ={k:g} | ρS(·,r²) κ={k:g} | R² κ={k:g}" for k in allk) + " |",
          "|---|---|" + "---|" * (3 * len(allk))]
    for cell, e in res["mechanism"].items():
        L.append(f"| {cell} | {e['median_sx']:.3f} | " + " | ".join(
            f"{e[f'spearman_hsic_corrsq2|{k}']:.2f} | {e[f'spearman_hsic_r2|{k}']:.2f} | {e[f'r2_fit|{k}']:.2f}"
            for k in allk) + " |")
    open(os.path.join(HERE, f"tables_{tag}.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
