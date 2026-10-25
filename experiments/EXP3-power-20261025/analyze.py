"""EXP3 analysis: raw and size-adjusted power, Wilson and bootstrap intervals, reading rule.

Usage: python -I analyze.py --in main --tag main
  -> summary_<tag>.csv, tables_<tag>.md, results_<tag>.json, null_regression_<tag>.csv
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

import common as C  # noqa: E402
import run_exp3 as R  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXP2_MAIN = os.path.join(HERE, "..", "EXP2-20261023", "main")
ALPHA = 0.05
NBOOT = 2000
SEED = 42

# rule -> (ratio column, raw-rejection function, p-value column or None)
RULES = {
    "ind": {
        "hsic_wild": ("r_hsic_wb", lambda r: float(r["p_hsic_wb"]) <= ALPHA, "p_hsic_wb"),
        "hsic_shift": ("r_hsic_shift", lambda r: float(r["p_hsic_shift"]) <= ALPHA, "p_hsic_shift"),
        "hsic_sub": ("r_hsic_sub", lambda r: int(float(r["rej_hsic_sub"])) == 1, "p_hsic_sub"),
        "dcov_wild": ("r_dcov_wb", lambda r: float(r["p_dcov_wb"]) <= ALPHA, "p_dcov_wb"),
        "dcov_sub": ("r_dcov_sub", lambda r: int(float(r["rej_dcov_sub"])) == 1, "p_dcov_sub"),
        "pearson_nw": ("r_nw", lambda r: int(float(r["rej_nw"])) == 1, "p_nw"),
        "pearson_ewc": ("r_ewc", lambda r: int(float(r["rej_ewc"])) == 1, None),
    },
    "mmd": {
        "mmd_unc_wild": ("r_unc_05", lambda r: float(r["p_unc_05"]) <= ALPHA, "p_unc_05"),
        "mmd_sub_sqrt": ("r_sub_sqrt", lambda r: float(r["t_n"]) > float(r["q95_sqrt"]), None),
        "mmd_sn_block": ("r_sn_05", lambda r: float(r["sn_T"]) > float(r["sn_q95_05"]), None),
    },
}
# bare statistics, size-adjusted directly (oracle reference, no critical-value rule)
ORACLE = {"ind": {"oracle_hsic": lambda r: float(r["hsic_stat"]),
                  "oracle_dcov": lambda r: float(r["dcov_stat"]),
                  "oracle_pearson": lambda r: abs(float(r["pearson_r"]))},
          "mmd": {"oracle_mmd": lambda r: float(r["mmd2"])}}
EXP2_RAW = {"hsic_wild": "p_hsic_wb", "hsic_shift": "p_hsic_shift", "dcov_wild": "p_dcov_wb",
            "mmd_unc_wild": "p_unc_05"}


def load(path):
    rows = list(csv.DictReader(open(path)))
    ok = [r for r in rows if r["status"] == "ok"]
    return ok, len(rows) - len(ok), sum(float(r["cpu_s"]) for r in rows)


def cstar(null_vals: np.ndarray) -> float:
    """Cut-off c* with #(null > c*) = floor(0.05 m) when there are no ties: the k-th order stat."""
    v = np.sort(null_vals)
    return float(v[v.size - int(math.floor(ALPHA * v.size)) - 1])


def adj_power(null_vals, alt_vals):
    c = cstar(null_vals)
    return int(np.sum(alt_vals > c)), alt_vals.size, c, float(np.mean(null_vals > c))


def p_randomized(null_p, alt_p):
    """Size-adjusted power on p-values with a randomized boundary: null rate exactly 5%."""
    v = np.sort(null_p)
    k = int(math.floor(ALPHA * v.size))
    c = v[k - 1]
    less, eq = np.sum(null_p < c), np.sum(null_p == c)
    gam = (k - less) / eq
    return float(np.mean(alt_p < c) + gam * np.mean(alt_p == c)), float(c), int(eq)


def boot_idx(m, rng):
    return rng.integers(0, m, size=(NBOOT, m))


def boot_power(null_vals, alt_vals, bi_null, bi_alt):
    vn = np.sort(null_vals[bi_null], axis=1)
    kk = vn.shape[1] - int(math.floor(ALPHA * vn.shape[1])) - 1
    c = vn[:, kk]
    return np.mean(alt_vals[bi_alt] > c[:, None], axis=1)


def fmt(p, lo, hi):
    return f"{100 * p:.1f} [{100 * lo:.1f}, {100 * hi:.1f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    cells = R.main_cells()
    for c in cells:
        c["name"] = R.name_of(c)
    by = {(c["family"], c["gen"], c["alt"], c["H"], c["n"]): c for c in cells}
    data = {}
    for key, c in by.items():
        p = os.path.join(HERE, a.indir, c["name"] + ".csv")
        if os.path.exists(p):
            data[key] = load(p)
    rng = np.random.default_rng(SEED)
    out_rows, res, cpu_total, err_total = [], {}, 0.0, 0
    for k, (rows, err, cpu) in data.items():
        cpu_total += cpu
        err_total += err
    # ------------------------------------------------ null regression against EXP2
    reg = []
    for (fam, g, alt, h, n), (rows, err, cpu) in data.items():
        if alt != "null":
            continue
        p2 = os.path.join(EXP2_MAIN, f"{fam}__{g}__H{h}__n{n}.csv")
        old = [r for r in csv.DictReader(open(p2)) if r["status"] == "ok"] if os.path.exists(p2) else []
        for rule, (rc, raw, pc) in RULES[fam].items():
            if rule == "pearson_nw":
                rows_ = [r for r in rows if r["rej_nw"] != "-1"]
            else:
                rows_ = rows
            k1 = sum(raw(r) for r in rows_)
            m1 = len(rows_)
            if old:
                old_ = [r for r in old if not (rule == "pearson_nw" and r["rej_nw"] == "-1")]
                k2, m2 = sum(raw(r) for r in old_), len(old_)
                pp = (k1 + k2) / (m1 + m2)
                se = math.sqrt(pp * (1 - pp) * (1 / m1 + 1 / m2)) if 0 < pp < 1 else float("nan")
                z = (k1 / m1 - k2 / m2) / se if se and se == se else float("nan")
            else:
                k2 = m2 = z = float("nan")
            reg.append({"family": fam, "gen": g, "H": h, "n": n, "rule": rule, "exp3_rate": k1 / m1, "exp3_m": m1,
                        "exp2_rate": (k2 / m2) if old else float("nan"), "exp2_m": m2, "z": z})
    with open(os.path.join(HERE, f"null_regression_{a.tag}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(reg[0].keys()) if reg else ["none"])
        w.writeheader()
        w.writerows(reg)
    # ------------------------------------------------ power
    for (fam, g, alt, h, n), (rows, err, cpu) in sorted(data.items(), key=lambda t: str(t[0])):
        if alt == "null":
            continue
        nk = (fam, g, "null", h, n)
        if nk not in data:
            continue
        nrows = data[nk][0]
        bn = boot_idx(len(nrows), rng)
        ba = boot_idx(len(rows), rng)
        cellres = {}
        boots = {}
        allrules = {**{r: v for r, v in RULES[fam].items()}, **ORACLE[fam]}
        for rule, spec in allrules.items():
            if rule.startswith("oracle"):
                fn, raw, pc = spec, None, None
                nv = np.array([fn(r) for r in nrows])
                av = np.array([fn(r) for r in rows])
            else:
                rc, raw, pc = spec
                nv = np.array([float(r[rc]) for r in nrows])
                av = np.array([float(r[rc]) for r in rows])
            # NaN ratios (NW long-run variance <= 0) never reject and are dropped from c*
            nvf = nv[np.isfinite(nv)]
            avc = np.where(np.isfinite(av), av, -np.inf)
            nvc = np.where(np.isfinite(nv), nv, -np.inf)
            k_adj, m_adj, c, null_rate = adj_power(nvf, avc)
            lo, hi = C.wilson(k_adj, m_adj)
            bp = boot_power(nvc, avc, bn, ba)
            boots[rule] = bp
            rec = {"family": fam, "gen": g, "alt": alt, "param": by[(fam, g, alt, h, n)]["param"], "H": h, "n": n,
                   "rule": rule, "adj_k": k_adj, "adj_m": m_adj, "adj_power": k_adj / m_adj,
                   "adj_wilson_lo": lo, "adj_wilson_hi": hi,
                   "adj_boot_lo": float(np.quantile(bp, 0.025)), "adj_boot_hi": float(np.quantile(bp, 0.975)),
                   "cstar_ratio": c, "null_rate_at_cstar": null_rate}
            if raw is not None:
                rr = [r for r in rows if not (rule == "pearson_nw" and r["rej_nw"] == "-1")]
                nr = [r for r in nrows if not (rule == "pearson_nw" and r["rej_nw"] == "-1")]
                kr = sum(raw(r) for r in rr)
                ks = sum(raw(r) for r in nr)
                lo_r, hi_r = C.wilson(kr, len(rr))
                lo_s, hi_s = C.wilson(ks, len(nr))
                rec.update({"raw_k": kr, "raw_m": len(rr), "raw_power": kr / len(rr), "raw_wilson_lo": lo_r,
                            "raw_wilson_hi": hi_r, "null_raw_size": ks / len(nr), "null_size_lo": lo_s,
                            "null_size_hi": hi_s})
                if pc is not None:
                    npv = np.array([float(r[pc]) for r in nr])
                    apv = np.array([float(r[pc]) for r in rr])
                    pw, pcut, ties = p_randomized(npv, apv)
                    rec.update({"adj_power_pvalue": pw, "p_cut": pcut, "p_cut_ties_null": ties})
            rec.update({"alt_reps": len(rows), "null_reps": len(nrows), "alt_errors": err,
                        "alt_cpu_s": round(cpu, 1)})
            out_rows.append(rec)
            cellres[rule] = rec
        # reading rule: shift vs best other HSIC rule
        sh = cellres["hsic_shift"]["adj_power"] if fam == "ind" else None
        if fam == "ind":
            best = max(("hsic_wild", "hsic_sub"), key=lambda r: cellres[r]["adj_power"])
            gap = sh - cellres[best]["adj_power"]
            gb = boots["hsic_shift"] - boots[best]
            glo, ghi = float(np.quantile(gb, 0.025)), float(np.quantile(gb, 0.975))
            if gap >= -0.05:
                v = "competitive"
            elif ghi < 0:
                v = "loses"
            else:
                v = "unclear"
            og = sh - cellres["oracle_hsic"]["adj_power"]
            ob = boots["hsic_shift"] - boots["oracle_hsic"]
            res[f"{g}|{alt}|H{h}|n{n}"] = {
                "best_other_hsic": best, "gap_pp": round(100 * gap, 1),
                "gap_boot_pp": [round(100 * glo, 1), round(100 * ghi, 1)], "verdict": v,
                "shift_minus_oracle_hsic_pp": round(100 * og, 1),
                "shift_minus_oracle_boot_pp": [round(100 * float(np.quantile(ob, 0.025)), 1),
                                               round(100 * float(np.quantile(ob, 0.975)), 1)]}
    fields = []
    for r in out_rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(os.path.join(HERE, f"summary_{a.tag}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)
    with open(os.path.join(HERE, f"results_{a.tag}.json"), "w") as f:
        json.dump({"cells_loaded": len(data), "errors": err_total, "cpu_hours": round(cpu_total / 3600, 2),
                   "nboot": NBOOT, "seed": SEED, "verdicts_shift_vs_best_other_hsic": res}, f, indent=1)
    # ------------------------------------------------ tables
    L = [f"# EXP3 tables ({a.tag})", "",
         "Percent. `adj` = size-adjusted power [bootstrap 95% interval incl. threshold error];",
         "`raw` = nominal-5% power; `size` = raw null rejection on the matched fresh null cell.",
         "Oracle rows size-adjust the bare statistic (no critical-value rule).", ""]
    for (fam, g, alt) in [("ind", "G0", "lin"), ("ind", "G0", "sq"), ("ind", "SEAS-V", "lin"), ("mmd", "G0", "shift")]:
        for n in R.N_GRID:
            sub = [r for r in out_rows if (r["family"], r["gen"], r["alt"], r["n"]) == (fam, g, alt, n)]
            if not sub:
                continue
            L += [f"## {fam} {g} {alt} (param {sub[0]['param']:g}), n = {n}", "",
                  "| rule | " + " | ".join(f"H={h} adj | raw | size" for h in R.H_GRID) + " |",
                  "|---|" + "---|---|---|" * len(R.H_GRID)]
            rules = list(RULES[fam]) + list(ORACLE[fam])
            for rule in rules:
                cells_ = []
                for h in R.H_GRID:
                    rr = [r for r in sub if r["rule"] == rule and r["H"] == h]
                    if not rr:
                        cells_.append("– | – | –")
                        continue
                    r = rr[0]
                    adj = fmt(r["adj_power"], r["adj_boot_lo"], r["adj_boot_hi"])
                    raw = f"{100 * r['raw_power']:.1f}" if "raw_power" in r else "–"
                    sz = f"{100 * r['null_raw_size']:.1f}" if "null_raw_size" in r else "–"
                    cells_.append(f"{adj} | {raw} | {sz}")
                L.append(f"| {rule} | " + " | ".join(cells_) + " |")
            L.append("")
    L += ["## Reading rule: HSIC shift vs best other HSIC rule (size-adjusted)", "",
          "| cell | best other | gap pp [boot 95%] | verdict | shift − oracle HSIC pp [boot] |", "|---|---|---|---|---|"]
    for k, v in res.items():
        L.append(f"| {k} | {v['best_other_hsic']} | {v['gap_pp']} [{v['gap_boot_pp'][0]}, {v['gap_boot_pp'][1]}] | "
                 f"{v['verdict']} | {v['shift_minus_oracle_hsic_pp']} [{v['shift_minus_oracle_boot_pp'][0]}, "
                 f"{v['shift_minus_oracle_boot_pp'][1]}] |")
    L += ["", "## Null regression: fresh EXP3 null raw size vs EXP2 main (same generator, rule)", "",
          "| cell | rule | EXP3 % | EXP2 % | z |", "|---|---|---|---|---|"]
    for r in reg:
        L.append(f"| {r['family']} {r['gen']} H{r['H']} n{r['n']} | {r['rule']} | {100 * r['exp3_rate']:.1f} | "
                 f"{100 * r['exp2_rate']:.1f} | {r['z']:.2f} |")
    L += ["", f"Cells loaded: {len(data)}. Errors: {err_total}. CPU: {cpu_total / 3600:.2f} h."]
    with open(os.path.join(HERE, f"tables_{a.tag}.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"cells {len(data)}, errors {err_total}, cpu {cpu_total / 3600:.2f} h")


if __name__ == "__main__":
    main()
