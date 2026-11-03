"""EXP7 compact tables for RESULTS.md, read only from the saved summary and results files.

Usage: python -I results_tables.py <tag>  -> results_tables_<tag>.md. No number here is computed from replication rows.
"""
from __future__ import annotations

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

KERNEL_ORDER = ("raw", "inc1", "inc50", "herm1", "herm2", "sym", "datadiff", "hp10", "dcov_inc1", "frac_oracle", "hsplit1", "hsplit2", "dcov", "pearson")
SIZE_RULES = [("raw", r) for r in ("wild", "shift", "sub", "harris", "harris19", "harrisA")] + \
    [("inc1", r) for r in ("wild", "shift", "sub", "harris", "harris19", "harrisA")] + \
    [(k, r) for k in ("inc50", "herm1", "herm2", "sym", "datadiff", "hp10", "dcov_inc1", "frac_oracle") for r in ("wild", "shift")] + \
    [(k, r) for k in ("hsplit1", "hsplit2") for r in ("shift", "wild", "shift_blk", "wild_blk")] + \
    [("dcov", "wild"), ("dcov", "sub"), ("pearson", "nw"), ("pearson", "ewc")]
ARMS = ("HT-S", "HT-L", "NL-SV", "NL-TAR", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow")
POWER_RULES = [("raw", "shift"), ("raw", "wild"), ("raw", "harris"), ("raw", "oracle"), ("inc1", "shift"), ("inc1", "wild"), ("inc1", "oracle"),
               ("inc50", "shift"), ("herm1", "shift"), ("herm2", "shift"), ("sym", "shift"), ("datadiff", "shift"), ("hp10", "shift"), ("dcov_inc1", "shift"),
               ("frac_oracle", "shift"), ("hsplit1", "shift"), ("hsplit2", "shift"), ("hsplit1", "shift_blk"), ("hsplit2", "shift_blk"), ("dcov", "wild"), ("pearson", "nw")]
SHORT = {"loses": "L", "gains": "G", "no clear change": "=", "above Harris": "+", "below Harris": "−", "recovers": "R", "full recovery": "F", "does not recover": "n"}
CLASS = {"valid": "", "mild": "*", "over": "**"}


def pc(v: str | float) -> str:
    return f"{100 * float(v):.1f}"


def main(tag: str) -> None:
    rows = list(csv.DictReader(open(os.path.join(HERE, f"summary_{tag}.csv"))))
    res = json.load(open(os.path.join(HERE, f"results_{tag}.json")))
    L = [f"# EXP7 compact tables ({tag})", "", "Percent. Marks on rates: no mark valid, * mildly oversized, ** oversized.", ""]
    aud = res["audit"]
    L += ["## Audit", "", f"- Error rows: {aud['errors']}. CPU-hours: {aud['cpu_hours']:.1f}. Largest Cholesky rank: {aud['max_rank']}.",
          f"- Fast path against exact path: {aud['reps']} independence replications, largest relative T error {aud['max_rel_t']:.1e}, decisions that differ: wild {aud['wb_flip']}, "
          f"shift {aud['shift_flip']}, sub {aud['sub_flip']}.",
          f"- Split MMD: {aud['mmd_reps']} replications, largest relative error {aud['mmd_max_rel']:.1e}, decisions that differ {aud['mmd_flip']}.",
          f"- Largest |S(0) − T|: {aud['max_abs_check']:.1e}. Stop rule triggered: {aud['stop_rule_triggered']}.",
          f"- Harris at N = 19, whole-series against segment bandwidths: {aud['harris_reps']} replications, decisions that differ: raw {aud['harris_flip_raw']}, inc1 {aud['harris_flip_inc1']}, "
          f"largest relative difference in the ratio {aud['harris_max_rel']:.1e}.", ""]
    # verdicts
    L += ["## Verdicts on the 24 Gaussian memory cells", "", "| kernel, rule | valid | mild | over | removes memory | memory-flat | rises n=500→4000 | above H=0.5 baseline | worst rate (HX, HY, n) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for k, r in SIZE_RULES:
        v = res["verdicts_memory"].get(f"{k}|{r}")
        if v:
            L.append(f"| {k}, {r} | {v['n_valid']}/{v['cells']} | {v['n_mild']} | {v['n_over']} | {'yes' if v['removes_memory'] else 'no'} | {'yes' if v['memory_flat'] else 'no'} | "
                     f"{len(v['rises'])} | {len(v['above_baseline'])} | {pc(v['worst_rate'])} {tuple(v['worst_cell'])} |")
    # size at G0
    size = {(x["family"], x["gen"], x["H"], x["n"], x["kernel"], x["rule"]): x for x in rows if x["kind"] == "size"}
    for n in ("500", "2000", "4000"):
        hs = sorted({key[2] for key in size if key[0] == "ind" and key[1] == "G0" and key[3] == n}, key=lambda s: (("/" in s), s))
        if not hs:
            continue
        L += ["", f"## Size, Gaussian memory cells, n = {n}", "", "| kernel, rule | " + " | ".join(hs) + " |", "|---|" + "---|" * len(hs)]
        for k, r in SIZE_RULES:
            cells = [size.get(("ind", "G0", h, n, k, r)) for h in hs]
            if any(cells):
                L.append(f"| {k}, {r} | " + " | ".join(f"{pc(c['value'])}{CLASS[c['note']]}" if c else "–" for c in cells) + " |")
    # process arms
    for h in ("0.5", "0.8"):
        L += ["", f"## Size, process arms, n = 2000, H = {h}", "", "| kernel, rule | " + " | ".join(ARMS) + " |", "|---|" + "---|" * len(ARMS)]
        for k, r in SIZE_RULES:
            cells = [size.get(("ind", g, h, "2000", k, r)) for g in ARMS]
            if any(cells):
                L.append(f"| {k}, {r} | " + " | ".join(f"{pc(c['value'])}{CLASS[c['note']]}" if c else "–" for c in cells) + " |")
    # low-layer part
    L += ["", "## Low-layer part of hsplit alone, nominal 2.5%, n = 2000", "", "| cell | hsplit1 plug-in | hsplit1 block | hsplit2 plug-in | hsplit2 block |", "|---|---|---|---|---|"]
    for g, hs in (("G0", ("0.5", "0.6", "0.7", "0.8", "0.9", "0.95", "0.6/0.9", "0.7/0.9")), *((a, ("0.5", "0.8")) for a in ARMS)):
        for h in hs:
            vals = [size.get(("ind", g, h, "2000", f"hsplit{r}", rule)) for r in (1, 2) for rule in ("low", "low_blk")]
            if all(vals):
                L.append(f"| {g}, H = {h} | " + " | ".join(pc(c["value"]) for c in vals) + " |")
    # growth
    L += ["", "## Bare statistic: 95th percentile of T", ""]
    g = res["growth"]
    if g["q95"]:
        ns = sorted({int(k.split("|")[2]) for k in g["q95"]})
        hs = sorted({float(k.split("|")[1]) for k in g["q95"]})
        L += ["| kernel, H | " + " | ".join(f"n={n}" for n in ns) + " | ratio n=8000/n=500 | slope [95%] | reference | flat |", "|---|" + "---|" * (len(ns) + 4)]
        for kern in ("raw", "inc1"):
            for h in hs:
                s = g["slope"][f"{kern}|{h}"]
                vals = " | ".join(f"{g['q95'][f'{kern}|{h}|{n}']['q95']:.3g}" for n in ns)
                L.append(f"| {kern}, {h} | {vals} | {s['ratio_last_first']:.2f} | {s['slope']:.3f} [{s['lo']:.3f}, {s['hi']:.3f}] | {s['reference']:.2f} | {'yes' if s['flat'] else 'no'} |")
    # power
    adj = {(x["family"], x["alt"], x["param"], x["H"], x["kernel"], x["rule"]): x for x in rows if x["kind"] == "adj_power"}
    raw = {(x["family"], x["alt"], x["param"], x["H"], x["kernel"], x["rule"]): x for x in rows if x["kind"] == "raw_power"}
    gap = {(x["family"], x["alt"], x["param"], x["H"], x["kernel"], x["rule"]): x for x in rows if x["kind"] == "gap_vs_raw"}
    ghar = {(x["family"], x["alt"], x["param"], x["H"], x["kernel"], x["rule"]): x for x in rows if x["kind"] == "gap_vs_harris"}
    gher = {(x["family"], x["alt"], x["param"], x["H"], x["kernel"], x["rule"]): x for x in rows if x["kind"] == "gap_vs_herm"}
    L += ["", "## Independence power at n = 2000", "", "Each entry: size-adjusted power / raw power. Letters after it: against raw, same rule: L loses, G gains, = no clear change. "
          "Against (raw, harris): + above, − below. Against herm_r, same rule: R recovers, F full recovery, n does not recover.", ""]
    alts = []
    for x in rows:
        if x["kind"] == "adj_power" and x["family"] == "ind" and (x["alt"], x["param"]) not in alts:
            alts.append((x["alt"], x["param"]))

    def entry(key):
        a = adj.get(key)
        if not a:
            return "–"
        s = pc(a["value"])
        if key in raw:
            s += f" / {pc(raw[key]['value'])}"
        marks = ""
        if key in gap:
            marks += SHORT[gap[key]["note"]]
        if key in ghar:
            marks += SHORT[ghar[key]["note"]]
        if key in gher:
            marks += SHORT[gher[key]["note"]]
        return f"{s} {marks}".strip()

    for alt, par in alts:
        hs = sorted({x["H"] for x in rows if x["kind"] == "adj_power" and x["alt"] == alt and x["param"] == par and x["family"] == "ind"})
        L += [f"### {alt} (param {float(par):g})", "", "| kernel, rule | " + " | ".join(f"H={h}" for h in hs) + " |", "|---|" + "---|" * len(hs)]
        for k, r in POWER_RULES:
            L.append(f"| {k}, {r} | " + " | ".join(entry(("ind", alt, par, h, k, r)) for h in hs) + " |")
        L.append("")
    # split MMD
    L += ["## Split MMD: verdicts, wild", "", "| kernel, rule | G0 valid | G0 removes memory | G0 memory-flat | G0 worst | HT-S valid | HT-S removes memory | HT-S worst |", "|---|---|---|---|---|---|---|---|"]
    for rk, a in res["verdicts_mmd_G0"].items():
        b = res["verdicts_mmd_HT-S"].get(rk)
        k, r = rk.split("|")
        L.append(f"| {k}, {r} | {a['n_valid']}/{a['cells']} | {'yes' if a['removes_memory'] else 'no'} | {'yes' if a['memory_flat'] else 'no'} | {pc(a['worst_rate'])} {tuple(a['worst_cell'])[-2:]} | "
                 f"{b['n_valid']}/{b['cells']} | {'yes' if b['removes_memory'] else 'no'} | {pc(b['worst_rate'])} {tuple(b['worst_cell'])[-2:]} |")
    mmd_rules = [(k, "wild") for k in ("raw", "herm1", "herm2", "sym", "datadiff", "hsplit1", "hsplit2")] + [("hsplit1", "wild_blk"), ("hsplit2", "wild_blk"), ("raw", "oracle")]
    for gen in ("G0", "HT-S"):
        for n in ("500", "2000"):
            hs = sorted({key[2] for key in size if key[0] == "mmd" and key[1] == gen and key[3] == n})
            if hs:
                L += ["", f"### Split MMD size, {gen}, n = {n}", "", "| kernel, rule | " + " | ".join(hs) + " |", "|---|" + "---|" * len(hs)]
                for k, r in mmd_rules:
                    cells = [size.get(("mmd", gen, h, n, k, r)) for h in hs]
                    if any(cells):
                        L.append(f"| {k}, {r} | " + " | ".join(f"{pc(c['value'])}{CLASS[c['note']]}" if c else "–" for c in cells) + " |")
    L += ["", "## Split MMD: power at n = 2000", "", "Size-adjusted power / raw power. Letters as above.", ""]
    alts = []
    for x in rows:
        if x["kind"] == "adj_power" and x["family"] == "mmd" and (x["alt"], x["param"]) not in alts:
            alts.append((x["alt"], x["param"]))
    for alt, par in alts:
        hs = sorted({x["H"] for x in rows if x["kind"] == "adj_power" and x["alt"] == alt and x["param"] == par and x["family"] == "mmd"})
        L += [f"### {alt} (param {float(par):g})", "", "| kernel, rule | " + " | ".join(f"H={h}" for h in hs) + " |", "|---|" + "---|" * len(hs)]
        for k, r in mmd_rules:
            L.append(f"| {k}, {r} | " + " | ".join(entry(("mmd", alt, par, h, k, r)) for h in hs) + " |")
        L.append("")
    # predictions
    p = res["predictions"]
    L += ["## Prediction readouts", "",
          f"- P1: inc1 wild valid {p['P1']['inc1_wild_valid']}, shift valid {p['P1']['inc1_shift_valid']}, sub valid {p['P1']['inc1_sub_valid']}. T flat at all H: {p['P1']['T_flat_all_H']}. Holds: {p['P1']['holds']}.",
          f"- P2: datadiff wild valid {p['P2']['datadiff_wild_valid']}, shift valid {p['P2']['datadiff_shift_valid']}. Gaussian part holds: {p['P2']['gaussian_part_holds']}.", ""]
    L += ["| NL-SV cell, rule | datadiff | inc1 | McNemar p | only datadiff | only inc1 |", "|---|---|---|---|---|---|"]
    for key, v in p["P2"]["sv"].items():
        L.append(f"| {key} | {pc(v['datadiff']['rate'])}{CLASS[v['datadiff']['class']]} | {pc(v['inc1']['rate'])}{CLASS[v['inc1']['class']]} | {v['mcnemar_p']:.3g} | {v['only_datadiff']} | {v['only_inc1']} |")
    L += ["", f"- P3: helps in a monotone cell where raw is invalid: {p['P3']['helps_monotone_where_raw_invalid']} (of {p['P3']['n_monotone_cells_raw_invalid']} such cells). Helps at NL-SV H = 0.8: {p['P3']['helps_sv']}. Holds: {p['P3']['holds']}.", "",
          "| cell, rule | raw | herm1 | drop | McNemar p | helps |", "|---|---|---|---|---|---|"]
    for key, v in p["P3"]["cells"].items():
        L.append(f"| {key.replace('ind__', '').replace('__null0', '')} | {pc(v['raw_rate'])} ({v['raw_class']}) | {pc(v['herm1_rate'])} | {pc(v['drop'])} | {v['mcnemar_p']:.3g} | {'yes' if v['helps'] else 'no'} |")
    L += ["", "- P6, mean layer memory on NL-SV:"]
    for key, v in p["P6"].items():
        L.append(f"  - {key}: first layer {v['mean_H_layer1']:.3f}, second layer {v['mean_H_layer2']:.3f}.")
    with open(os.path.join(HERE, f"results_tables_{tag}.md"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
