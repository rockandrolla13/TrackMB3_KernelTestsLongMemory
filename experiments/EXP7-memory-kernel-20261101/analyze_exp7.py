"""EXP7 analysis: size classes, "removes memory" verdicts, bare-statistic growth, size-adjusted power.

Usage:
  python -I analyze_exp7.py --tune                 -> tune_report.txt
  python -I analyze_exp7.py --in pilot --extra pilot_extras --tag pilot -> summary_pilot.csv, tables_pilot.md, results_pilot.json
  python -I analyze_exp7.py --in main  --tag main  -> summary_main.csv, tables_main.md, results_main.json
  python plot_growth.py main                        -> growth_q95_main.png (matplotlib does not import under python -I here)

Size adjustment, c*, Wilson and bootstrap intervals follow EXP3 analyze.py. Reading rules are in DESIGN.md.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402

import exp7_core as K  # noqa: E402
import run_exp7 as R  # noqa: E402
from kmem.analysis import boot_power, cstar, mcnemar_p, size_class, wilson  # noqa: E402,F401
from kmem.calibration import registry as cal  # noqa: E402

ALPHA = 0.05
NBOOT = 2000
SEED = 42
GAP_PP = 0.05
HELP_PP = 0.02
FLAT_SLOPE = 0.05
Z = 1.959963984540054


# ---------------------------------------------------------------- loading

def load(path: str) -> dict:
    """Column arrays for the complete status-ok rows, one per replication, sorted by rep."""
    with open(path) as f:
        rows = [r for r in csv.DictReader(f)]
    err = sum(1 for r in rows if r["status"] != "ok")
    cpu = sum(float(r["cpu_s"]) for r in rows)
    ok = {}
    for r in rows:
        if r["status"] == "ok" and None not in r.values():
            ok[int(r["rep"])] = r
    reps = sorted(ok)
    out = {"_errors": err, "_cpu": cpu, "_m": len(reps), "rep": np.array(reps)}
    if not reps:
        return out
    for c in rows[0]:
        if c in ("rep", "status", "error"):
            continue
        out[c] = np.array([float(ok[i][c]) if ok[i][c] != "" else np.nan for i in reps])
    return out


def derive(c: dict, d: dict) -> None:
    """hsplit ratios. The rule rejects iff the remainder has p <= alpha/2 or the low-layer ratio exceeds 1."""
    if d["_m"] == 0 or "rl1" not in d:
        return
    for r in K.HSPLIT_R:
        if c["family"] == "ind":
            d[f"R_shift_hsplit{r}"] = np.maximum(d[f"r975_shift_herm{r}"], d[f"rl{r}"])
            d[f"R_shift_blk_hsplit{r}"] = np.maximum(d[f"r975_shift_herm{r}"], d[f"rb{r}"])
        d[f"R_wild_hsplit{r}"] = np.maximum(d[f"r975_wb_herm{r}"], d[f"rl{r}"])
        d[f"R_wild_blk_hsplit{r}"] = np.maximum(d[f"r975_wb_herm{r}"], d[f"rb{r}"])


def load_all(indir: str, cells: list[dict], extra: str | None = None) -> dict:
    """extra: a folder from --extras-only, same seeds. Its columns join by replication, and shared columns must agree."""
    data = {}
    for c in cells:
        p = os.path.join(HERE, indir, c["name"] + ".csv")
        if not os.path.exists(p):
            continue
        d = load(p)
        d["_join_disagree"] = 0
        if extra and os.path.exists(os.path.join(HERE, extra, c["name"] + ".csv")):
            e = load(os.path.join(HERE, extra, c["name"] + ".csv"))
            common = np.intersect1d(d["rep"], e["rep"])
            ia, ib = np.searchsorted(d["rep"], common), np.searchsorted(e["rep"], common)
            for key in list(d):
                if not key.startswith("_") and isinstance(d[key], np.ndarray):
                    d[key] = d[key][ia]
            d["_m"] = common.size
            for key in e:
                if key.startswith("_") or key in ("rep", "cpu_s", "exact") or not isinstance(e[key], np.ndarray):
                    continue
                v = e[key][ib]
                if key in d and np.isfinite(d[key]).any():
                    d["_join_disagree"] += int(np.sum(np.abs(d[key] - v) > 1e-9 * np.maximum(np.abs(v), 1.0)))
                else:
                    d[key] = v
        derive(c, d)
        data[c["name"]] = (c, d)
    return data


# ---------------------------------------------------------------- statistics


# ---------------------------------------------------------------- rules

def ind_rules(kernels=K.IND_KERNELS) -> list[tuple]:
    """(kernel, rule, ratio column, function from data dict to a rejection array or None)."""
    out = []
    for k in kernels:
        out.append((k, "wild", f"r_wb_{k}", lambda d, k=k: cal.WILD.reject(d, suffix=f"_{k}")))
        out.append((k, "shift", f"r_shift_{k}", lambda d, k=k: cal.SHIFT.reject(d, suffix=f"_{k}")))
        if k in K.SUB_KERNELS:
            out.append((k, "sub", f"r_sub_{k}", lambda d, k=k: cal.SUBSAMPLE.reject(d, suffix=f"_{k}")))
        if k in K.HARRIS_KERNELS:
            out.append((k, "harris", f"rhc_{k}", lambda d, k=k: cal.HARRIS.reject(d, suffix=f"_{k}")))
            out.append((k, "harris19", f"rh19_{k}", lambda d, k=k: cal.HARRIS19.reject(d, suffix=f"_{k}")))
            out.append((k, "harrisA", f"rha_{k}", lambda d, k=k: cal.HARRIS_A.reject(d, suffix=f"_{k}")))
        out.append((k, "oracle", f"t_{k}", None))
    out += [("dcov", "wild", "r_dcov_wb", lambda d: cal.WILD.reject(d, "p_dcov_wb")),
            ("dcov", "sub", "r_dcov_sub", lambda d: cal.SUBSAMPLE.reject(d, "rej_dcov_sub")),
            ("dcov", "oracle", "dcov_stat", None),
            ("pearson", "nw", "r_nw", lambda d: cal.NEWEY_WEST.reject(d, "rej_nw")),
            ("pearson", "ewc", "r_ewc", lambda d: cal.EWC.reject(d, "rej_ewc")),
            ("pearson", "oracle", "abs_r", None)]
    out += hsplit_rules(("shift", "wild"))
    return out


def hsplit_rules(rules) -> list[tuple]:
    """Primary rules use the plug-in table (rl), *_blk use windows (rb). low, low_blk are the low-layer part alone."""
    out = []
    for r in K.HSPLIT_R:
        for rule in rules:
            pc = "p_shift" if rule == "shift" else "p_wb"
            out.append((f"hsplit{r}", rule, f"R_{rule}_hsplit{r}",
                        lambda d, r=r, pc=pc: (d[f"{pc}_herm{r}"] <= ALPHA / 2) | (d[f"rl{r}"] > 1)))
            out.append((f"hsplit{r}", f"{rule}_blk", f"R_{rule}_blk_hsplit{r}",
                        lambda d, r=r, pc=pc: (d[f"{pc}_herm{r}"] <= ALPHA / 2) | (d[f"rb{r}"] > 1)))
        out.append((f"hsplit{r}", "low", f"rl{r}", lambda d, r=r: cal.SN_LOW.reject(d, f"rl{r}")))
        out.append((f"hsplit{r}", "low_blk", f"rb{r}", lambda d, r=r: cal.SN_LOW_BLOCK.reject(d, f"rb{r}")))
    return out


def mmd_rules() -> list[tuple]:
    out = []
    for k in K.MMD_KERNELS:
        out.append((k, "wild", f"r_wb_{k}", lambda d, k=k: cal.WILD.reject(d, suffix=f"_{k}")))
        out.append((k, "oracle", f"mmd2_{k}", None))
    return out + hsplit_rules(("wild",))


def hlabel(c: dict) -> str:
    return f"{c['HX']}" if c["HX"] == c["HY"] else f"{c['HX']}/{c['HY']}"


def null_of(c: dict, data: dict) -> dict | None:
    """The matched null cell of an alternative cell."""
    if c["family"] == "mmd":
        nc = {**c, "alt": "null", "param": 0.0}
    elif c["alt"] == "factor":
        nc = {**c, "gen": "FACTOR", "alt": "null", "param": 0.0}
    elif c["alt"] == "scale":
        nc = {**c, "gen": "SCALE", "alt": "null"}
    else:
        nc = {**c, "alt": "null", "param": 0.0}
    name = R.name_of(nc)
    return data[name] if name in data else None


# ---------------------------------------------------------------- analysis blocks

def audit(data: dict) -> dict:
    """Fast path against the exact path on the exact replications; shift check on every row."""
    res = {"reps": 0, "max_rel_t": 0.0, "wb_flip": 0, "shift_flip": 0, "sub_flip": 0, "worst_cell": None,
           "mmd_reps": 0, "mmd_max_rel": 0.0, "mmd_flip": 0, "max_abs_check": 0.0, "max_rank": 0, "errors": 0, "cpu_hours": 0.0,
           "max_rel_r": 0.0, "join_disagree": 0, "harris_reps": 0, "harris_flip_raw": 0, "harris_flip_inc1": 0, "harris_max_rel": 0.0}
    for name, (c, d) in data.items():
        res["join_disagree"] += d.get("_join_disagree", 0)
        res["errors"] += d["_errors"]
        res["cpu_hours"] += d["_cpu"] / 3600
        if c["family"] == "growth" or d["_m"] == 0:
            res["max_rank"] = max(res["max_rank"], int(np.max(d["rx_rank"]))) if d["_m"] else res["max_rank"]
            continue
        ex = d["exact"] == 1
        if c["family"] == "ind":
            res["max_rank"] = max(res["max_rank"], int(np.max(d["rx_rank"])), int(np.max(d["ry_rank"])))
            if "xh19_raw" in d and np.isfinite(d["xh19_raw"]).any():
                ok = np.isfinite(d["xh19_raw"])
                res["harris_reps"] += int(ok.sum())
                for k in K.HARRIS_KERNELS:
                    res[f"harris_flip_{k}"] += int(np.sum((d[f"xh19_{k}"][ok] > 1) != (d[f"rh19_{k}"][ok] > 1)))
                    res["harris_max_rel"] = max(res["harris_max_rel"], float(np.max(np.abs(d[f"rh19_{k}"][ok] / d[f"xh19_{k}"][ok] - 1.0))))
            for k in K.IND_KERNELS:
                if f"check_{k}" in d and np.isfinite(d[f"check_{k}"]).any():
                    res["max_abs_check"] = max(res["max_abs_check"], float(np.nanmax(np.abs(d[f"check_{k}"]))))
            if ex.any():
                e = float(np.max(np.abs(d["t_raw"][ex] / d["xraw_t"][ex] - 1.0)))
                er = float(np.max(np.abs(d["r_shift_raw"][ex] / d["xraw_r_shift"][ex] - 1.0)))
                res["reps"] += int(ex.sum())
                if e > res["max_rel_t"]:
                    res["max_rel_t"], res["worst_cell"] = e, name
                res["max_rel_r"] = max(res["max_rel_r"], er)
                res["wb_flip"] += int(np.sum((d["p_wb_raw"][ex] <= ALPHA) != (d["xraw_p_wb"][ex] <= ALPHA)))
                res["shift_flip"] += int(np.sum((d["p_shift_raw"][ex] <= ALPHA) != (d["xraw_p_shift"][ex] <= ALPHA)))
                res["sub_flip"] += int(np.sum(d["rej_sub_raw"][ex] != d["xraw_rej_sub"][ex]))
        else:
            res["max_rank"] = max(res["max_rank"], int(np.max(d["rank"])))
            if ex.any():
                res["mmd_reps"] += int(ex.sum())
                res["mmd_max_rel"] = max(res["mmd_max_rel"], float(np.max(np.abs(d["mmd2_raw"][ex] / d["xraw_mmd2"][ex] - 1.0))))
                res["mmd_flip"] += int(np.sum((d["p_wb_raw"][ex] <= ALPHA) != (d["xraw_p_wb"][ex] <= ALPHA)))
    res["stop_rule_triggered"] = bool(res["max_rel_t"] > 1e-4 or res["mmd_max_rel"] > 1e-4 or res["join_disagree"] > 0
                                      or max(res["wb_flip"], res["shift_flip"], res["sub_flip"]) > 0.001 * max(res["reps"], 1)
                                      or res["mmd_flip"] > 0.001 * max(res["mmd_reps"], 1))
    return res


def size_block(data: dict, summary: list) -> dict:
    res = {}
    for name, (c, d) in data.items():
        if c["family"] == "growth" or c["alt"] != "null" or d["_m"] == 0:
            continue
        if c["family"] == "ind":
            d["abs_r"] = np.abs(d["pearson_r"])
        cell = {}
        for kern, rule, rc, raw in (ind_rules() if c["family"] == "ind" else mmd_rules()):
            if raw is None or rc not in d or not np.isfinite(d[rc]).any():
                continue
            # a rule counts on the rows where it ran. NaN also marks a non-positive NW long-run variance, as in EXP3.
            ok = np.isfinite(d[rc])
            v = raw(d)[ok]
            k, m = int(v.sum()), int(v.size)
            lo, hi = wilson(k, m)
            e = {"k": k, "m": m, "rate": k / m, "lo": lo, "hi": hi, "class": size_class(k, m)}
            cell[f"{kern}|{rule}"] = e
            summary.append({"family": c["family"], "gen": c["gen"], "alt": "null", "H": hlabel(c), "n": c["n"], "param": c["param"],
                            "kernel": kern, "rule": rule, "kind": "size", "value": k / m, "lo": lo, "hi": hi, "m": m, "note": e["class"]})
        res[name] = cell
    return res


def memory_verdicts(data: dict, size: dict, family: str = "ind", gen: str = "G0") -> dict:
    """The "removes memory" rule and the weaker "memory-flat" reading, per (kernel, rule)."""
    cells = {(c["HX"], c["HY"], c["n"]): name for name, (c, d) in data.items()
             if c["family"] == family and c["gen"] == gen and c["alt"] == "null" and name in size}
    if not cells:
        return {}
    ns = sorted({k[2] for k in cells})
    pairs = sorted({k[:2] for k in cells})
    out = {}
    for rk in size[next(iter(cells.values()))]:
        if rk.split("|")[1] in ("low", "low_blk") or any(rk not in size[n] for n in cells.values()):
            continue
        classes = {kc: size[n][rk]["class"] for kc, n in cells.items()}
        n_valid = sum(v == "valid" for v in classes.values())
        rises = []
        for p in pairs:
            a, b = cells.get((*p, ns[0])), cells.get((*p, ns[-1]))
            if a and b:
                ea, eb = size[a][rk], size[b][rk]
                if eb["rate"] > ea["rate"] and eb["lo"] > ea["hi"]:
                    rises.append({"H": p, "from": ea["rate"], "to": eb["rate"]})
        above = []
        for n in ns:
            base = cells.get((0.5, 0.5, n))
            if base is None:
                continue
            e0 = size[base][rk]
            for p in pairs:
                e = size[cells[(*p, n)]][rk] if (*p, n) in cells else None
                if e and e["rate"] > e0["rate"] and e["lo"] > e0["hi"]:
                    above.append({"H": p, "n": n, "rate": e["rate"], "base": e0["rate"]})
        worst = max(((size[cells[kc]][rk]["rate"], kc) for kc in cells), key=lambda t: t[0])
        out[rk] = {"cells": len(cells), "n_valid": n_valid, "removes_memory": n_valid == len(cells) and not rises,
                   "memory_flat": not rises and not above, "rises": rises, "above_baseline": above,
                   "n_mild": sum(v == "mild" for v in classes.values()), "n_over": sum(v == "over" for v in classes.values()),
                   "worst_rate": worst[0], "worst_cell": list(worst[1])}
    return out


def power_block(data: dict, summary: list, rng: np.random.Generator) -> dict:
    res = {}
    for name, (c, d) in sorted(data.items(), key=lambda t: t[1][0]["idx"]):
        if c["family"] == "growth" or c["alt"] == "null" or d["_m"] == 0:
            continue
        nc = null_of(c, data)
        if nc is None or nc[1]["_m"] == 0:
            continue
        nd = nc[1]
        if c["family"] == "ind":
            d["abs_r"] = np.abs(d["pearson_r"])
            nd["abs_r"] = np.abs(nd["pearson_r"])
        bi_n = rng.integers(0, nd["_m"], size=(NBOOT, nd["_m"]))
        bi_a = rng.integers(0, d["_m"], size=(NBOOT, d["_m"]))
        boots, cell = {}, {}
        for kern, rule, rc, raw in (ind_rules() if c["family"] == "ind" else mmd_rules()):
            if rc not in nd or rc not in d:
                continue
            nv, av = nd[rc], d[rc]
            if not np.isfinite(nv).any() or not np.isfinite(av).any():
                continue
            if kern == "pearson":  # NaN ratios (NW long-run variance <= 0) never reject, as in EXP3
                nvf = nv[np.isfinite(nv)]
                nvc = np.where(np.isfinite(nv), nv, -np.inf)
                avc = np.where(np.isfinite(av), av, -np.inf)
                ok_n, ok_a = np.ones(nv.size, bool), np.ones(av.size, bool)
            else:  # NaN means this rule did not run on that row (DESIGN.md cuts)
                ok_n, ok_a = np.isfinite(nv), np.isfinite(av)
                nvf = nvc = nv[ok_n]
                avc = av[ok_a]
            cs = cstar(nvf)
            k_adj = int(np.sum(avc > cs))
            lo, hi = wilson(k_adj, avc.size)
            b = boot_power(nvc, avc, bi_n[:, :nvc.size] % nvc.size, bi_a[:, :avc.size] % avc.size)
            boots[(kern, rule)] = b
            e = {"adj": k_adj / avc.size, "w_lo": lo, "w_hi": hi, "b_lo": float(np.quantile(b, 0.025)), "b_hi": float(np.quantile(b, 0.975)),
                 "cstar": cs, "m_alt": int(avc.size), "m_null": int(nvf.size), "kern": kern, "rule": rule}
            if raw is not None:
                e["raw"] = float(np.mean(raw(d)[ok_a]))
                e["null_size"] = float(np.mean(raw(nd)[ok_n]))
            cell[f"{kern}|{rule}"] = e
        for e_key, e in cell.items():
            kern, rule = e["kern"], e["rule"]
            base_rule = {"shift_blk": "shift", "wild_blk": "wild"}.get(rule, rule)
            if c["family"] == "ind" and rule in ("shift", "wild", "shift_blk", "wild_blk", "sub") and ("raw", "harris") in boots:
                g = boots[(kern, rule)] - boots[("raw", "harris")]
                gh = e["adj"] - cell["raw|harris"]["adj"]
                hlo, hhi = float(np.quantile(g, 0.025)), float(np.quantile(g, 0.975))
                hv = "above Harris" if gh > GAP_PP and hlo > 0 else "below Harris" if gh < -GAP_PP and hhi < 0 else "no clear change"
                e.update({"gh": gh, "gh_lo": hlo, "gh_hi": hhi, "gh_verdict": hv})
                summary.append({"family": c["family"], "gen": c["gen"], "alt": c["alt"], "H": hlabel(c), "n": c["n"], "param": c["param"],
                                "kernel": kern, "rule": rule, "kind": "gap_vs_harris", "value": gh, "lo": hlo, "hi": hhi, "m": e["m_alt"], "note": hv})
            if kern in ("raw", "dcov", "pearson") or rule in ("low", "low_blk") or ("raw", base_rule) not in boots:
                e["gap"] = None
            else:
                base = cell[f"raw|{base_rule}"]
                g = boots[(kern, rule)] - boots[("raw", base_rule)]
                gap = e["adj"] - base["adj"]
                glo, ghi = float(np.quantile(g, 0.025)), float(np.quantile(g, 0.975))
                v = "loses" if gap < -GAP_PP and ghi < 0 else "gains" if gap > GAP_PP and glo > 0 else "no clear change"
                e.update({"gap": gap, "gap_lo": glo, "gap_hi": ghi, "verdict": v})
                if kern.startswith("hsplit"):
                    herm = f"herm{kern[-1]}"
                    g = boots[(kern, rule)] - boots[(herm, base_rule)]
                    rec = e["adj"] - cell[f"{herm}|{base_rule}"]["adj"]
                    rlo, rhi = float(np.quantile(g, 0.025)), float(np.quantile(g, 0.975))
                    rv = "recovers" if rec > GAP_PP and rlo > 0 else "does not recover"
                    if rv == "recovers" and gap >= -GAP_PP:
                        rv = "full recovery"
                    e.update({"rec": rec, "rec_lo": rlo, "rec_hi": rhi, "rec_verdict": rv})
                    summary.append({"family": c["family"], "gen": c["gen"], "alt": c["alt"], "H": hlabel(c), "n": c["n"], "param": c["param"],
                                    "kernel": kern, "rule": rule, "kind": "gap_vs_herm", "value": rec, "lo": rlo, "hi": rhi, "m": d["_m"], "note": rv})
            summary.append({"family": c["family"], "gen": c["gen"], "alt": c["alt"], "H": hlabel(c), "n": c["n"], "param": c["param"],
                            "kernel": kern, "rule": rule, "kind": "adj_power", "value": e["adj"], "lo": e["b_lo"], "hi": e["b_hi"], "m": e["m_alt"],
                            "note": f"cstar {e['cstar']:.5g}; m_null {e['m_null']}; Wilson [{e['w_lo']:.3f}, {e['w_hi']:.3f}]"})
            if "raw" in e:
                rk = int(round(e["raw"] * e["m_alt"]))
                rlo, rhi = wilson(rk, e["m_alt"])
                summary.append({"family": c["family"], "gen": c["gen"], "alt": c["alt"], "H": hlabel(c), "n": c["n"], "param": c["param"],
                                "kernel": kern, "rule": rule, "kind": "raw_power", "value": e["raw"], "lo": rlo, "hi": rhi, "m": e["m_alt"],
                                "note": f"null size {e['null_size']:.3f}"})
            if e.get("gap") is not None:
                summary.append({"family": c["family"], "gen": c["gen"], "alt": c["alt"], "H": hlabel(c), "n": c["n"], "param": c["param"],
                                "kernel": kern, "rule": rule, "kind": "gap_vs_raw", "value": e["gap"], "lo": e["gap_lo"], "hi": e["gap_hi"], "m": e["m_alt"],
                                "note": e["verdict"]})
        res[name] = cell
    return res


def growth_block(data: dict, summary: list, rng: np.random.Generator) -> dict:
    """95th percentile of T against n, with a bootstrap interval, and the OLS slope of log q95 on log n."""
    res = {"q95": {}, "slope": {}}
    hs = sorted({c["HX"] for c, d in data.values() if c["family"] == "growth"})
    ns = sorted({c["n"] for c, d in data.values() if c["family"] == "growth"})
    for h in hs:
        for kern in ("raw", "inc1"):
            boots, q = [], []
            have = []
            for n in ns:
                name = f"growth__G0__null0__H{h}__n{n}"
                if name not in data or data[name][1]["_m"] == 0:
                    continue
                d = data[name][1]
                t = d[f"t_{kern}"]
                bi = rng.integers(0, t.size, size=(NBOOT, t.size))
                b = np.quantile(t[bi], 0.95, axis=1)
                val = float(np.quantile(t, 0.95))
                lo, hi = float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))
                boots.append(b)
                q.append(val)
                have.append(n)
                res["q95"][f"{kern}|{h}|{n}"] = {"q95": val, "lo": lo, "hi": hi, "m": int(t.size), "median": float(np.median(t))}
                summary.append({"family": "growth", "gen": "G0", "alt": "null", "H": str(h), "n": n, "param": 0.0, "kernel": kern, "rule": "bare",
                                "kind": "growth_q95", "value": val, "lo": lo, "hi": hi, "m": int(t.size), "note": ""})
            if len(have) >= 2:
                lx = np.log(np.array(have, float))
                lx = lx - lx.mean()
                w = lx / np.dot(lx, lx)
                slope = float(np.dot(w, np.log(np.array(q))))
                bs = np.log(np.array(boots)).T @ w
                lo, hi = float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))
                ref = max(0.0, 4 * h - 3) if kern == "raw" else 0.0
                flat = bool(lo <= 0.0 <= hi or abs(slope) <= FLAT_SLOPE)
                res["slope"][f"{kern}|{h}"] = {"slope": slope, "lo": lo, "hi": hi, "reference": ref, "flat": flat,
                                               "ratio_last_first": q[-1] / q[0]}
                summary.append({"family": "growth", "gen": "G0", "alt": "null", "H": str(h), "n": "all", "param": 0.0, "kernel": kern, "rule": "bare",
                                "kind": "growth_slope", "value": slope, "lo": lo, "hi": hi, "m": int(len(have)),
                                "note": f"reference {ref:.2f}; {'flat' if flat else 'not flat'}"})
    return res


def rej_of(data: dict, name: str, kern: str, rule: str) -> np.ndarray:
    d = data[name][1]
    for k, r, rc, raw in ind_rules():
        if (k, r) == (kern, rule):
            return raw(d)
    raise KeyError((kern, rule))


def predictions(data: dict, size: dict, verdicts: dict, growth: dict) -> dict:
    out = {}
    def nm(gen, h, n=2000):
        return f"ind__{gen}__null0__H{h}__n{n}"

    # P1
    v_w, v_s, v_b = verdicts.get("inc1|wild", {}), verdicts.get("inc1|shift", {}), verdicts.get("inc1|sub", {})
    slopes = {k.split("|")[1]: v for k, v in growth.get("slope", {}).items() if k.startswith("inc1|")}
    out["P1"] = {"inc1_wild_valid": f"{v_w.get('n_valid')}/{v_w.get('cells')}", "inc1_shift_valid": f"{v_s.get('n_valid')}/{v_s.get('cells')}",
                 "inc1_sub_valid": f"{v_b.get('n_valid')}/{v_b.get('cells')}", "T_flat_all_H": all(v["flat"] for v in slopes.values()) if slopes else None,
                 "holds": bool(v_w.get("n_valid") == v_w.get("cells") and v_s.get("n_valid") == v_s.get("cells")
                               and slopes and all(v["flat"] for v in slopes.values()))}
    # P2
    dw, ds = verdicts.get("datadiff|wild", {}), verdicts.get("datadiff|shift", {})
    sv = {}
    for h in (0.5, 0.8):
        name = nm("NL-SV", h)
        if name not in data:
            continue
        for rule in ("wild", "shift"):
            a, b = rej_of(data, name, "datadiff", rule), rej_of(data, name, "inc1", rule)
            p, n10, n01 = mcnemar_p(a, b)
            sv[f"H{h}|{rule}"] = {"datadiff": size[name][f"datadiff|{rule}"], "inc1": size[name][f"inc1|{rule}"],
                                  "mcnemar_p": p, "only_datadiff": n10, "only_inc1": n01}
    out["P2"] = {"datadiff_wild_valid": f"{dw.get('n_valid')}/{dw.get('cells')}", "datadiff_shift_valid": f"{ds.get('n_valid')}/{ds.get('cells')}",
                 "gaussian_part_holds": bool(dw.get("n_valid") == dw.get("cells") and ds.get("n_valid") == ds.get("cells")), "sv": sv}
    # P3
    p3 = {}
    tests = [(nm("G0", h, n), "monotone") for n in (2000, 4000) for h in (0.8, 0.9, 0.95)]
    tests += [(nm("HT-S", 0.8), "monotone"), (nm("NL-SV", 0.8), "sv"), (nm("NL-SV", 0.5), "sv_h05")]
    for name, tag in tests:
        if name not in data:
            continue
        for rule in ("wild", "shift", "sub"):
            if rule == "sub":
                continue  # herm1 has no subsampling rule
            a, b = rej_of(data, name, "herm1", rule), rej_of(data, name, "raw", rule)
            p, n10, n01 = mcnemar_p(a, b)
            drop = float(b.mean() - a.mean())
            p3[f"{name}|{rule}"] = {"tag": tag, "raw_rate": float(b.mean()), "herm1_rate": float(a.mean()), "raw_class": size[name][f"raw|{rule}"]["class"],
                                    "drop": drop, "mcnemar_p": p, "helps": bool(drop > HELP_PP and p < 0.05)}
    mono = [v for v in p3.values() if v["tag"] == "monotone" and v["raw_class"] != "valid"]
    out["P3"] = {"cells": p3, "helps_monotone_where_raw_invalid": any(v["helps"] for v in mono) if mono else None,
                 "n_monotone_cells_raw_invalid": len(mono),
                 "helps_sv": any(v["helps"] for v in p3.values() if v["tag"] == "sv")}
    out["P3"]["holds"] = bool(out["P3"]["helps_monotone_where_raw_invalid"] and not out["P3"]["helps_sv"])
    # P4
    out["P4"] = {rk: {"valid": f"{v['n_valid']}/{v['cells']}", "removes_memory": v["removes_memory"], "worst": v["worst_rate"], "worst_cell": v["worst_cell"]}
                 for rk, v in verdicts.items() if rk.startswith("hsplit")}
    low = {}
    for name, (c, d) in data.items():
        if name in size and c["family"] == "ind" and c["gen"] == "G0":
            low[f"H{hlabel(c)}|n{c['n']}"] = {k: size[name][k]["rate"] for k in ("hsplit1|low", "hsplit2|low", "hsplit1|low_blk", "hsplit2|low_blk") if k in size[name]}
    out["P4"]["low_layer_rates_nominal_2.5pc"] = low
    # P6
    p6 = {}
    for h in (0.5, 0.8):
        name = nm("NL-SV", h)
        if name in data and "hl_x1" in data[name][1]:
            d = data[name][1]
            p6[f"H{h}"] = {"mean_H_layer1": float(np.mean(np.r_[d["hl_x1"], d["hl_y1"]])), "mean_H_layer2": float(np.mean(np.r_[d["hl_x2"], d["hl_y2"]])),
                           "rates": {k: size[name][k] for k in ("raw|wild", "raw|shift", "herm1|wild", "herm1|shift", "herm2|wild", "herm2|shift", "hsplit1|wild", "hsplit1|shift",
                                                                 "hsplit2|wild", "hsplit2|shift") if k in size[name]}}
    out["P6"] = p6
    # P7
    p7 = {}
    for g in ("SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow"):
        for h in (0.5, 0.8):
            name = nm(g, h)
            if name in size:
                p7[f"{g}|H{h}"] = {k: size[name][k]["rate"] for k in ("raw|shift", "herm1|shift", "hsplit1|shift", "hsplit2|shift", "hsplit1|low", "hsplit2|low") if k in size[name]}
    out["P7"] = p7
    return out


# ---------------------------------------------------------------- tables

def pct(x: float) -> str:
    return f"{100 * x:.1f}"


CLS = {"valid": "", "mild": "*", "over": "**"}


def tables(tag: str, data: dict, size: dict, verdicts: dict, mmd_verdicts: dict, mmd_ht: dict, growth: dict, power: dict, pred: dict, aud: dict) -> str:
    L = [f"# EXP7 tables ({tag})", "",
         "Rates and power in percent. Size class marks: no mark valid, * mildly oversized, ** oversized (DESIGN.md).", "",
         f"Errors: {aud['errors']}. CPU-hours: {aud['cpu_hours']:.2f}. Largest pivoted-Cholesky rank: {aud['max_rank']}.", "",
         "## Fast path against the exact path", "",
         f"Independence: {aud['reps']} exact replications. Largest relative error in T: {aud['max_rel_t']:.2e} ({aud['worst_cell']}). "
         f"Largest in the shift ratio: {aud['max_rel_r']:.2e}. Decisions that differ: wild {aud['wb_flip']}, shift {aud['shift_flip']}, sub {aud['sub_flip']}.", "",
         f"Split MMD: {aud['mmd_reps']} exact replications. Largest relative error in MMD²: {aud['mmd_max_rel']:.2e}. Decisions that differ: {aud['mmd_flip']}.", "",
         f"Largest |S(0) − T| over every kernel and row: {aud['max_abs_check']:.2e}. Stop rule triggered: {aud['stop_rule_triggered']}.", ""]
    rules_all = [(k, r) for k, r, rc, raw in ind_rules() if raw is not None]
    L += [f"Join of the extras-only pass, shared columns that disagree: {aud['join_disagree']}.", "",
          f"Harris at N = 19, whole-series bandwidths against segment bandwidths: {aud['harris_reps']} replications. Decisions that differ: raw {aud['harris_flip_raw']}, inc1 "
          f"{aud['harris_flip_inc1']}. Largest relative difference in the ratio: {aud['harris_max_rel']:.2e}.", ""]

    def size_table(cells: list[tuple], rule_list, head: str):
        L.extend([f"### {head}", "", "| kernel, rule | " + " | ".join(c[0] for c in cells) + " |", "|---|" + "---|" * len(cells)])
        for k, r in rule_list:
            row = []
            for _, name in cells:
                e = size.get(name, {}).get(f"{k}|{r}")
                mark = "" if not e or r in ("low", "low_blk") else CLS[e["class"]]
                row.append(f"{pct(e['rate'])}{mark}" if e else "–")
            L.append(f"| {k}, {r} | " + " | ".join(row) + " |")
        L.append("")

    L += ["## Size: Gaussian memory cells", ""]
    ns = sorted({c["n"] for c, d in data.values() if c["family"] == "ind" and c["gen"] == "G0" and c["alt"] == "null"})
    for n in ns:
        cells = [(hlabel(c), name) for name, (c, d) in data.items() if c["family"] == "ind" and c["gen"] == "G0" and c["alt"] == "null" and c["n"] == n]
        ms = sorted({data[nm][1]["_m"] for _, nm in cells})
        size_table(cells, rules_all, f"n = {n}, replications {ms}")
    L += ["## Verdicts: Gaussian memory cells", "",
          "| kernel, rule | valid cells | removes memory | memory-flat | worst rate (HX, HY, n) | rises n=500→4000 |", "|---|---|---|---|---|---|"]
    for rk, v in verdicts.items():
        k, r = rk.split("|")
        L.append(f"| {k}, {r} | {v['n_valid']}/{v['cells']} | {'yes' if v['removes_memory'] else 'no'} | {'yes' if v['memory_flat'] else 'no'} | "
                 f"{pct(v['worst_rate'])} {tuple(v['worst_cell'])} | {len(v['rises'])} |")
    L += ["", "## Size: process arms, n = 2000", ""]
    for h in (0.5, 0.8):
        cells = [(c["gen"], name) for name, (c, d) in data.items() if c["family"] == "ind" and c["alt"] == "null" and c["n"] == 2000 and c["HX"] == h
                 and c["HY"] == h and c["gen"] not in ("FACTOR", "SCALE")]
        size_table(cells, rules_all, f"H = {h}")
    L += ["## Size: matched nulls of the stress alternatives", ""]
    cells = [(f"{c['gen']} H{c['HX']}", name) for name, (c, d) in data.items() if c["family"] == "ind" and c["gen"] in ("FACTOR", "SCALE")]
    if cells:
        size_table(cells, rules_all, "n = 2000")
    L += ["## Bare-statistic growth: 95th percentile of T", ""]
    if growth["q95"]:
        hs = sorted({float(k.split('|')[1]) for k in growth["q95"]})
        ng = sorted({int(k.split('|')[2]) for k in growth["q95"]})
        for kern in ("raw", "inc1"):
            L += [f"### {kern}", "", "| H | " + " | ".join(f"n={n}" for n in ng) + " | slope [95%] | reference | flat? |", "|---|" + "---|" * (len(ng) + 3)]
            for h in hs:
                row = []
                for n in ng:
                    e = growth["q95"].get(f"{kern}|{h}|{n}")
                    row.append(f"{e['q95']:.4g} [{e['lo']:.4g}, {e['hi']:.4g}]" if e else "–")
                s = growth["slope"].get(f"{kern}|{h}")
                ss = f"{s['slope']:.3f} [{s['lo']:.3f}, {s['hi']:.3f}] | {s['reference']:.2f} | {'yes' if s['flat'] else 'no'}" if s else "– | – | –"
                L.append(f"| {h} | " + " | ".join(row) + f" | {ss} |")
            L.append("")
    L += ["## Power (independence): size-adjusted [bootstrap 95%] / raw / null size, n = 2000", ""]
    def power_tables(fam: str, kernels_rules):
        alts = []
        for name, (c, d) in sorted(data.items(), key=lambda t: t[1][0]["idx"]):
            if c["family"] == fam and c["alt"] != "null" and (c["alt"], c["param"]) not in alts:
                alts.append((c["alt"], c["param"]))
        for alt, par in alts:
            cs = [(c["HX"], name) for name, (c, d) in data.items() if c["family"] == fam and c["alt"] == alt and c["param"] == par and name in power]
            cs.sort()
            if not cs:
                continue
            L.extend([f"### {alt} (param {par:g})", "", "| kernel, rule | " + " | ".join(f"H={h}" for h, _ in cs) + " |", "|---|" + "---|" * len(cs)])
            for k, r in kernels_rules:
                if not any(f"{k}|{r}" in power[name] for _, name in cs):
                    continue
                row = []
                for _, name in cs:
                    e = power[name].get(f"{k}|{r}")
                    if not e:
                        row.append("–")
                        continue
                    s = f"{pct(e['adj'])} [{pct(e['b_lo'])}, {pct(e['b_hi'])}]"
                    if "raw" in e:
                        s += f" / {pct(e['raw'])} / {pct(e['null_size'])}"
                    if e.get("gap") is not None:
                        s += f"; gap {100 * e['gap']:+.1f} {e['verdict']}"
                    if "rec" in e:
                        s += f"; vs herm {100 * e['rec']:+.1f} {e['rec_verdict']}"
                    if "gh" in e:
                        s += f"; vs Harris {100 * e['gh']:+.1f}"
                    row.append(s)
                L.append(f"| {k}, {r} | " + " | ".join(row) + " |")
            L.append("")
    power_tables("ind", [(k, r) for k, r, _, _ in ind_rules()])
    L += ["## Split MMD: size, uncentred wild", ""]
    for gen in ("G0", "HT-S"):
        for n in (500, 2000):
            cells = sorted([(c["HX"], name) for name, (c, d) in data.items() if c["family"] == "mmd" and c["gen"] == gen and c["alt"] == "null" and c["n"] == n])
            if cells:
                size_table([(f"H={h}", nm) for h, nm in cells], [(k, r) for k, r, _, raw in mmd_rules() if raw is not None], f"{gen}, n = {n}")
    L += ["### Verdicts, split MMD, wild", "", "| kernel, rule | G0 valid cells | G0 removes memory | G0 memory-flat | HT-S valid cells | HT-S removes memory |", "|---|---|---|---|---|---|"]
    for rk in mmd_verdicts:
        a, b = mmd_verdicts.get(rk), mmd_ht.get(rk)
        if a and b:
            k, r = rk.split("|")
            L.append(f"| {k}, {r} | {a['n_valid']}/{a['cells']} | {'yes' if a['removes_memory'] else 'no'} | {'yes' if a['memory_flat'] else 'no'} | "
                     f"{b['n_valid']}/{b['cells']} | {'yes' if b['removes_memory'] else 'no'} |")
    L += ["", "## Split MMD: power, size-adjusted [bootstrap 95%] / raw / null size, n = 2000", ""]
    power_tables("mmd", [(k, r) for k, r, _, _ in mmd_rules()])
    L += ["## Prediction readouts", "", "```json", json.dumps(pred, indent=1, default=float), "```", ""]
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------- tuning report

def tune_report() -> None:
    cells = R.tune_cells()
    data = load_all("tune", cells)
    rng = np.random.default_rng(SEED)
    L = ["EXP7 tuning report: size-adjusted power at G0, H = 0.5, n = 2000, 200 replications per candidate.", ""]
    def adj(nv, av):
        return float(np.mean(av > cstar(nv)))
    g0 = data["ind__G0__null0__H0.5__n2000"][1]
    L.append(f"errors: {sum(d['_errors'] for _, d in data.values())}")
    L.append("scale coupling (matched SCALE-null with the same gamma):")
    for g in (0.05, 0.1, 0.15, 0.2, 0.3):
        a = data[f"ind__G0__scale{g:g}__H0.5__n2000"][1]
        n = data[f"ind__SCALE__null{g:g}__H0.5__n2000"][1]
        L.append(f"  gamma {g}: oracle raw {adj(n['t_raw'], a['t_raw']):.3f}  oracle inc1 {adj(n['t_inc1'], a['t_inc1']):.3f}")
    for alt in ("lag", "flat"):
        a = data[f"ind__G0__{alt}0.13__H0.5__n2000"][1]
        L.append(f"{alt} rho 0.13: oracle raw {adj(g0['t_raw'], a['t_raw']):.3f}  oracle inc1 {adj(g0['t_inc1'], a['t_inc1']):.3f}")
    m0 = data["mmd__G0__null0__H0.5__n2000"][1]
    for alt, grid in (("scale", (1.05, 1.1, 1.15, 1.2, 1.3)), ("shape", (0.1, 0.15, 0.2, 0.3, 0.4))):
        L.append(f"split MMD {alt}:")
        for v in grid:
            a = data[f"mmd__G0__{alt}{v:g}__H0.5__n2000"][1]
            s = "  ".join(f"{k} {adj(m0[f'mmd2_{k}'], a[f'mmd2_{k}']):.3f}" for k in K.MMD_KERNELS)
            L.append(f"  {v}: {s}")
    open(os.path.join(HERE, "tune_report.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", default=None)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--tune", action="store_true")
    ap.add_argument("--extra", default=None, help="folder from an --extras-only pass, joined by replication")
    a = ap.parse_args()
    if a.tune:
        tune_report()
        return
    cells = R.main_cells()
    data = load_all(a.indir, cells, a.extra)
    print(f"cells with data: {len(data)} of {len(cells)}", flush=True)
    rng = np.random.default_rng(SEED)
    summary: list = []
    aud = audit(data)
    size = size_block(data, summary)
    ver = memory_verdicts(data, size)
    mmd_ver = memory_verdicts(data, size, family="mmd", gen="G0")
    mmd_ht = memory_verdicts(data, size, family="mmd", gen="HT-S")
    growth = growth_block(data, summary, rng)
    power = power_block(data, summary, rng)
    pred = predictions(data, size, ver, growth)
    fields = ["family", "gen", "alt", "H", "n", "param", "kernel", "rule", "kind", "value", "lo", "hi", "m", "note"]
    with open(os.path.join(HERE, f"summary_{a.tag}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(summary)
    res = {"audit": aud, "verdicts_memory": ver, "verdicts_mmd_G0": mmd_ver, "verdicts_mmd_HT-S": mmd_ht, "growth": growth,
           "predictions": pred, "size": size, "power": power, "reps": {n: d["_m"] for n, (c, d) in data.items()}}
    with open(os.path.join(HERE, f"results_{a.tag}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    with open(os.path.join(HERE, f"tables_{a.tag}.md"), "w") as f:
        f.write(tables(a.tag, data, size, ver, mmd_ver, mmd_ht, growth, power, pred, aud))
    print(f"errors {aud['errors']}, CPU-h {aud['cpu_hours']:.2f}, stop rule triggered: {aud['stop_rule_triggered']}", flush=True)


if __name__ == "__main__":
    main()
