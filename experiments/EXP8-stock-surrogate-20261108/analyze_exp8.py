"""EXP8 analysis: rejection rates per cell and the pass/fail verdicts of P0-P5 (DESIGN.md, amendments A1, A2).

  python experiments/EXP8-stock-surrogate-20261108/analyze_exp8.py  -> results.json, summary.csv, tables.md
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "EXP7-memory-kernel-20261101"))
import numpy as np  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import analyze_exp7 as A7  # noqa: E402
import run_exp8 as R  # noqa: E402
from kmem.analysis import wilson  # noqa: E402

MAIN = os.path.join(HERE, "main")
STOCK = os.path.join(ROOT, "experiments", "STOCK-20261023", "results")
BAND = (3.4, 6.8)  # P1 band, percent
TEXT = {"status", "error", "day_a", "day_b"}
RULES = [r for r in A7.ind_rules(R.KERNELS) if r[3] is not None]


def load(name: str) -> dict:
    with open(os.path.join(MAIN, name + ".csv"), newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["status"] == "ok"]
    rows.sort(key=lambda r: int(r["rep"]))
    return {k: (np.array([r[k] for r in rows]) if k in TEXT else np.array([float(r[k]) if r[k] != "" else np.nan for r in rows]))
            for k in rows[0]}


def rej(d: dict, kern: str, rule: str) -> np.ndarray:
    for k, r, _, fn in RULES:
        if k == kern and r == rule:
            return np.asarray(fn(d), dtype=bool)
    raise KeyError((kern, rule))


def pct(k: int, m: int) -> dict:
    lo, hi = wilson(k, m)
    return {"k": int(k), "m": int(m), "rate": 100 * k / m, "lo": 100 * lo, "hi": 100 * hi}


def newcombe(k1, m1, k2, m2) -> tuple[float, float, float]:
    """Difference p1 - p2 in percent with the Newcombe hybrid-score 95% interval."""
    p1, p2 = k1 / m1, k2 / m2
    l1, u1 = wilson(k1, m1)
    l2, u2 = wilson(k2, m2)
    d = p1 - p2
    lo = d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    hi = d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return 100 * d, 100 * lo, 100 * hi


def pooled(data: dict, arm: str, n: int, kern: str, rule: str) -> np.ndarray:
    return np.concatenate([rej(data[f"{arm}__{a}-{b}__n{n}"], kern, rule) for a, b in R.PAIRS])


def trend_flags(d: dict) -> tuple[np.ndarray, np.ndarray]:
    cut = np.quantile(np.abs(np.concatenate([d["slope_perm_x"], d["slope_perm_y"]])), 0.95)
    return np.abs(d["slope_x"]) > cut, np.abs(d["slope_y"]) > cut


def p0(data: dict) -> dict:
    out, ok = {}, True
    for a, b in R.PAIRS:
        for n in R.NS:
            d = data[f"A__{a}-{b}__n{n}"]
            with open(os.path.join(STOCK, f"{a}-{b}_n{n}.csv"), newline="") as f:
                s = {int(r["pair"]): r for r in csv.DictReader(f)}
            bad = 0
            for i, rep in enumerate(d["rep"].astype(int)):
                r = s[rep]
                same = (d["day_a"][i] == r["day_a"] and d["day_b"][i] == r["day_b"]
                        and int(d["start_a"][i]) == int(r["start_a"]) and int(d["start_b"][i]) == int(r["start_b"])
                        and f"{d['p_shift_raw'][i]:.12g}" == f"{float(r['p_hsic_shift']):.12g}"
                        and abs(d["pearson_r"][i] - float(r["pearson_r"])) <= 1e-12)
                bad += not same
            out[f"{a}-{b}_n{n}"] = {"pairs": int(d["rep"].size), "mismatches": bad}
            ok &= bad == 0 and d["rep"].size == 500
    return {"pass": bool(ok), "cells": out}


def p1(data: dict, arm: str = "S1") -> dict:
    cells, ok = {}, True
    for a, b in R.PAIRS:
        d = data[f"{arm}__{a}-{b}__n2000"]
        for k in ("raw", "inc1"):
            v = rej(d, k, "shift")
            c = pct(v.sum(), v.size)
            cells[f"{a}-{b}_{k}"] = c
            ok &= BAND[0] <= c["rate"] <= BAND[1]
    return {"pass": bool(ok), "cells": cells}


def p2(data: dict) -> dict:
    parts = {}
    for rule in ("wild", "sub"):
        v = pooled(data, "S1", 2000, "raw", rule)
        c = pct(v.sum(), v.size)
        parts[rule] = {**c, "pass": c["lo"] > 5}
    hbar = np.concatenate([(data[f"S1__{a}-{b}__n2000"]["H_x"] + data[f"S1__{a}-{b}__n2000"]["H_y"]) / 2
                           for a, b in R.PAIRS])
    v = pooled(data, "S1", 2000, "raw", "wild").astype(float)
    fit = sm.Logit(v, sm.add_constant(hbar)).fit(disp=0)
    slope, p_two = float(fit.params[1]), float(fit.pvalues[1])
    p_one = p_two / 2 if slope > 0 else 1 - p_two / 2
    parts["H_trend_wild"] = {"slope": slope, "p_one_sided": p_one, "pass": slope > 0 and p_one < 0.05,
                             "mean_H": float(hbar.mean()), "sd_H": float(hbar.std())}
    return {"pass": all(p["pass"] for p in parts.values()), "parts": parts}


def p3(data: dict) -> dict:
    parts = {}
    for rule in ("wild", "sub"):
        s2, s1 = pooled(data, "S2", 2000, "raw", rule), pooled(data, "S1", 2000, "raw", rule)
        d, lo, hi = newcombe(s2.sum(), s2.size, s1.sum(), s1.size)
        parts[rule] = {"S2": pct(s2.sum(), s2.size), "S1": pct(s1.sum(), s1.size), "diff": d, "lo": lo, "hi": hi,
                       "pass": lo > 0}
    shift = p1(data, "S2")
    parts["shift_S2"] = shift
    return {"pass": all(p["pass"] for p in parts.values()), "parts": parts}


def p4(data: dict) -> dict:
    raw = pooled(data, "A", 2000, "raw", "wild")
    inc = pooled(data, "A", 2000, "inc1", "wild")
    n10, n01 = int(np.sum(raw & ~inc)), int(np.sum(~raw & inc))
    p_one = binomtest(n10, n10 + n01, 0.5, alternative="greater").pvalue if n10 + n01 else 1.0
    wild = {"raw": pct(raw.sum(), raw.size), "inc1": pct(inc.sum(), inc.size), "raw_only": n10, "inc1_only": n01,
            "p_one_sided": float(p_one), "pass": bool(raw.mean() > inc.mean() and p_one < 0.05)}
    gaps = {}
    for k in ("raw", "inc1"):
        tt, nn = [], []
        for a, b in R.PAIRS:
            d = data[f"A__{a}-{b}__n2000"]
            fx, fy = trend_flags(d)
            v = rej(d, k, "sub")
            tt.append(v[fx & fy])
            nn.append(v[~fx & ~fy])
        tt, nn = np.concatenate(tt), np.concatenate(nn)
        gaps[k] = {"trend_both": pct(tt.sum(), tt.size), "trend_neither": pct(nn.sum(), nn.size),
                   "gap": 100 * (tt.mean() - nn.mean())}
    trend = {**gaps, "paper_raw_gap": 3.5, "pass": gaps["inc1"]["gap"] < gaps["raw"]["gap"]}
    return {"pass": wild["pass"] and trend["pass"], "wild": wild, "trend": trend}


def p5(data: dict) -> dict:
    cells, ok = {}, True
    for arm in ("S1", "S2"):
        for a, b in R.PAIRS:
            for n in R.NS:
                d = data[f"{arm}__{a}-{b}__n{n}"]
                for k in ("raw", "inc1"):
                    v = rej(d, k, "harris")
                    c = pct(v.sum(), v.size)
                    cells[f"{arm}_{a}-{b}_n{n}_{k}"] = c
                    ok &= c["lo"] <= 5
    return {"pass": bool(ok), "cells": cells}


def summary(data: dict) -> list[dict]:
    out = []
    for name, d in data.items():
        for k, r, _, fn in RULES:
            v = np.asarray(fn(d), dtype=bool)
            out.append({"cell": name, "kernel": k, "rule": r, **pct(v.sum(), v.size)})
    return out


def pooled_table(data: dict) -> list[dict]:
    out = []
    for arm in ("A", "S1", "S2"):
        for n in R.NS:
            for k, r, _, fn in RULES:
                v = np.concatenate([np.asarray(fn(data[f"{arm}__{a}-{b}__n{n}"]), dtype=bool) for a, b in R.PAIRS])
                out.append({"arm": arm, "n": n, "kernel": k, "rule": r, **pct(v.sum(), v.size)})
    return out


def main():
    data = {c["name"]: load(c["name"]) for c in R.cells()}
    res = {"P0": p0(data), "P1": p1(data), "P2": p2(data), "P3": p3(data), "P4": p4(data), "P5": p5(data)}
    pt = pooled_table(data)
    res["pooled"] = pt
    res["window_H"] = {arm: {"mean_H_x": float(np.mean(np.concatenate([data[f"{arm}__{a}-{b}__n2000"]["H_x"] for a, b in R.PAIRS]))),
                             "mean_H_y": float(np.mean(np.concatenate([data[f"{arm}__{a}-{b}__n2000"]["H_y"] for a, b in R.PAIRS])))}
                       for arm in ("A", "S1", "S2")}
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    rows = summary(data)
    with open(os.path.join(HERE, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = ["# EXP8 pooled rejection rates (percent, Wilson 95%)", "",
             "| arm | n | kernel | rule | rate | interval |", "|---|---|---|---|---:|---|"]
    lines += [f"| {r['arm']} | {r['n']} | {r['kernel']} | {r['rule']} | {r['rate']:.1f} | [{r['lo']:.1f}, {r['hi']:.1f}] |"
              for r in pt]
    with open(os.path.join(HERE, "tables.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    for p in ("P0", "P1", "P2", "P3", "P4", "P5"):
        print(p, "PASS" if res[p]["pass"] else "FAIL")


if __name__ == "__main__":
    main()
