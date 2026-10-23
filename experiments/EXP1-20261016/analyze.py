"""Analysis for EXP1-20261016: H5, E1, E2, E3. Exploratory; nothing here is evidence for a claim.

Usage: python -I analyze.py [h5] [e1] [e2] [e3]     (default: all parts whose output folders exist)
Writes <part>_summary.csv and <part>_tables.md next to this script.
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KL = os.path.join(os.path.dirname(HERE), "JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear", "cells", "cell_summary.csv")
ALPHA = 0.05
Z = 1.959963984540054


def wilson(k, r):
    if r == 0:
        return float("nan"), float("nan")
    p = k / r
    den = 1 + Z * Z / r
    c = (p + Z * Z / (2 * r)) / den
    h = Z * math.sqrt(p * (1 - p) / r + Z * Z / (4 * r * r)) / den
    return c - h, c + h


def rows_of(path):
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


def pct(v):
    return "—" if v is None or not np.isfinite(v) else f"{100 * v:.1f}"


def fmt_ci(s, lo, hi):
    return f"{pct(s)} [{pct(lo)}, {pct(hi)}]"


def write_csv(path, recs):
    if not recs:
        return
    keys = list(recs[0].keys())
    for r in recs:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(recs)


def rate(flags):
    flags = [bool(f) for f in flags]
    r = len(flags)
    k = int(sum(flags))
    lo, hi = wilson(k, r)
    return {"R": r, "rej": k, "size": k / r if r else float("nan"), "ci_lo": lo, "ci_hi": hi}


# ======================================================================= H5

def h5():
    d = os.path.join(HERE, "h5")
    man = json.load(open(os.path.join(d, "manifest.json")))
    recs, L = [], ["# H5 — wild-bootstrap variants at H = 0.5 (MMD homogeneity, l = floor(sqrt n), B = 500)\n",
                   "Size at nominal 5% (%, Wilson 95% CI). Exploratory simulation.\n",
                   "| cell | reps ok / errors | old (mmd-dml) | csg (published) | iid | unc | ustat |", "|---|---|---|---|---|---|---|"]
    for c in man["cells"]:
        rs = rows_of(os.path.join(d, c["name"] + ".csv"))
        ok = [r for r in rs if r["status"] == "ok"]
        line = [c["name"], f"{len(ok)} / {len(rs) - len(ok)}"]
        for v in ("old", "csg", "iid", "unc", "ustat"):
            s = rate(float(r[f"p_{v}"]) <= ALPHA for r in ok)
            recs.append({"cell": c["name"], "process": c["process"], "n": c["n"], "variant": v, **s,
                         "errors": len(rs) - len(ok), "cpu_s": sum(float(r["cpu_s"]) for r in rs)})
            line.append(fmt_ci(s["size"], s["ci_lo"], s["ci_hi"]))
        L.append("| " + " | ".join(line) + " |")
    on = {}
    for v in ("old", "csg", "iid", "unc", "ustat"):
        on[v] = all(abs(r["size"] - ALPHA) <= 0.01 for r in recs if r["variant"] == v)
    L.append("\nRule (fixed before the run): a variant is on size if it is within 5 ± 1 pp in every H5 cell.\n")
    for v, ok in on.items():
        L.append(f"- {v}: {'on size' if ok else 'not on size'}")
    write_csv(os.path.join(HERE, "h5_summary.csv"), recs)
    open(os.path.join(HERE, "h5_tables.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


# ======================================================================= E2

E2_WILD = None  # discovered from CSV columns


def e2():
    d = os.path.join(HERE, "e2")
    man = json.load(open(os.path.join(d, "manifest.json")))
    data = {}
    for c in man["cells"]:
        rs = rows_of(os.path.join(d, c["name"] + ".csv"))
        data[c["name"]] = (c, rs)
    # i.i.d. SN critical value: 95% quantile of sn_T in G1 H=0.5 at the same n
    sn_crit = {}
    for name, (c, rs) in data.items():
        if c["process"] == "G1" and c["H"] == 0.5 and not c["shift"]:
            v = np.sort([float(r["sn_T"]) for r in rs if r["status"] == "ok"])
            sn_crit[c["n"]] = float(v[int(math.ceil(0.95 * v.size)) - 1]) if v.size else float("nan")
    recs = []
    cal_names = None
    for name, (c, rs) in data.items():
        ok = [r for r in rs if r["status"] == "ok"]
        if not rs:
            continue
        n, h = c["n"], c["H"]
        wild = sorted(k for k in rs[0] if k.startswith("p_") and rs[0][k] != "")
        cals = {}
        for k in wild:
            cals["wild_" + k[2:]] = [float(r[k]) <= ALPHA for r in ok if r[k] != ""]
        tn = np.array([float(r["t_n"]) for r in ok])
        q_sq = np.array([float(r["q95_sqrt"]) for r in ok])
        l_sq = np.array([float(r["l_sqrt"]) for r in ok])
        q_bd = np.array([float(r["q95_bd"]) for r in ok])
        hl = np.array([float(r["h_lw"]) for r in ok])
        cals["sub_sqrt"] = tn > q_sq
        cals["sub_bd"] = tn > q_bd
        cals["sub_rc_hat"] = tn > q_sq * (n / l_sq) ** (2 * np.maximum(hl, 0.5) - 1)
        cals["sub_rc_oracle_rank1"] = tn > q_sq * (n / l_sq) ** (2 * h - 1)
        cals["sub_rc_oracle_rank2"] = tn > q_sq * (n / l_sq) ** max(0.0, 4 * h - 3)
        snT = np.array([float(r["sn_T"]) for r in ok])
        cals["sn_block05"] = snT > np.array([float(r["sn_q95_05"]) for r in ok])
        cals["sn_block07"] = snT > np.array([float(r["sn_q95_07"]) for r in ok])
        cals["sn_iid"] = snT > sn_crit.get(n, float("nan"))
        cal_names = list(cals)
        for k, v in cals.items():
            recs.append({"cell": name, "process": c["process"], "H": h, "n": n, "power": bool(c["shift"]),
                         "calibration": k, **rate(v), "errors": len(rs) - len(ok),
                         "h_lw_mean": float(hl.mean()) if hl.size else float("nan"),
                         "h_gph_mean": float(np.mean([float(r["h_gph"]) for r in ok])) if ok else float("nan"),
                         "l_bd_median": float(np.median([float(r["l_bd"]) for r in ok])) if ok else float("nan"),
                         "cpu_s": sum(float(r["cpu_s"]) for r in rs)})
    write_csv(os.path.join(HERE, "e2_summary.csv"), recs)

    def get(p, h, n, cal, power=False):
        for r in recs:
            if r["process"] == p and r["H"] == h and r["n"] == n and r["calibration"] == cal and r["power"] == power:
                return r
        return None

    L = ["# E2 — MMD homogeneity: rank-2 processes and memory-aware calibrations\n",
         "Size at nominal 5% (%, Wilson 95% CI). Exploratory simulation. "
         "sn_iid uses the G1 H = 0.5 cell's own 95% quantile as critical value, so it is 5% there by construction.\n"]
    hs = sorted({r["H"] for r in recs})
    ns = sorted({r["n"] for r in recs})
    procs = [p for p in ("G1", "G3", "G5", "G6") if any(r["process"] == p for r in recs)]
    for p in procs:
        for n in ns:
            L.append(f"\n## {p}, n = {n}\n")
            L.append("| calibration | " + " | ".join(f"H={h}" for h in hs) + " |")
            L.append("|---|" + "---|" * len(hs))
            for cal in cal_names:
                cells = [get(p, h, n, cal) for h in hs]
                L.append(f"| {cal} | " + " | ".join("—" if r is None else fmt_ci(r["size"], r["ci_lo"], r["ci_hi"]) for r in cells) + " |")
            hl = [get(p, h, n, "sub_sqrt") for h in hs]
            L.append("| mean Ĥ (LW, m=n^0.65) | " + " | ".join("—" if r is None else f"{r['h_lw_mean']:.3f}" for r in hl) + " |")
            L.append("| mean Ĥ (GPH) | " + " | ".join("—" if r is None else f"{r['h_gph_mean']:.3f}" for r in hl) + " |")
            L.append("| median l_bd | " + " | ".join("—" if r is None else f"{r['l_bd_median']:.0f}" for r in hl) + " |")
    pw = [r for r in recs if r["power"]]
    if pw:
        L.append("\n## Power cells (H = 0.7, second half shifted by 0.3 sd)\n")
        L.append("| process | n | " + " | ".join(cal_names) + " |")
        L.append("|---|---|" + "---|" * len(cal_names))
        for p in ("G1", "G5"):
            for n in ns:
                L.append(f"| {p} | {n} | " + " | ".join(pct(get(p, 0.7, n, cal, True)["size"]) if get(p, 0.7, n, cal, True) else "—" for cal in cal_names) + " |")

    # ---- verdicts (rules fixed in parameters_fixed)
    nmax = max(ns)
    L.append(f"\n## Verdicts at n = {nmax}\n")
    L.append("**H1c rule:** G5/G6 hold size within ±2 pp of their own H = 0.5 level for H ≤ 0.7 and break (excess > 2 pp) for "
             "H ≥ 0.8, for sub_sqrt and wild unc (l = n^0.5; amended from csg after H5), old shown for reference.\n")
    verdict = {}
    for p in ("G5", "G6", "G1", "G3"):
        for cal in ("sub_sqrt", "wild_unc_05", "wild_old_05"):
            base = get(p, 0.5, nmax, cal)
            if base is None:
                continue
            low = [get(p, h, nmax, cal) for h in hs if h <= 0.7]
            high = [get(p, h, nmax, cal) for h in hs if h >= 0.8]
            hold = all(r is not None and abs(r["size"] - base["size"]) <= 0.02 for r in low)
            brk = all(r is not None and r["size"] - base["size"] > 0.02 for r in high)
            verdict[(p, cal)] = (hold, brk)
            L.append(f"- {p} / {cal}: holds to 0.7 = {hold}; breaks at ≥ 0.8 = {brk} "
                     f"(H=0.5 {pct(base['size'])}%; " + ", ".join(f"H={r['H']}: {pct(r['size'])}%" for r in low + high if r) + ")")
    L.append("\n**Fix rule:** a calibration is a fix if its size is within 5 ± 2 pp in every G1 and G3 cell with H ≤ 0.8 at n = "
             f"{nmax}.\n")
    for cal in cal_names:
        cells = [get(p, h, nmax, cal) for p in ("G1", "G3") for h in hs if h <= 0.8]
        cells = [c for c in cells if c]
        if not cells:
            continue
        ok = all(abs(c["size"] - ALPHA) <= 0.02 for c in cells)
        worst = max(cells, key=lambda c: abs(c["size"] - ALPHA))
        L.append(f"- {cal}: fix = {ok}; worst cell {worst['process']} H={worst['H']} at {pct(worst['size'])}%")
    open(os.path.join(HERE, "e2_tables.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[-40:]))


# ======================================================================= E1

E1_TESTS = {"hsic_wb": ("p_hsic_wb", "p"), "hsic_shift": ("p_hsic_shift", "p"), "hsic_sub": ("rej_hsic_sub", "rej"),
            "dcov_wb": ("p_dcov_wb", "p"), "dcov_sub": ("rej_dcov_sub", "rej"),
            "pearson_nw": ("rej_nw", "rej"), "pearson_ewc": ("rej_ewc", "rej")}


def cluster_ci(flags, groups):
    f = np.asarray(flags, float)
    g = np.asarray(groups)
    nobs = f.size
    m = f.mean()
    s = 0.0
    for v in np.unique(g):
        s += (f[g == v] - m).sum() ** 2
    se = math.sqrt(s) / nobs
    return m - Z * se, m + Z * se, se


def e1():
    d = os.path.join(HERE, "e1", "cells")
    man = json.load(open(os.path.join(d, "manifest.json")))
    recs = []
    L = ["# E1 — real-data placebo pairs (LQD/HYG/JNK, Arca, 2023-05-30..2025-05-29)\n",
         "Windows start at least 20 trading days apart; process independence is not established. Placebo rejection at nominal 5% (%). "
         "Wilson 95% CI, then the day-clustered 95% CI (wider of clustering by A-day and by B-day). Exploratory.\n"]
    for c in man["cells"]:
        rs = rows_of(os.path.join(d, c["name"] + ".csv"))
        ok = [r for r in rs if r["status"] == "ok"]
        if not ok:
            continue
        da = [int(r["day_a"]) for r in ok]
        db = [int(r["day_b"]) for r in ok]
        mem = {}
        for side in ("x", "y"):
            for est in ("lw", "gph"):
                for bw in ("050", "065", "080"):
                    v = np.array([float(r[f"{side}_{est}{bw}"]) for r in ok])
                    mem[f"{side}_{est}{bw}_median"] = float(np.median(v))
                    mem[f"{side}_{est}{bw}_p10"] = float(np.quantile(v, 0.1))
                    mem[f"{side}_{est}{bw}_p90"] = float(np.quantile(v, 0.9))
        sumh = np.array([float(r["x_lw065"]) + float(r["y_lw065"]) for r in ok])
        for t, (col, kind) in E1_TESTS.items():
            fl = [(float(r[col]) <= ALPHA) if kind == "p" else (int(float(r[col])) == 1) for r in ok if r[col] not in ("", "nan", "-1")]
            s = rate(fl)
            if len(fl) == len(ok):
                a = cluster_ci(fl, da)
                b = cluster_ci(fl, db)
                clo, chi = (a[0], a[1]) if a[2] >= b[2] else (b[0], b[1])
                fa = np.array(fl, float)
                q = np.quantile(sumh, [1 / 3, 2 / 3])
                terc = [float(fa[sumh <= q[0]].mean()), float(fa[(sumh > q[0]) & (sumh <= q[1])].mean()), float(fa[sumh > q[1]].mean())]
            else:
                clo = chi = float("nan")
                terc = [float("nan")] * 3
            recs.append({"cell": c["name"], "A": c["A"], "B": c["B"], "series": c["series"], "n": c["n"], "test": t, **s,
                         "clu_lo": clo, "clu_hi": chi, "errors": len(rs) - len(ok),
                         "size_by_sumH_tercile_low": terc[0], "size_by_sumH_tercile_mid": terc[1], "size_by_sumH_tercile_high": terc[2],
                         "sumH_lw065_median": float(np.median(sumh)), **mem,
                         "cpu_s": sum(float(r["cpu_s"]) for r in rs)})
    write_csv(os.path.join(HERE, "e1_summary.csv"), recs)

    def get(a, b, s, n, t):
        return next((r for r in recs if r["A"] == a and r["B"] == b and r["series"] == s and r["n"] == n and r["test"] == t), None)

    pairings = [(c["A"], c["B"]) for c in man["cells"]]
    pairings = list(dict.fromkeys(pairings))
    for s in ("sign", "dur"):
        L.append(f"\n## Series: {'trade signs' if s == 'sign' else 'adjusted log durations'}\n")
        L.append("| pairing | n | errors | median Ĥ_X+Ĥ_Y (LW n^0.65) | " + " | ".join(E1_TESTS) + " |")
        L.append("|---|---|---|---|" + "---|" * len(E1_TESTS))
        for a, b in pairings:
            for n in (500, 1000, 2000):
                rr = [get(a, b, s, n, t) for t in E1_TESTS]
                if rr[0] is None:
                    continue
                L.append(f"| {a}–{b} | {n} | {rr[0]['errors']} | {rr[0]['sumH_lw065_median']:.2f} | " +
                         " | ".join(f"{pct(r['size'])} [{pct(r['ci_lo'])}, {pct(r['ci_hi'])}] / [{pct(r['clu_lo'])}, {pct(r['clu_hi'])}]" for r in rr) + " |")
        L.append("\nSize by tercile of Ĥ_X + Ĥ_Y (LW, m = n^0.65) within the cell, at n = 2000 (low / mid / high, %):\n")
        L.append("| pairing | " + " | ".join(E1_TESTS) + " |")
        L.append("|---|" + "---|" * len(E1_TESTS))
        for a, b in pairings:
            rr = [get(a, b, s, 2000, t) for t in E1_TESTS]
            if rr[0] is None:
                continue
            L.append(f"| {a}–{b} | " + " | ".join(f"{pct(r['size_by_sumH_tercile_low'])} / {pct(r['size_by_sumH_tercile_mid'])} / {pct(r['size_by_sumH_tercile_high'])}" for r in rr) + " |")
        L.append("\nPer-window memory estimates (median [p10, p90]) of the X window, by estimator and bandwidth m = n^e:\n")
        L.append("| pairing | n | LW e=0.5 | LW 0.65 | LW 0.8 | GPH 0.5 | GPH 0.65 | GPH 0.8 |")
        L.append("|---|---|---|---|---|---|---|---|")
        for a, b in pairings:
            for n in (500, 1000, 2000):
                r = get(a, b, s, n, "hsic_wb")
                if r is None:
                    continue
                L.append(f"| {a} (vs {b}) | {n} | " + " | ".join(
                    f"{r[f'x_{e}{bw}_median']:.2f} [{r[f'x_{e}{bw}_p10']:.2f}, {r[f'x_{e}{bw}_p90']:.2f}]"
                    for e in ("lw", "gph") for bw in ("050", "065", "080")) + " |")
    ml = os.path.join(HERE, "e1", "memory_long.json")
    if os.path.exists(ml):
        mj = json.load(open(ml))
        L.append("\n## Long-sample memory per fund (segments of 32,768 events; median across segments [p10, p90])\n")
        L.append("| fund | series | segments | LW 0.5 | LW 0.65 | LW 0.8 | GPH 0.5 | GPH 0.65 | GPH 0.8 | band (range of medians) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for f, v in mj.items():
            for s, r in v.items():
                L.append(f"| {f} | {s} | {r['segments']} | " + " | ".join(
                    f"{r[k]['median']:.2f} [{r[k]['p10']:.2f}, {r[k]['p90']:.2f}]" for k in
                    ("lw_0.5", "lw_0.65", "lw_0.8", "gph_0.5", "gph_0.65", "gph_0.8")) + f" | {r['band'][0]:.2f}–{r['band'][1]:.2f} |")
    open(os.path.join(HERE, "e1_tables.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:60]))


# ======================================================================= E3

E3_TESTS = E1_TESTS


def e3():
    d = os.path.join(HERE, "e3")
    man = json.load(open(os.path.join(d, "manifest.json")))
    recs = []
    for c in man["cells"]:
        rs = rows_of(os.path.join(d, c["name"] + ".csv"))
        ok = [r for r in rs if r["status"] == "ok"]
        for t, (col, kind) in E3_TESTS.items():
            fl = [(float(r[col]) <= ALPHA) if kind == "p" else (int(float(r[col])) == 1) for r in ok if r[col] not in ("", "nan", "-1")]
            recs.append({"cell": c["name"], "HX": c["HX"], "HY": c["HY"], "sumH": round(c["HX"] + c["HY"], 4), "n": c["n"],
                         "test": t, **rate(fl), "reps_target": c["reps"], "reps_done": len(rs), "errors": len(rs) - len(ok),
                         "cpu_s": sum(float(r["cpu_s"]) for r in rs)})
    write_csv(os.path.join(HERE, "e3_summary.csv"), recs)
    kl = list(csv.DictReader(open(KL))) if os.path.exists(KL) else []

    def kl_size(proc, hx, hy, n, t):
        for r in kl:
            if r["proc"] == proc and float(r["HX"]) == hx and float(r["HY"]) == hy and int(r["n"]) == n and r["test"] == t and float(r["r"] or 0) == 0:
                return float(r["size"])
        return None

    L = ["# E3 — fine grid at the cross-rule line (independent fGn pairs)\n",
         "Size at nominal 5% (%, Wilson 95% CI). Exploratory simulation. Same code as the kernel-linear arm.\n"]
    pairs = list(dict.fromkeys((c["HX"], c["HY"]) for c in man["cells"]))
    for n in (500, 1000, 2000):
        L.append(f"\n## n = {n}\n")
        L.append("| (H_X, H_Y) | sum | " + " | ".join(E3_TESTS) + " |")
        L.append("|---|---|" + "---|" * len(E3_TESTS))
        for hx, hy in pairs:
            rr = [next((r for r in recs if r["HX"] == hx and r["HY"] == hy and r["n"] == n and r["test"] == t), None) for t in E3_TESTS]
            L.append(f"| ({hx}, {hy}) | {hx + hy:.3f} | " + " | ".join("—" if r is None else fmt_ci(r["size"], r["ci_lo"], r["ci_hi"]) for r in rr) + " |")
        L.append("| KL reference H = 0.5 (G1) | 1.0 | " + " | ".join(pct(kl_size("G1", 0.5, 0.5, n, t)) if kl_size("G1", 0.5, 0.5, n, t) is not None else "—" for t in E3_TESTS) + " |")
    # reads
    L.append("\n## Reads (n = 2000)\n")
    L.append("KC2 read: excess over 5% greater than 2 pp at sums ≤ 1.4. Also shown: excess over the test's own KL H = 0.5 size "
             "(removes the memory-free wild-bootstrap oversize).\n")
    for t in E3_TESTS:
        low = [r for r in recs if r["test"] == t and r["n"] == 2000 and r["sumH"] <= 1.4 + 1e-9]
        base = kl_size("G1", 0.5, 0.5, 2000, t)
        L.append(f"- {t}: " + "; ".join(f"({r['HX']}, {r['HY']}) {pct(r['size'])}% (excess {100 * (r['size'] - ALPHA):+.1f} pp"
                                          + (f", vs own H=0.5 {100 * (r['size'] - base):+.1f} pp" if base is not None else "") + ")" for r in low))
    L.append("\nUnequal vs equal pairs at similar sums (n = 2000): (0.55, 0.85) sum 1.40 vs (0.7, 0.7); (0.6, 0.9) sum 1.50 vs (0.75, 0.75); "
             "(0.65, 0.9) sum 1.55 vs (0.775, 0.775).\n")
    for t in E3_TESTS:
        def sz(hx, hy):
            r = next((r for r in recs if r["HX"] == hx and r["HY"] == hy and r["n"] == 2000 and r["test"] == t), None)
            return r
        comps = [((0.55, 0.85), (0.7, 0.7)), ((0.6, 0.9), (0.75, 0.75)), ((0.65, 0.9), (0.775, 0.775))]
        L.append(f"- {t}: " + "; ".join(
            f"{u} {fmt_ci(sz(*u)['size'], sz(*u)['ci_lo'], sz(*u)['ci_hi'])} vs {e} {fmt_ci(sz(*e)['size'], sz(*e)['ci_lo'], sz(*e)['ci_hi'])}"
            for u, e in comps if sz(*u) and sz(*e)))
    open(os.path.join(HERE, "e3_tables.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    parts = sys.argv[1:] or [p for p in ("h5", "e1", "e2", "e3") if os.path.isdir(os.path.join(HERE, p if p != "e1" else "e1/cells"))]
    for p in parts:
        {"h5": h5, "e1": e1, "e2": e2, "e3": e3}[p]()
