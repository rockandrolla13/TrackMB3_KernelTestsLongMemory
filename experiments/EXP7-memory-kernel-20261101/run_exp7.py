"""EXP7 driver: memory kernels, size, growth and power (see DESIGN.md).

Usage:
  python -I run_exp7.py --mode tune  --out tune  --workers 12
  python -I run_exp7.py --mode pilot --out pilot --workers 12
  python -I run_exp7.py --mode pilot --out pilot_hsplit --extras-only --workers 12
  python -I run_exp7.py --mode main  --out main  --workers 16

One CSV per cell in --out, one row per replication, appended as chunks finish. Re-running resumes:
a replication counts as done only with a complete status-ok row, so error rows are retried.
Seeds: SeedSequence(entropy=MASTER[mode], spawn_key=(cell_idx, rep)); one default_rng per replication.
The cell list is append-only, so adding cells keeps old seeds.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402

import exp7_core as K  # noqa: E402

MASTER = {"pilot": 2026110171, "main": 2026110172, "tune": 2026110173}
MAX_WORKERS = 16
RHO_LIN, RHO_SQ, RHO_STRESS = 0.065, 0.075, 0.13
GAMMA = 0.1  # tuning run, DESIGN.md amendment A1
SCALE_C = 1.1  # tuning run, DESIGN.md amendment A1
SHAPE_A = 0.1  # tuning run, DESIGN.md amendment A1
DELTAS = (0.1, 0.5)
H_PAIRS = ((0.5, 0.5), (0.6, 0.6), (0.7, 0.7), (0.8, 0.8), (0.9, 0.9), (0.95, 0.95), (0.6, 0.9), (0.7, 0.9))
PROCESS_ARMS = ("HT-S", "HT-L", "NL-SV", "NL-TAR", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow")
POWER_H = (0.5, 0.7, 0.8, 0.9)
MMD_H = (0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9)
MMD_POWER_H = (0.5, 0.7, 0.8)
GROWTH_N = (500, 1000, 2000, 4000, 8000)
GROWTH_H = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
CHUNK = {"ind": {500: 25, 2000: 5, 4000: 2}, "mmd": {500: 100, 2000: 25},
         "growth": {500: 200, 1000: 100, 2000: 100, 4000: 50, 8000: 25}}
EXACT_REPS = {"main": {500: 100, 2000: 100, 4000: 25}, "pilot": {500: 20, 2000: 20, 4000: 5}, "tune": {500: 0, 2000: 0, 4000: 0}}


def ind_cell(gen, alt, HX, HY, n, par=0.0, target=None):
    return {"family": "ind", "gen": gen, "alt": alt, "HX": HX, "HY": HY, "n": n, "param": par, "target": target}


def main_cells() -> list[dict]:
    """Append-only. The position in this list is the cell index in the seed."""
    cells = []
    for n, target in ((2000, 1000), (500, 1000), (4000, 300)):  # Gaussian memory cells
        for hx, hy in H_PAIRS:
            cells.append(ind_cell("G0", "null", hx, hy, n, target=target))
    for g in PROCESS_ARMS:  # process arms
        for h in (0.5, 0.8):
            cells.append(ind_cell(g, "null", h, h, 2000, target=1000))
    for g, par in (("FACTOR", 0.0), ("SCALE", GAMMA)):  # matched nulls of stress alternatives
        for h in POWER_H:
            cells.append(ind_cell(g, "null", h, h, 2000, par, 1000))
    for alt, par in (("lin", RHO_LIN), ("sq", RHO_SQ), ("lag", RHO_STRESS), ("flat", RHO_STRESS), ("factor", 0.0), ("scale", GAMMA)):
        for h in POWER_H:
            cells.append(ind_cell("G0", alt, h, h, 2000, par, 1000))
    for n in GROWTH_N:  # bare-statistic growth
        for h in GROWTH_H:
            cells.append({"family": "growth", "gen": "G0", "alt": "null", "HX": h, "HY": h, "n": n, "param": 0.0, "target": 2000})
    for g in ("G0", "HT-S"):  # split MMD null
        for n in (2000, 500):
            for h in MMD_H:
                cells.append({"family": "mmd", "gen": g, "alt": "null", "HX": h, "HY": h, "n": n, "param": 0.0, "target": 1000})
    for alt, par in (("mean", DELTAS[0]), ("mean", DELTAS[1]), ("scale", SCALE_C), ("shape", SHAPE_A)):  # split MMD power
        for h in MMD_POWER_H:
            cells.append({"family": "mmd", "gen": "G0", "alt": alt, "HX": h, "HY": h, "n": 2000, "param": par, "target": 1000})
    return finish(cells)


def tune_cells() -> list[dict]:
    cells = []
    h = 0.5
    for g in (0.05, 0.1, 0.15, 0.2, 0.3):
        cells.append(ind_cell("SCALE", "null", h, h, 2000, g, 200))
        cells.append(ind_cell("G0", "scale", h, h, 2000, g, 200))
    cells.append(ind_cell("G0", "null", h, h, 2000, 0.0, 200))
    cells.append(ind_cell("G0", "lag", h, h, 2000, RHO_STRESS, 200))
    cells.append(ind_cell("G0", "flat", h, h, 2000, RHO_STRESS, 200))
    cells.append({"family": "mmd", "gen": "G0", "alt": "null", "HX": h, "HY": h, "n": 2000, "param": 0.0, "target": 200})
    for c in (1.05, 1.1, 1.15, 1.2, 1.3):
        cells.append({"family": "mmd", "gen": "G0", "alt": "scale", "HX": h, "HY": h, "n": 2000, "param": c, "target": 200})
    for a in (0.1, 0.15, 0.2, 0.3, 0.4):
        cells.append({"family": "mmd", "gen": "G0", "alt": "shape", "HX": h, "HY": h, "n": 2000, "param": a, "target": 200})
    return finish(cells)


def finish(cells: list[dict]) -> list[dict]:
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = name_of(c)
    if len({c["name"] for c in cells}) != len(cells):
        raise ValueError("cell names are not unique")
    return cells


def name_of(c: dict) -> str:
    h = f"H{c['HX']}" if c["HX"] == c["HY"] else f"HX{c['HX']}_HY{c['HY']}"
    return f"{c['family']}__{c['gen']}__{c['alt']}{c['param']:g}__{h}__n{c['n']}"


def columns(c: dict) -> list[str]:
    head = ["rep", "status", "error", "cpu_s", "exact"]
    if c["family"] == "growth":
        return head + K.GROWTH_COLS
    if c["family"] == "mmd":
        return head + K.mmd_columns(K.MMD_KERNELS)
    return head + K.ind_columns(K.IND_KERNELS)


def run_chunk(args):
    cell, reps, master, tune, exact_reps, extras_only = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        exact = False
        try:
            if cell["family"] == "growth":
                x, y = K.ind_pair("G0", "null", cell["HX"], cell["HY"], cell["n"], 0.0, rng)
                out = K.growth_bundle(x, y)
            elif cell["family"] == "mmd":
                x = K.mmd_series(cell["gen"], cell["alt"], cell["HX"], cell["n"], cell["param"], rng)
                exact = r < exact_reps[cell["n"]] and not extras_only
                out = K.mmd_bundle(x, rng, kernels=("raw", "herm1", "herm2") if extras_only else K.MMD_KERNELS, exact=exact,
                                   hsplit=not tune)
            else:
                x, y = K.ind_pair(cell["gen"], cell["alt"], cell["HX"], cell["HY"], cell["n"], cell["param"], rng)
                exact = r < exact_reps[cell["n"]] and not extras_only
                kernels = ("raw", "inc1") if tune else ("raw", "inc1", "herm1", "herm2") if extras_only else K.IND_KERNELS
                out = K.ind_bundle(x, y, rng, kernels=kernels, HX=cell["HX"], HY=cell["HY"], exact=exact,
                                   full=not (tune or extras_only), hsplit=not tune, strict=(not tune) and r < K.HARRIS_STRICT_REPS[cell["n"]])
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4),
                     "exact": int(exact), **out})
    return cell["name"], rows


def fmt(key, v):
    """Twelve significant digits keep files small. Exact-path columns stay at full precision for the check."""
    if key.startswith("xraw_") or isinstance(v, (str, int)):
        return v
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.12g}"
    return v


def completed(path: str) -> set[int]:
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {int(row["rep"]) for row in csv.DictReader(f) if None not in row.values() and row["status"] == "ok"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=list(MASTER), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--reps", type=int, default=None, help="pilot: reps per cell (default 100)")
    ap.add_argument("--reps4000", type=int, default=None, help="main: override the n = 4000 target")
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    ap.add_argument("--extras-only", action="store_true", help="raw, inc1, herm1, herm2, low-layer tests and Harris rules, same seeds")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes (DESIGN.md)")
    master = MASTER[a.mode]
    cells = tune_cells() if a.mode == "tune" else main_cells()
    if a.extras_only:
        cells = [c for c in cells if c["family"] != "growth"]
    if a.only:
        cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
    outdir = os.path.join(HERE, a.out)
    os.makedirs(outdir, exist_ok=True)
    tune = a.mode == "tune"
    tasks, mcells = [], []
    for c in cells:
        target = c["target"]
        if a.mode == "pilot":
            target = a.reps if a.reps else 100
        elif a.mode == "main" and c["family"] == "ind" and c["n"] == 4000 and a.reps4000:
            target = a.reps4000
        path = os.path.join(outdir, c["name"] + ".csv")
        done = completed(path)
        todo = [r for r in range(target) if r not in done]
        mcells.append({**c, "target": target, "already_done": len(done),
                       "seed": {"entropy": master, "spawn_key": [c["idx"], "rep"]}})
        cs = CHUNK[c["family"]][c["n"]]
        for i in range(0, len(todo), cs):
            tasks.append((c, todo[i:i + cs], master, tune, todo[i] / target))
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    with open(os.path.join(outdir, f"manifest_{stamp}.json"), "w") as f:
        json.dump({"mode": a.mode, "master": master, "gamma": GAMMA, "scale_c": SCALE_C, "shape_a": SHAPE_A,
                   "seed_rule": "np.random.SeedSequence(entropy=master, spawn_key=(cell_idx, rep)); one default_rng per replication",
                   "B_wild": K.B_WILD, "tol_chol": K.TOL_CHOL, "tol_pair": K.TOL_PAIR, "exact_reps": EXACT_REPS[a.mode],
                   "numpy": np.__version__, "python": sys.version.split()[0], "workers": a.workers, "argv": sys.argv,
                   "cells": mcells}, f, indent=1)
    # every cell advances in proportion, so a stopped run is balanced
    tasks.sort(key=lambda t: (t[4], t[0]["idx"]))
    head = ["rep", "status", "error", "cpu_s", "exact"]

    def fields_of(c):
        if a.extras_only:
            if c["family"] == "mmd":
                return head + ["sigma", "rank"] + K.HSPLIT_MMD_COLS + [f"{key}_{k}" for k in ("raw", "herm1", "herm2") for key in ("mmd2", "p_wb", "r_wb")] \
                    + ["r975_wb_herm1", "r975_wb_herm2"]
            return head + ["sx", "sy", "rx_rank", "ry_rank"] + K.HSPLIT_IND_COLS + K.kernel_columns(("raw", "inc1", "herm1", "herm2"))
        if tune and c["family"] == "ind":
            return head + ["sx", "sy", "rx_rank", "ry_rank"] + K.kernel_columns(("raw", "inc1"))
        return columns(c)

    fields = {c["name"]: fields_of(c) for c in cells}
    print(f"{stamp} {len(tasks)} chunks, {sum(len(t[1]) for t in tasks)} replications, {len(cells)} cells", flush=True)
    t0 = time.time()
    nd = 0
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=200) as pool:
        for name, rows in pool.imap_unordered(run_chunk, [(*t[:4], EXACT_REPS[a.mode], a.extras_only) for t in tasks]):
            path = os.path.join(outdir, name + ".csv")
            new = not os.path.exists(path)
            with open(path, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields[name], extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows([{k: fmt(k, v) for k, v in r.items()} for r in rows])
            nd += 1
            if nd % 250 == 0 or nd == len(tasks):
                print(f"{nd}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print(f"DONE {time.time() - t0:.0f}s wall", flush=True)


if __name__ == "__main__":
    main()
