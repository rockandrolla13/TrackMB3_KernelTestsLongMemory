"""E3: finer (H_X, H_Y) grid at the cross-rule line H_X + H_Y = 3/2 (fGn pairs, independent).

Same tests and code as the kernel-linear arm (mbcore.kind_tests, mbcore.lind_tests).

Usage:
  python -I run_e3.py --mode pilot --reps 5 --out e3_pilot
  python -I run_e3.py --mode main  --reps 2000 --out e3
Seeds: replication r of cell c uses SeedSequence(entropy=MASTER[mode], spawn_key=(c, r)).
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import mbcore as mb  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASTER = {"pilot": 2026101631, "main": 2026101632}
PAIRS = ((0.65, 0.65), (0.7, 0.7), (0.725, 0.725), (0.75, 0.75), (0.775, 0.775), (0.8, 0.8), (0.85, 0.85),
         (0.55, 0.85), (0.6, 0.9), (0.65, 0.9))
N_GRID = (500, 1000, 2000)
B_WB = 500
CHUNK = {500: 25, 1000: 10, 2000: 5}
KIND = ["hsic_stat", "dcov_stat", "p_hsic_wb", "p_hsic_shift", "p_hsic_sub", "rej_hsic_sub",
        "p_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "shift_check"]
LIND = ["pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw", "ewc_nu", "t_ewc", "rej_ewc", "nw_lrv_nonpos"]


def cell_list():
    cells = []
    for hx, hy in PAIRS:
        for n in N_GRID:
            cells.append({"proc": "G4", "HX": hx, "HY": hy, "lam": None, "n": n})
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = f"e3__HX{c['HX']}_HY{c['HY']}__n{c['n']}"
    return cells


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            x, y = mb.gen_series(cell, cell["n"], rng)
            out = mb.kind_tests(x, y, rng, B=B_WB)
            out.update(mb.lind_tests(x, y))
        except Exception as e:  # recorded in the row
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    if a.workers > 6:
        sys.exit("at most 6 worker processes")
    outdir = os.path.join(HERE, a.out)
    os.makedirs(outdir, exist_ok=True)
    master = MASTER[a.mode]
    cells = cell_list()
    if a.only:
        cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
    tasks, mcells = [], []
    for c in cells:
        path = os.path.join(outdir, c["name"] + ".csv")
        done = set()
        if os.path.exists(path):
            with open(path) as f:
                done = {int(r["rep"]) for r in csv.DictReader(f)}
        todo = [r for r in range(a.reps) if r not in done]
        mcells.append({**c, "reps": a.reps, "seed": {"entropy": master, "spawn_key": [c["idx"], "rep"]}})
        for i in range(0, len(todo), CHUNK[c["n"]]):
            tasks.append((c, todo[i:i + CHUNK[c["n"]]], master))
    mname = "manifest.json" if not a.only else "manifest_partial_" + time.strftime("%Y%m%dT%H%M%S", time.gmtime()) + ".json"
    with open(os.path.join(outdir, mname), "w") as f:
        json.dump({"mode": a.mode, "master_entropy": master, "B_wild": B_WB, "numpy": np.__version__,
                   "python": sys.version.split()[0], "workers": a.workers, "cells": mcells}, f, indent=1)
    tasks.sort(key=lambda t: -t[0]["n"])
    print(f"{len(tasks)} chunks, {sum(len(t[1]) for t in tasks)} replications", flush=True)
    fields = ["rep", "status", "error", "cpu_s"] + KIND + LIND
    t0 = time.time()
    nd = 0
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=50) as pool:
        for name, rows in pool.imap_unordered(run_chunk, tasks):
            path = os.path.join(outdir, name + ".csv")
            new = not os.path.exists(path)
            with open(path, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows(rows)
            nd += 1
            if nd % 20 == 0 or nd == len(tasks):
                print(f"{nd}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
