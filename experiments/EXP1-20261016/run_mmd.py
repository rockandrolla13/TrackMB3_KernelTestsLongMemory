"""Driver for H5 (wild-bootstrap hygiene) and E2 (rank-2 MMD and memory-aware calibration).

Usage:
  python -I run_mmd.py --exp h5 --mode pilot --reps 20 --out h5_pilot
  python -I run_mmd.py --exp h5 --mode main  --reps 2000 --out h5
  python -I run_mmd.py --exp e2 --mode pilot --reps 20 --out e2_pilot
  python -I run_mmd.py --exp e2 --mode main  --reps 2000 --reps-n2000 1000 --power-reps 500 --out e2

One CSV per cell in --out, one row per replication, appended as chunks finish. Re-running resumes.
Seeds: replication r of cell c uses SeedSequence(entropy=MASTER[(exp, mode)], spawn_key=(c, r)).
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

import mmdcal as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASTER = {("h5", "pilot"): 2026101651, ("h5", "main"): 2026101652,
          ("e2", "pilot"): 2026101621, ("e2", "main"): 2026101622}
H_E2 = (0.5, 0.55, 0.6, 0.7, 0.75, 0.8, 0.9)
N_GRID = (500, 1000, 2000)
POWER_SHIFT_SD = 0.3
SD = {"G1": 1.0, "G3": 1.0, "G5": 2 ** 0.5, "G6": (1 - 2 / np.pi) ** 0.5}
CHUNK = {500: 50, 1000: 20, 2000: 10}


def cell_list(exp: str) -> list[dict]:
    cells = []
    if exp == "h5":
        for proc in ("G1", "G3"):
            for n in (500, 2000):
                cells.append({"process": proc, "H": 0.5, "n": n, "shift": 0.0})
    else:
        for proc in ("G1", "G3", "G5", "G6"):
            for h in H_E2:
                for n in N_GRID:
                    cells.append({"process": proc, "H": h, "n": n, "shift": 0.0})
        for proc in ("G1", "G5"):
            for n in N_GRID:
                cells.append({"process": proc, "H": 0.7, "n": n, "shift": POWER_SHIFT_SD * SD[proc]})
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = f"{exp}__{c['process']}_H{c['H']}__n{c['n']}" + ("__power" if c["shift"] else "")
    return cells


def run_chunk(args):
    exp, cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            x = M.generate(cell["process"], cell["H"], cell["n"], rng)
            x[cell["n"] // 2:] += cell["shift"]
            out = M.rep_h5(x, rng) if exp == "h5" else M.rep_e2(x, rng)
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


FIELDS = {
    "h5": ["mmd2", "p_old", "p_csg", "p_iid", "p_unc", "p_ustat"],
    "e2": ["mmd2", "sigma", "t_n"] + [f"p_{f}_{int(e * 10):02d}" for f, e in M.WILD_E2]
          + ["h_lw", "h_gph", "l_sqrt", "q95_sqrt", "l_bd", "q95_bd", "sn_T", "sn_q95_05", "sn_q95_07"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", choices=["h5", "e2"], required=True)
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--power-reps", type=int, default=None)
    ap.add_argument("--reps-n2000", type=int, default=None, help="null reps for n = 2000 cells (fallback rule)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    if not os.path.exists(os.path.join(HERE, "validate_mmdcal.PASS")):
        sys.exit("validate_mmdcal.PASS missing")
    if a.workers > 6:
        sys.exit("at most 6 worker processes")
    outdir = os.path.join(HERE, a.out)
    os.makedirs(outdir, exist_ok=True)
    master = MASTER[(a.exp, a.mode)]
    cells = cell_list(a.exp)
    if a.only:
        cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
    tasks, mcells = [], []
    for c in cells:
        target = a.power_reps if (c["shift"] and a.power_reps is not None) else a.reps
        if not c["shift"] and c["n"] == 2000 and a.reps_n2000 is not None:
            target = a.reps_n2000
        path = os.path.join(outdir, c["name"] + ".csv")
        done = set()
        if os.path.exists(path):
            with open(path) as f:
                done = {int(row["rep"]) for row in csv.DictReader(f)}
        todo = [r for r in range(target) if r not in done]
        mcells.append({**c, "reps": target, "seed": {"entropy": master, "spawn_key": [c["idx"], "rep"]}})
        cs = CHUNK[c["n"]]
        for i in range(0, len(todo), cs):
            tasks.append((a.exp, c, todo[i:i + cs], master))
    mname = "manifest.json" if not a.only else "manifest_partial_" + time.strftime("%Y%m%dT%H%M%S", time.gmtime()) + ".json"
    with open(os.path.join(outdir, mname), "w") as f:
        json.dump({"exp": a.exp, "mode": a.mode, "master_entropy": master,
                   "seed_rule": "np.random.SeedSequence(entropy=master_entropy, spawn_key=(cell_idx, rep)); one default_rng per replication",
                   "B_wild": M.B_WILD, "numpy": np.__version__, "python": sys.version.split()[0],
                   "workers": a.workers, "cells": mcells}, f, indent=1)
    tasks.sort(key=lambda t: -t[1]["n"])
    print(f"{len(tasks)} chunks, {sum(len(t[2]) for t in tasks)} replications", flush=True)
    fields = ["rep", "status", "error", "cpu_s"] + FIELDS[a.exp]
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
            if nd % 10 == 0 or nd == len(tasks):
                print(f"{nd}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
