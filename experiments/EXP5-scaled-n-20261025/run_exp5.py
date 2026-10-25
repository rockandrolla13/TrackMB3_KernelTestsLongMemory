"""EXP5 driver (see DESIGN.md). Part 1: null-SD-scaled split-MMD shift. Part 2: independence
null at larger n for memory sums 1.0, 1.4, 1.6.

Usage:
  python -I run_exp5.py --mode pilot --reps 50 --out pilot --workers 3 --only mmd__
  python -I run_exp5.py --mode main  --reps 1000 --out main --workers 3 --only ind__G4__H0.8__n4000

One CSV per cell in --out, one row per replication, appended as chunks finish. Re-running resumes
(rows already on disk are skipped). Seeds: replication r of cell c uses
SeedSequence(entropy=MASTER[mode], spawn_key=(cell_idx, rep)); one default_rng per replication.
Cell indices are fixed by the cell list below (Part 2 starts at index 100).
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import exp5_core as P  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASTER = {"pilot": 2026102551, "main": 2026102552}
H_GRID = (0.5, 0.7, 0.8)
N_MMD = (2000, 500)
N_IND = (4000, 8000)
C_GRID = (5 ** 0.5, 1.5)  # c = sqrt(5) gives delta = 0.1 at H = 0.5, n = 2000 (EXP3's delta)
B_WB = 500
CHUNK = {500: 50, 2000: 10, 4000: 2, 8000: 1}
MAX_WORKERS = 6

MMD = ["mmd2", "sigma", "t_n", "hm_diff", "t_naive", "p_unc_05", "r_unc_05", "l_sqrt", "q95_sqrt",
       "r_sub_sqrt", "sn_T", "sn_q95_05", "r_sn_05"]
IND = ["hsic_stat", "dcov_stat", "p_hsic_wb", "p_hsic_shift", "p_hsic_sub", "rej_hsic_sub",
       "p_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "shift_check", "pearson_r", "nw_bw", "t_nw", "p_nw",
       "rej_nw", "ewc_nu", "t_ewc", "rej_ewc", "nw_lrv_nonpos"]
FIELDS = {"mmd": MMD, "ind": IND}


def cell_list() -> list[dict]:
    cells = []
    for n in N_MMD:
        for h in H_GRID:
            cells.append({"family": "mmd", "alt": "null", "c": 0.0, "H": h, "n": n, "delta": 0.0})
    for c in C_GRID:
        for n in N_MMD:
            for h in H_GRID:
                cells.append({"family": "mmd", "alt": "shift", "c": round(c, 4), "H": h, "n": n,
                              "delta": round(c * P.sd_half(h, n), 6)})
    for i, c in enumerate(cells):
        c["idx"] = i
    k = 100
    for n in N_IND:
        for h in H_GRID:
            cells.append({"family": "ind", "alt": "null", "H": h, "n": n, "sumH": round(2 * h, 3), "idx": k})
            k += 1
    for c in cells:
        c["name"] = name_of(c)
    return cells


def name_of(c: dict) -> str:
    if c["family"] == "mmd":
        tag = "null" if c["alt"] == "null" else f"shiftc{c['c']:g}"
        return f"mmd__G0__{tag}__H{c['H']}__n{c['n']}"
    return f"ind__G4__H{c['H']}__n{c['n']}"


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        w0 = time.time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            if cell["family"] == "mmd":
                x = P.mmd_series(cell["H"], cell["n"], cell["delta"], rng)
                out = P.mmd_ext(x, rng)
            else:
                out = P.ind_rep(cell["H"], cell["n"], rng, B=B_WB)
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4),
                     "wall_s": round(time.time() - w0, 3),
                     "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1), **out})
    return cell["name"], cell["family"], rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes (shared machine)")
    master = MASTER[a.mode]
    cells = cell_list()
    if a.only:
        cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
    outdir = os.path.join(HERE, a.out)
    os.makedirs(outdir, exist_ok=True)
    tasks, mcells = [], []
    for c in cells:
        path = os.path.join(outdir, c["name"] + ".csv")
        done = set()
        if os.path.exists(path):
            with open(path) as f:
                done = {int(row["rep"]) for row in csv.DictReader(f)}
        todo = [r for r in range(a.reps) if r not in done]
        mcells.append({**c, "reps": a.reps, "already_done": len(done),
                       "seed": {"entropy": master, "spawn_key": [c["idx"], "rep"]}})
        cs = CHUNK[c["n"]]
        for i in range(0, len(todo), cs):
            tasks.append((c, todo[i:i + cs], master))
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    with open(os.path.join(outdir, f"manifest_{stamp}.json"), "w") as f:
        json.dump({"mode": a.mode, "master": master,
                   "seed_rule": "np.random.SeedSequence(entropy=master, spawn_key=(cell_idx, rep)); one default_rng per replication",
                   "B_wild": B_WB, "numpy": np.__version__, "python": sys.version.split()[0],
                   "workers": a.workers, "argv": sys.argv, "cells": mcells}, f, indent=1)
    tasks.sort(key=lambda t: -t[0]["n"])
    print(f"{stamp} {len(tasks)} chunks, {sum(len(t[1]) for t in tasks)} replications", flush=True)
    t0 = time.time()
    nd = 0
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=20) as pool:
        for name, fam, rows in pool.imap_unordered(run_chunk, tasks):
            fields = ["rep", "status", "error", "cpu_s", "wall_s", "maxrss_mb"] + FIELDS[fam]
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
    print(f"DONE {time.time() - t0:.0f}s wall", flush=True)


if __name__ == "__main__":
    main()
