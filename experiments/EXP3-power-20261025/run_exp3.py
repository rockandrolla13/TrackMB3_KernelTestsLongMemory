"""EXP3 driver: null and alternative cells for size-adjusted power (see DESIGN.md).

Usage:
  python -I run_exp3.py --mode pilot --reps 100  --out pilot --workers 10
  python -I run_exp3.py --mode main  --reps 1000 --out main  --workers 10

One CSV per cell in --out, one row per replication, appended as chunks finish. Re-running resumes
(rows already on disk are skipped). Seeds: replication r of cell c uses
SeedSequence(entropy=MASTER[mode], spawn_key=(c, r)); one default_rng per replication.
Cell indices are fixed by the cell list below, so adding cells at the end keeps old seeds.
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
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))  # kmem (D7: not installed)
import numpy as np  # noqa: E402

import power_core as P  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASTER = {"pilot": 2026102541, "main": 2026102542}
H_GRID = (0.5, 0.7, 0.8)
N_GRID = (2000, 500)
B_WB = 500
CHUNK = {500: 25, 2000: 5}
MAX_WORKERS = 10

# Effect sizes at n = 2000, fixed by the pilot (DESIGN.md amendment A1). n = 500 uses 2x
# (local-alternative scaling, effect proportional to n^{-1/2}).
RHO_LIN = 0.065
RHO_SQ = 0.075
DELTA = 0.10
PILOT_GRID = {"lin": (0.04, 0.06, 0.08), "sq": (0.06, 0.09, 0.12), "shift": (0.08, 0.12, 0.16)}

KIND = ["hsic_stat", "dcov_stat", "p_hsic_wb", "r_hsic_wb", "p_hsic_shift", "r_hsic_shift",
        "p_hsic_sub", "rej_hsic_sub", "r_hsic_sub", "p_dcov_wb", "r_dcov_wb", "p_dcov_sub",
        "rej_dcov_sub", "r_dcov_sub", "shift_check", "pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw",
        "r_nw", "ewc_nu", "t_ewc", "rej_ewc", "r_ewc", "nw_lrv_nonpos"]
MMD = ["mmd2", "sigma", "t_n", "p_unc_05", "r_unc_05", "l_sqrt", "q95_sqrt", "r_sub_sqrt",
       "sn_T", "sn_q95_05", "r_sn_05"]
FIELDS = {"mmd": MMD, "ind": KIND}


def scale(param: float, n: int) -> float:
    return param * (2000 / n) ** 0.5


def main_cells() -> list[dict]:
    spec = [("ind", "G0", "null", 0.0), ("ind", "SEAS-V", "null", 0.0), ("mmd", "G0", "null", 0.0),
            ("ind", "G0", "lin", RHO_LIN), ("ind", "G0", "sq", RHO_SQ), ("ind", "SEAS-V", "lin", RHO_LIN),
            ("mmd", "G0", "shift", DELTA)]
    cells = []
    for fam, g, alt, par in spec:
        for n in N_GRID:
            for h in H_GRID:
                cells.append({"family": fam, "gen": g, "alt": alt, "H": h, "n": n,
                              "param": round(scale(par, n), 6), "param_n2000": par})
    return cells


def pilot_cells() -> list[dict]:
    cells = []
    for alt, grid in PILOT_GRID.items():
        fam = "mmd" if alt == "shift" else "ind"
        for par in grid:
            cells.append({"family": fam, "gen": "G0", "alt": alt, "H": 0.5, "n": 2000,
                          "param": par, "param_n2000": par})
    # timing and sanity for the seasonal alternative and for n = 500
    cells.append({"family": "ind", "gen": "SEAS-V", "alt": "lin", "H": 0.5, "n": 2000, "param": 0.06,
                  "param_n2000": 0.06})
    cells.append({"family": "ind", "gen": "G0", "alt": "lin", "H": 0.5, "n": 500, "param": 0.12,
                  "param_n2000": 0.06})
    return cells


def name_of(c: dict) -> str:
    return f"{c['family']}__{c['gen']}__{c['alt']}{c['param']:g}__H{c['H']}__n{c['n']}"


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            if cell["family"] == "mmd":
                x = P.mmd_series(cell["gen"], cell["alt"], cell["H"], cell["n"], cell["param"], rng)
                out = P.mmd_ext(x, rng)
            else:
                x, y = P.ind_pair(cell["gen"], cell["alt"], cell["H"], cell["n"], cell["param"], rng)
                out = P.kind_ext(x, y, rng, B=B_WB)
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], cell["family"], rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes (shared machine)")
    master = MASTER[a.mode]
    cells = pilot_cells() if a.mode == "pilot" else main_cells()
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = name_of(c)
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
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=50) as pool:
        for name, fam, rows in pool.imap_unordered(run_chunk, tasks):
            fields = ["rep", "status", "error", "cpu_s"] + FIELDS[fam]
            path = os.path.join(outdir, name + ".csv")
            new = not os.path.exists(path)
            with open(path, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows(rows)
            nd += 1
            if nd % 25 == 0 or nd == len(tasks):
                print(f"{nd}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print(f"DONE {time.time() - t0:.0f}s wall", flush=True)


if __name__ == "__main__":
    main()
