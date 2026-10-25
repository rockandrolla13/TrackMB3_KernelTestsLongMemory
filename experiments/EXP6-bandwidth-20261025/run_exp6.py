"""EXP6 driver: bandwidth sensitivity of the HSIC rules (see DESIGN.md).

Usage:
  python -I run_exp6.py --mode pilot   --reps 40   --out pilot   --workers 2
  python -I run_exp6.py --mode main    --reps 1000 --out main    --workers 5
  python -I run_exp6.py --mode regress --reps 5    --out regress --workers 2

One CSV per cell in --out, one row per replication, appended as chunks finish. Re-running resumes.
Seeds: SeedSequence(entropy=MASTER[mode], spawn_key=(cell_idx, rep)); one default_rng per rep.
Regress mode uses the EXP3 main master and EXP3 cell indices, for the kappa = 1 check only.
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402

import bw_core as K  # noqa: E402

MASTER = {"pilot": 2026102561, "main": 2026102562, "regress": 2026102542}
B_WB = 500
CHUNK = {500: 20, 2000: 4}
MAX_WORKERS = 5
RHO_LIN, RHO_SQ = 0.065, 0.075  # EXP3 values at n = 2000; doubled at n = 500

# (gen, alt, rho at n = 2000, H); index in this list is the cell index (main, pilot)
SPEC = [("G0", "null", 0.0, 0.5), ("G0", "null", 0.0, 0.8), ("SEAS-V", "null", 0.0, 0.5),
        ("G0", "lin", RHO_LIN, 0.5), ("G0", "sq", RHO_SQ, 0.5), ("SEAS-V", "lin", RHO_LIN, 0.5)]
EXP3_IDX = [0, 2, 6, 18, 24, 30]  # same cells in EXP3 main_cells() at n = 2000


def cells_for(mode: str, n_grid) -> list[dict]:
    cells = []
    for j, n in enumerate(n_grid):
        for i, (g, alt, rho, h) in enumerate(SPEC):
            idx = EXP3_IDX[i] if mode == "regress" else i + 100 * (0 if n == 2000 else 1)
            par = round(rho * (2000 / n) ** 0.5, 6)
            c = {"gen": g, "alt": alt, "H": h, "n": n, "param": par, "idx": idx}
            c["name"] = f"ind__{g}__{alt}{par:g}__H{h}__n{n}"
            cells.append(c)
    return cells


def run_chunk(args):
    cell, reps, master, kappas, okap = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            x, y = K.ind_pair(cell["gen"], cell["alt"], cell["H"], cell["n"], cell["param"], rng)
            out = K.kind_bw(x, y, rng, kappas=kappas, B=B_WB, oracle_kappas=okap)
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=list(MASTER), required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--kappas", default="0.5,1,2,4")
    ap.add_argument("--oracle-kappas", default="", help="kappas with the bare statistic only")
    ap.add_argument("--n", default="2000", help="comma list from {2000,500}")
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes (shared machine)")
    kappas = tuple(float(k) for k in a.kappas.split(","))
    okap = tuple(float(k) for k in a.oracle_kappas.split(",") if k)
    n_grid = tuple(int(v) for v in a.n.split(","))
    if a.mode == "regress" and n_grid != (2000,):
        sys.exit("regress mode is defined for n = 2000 only")
    master = MASTER[a.mode]
    cells = cells_for(a.mode, n_grid)
    if a.only:
        cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
    outdir = os.path.join(HERE, a.out)
    os.makedirs(outdir, exist_ok=True)
    fields = ["rep", "status", "error", "cpu_s"] + K.fields(kappas, okap)
    tasks, mcells = [], []
    for c in cells:
        path = os.path.join(outdir, c["name"] + ".csv")
        done = set()
        if os.path.exists(path):
            with open(path) as f:
                rd = csv.DictReader(f)
                if rd.fieldnames != fields:
                    sys.exit(f"{path}: header differs from the requested kappa set")
                done = {int(row["rep"]) for row in rd}
        todo = [r for r in range(a.reps) if r not in done]
        mcells.append({**c, "reps": a.reps, "already_done": len(done),
                       "seed": {"entropy": master, "spawn_key": [c["idx"], "rep"]}})
        cs = CHUNK[c["n"]]
        for i in range(0, len(todo), cs):
            tasks.append((c, todo[i:i + cs], master, kappas, okap))
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    with open(os.path.join(outdir, f"manifest_{stamp}.json"), "w") as f:
        json.dump({"mode": a.mode, "master": master, "kappas": kappas, "oracle_kappas": okap,
                   "seed_rule": "np.random.SeedSequence(entropy=master, spawn_key=(cell_idx, rep)); one default_rng per replication",
                   "B_wild": B_WB, "numpy": np.__version__, "python": sys.version.split()[0],
                   "workers": a.workers, "argv": sys.argv, "cells": mcells}, f, indent=1)
    # interleave cells so a partial run has balanced reps per cell
    tasks.sort(key=lambda t: (-t[0]["n"], t[1][0], t[0]["idx"]))
    print(f"{stamp} {len(tasks)} chunks, {sum(len(t[1]) for t in tasks)} replications, kappas {kappas}", flush=True)
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
            if nd % 25 == 0 or nd == len(tasks):
                print(f"{nd}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print(f"DONE {time.time() - t0:.0f}s wall", flush=True)


if __name__ == "__main__":
    main()
