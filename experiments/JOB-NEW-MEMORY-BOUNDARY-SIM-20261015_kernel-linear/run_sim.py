"""Simulation driver for the kernel-linear arm (K-ind, L-ind).

Usage:
  python -I run_sim.py --mode pilot --reps 20 --out pilot
  python -I run_sim.py --mode main --reps 2000 [--reps-largest 1000] --out cells

Writes one CSV per cell into --out, one row per replication, appended as results arrive.
Re-running the same command resumes: replications already in a CSV are skipped.
Seeds: replication r of cell c uses SeedSequence(entropy=MASTER, spawn_key=(c, r)), where c is
the cell's index in the fixed cell list below, so results do not depend on chunking or worker count.
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
MASTER = {"pilot": 20261015001, "main": 20261015002}
KIND_N = (500, 1000, 2000, 4000)
LIND_T = (500, 2000, 8000)
POWER_R = 0.25
POWER_SID = "G2_H0.7"
B_WB = 500
CHUNK = {500: 25, 1000: 20, 2000: 10, 4000: 5, 8000: 50}

KIND_FIELDS = ["hsic_stat", "dcov_stat", "p_hsic_wb", "p_hsic_shift", "p_hsic_sub", "rej_hsic_sub",
               "p_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "shift_check"]
LIND_FIELDS = ["pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw", "ewc_nu", "t_ewc", "rej_ewc", "nw_lrv_nonpos"]


def cell_list() -> list[dict]:
    """Fixed, ordered cell list. Index in this list is the seed spawn key's first element."""
    cells = []
    settings = mb.make_settings()
    for st in settings:
        for n in KIND_N:
            cells.append({"family": "kind", "n": n, "r": 0.0, **st})
    for st in settings:
        for T in LIND_T:
            cells.append({"family": "lind", "n": T, "r": 0.0, **st})
    pst = next(s for s in settings if s["sid"] == POWER_SID)
    for n in KIND_N:
        cells.append({"family": "kind", "n": n, "r": POWER_R, **pst})
    for T in LIND_T:
        cells.append({"family": "lind", "n": T, "r": POWER_R, **pst})
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = f"{c['family']}__{c['sid']}__n{c['n']}" + (f"__r{c['r']}" if c["r"] else "")
    return cells


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        ss = np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r))
        rng = np.random.default_rng(ss)
        status, err = "ok", ""
        out: dict = {}
        try:
            x, y = mb.gen_series(cell, cell["n"], rng, r=cell["r"])
            if cell["family"] == "kind":
                out = mb.kind_tests(x, y, rng, B=B_WB)
            else:
                out = mb.lind_tests(x, y)
        except Exception as e:  # recorded, never swallowed: the row carries the failure
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err,
                     "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--reps-largest", type=int, default=None,
                    help="reps for kind cells at n=4000 (fallback rule)")
    ap.add_argument("--power-reps", type=int, default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--exclude", default=None,
                    help="comma-separated substrings; matching cells are dropped (listed in manifest with reps 0)")
    ap.add_argument("--only", default=None, help="comma-separated substrings; a cell runs if any matches")
    a = ap.parse_args()
    if not os.path.exists(os.path.join(HERE, "validation.PASS")):
        sys.exit("validation.PASS missing: run validate.py first")
    if a.workers > 8:
        sys.exit("at most 8 worker processes")
    outdir = os.path.join(HERE, a.out)
    os.makedirs(outdir, exist_ok=True)
    master = MASTER[a.mode]
    cells = cell_list()
    if a.only:
        pats = a.only.split(",")
        cells = [c for c in cells if any(p in c["name"] for p in pats)]

    excl = a.exclude.split(",") if a.exclude else []

    def target_reps(c):
        if any(p in c["name"] for p in excl):
            return 0
        if c["r"] and a.power_reps is not None:
            return a.power_reps
        if c["family"] == "kind" and c["n"] == 4000 and a.reps_largest is not None:
            return a.reps_largest
        return a.reps

    tasks = []
    manifest_cells = []
    for c in cells:
        path = os.path.join(outdir, c["name"] + ".csv")
        done = set()
        if os.path.exists(path):
            with open(path) as f:
                done = {int(row["rep"]) for row in csv.DictReader(f)}
        todo = [r for r in range(target_reps(c)) if r not in done]
        manifest_cells.append({k: c[k] for k in ("idx", "name", "family", "proc", "sid", "HX", "HY", "lam", "n", "r")}
                              | {"reps": target_reps(c), "dropped": target_reps(c) == 0, "seed": {"entropy": master, "spawn_key": [c["idx"], "rep"]}})
        cs = CHUNK[c["n"]]
        for i in range(0, len(todo), cs):
            tasks.append((c, todo[i:i + cs], master))
    mname = "manifest.json" if not a.only else "manifest_partial_" + time.strftime("%Y%m%dT%H%M%S", time.gmtime()) + ".json"
    with open(os.path.join(outdir, mname), "w") as f:
        json.dump({"mode": a.mode, "exclude": excl, "master_entropy": master,
                   "seed_rule": "np.random.SeedSequence(entropy=master_entropy, spawn_key=(cell_idx, rep)); one default_rng per replication",
                   "B_wild_bootstrap": B_WB, "numpy": np.__version__, "python": sys.version.split()[0],
                   "workers": a.workers, "cells": manifest_cells}, f, indent=1)
    # biggest jobs first for load balance
    tasks.sort(key=lambda t: -(t[0]["n"] if t[0]["family"] == "kind" else t[0]["n"] / 16))
    print(f"{len(tasks)} chunks, {sum(len(t[1]) for t in tasks)} replications", flush=True)
    fields = {"kind": ["rep", "status", "error", "cpu_s"] + KIND_FIELDS,
              "lind": ["rep", "status", "error", "cpu_s"] + LIND_FIELDS}
    fam = {c["name"]: c["family"] for c in cells}
    t0 = time.time()
    ndone = 0
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=50) as pool:
        for name, rows in pool.imap_unordered(run_chunk, tasks):
            path = os.path.join(outdir, name + ".csv")
            new = not os.path.exists(path)
            with open(path, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields[fam[name]], extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows(rows)
            ndone += 1
            if ndone % 20 == 0 or ndone == len(tasks):
                print(f"{ndone}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
