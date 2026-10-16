"""Driver: run every cell of one family, write per-cell CSVs incrementally, resume on restart.

usage: python -I run.py {khom,dml} {pilot,main} --reps R [--reps-largest R2] [--procs 6]

Seeds: rep r of cell c in family f uses SeedSequence(entropy=seed, spawn_key=(f, c, r)), i.e. the
r-th child of the c-th child of the f-th child of SeedSequence(seed). Results do not depend on
chunking or process count. f = 0 for khom, 1 for dml. seed = MASTER_SEED (main) or PILOT_SEED.
"""
import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per worker process; set before numpy import

import numpy as np  # noqa: E402

import common as C  # noqa: E402

FAMILY_IDX = {"khom": 0, "dml": 1}
_MOD = None


def _init(family):
    global _MOD
    _MOD = __import__(family)


def _work(task):
    fam_idx, seed, cell, r0, r1 = task
    rows = []
    for r in range(r0, r1):
        rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(fam_idx, cell["cell_idx"], r)))
        t0 = time.process_time()
        try:
            out = _MOD.run_rep(cell, rng)
            status, err = "ok", ""
        except Exception as exc:  # failures stay visible (invariant 4)
            out, status, err = {}, "error", f"{type(exc).__name__}: {exc}"
        out.update({"rep": r, "status": status, "error": err, "cpu_s": time.process_time() - t0})
        rows.append(out)
    return cell["cell_id"], rows


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("family", choices=("khom", "dml"))
    ap.add_argument("tag", choices=("pilot", "main"))
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--reps-largest", type=int, default=None)
    ap.add_argument("--procs", type=int, default=6)
    ap.add_argument("--chunk", type=int, default=25)
    a = ap.parse_args()
    if a.procs > 6:
        sys.exit("at most 6 processes")
    mod = __import__(a.family)
    seed = C.PILOT_SEED if a.tag == "pilot" else C.MASTER_SEED
    fam_idx = FAMILY_IDX[a.family]
    cells = mod.cell_list()
    largest = max(c[mod.SIZE_KEY] for c in cells)
    outdir = os.path.join(HERE, a.tag, a.family)
    os.makedirs(outdir, exist_ok=True)
    cols = ("rep", "status", "error", "cpu_s") + tuple(mod.COLUMNS)

    tasks = []
    reps_of = {}
    for c in sorted(cells, key=lambda c: -c[mod.SIZE_KEY]):
        reps = a.reps_largest if (a.reps_largest and c[mod.SIZE_KEY] == largest) else a.reps
        reps_of[c["cell_id"]] = reps
        path = os.path.join(outdir, c["cell_id"] + ".csv")
        done = set()
        if os.path.exists(path):
            with open(path) as f:
                done = {int(row["rep"]) for row in csv.DictReader(f)}
        else:
            with open(path, "w", newline="") as f:
                csv.writer(f).writerow(cols)
        todo = [r for r in range(reps) if r not in done]
        for i in range(0, len(todo), a.chunk):
            blk = todo[i:i + a.chunk]
            if blk == list(range(blk[0], blk[-1] + 1)):
                tasks.append((fam_idx, seed, c, blk[0], blk[-1] + 1))
            else:
                tasks.extend((fam_idx, seed, c, r, r + 1) for r in blk)

    manifest_path = os.path.join(outdir, "manifest.json")
    manifest = {
        "job": "JOB-NEW-MEMORY-BOUNDARY-SIM-20261015",
        "worker": "claude-mb-mmddml-20261015",
        "family": a.family,
        "tag": a.tag,
        "evidence_class": "exploratory simulation" + (" (timing pilot; not evidence)" if a.tag == "pilot" else ""),
        "master_seed": seed,
        "seed_scheme": f"SeedSequence(entropy={seed}, spawn_key=({fam_idx}, cell_idx, rep))",
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "procs": a.procs,
        "script_sha256": {f: sha256(os.path.join(HERE, f)) for f in
                          ("common.py", "khom.py", "dml.py", "run.py")},
        "cells": [dict(c, reps=reps_of[c["cell_id"]],
                       seed_spawn_key=[fam_idx, c["cell_idx"], "rep"]) for c in cells],
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "tasks_pending": len(tasks),
    }
    with open(manifest_path + ".tmp", "w") as f:
        json.dump(manifest, f, indent=1)
    os.replace(manifest_path + ".tmp", manifest_path)

    t_wall = time.time()
    cpu = 0.0
    n_err = 0
    with mp.get_context("fork").Pool(a.procs, initializer=_init, initargs=(a.family,)) as pool:
        for k, (cell_id, rows) in enumerate(pool.imap_unordered(_work, tasks), 1):
            with open(os.path.join(outdir, cell_id + ".csv"), "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols, extrasaction="raise")
                for row in rows:
                    w.writerow({c: row.get(c, "") for c in cols})
            cpu += sum(r["cpu_s"] for r in rows)
            n_err += sum(r["status"] != "ok" for r in rows)
            if k % 20 == 0 or k == len(tasks):
                print(f"{time.strftime('%H:%M:%S')} {k}/{len(tasks)} tasks, cpu {cpu / 3600:.2f} h, "
                      f"wall {(time.time() - t_wall) / 60:.1f} min, errors {n_err}", flush=True)

    manifest["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["wall_s_this_invocation"] = time.time() - t_wall
    manifest["cpu_s_this_invocation"] = cpu
    manifest["errors_this_invocation"] = n_err
    with open(manifest_path + ".tmp", "w") as f:
        json.dump(manifest, f, indent=1)
    os.replace(manifest_path + ".tmp", manifest_path)


if __name__ == "__main__":
    main()
