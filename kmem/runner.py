"""One runner for every experiment (plan step 4.1, AR-DRY-004, AR-PAR-001, CR-002, CR-003, CR-014, CR-021).

An experiment supplies its row function (the driver's run_chunk) and one Job per cell. The runner does the rest,
identically for all experiments:

- resume: a replication counts as done only with a complete status-ok row (the EXP4/EXP7 rule), so error rows
  and rows cut by a crash are recomputed; an existing CSV must have exactly this run's header;
- tasks: chunks of missing replications, ordered by progress through each cell, so a stopped run is balanced;
- manifest: masters, seed rule, options, writer format, git commit, versions, argv and every cell's seed;
- writing: one CSV per cell, appended as chunks finish, in the experiment's float format.

Values depend only on (master, cell seed index, replication), never on scheduling, chunking or resumes.
"""
from __future__ import annotations

import csv
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass(frozen=True)
class Job:
    """One cell to run. args for run_chunk: (fields, reps, master, *extra[, missing subkeys])."""
    name: str
    fields: dict
    seed_index: int
    target: int
    master: int
    columns: tuple
    chunk: int
    extra: tuple = ()
    subkeys: tuple | None = None  # EXP4: the adjustments; a replication is done when each has its row
    subkey_column: str = "adj"
    seed_doc: dict | None = None
    manifest_cell: dict = field(default_factory=dict)


def fmt12g(key: str, v):
    """EXP7's CSV format: twelve significant digits, exact-path (xraw_) columns at full precision."""
    if key.startswith("xraw_") or isinstance(v, (str, int)):
        return v
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.12g}"
    return v


FORMATTERS: dict[str, Callable | None] = {"str": None, "fmt12g": fmt12g}


def completed(path: str, columns, subkey_column: str | None = None) -> set:
    """(rep, subkey) pairs with a complete status-ok row. A different header is an error, not a resume."""
    if not os.path.exists(path):
        return set()
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        if rd.fieldnames != list(columns):
            raise SystemExit(f"{path}: header differs from this run's columns; use another --out")
        return {(int(r["rep"]), r[subkey_column] if subkey_column else None) for r in rd
                if None not in r.values() and r["status"] == "ok"}


def plan(jobs: list[Job], out_dir: str) -> tuple[list[tuple], list[dict]]:
    tasks, mcells = [], []
    for job in jobs:
        done = completed(os.path.join(out_dir, job.name + ".csv"), job.columns,
                         job.subkey_column if job.subkeys is not None else None)
        subs = job.subkeys if job.subkeys is not None else (None,)
        groups: dict[tuple, list[int]] = {}
        for r in range(job.target):
            missing = tuple(s for s in subs if (r, s) not in done)
            if missing:
                groups.setdefault(missing, []).append(r)
        seed = job.seed_doc or {"entropy": job.master, "spawn_key": [job.seed_index, "rep"]}
        mcells.append({**job.fields, "target": job.target, "reps": job.target, "already_done": len(done),
                       "seed": seed, **job.manifest_cell})
        for missing, reps in groups.items():
            for i in range(0, len(reps), job.chunk):
                part = reps[i:i + job.chunk]
                args = (job.fields, part, job.master) + tuple(job.extra) + ((missing,) if job.subkeys is not None else ())
                tasks.append((part[0] / max(job.target, 1), job.seed_index, args))
    tasks.sort(key=lambda t: (t[0], t[1]))
    return [t[2] for t in tasks], mcells


def git_state() -> dict:
    def git(*a):
        return subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True)
    head = git("rev-parse", "HEAD")
    dirty = git("diff", "--name-only", "HEAD", "--", "kmem", ":(glob)experiments/**/*.py")
    return {"commit": head.stdout.strip() if head.returncode == 0 else "unknown",
            "code_modified": sorted(dirty.stdout.split()) if dirty.returncode == 0 else None}


def run_jobs(jobs: list[Job], row_fn: Callable, out_dir: str, workers: int, manifest: dict, writer: str = "str",
             max_tasks_per_child: int = 50, progress_every: int = 25) -> str:
    """Run every missing replication of every job; return the manifest path."""
    os.makedirs(out_dir, exist_ok=True)
    tasks, mcells = plan(jobs, out_dir)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    mpath = os.path.join(out_dir, f"manifest_{stamp}.json")
    from kmem import __version__
    with open(mpath, "w") as f:
        json.dump({**manifest, "writer": writer, "kmem": __version__, "git": git_state(), "numpy": np.__version__,
                   "python": sys.version.split()[0], "workers": workers, "argv": sys.argv, "cells": mcells}, f, indent=1)
    by_name = {j.name: j for j in jobs}
    fmt = FORMATTERS[writer]
    nrep = sum(len(t[1]) for t in tasks)
    print(f"{stamp} {len(tasks)} chunks, {nrep} replications, {len(jobs)} cells", flush=True)
    t0 = time.time()
    if tasks:
        nd = 0
        with mp.get_context("fork").Pool(workers, maxtasksperchild=max_tasks_per_child) as pool:
            for result in pool.imap_unordered(row_fn, tasks):
                name, rows = result[0], result[-1]
                job = by_name[name]
                path = os.path.join(out_dir, name + ".csv")
                new = not os.path.exists(path)
                with open(path, "a", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=list(job.columns), extrasaction="ignore")
                    if new:
                        w.writeheader()
                    w.writerows(rows if fmt is None else [{k: fmt(k, v) for k, v in r.items()} for r in rows])
                nd += 1
                if nd % progress_every == 0 or nd == len(tasks):
                    print(f"{nd}/{len(tasks)} chunks, {time.time() - t0:.0f}s wall", flush=True)
    print(f"DONE {time.time() - t0:.0f}s wall", flush=True)
    return mpath
