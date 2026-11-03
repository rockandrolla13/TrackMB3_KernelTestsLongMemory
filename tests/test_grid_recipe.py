"""The "add an experiment grid" recipe of kmem/EXTENDING.md, as a toy experiment (plan step 5.2).

Two cells, built only from kmem pieces (generators, kernels, features, rules), run through kmem.runner.
A second run, stopped and resumed, must give the same values.
"""
from __future__ import annotations

import csv
import math
import time

import numpy as np

from kmem import kernels
from kmem.bandwidth import median_bw_fast
from kmem.calibration.wild import Multipliers
from kmem.features import compress, gauss_features
from kmem.generators.independence import pair
from kmem.rules import rules_fast
from kmem.runner import Job, run_jobs

MASTER = 999  # a toy master; a real experiment registers its masters in kmem/seeds.py
KERNELS = ("raw", "inc1")
COLUMNS = ("rep", "status", "error", "cpu_s") + tuple(f"{k}_{kn}" for kn in KERNELS for k in ("t", "p_wb", "r_wb"))


def cell_list() -> list[dict]:
    cells = [{"gen": "G0", "alt": "null", "H": 0.8, "n": 200}, {"gen": "G0", "alt": "lin", "H": 0.5, "n": 200}]
    for i, c in enumerate(cells):  # explicit, append-only seed indices
        c["idx"], c["name"] = i, f"ind__{c['gen']}__{c['alt']}__H{c['H']}__n{c['n']}"
    return cells


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        x, y = pair(cell["gen"], cell["alt"], cell["H"], cell["H"], cell["n"], 0.3, rng)
        mult = Multipliers(cell["n"], 100, rng)
        fx, fy = gauss_features(x, median_bw_fast(x)), gauss_features(y, median_bw_fast(y))
        out = {}
        for kn in KERNELS:
            k = kernels.get(kn)
            gx, gy = k.features(fx, x), k.features(fy, y)
            res = rules_fast(compress(gx), compress(gy), mult.get(gx.shape[0]), rules=("wild",))
            out.update({f"{key}_{kn}": res[key] for key in ("t", "p_wb", "r_wb")})
        rows.append({"rep": r, "status": "ok", "error": "", "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


def job_for(cell: dict, target: int, master: int, opts: dict | None = None) -> Job:
    return Job(name=cell["name"], fields=cell, seed_index=cell["idx"], target=target, master=master,
               columns=COLUMNS, chunk=2)


def _values(path) -> dict:
    with open(path, newline="") as f:
        return {r["rep"]: {k: v for k, v in r.items() if k != "cpu_s"} for r in csv.DictReader(f)}


def test_toy_grid_runs_and_reruns_identically(tmp_path) -> None:
    jobs = [job_for(c, 4, MASTER) for c in cell_list()]
    run_jobs(jobs, run_chunk, str(tmp_path / "a"), 2, {"mode": "toy", "master": MASTER})
    # a second run, stopped after two replications per cell and then resumed
    run_jobs([job_for(c, 2, MASTER) for c in cell_list()], run_chunk, str(tmp_path / "b"), 2, {"mode": "toy"})
    run_jobs(jobs, run_chunk, str(tmp_path / "b"), 2, {"mode": "toy"})
    for c in cell_list():
        a, b = _values(tmp_path / "a" / f"{c['name']}.csv"), _values(tmp_path / "b" / f"{c['name']}.csv")
        assert sorted(a) == ["0", "1", "2", "3"] and a == b
    lin = _values(tmp_path / "a" / "ind__G0__lin__H0.5__n200.csv")
    assert all(math.isfinite(float(r["t_raw"])) for r in lin.values())
