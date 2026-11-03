"""EXP4 driver: seasonal adjustment (deseason.adjust) applied to each series, then the EXP2 independence tests.

Usage:
  python -I run_exp4.py --mode pilot --reps 100  --out pilot --workers 5
  python -I run_exp4.py --mode main  --reps 1000 --out main  --workers 5

One CSV per (generator, H, n) cell in --out, one row per (replication, adjustment), appended as
chunks finish. Re-running resumes per (replication, adjustment): a pair counts as done only if a
complete row with status ok exists, so error rows are retried and a new adjustment is run for
replications already on disk. Seeds: replication r
of cell c draws the pair from SeedSequence(entropy=MASTER[mode], spawn_key=(c, r, 0)). Every
adjustment's test multipliers come from a fresh rng on spawn_key=(c, r, 1), so all adjustments of
one replication share the same data and the same multipliers (common random numbers).
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))  # kmem (D7: not installed)
import numpy as np  # noqa: E402

import deseason as D  # noqa: E402
import gens_exp4  # noqa: E402
from kmem.legacy import mbcore as mb  # noqa: E402
from kmem.runner import Job, run_jobs  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASTER = {"pilot": 2026102541, "main": 2026102542}
GENS = ("G0", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow")  # append only: order fixes cell indices
H_GRID = (0.5, 0.7, 0.8)
N_GRID = (500, 2000)
ADJ = ("A0", "A1", "A2", "A3")
B_WB = 500
CHUNK = {500: 10, 2000: 2}
MAX_WORKERS = 16
KIND = ["hsic_stat", "dcov_stat", "p_hsic_wb", "p_hsic_shift", "p_hsic_sub", "rej_hsic_sub",
        "p_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "shift_check"]
LIND = ["pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw", "ewc_nu", "t_ewc", "rej_ewc", "nw_lrv_nonpos"]
FIELDS = ["rep", "adj", "status", "error", "cpu_s", "px", "py", "srcx", "srcy", "reused"] + KIND + LIND


def cell_list() -> list[dict]:
    cells = []
    for g in GENS:
        for h in H_GRID:
            for n in N_GRID:
                cells.append({"gen": g, "H": h, "n": n})
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = f"ind__{c['gen']}__H{c['H']}__n{c['n']}"
    return cells


def run_chunk(args):
    cell, reps, master, adjs = args
    rows = []
    for r in reps:
        try:
            rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r, 0)))
            x, y = gens_exp4.pair(cell["gen"], cell["H"], cell["n"], rng)
            src = (D.estimate_period(x)[1], D.estimate_period(y)[1])
        except Exception as e:  # recorded in the row, never swallowed
            rows += [{"rep": r, "adj": a, "status": "error", "error": repr(e)[:300], "cpu_s": 0} for a in adjs]
            continue
        done = {}
        for a in adjs:
            t0 = time.process_time()
            status, err, out, reused = "ok", "", {}, 0
            try:
                xa, px = D.adjust(x, a)
                ya, py = D.adjust(y, a)
                if a == "A3" and px == py == D.P_TRUE and "A2" in done:
                    out, reused = dict(done["A2"]), 1
                else:
                    trng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r, 1)))
                    out = mb.kind_tests(xa, ya, trng, B=B_WB)
                    out.update(mb.lind_tests(xa, ya))
                done[a] = out
            except Exception as e:  # recorded in the row, never swallowed
                status, err, px, py = "error", repr(e)[:300], -1, -1
            rows.append({"rep": r, "adj": a, "status": status, "error": err,
                         "cpu_s": round(time.process_time() - t0, 4), "px": px, "py": py,
                         "srcx": src[0] if a == "A3" else "", "srcy": src[1] if a == "A3" else "",
                         "reused": reused, **out})
    return cell["name"], rows


def job_for(cell: dict, target: int, master: int, opts: dict | None = None, adjs: tuple | None = None) -> Job:
    """The runner's job for one cell (kmem.runner). A replication is done when every adjustment has its row.
    adjs defaults to the cell's recorded adjustments (as in the EXP4 spec)."""
    adjs = tuple(adjs if adjs is not None else cell["adjustments"])
    return Job(name=cell["name"], fields=cell, seed_index=cell["idx"], target=target, master=master,
               columns=tuple(FIELDS), chunk=CHUNK[cell["n"]], subkeys=adjs, subkey_column="adj",
               seed_doc={"entropy": master, "spawn_key_data": [cell["idx"], "rep", 0],
                         "spawn_key_test": [cell["idx"], "rep", 1]},
               manifest_cell={"adjustments": list(adjs)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    ap.add_argument("--adj", default=",".join(ADJ), help="adjustments; n=2000 may drop A0 (see --adj2000)")
    ap.add_argument("--adj2000", default=None, help="adjustment list for n = 2000 cells, if different")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes (DESIGN.md amendment 5)")
    master = MASTER[a.mode]
    cells = cell_list()
    if a.only:
        cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
    jobs = [job_for(c, a.reps, master, adjs=tuple((a.adj2000 if (c["n"] == 2000 and a.adj2000) else a.adj).split(",")))
            for c in cells]
    run_jobs(jobs, run_chunk, os.path.join(HERE, a.out), a.workers,
             {"mode": a.mode, "master": master,
              "seed_rule": "data: SeedSequence(master, spawn_key=(cell_idx, rep, 0)); tests: fresh rng on (cell_idx, rep, 1) per adjustment",
              "B_wild": B_WB, "period_true": D.P_TRUE})


if __name__ == "__main__":
    main()
