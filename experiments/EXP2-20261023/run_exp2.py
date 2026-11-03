"""EXP2 driver: feature x memory factorial for split MMD (mmdcal.rep_e2) and HSIC/dCov (mbcore.kind_tests).

Usage:
  python -I run_exp2.py --mode pilot --reps 100  --out pilot --workers 16
  python -I run_exp2.py --mode main  --reps 1000 --out main  --workers 16

One CSV per cell in --out, one row per replication, appended as chunks finish. Re-running resumes
(rows already on disk are skipped). Seeds: replication r of cell c uses
SeedSequence(entropy=MASTER[(family, mode)], spawn_key=(c, r)); one default_rng per replication.
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

from kmem.legacy import gens  # noqa: E402
from kmem.legacy import mbcore as mb  # noqa: E402
from kmem.legacy import mmdcal as M  # noqa: E402
from kmem.runner import Job, run_jobs  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASTER = {("mmd", "pilot"): 2026102321, ("mmd", "main"): 2026102322,
          ("ind", "pilot"): 2026102331, ("ind", "main"): 2026102332}
H_GRID = (0.5, 0.7, 0.8)
N_GRID = (500, 2000)
B_WB = 500
CHUNK = {500: 25, 2000: 5}
MAX_WORKERS = 16
KIND = ["hsic_stat", "dcov_stat", "p_hsic_wb", "p_hsic_shift", "p_hsic_sub", "rej_hsic_sub",
        "p_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "shift_check"]
LIND = ["pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw", "ewc_nu", "t_ewc", "rej_ewc", "nw_lrv_nonpos"]
MMD = (["mmd2", "sigma", "t_n"] + [f"p_{f}_{int(e * 10):02d}" for f, e in M.WILD_E2]
       + ["h_lw", "h_gph", "l_sqrt", "q95_sqrt", "l_bd", "q95_bd", "sn_T", "sn_q95_05", "sn_q95_07"])
FIELDS = {"mmd": MMD, "ind": KIND + LIND}


def cell_list(family: str) -> list[dict]:
    gl = list(gens.MAIN) + (["SEAS-M-shared"] if family == "ind" else [])
    cells = []
    for g in gl:
        for h in H_GRID:
            for n in N_GRID:
                cells.append({"family": family, "gen": g, "H": h, "n": n})
    for i, c in enumerate(cells):
        c["idx"] = i
        c["name"] = f"{family}__{c['gen']}__H{c['H']}__n{c['n']}"
    return cells


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            if cell["family"] == "mmd":
                x = gens.generate(cell["gen"], cell["H"], cell["n"], rng)
                out = M.rep_e2(x, rng)
            else:
                x, y = gens.pair(cell["gen"], cell["H"], cell["n"], rng)
                out = mb.kind_tests(x, y, rng, B=B_WB)
                out.update(mb.lind_tests(x, y))
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


def job_for(cell: dict, target: int, master: int, opts: dict | None = None) -> Job:
    """The runner's job for one cell (kmem.runner); used by main() and by `python -m kmem run EXP2`."""
    return Job(name=cell["name"], fields=cell, seed_index=cell["idx"], target=target, master=master,
               columns=("rep", "status", "error", "cpu_s", *FIELDS[cell["family"]]), chunk=CHUNK[cell["n"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "main"], required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--family", default="mmd,ind")
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes")
    jobs = []
    for fam in a.family.split(","):
        master = MASTER[(fam, a.mode)]
        cells = cell_list(fam)
        if a.only:
            cells = [c for c in cells if any(p in c["name"] for p in a.only.split(","))]
        jobs += [job_for(c, a.reps, master) for c in cells]
    run_jobs(jobs, run_chunk, os.path.join(HERE, a.out), a.workers,
             {"mode": a.mode, "masters": {f"{k[0]}_{k[1]}": v for k, v in MASTER.items()},
              "seed_rule": "np.random.SeedSequence(entropy=master, spawn_key=(cell_idx, rep)); one default_rng per replication", "B_wild": B_WB})


if __name__ == "__main__":
    main()
