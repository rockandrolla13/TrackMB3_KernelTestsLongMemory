"""EXP8 driver: memory kernels on the stock placebo (arm A) and a surrogate null (arm B). See DESIGN.md.

Usage, from the repo root:
  python experiments/EXP8-stock-surrogate-20261108/run_exp8.py --mode sim  --out sim  --workers 16
  python experiments/EXP8-stock-surrogate-20261108/run_exp8.py --mode main --out main --workers 16

One CSV per cell in --out, one row per pair, appended as chunks finish; re-running resumes.
Arm A seeds: STOCK-20261023's SeedSequence(2026102301, spawn_key=(stock_cell, pair)) and its draw order, so the
windows are the saved STOCK windows. Arm B and sim seeds: SeedSequence(MASTER, spawn_key=(cell_idx, pair)).
Raw durations stay in data/raw/; rows hold derived statistics, days and window starts only.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "EXP7-memory-kernel-20261101"))
import numpy as np  # noqa: E402

import exp7_core as K  # noqa: E402
from kmem.generators.registry import get_process  # noqa: E402
from kmem.runner import Job, run_jobs  # noqa: E402
from kmem.surrogates import hurst_lw, iaaft, iaaft_drift, variance_trend  # noqa: E402

MASTER = {"main": 2026110801, "sim": 2026110802}
STOCK_SEED = 2026102301
STOCK_N = (500, 1000, 2000)  # STOCK-20261023 cell order: pairing-major, then n
SYMBOLS = ("AAPL", "MSFT", "NVDA")
PAIRS = (("AAPL", "MSFT"), ("AAPL", "NVDA"), ("MSFT", "NVDA"), ("AAPL", "AAPL"))
MIN_GAP_DAYS = 20
NS = (500, 2000)
MAX_WORKERS = 16
KERNELS = tuple(k for k in K.IND_KERNELS if k != "frac_oracle")  # frac_oracle needs the true H
SIM_KERNELS = ("raw", "inc1")
SIM_H = 0.8
CHUNK = {500: 25, 2000: 5}
DATA: dict = {}

META = ["day_a", "day_b", "start_a", "start_b", "sd_x", "sd_y", "H_x", "H_y", "slope_x", "slope_y",
        "slope_perm_x", "slope_perm_y", "slope_src_y"]


def load_data() -> None:
    cache = os.path.join(ROOT, "data/raw/stock_duration_cache_2024")
    for s in SYMBOLS:
        with np.load(os.path.join(cache, f"{s}.npz")) as z:
            x = z["dur"].astype(float)
            day = z["dur_day"].astype(int)
            days = z["days"].astype(str)
        first = np.searchsorted(day, np.arange(days.size + 1))
        if np.any(np.diff(first) < max(NS)):
            raise ValueError(f"A day has fewer than {max(NS)} durations: {s}")
        DATA[s] = {"x": x, "days": days, "first": first}
    if not all(np.array_equal(DATA["AAPL"]["days"], DATA[s]["days"]) for s in SYMBOLS[1:]):
        raise ValueError("Stock trading calendars differ")


def draw_windows(a: str, b: str, n: int, rng: np.random.Generator):
    """STOCK-20261023 run_chunk's draws, in its order: two separated days, then one window start per side."""
    ndays = DATA[a]["days"].size
    for _ in range(1000):
        i, j = (int(v) for v in rng.integers(0, ndays, size=2))
        if abs(i - j) >= MIN_GAP_DAYS:
            break
    else:
        raise RuntimeError("Could not sample separated days")
    out = []
    for s, d in ((a, i), (b, j)):
        lo, hi = DATA[s]["first"][d:d + 2]
        st = int(rng.integers(lo, hi - n + 1))
        out.append((DATA[s]["x"][st:st + n], st, str(DATA[s]["days"][d])))
    return out


def window_meta(x, y, y_src, rng) -> dict:
    sx_, _ = variance_trend(x)
    sy_, _ = variance_trend(y)
    return {"sd_x": float(x.std()), "sd_y": float(y.std()), "H_x": hurst_lw(x), "H_y": hurst_lw(y),
            "slope_x": sx_, "slope_y": sy_,
            "slope_perm_x": variance_trend(x[rng.permutation(x.size)])[0],
            "slope_perm_y": variance_trend(y[rng.permutation(y.size)])[0],
            "slope_src_y": variance_trend(y_src)[0]}


def cells() -> list[dict]:
    out = []
    for arm in ("A", "S1", "S2"):
        for pi, (a, b) in enumerate(PAIRS):
            for n in NS:
                out.append({"arm": arm, "a": a, "b": b, "n": n, "stock_cell": pi * len(STOCK_N) + STOCK_N.index(n),
                            "target": 500 if arm == "A" else 1000})
    for i, c in enumerate(out):
        c["idx"] = i
        c["name"] = f"{c['arm']}__{c['a']}-{c['b']}__n{c['n']}"
    return out


def sim_cells() -> list[dict]:
    return [{"arm": "sim", "a": "fGn", "b": "fGn", "n": 2000, "H": SIM_H, "idx": 0, "target": 1000,
             "name": f"sim__S1__H{SIM_H}__n2000"},
            {"arm": "simraw", "a": "fGn", "b": "fGn", "n": 2000, "H": SIM_H, "idx": 1, "target": 1000,
             "name": f"sim__none__H{SIM_H}__n2000"}]


def columns(c: dict) -> list[str]:
    head = ["rep", "status", "error", "cpu_s"]
    if c["arm"].startswith("sim"):
        return head + ["H_x", "H_y", "sx", "sy", "rx_rank", "ry_rank"] + K.kernel_columns(SIM_KERNELS)
    return head + META + K.ind_columns(KERNELS)


def run_chunk(args):
    cell, reps, master = args
    rows = []
    for r in reps:
        t0 = time.process_time()
        if cell["arm"] == "A":
            rng = np.random.default_rng(np.random.SeedSequence(STOCK_SEED, spawn_key=(cell["stock_cell"], r)))
        else:
            rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], r)))
        status, err, out = "ok", "", {}
        try:
            if cell["arm"].startswith("sim"):
                g = get_process("G0")
                x = g(cell["H"], cell["n"], rng)
                y_src = g(cell["H"], cell["n"], rng)
                y = iaaft(y_src, rng) if cell["arm"] == "sim" else y_src
                out = {"H_x": hurst_lw(x), "H_y": hurst_lw(y)}
                out.update(K.ind_bundle(x, y, rng, kernels=SIM_KERNELS, full=False, hsplit=False))
            else:
                (x, sa, da), (y_src, sb, db) = draw_windows(cell["a"], cell["b"], cell["n"], rng)
                if x.std() == 0 or y_src.std() == 0:
                    raise ValueError("Constant duration window")
                y = {"A": lambda: y_src, "S1": lambda: iaaft(y_src, rng), "S2": lambda: iaaft_drift(y_src, rng)}[cell["arm"]]()
                out = {"day_a": da, "day_b": db, "start_a": sa, "start_b": sb}
                out.update(K.ind_bundle(x, y, rng, kernels=KERNELS, full=True, hsplit=True, strict=False))
                out.update(window_meta(x, y, y_src, rng))
        except Exception as e:  # recorded in the row, never swallowed
            status, err = "error", repr(e)[:300]
        rows.append({"rep": r, "status": status, "error": err, "cpu_s": round(time.process_time() - t0, 4), **out})
    return cell["name"], rows


def job_for(cell: dict, target: int, master: int, opts: dict | None = None) -> Job:
    seed = ({"entropy": STOCK_SEED, "spawn_key": [cell["stock_cell"], "rep"]} if cell["arm"] == "A" else None)
    return Job(name=cell["name"], fields=cell, seed_index=cell["idx"], target=target, master=master,
               columns=tuple(columns(cell)), chunk=CHUNK[cell["n"]], seed_doc=seed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=list(MASTER), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--reps", type=int, default=None, help="override every target (smoke runs)")
    ap.add_argument("--only", default=None, help="comma list of substrings of cell names")
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} worker processes")
    cs = sim_cells() if a.mode == "sim" else cells()
    if a.mode == "main":
        load_data()
    if a.only:
        cs = [c for c in cs if any(p in c["name"] for p in a.only.split(","))]
    jobs = [job_for(c, a.reps or c["target"], MASTER[a.mode]) for c in cs]
    run_jobs(jobs, run_chunk, os.path.join(HERE, a.out), a.workers,
             {"experiment": "EXP8", "mode": a.mode, "master": MASTER[a.mode], "stock_seed": STOCK_SEED,
              "kernels": list(KERNELS if a.mode == "main" else SIM_KERNELS), "B_wild": K.B_WILD,
              "min_gap_days": MIN_GAP_DAYS, "series": "dur (log duration minus pooled clock-bin mean)",
              "seed_rule": "arm A: STOCK-20261023 seeds and draw order; others SeedSequence(master, (cell_idx, rep))"},
             writer="fmt12g", max_tasks_per_child=200, progress_every=50)


if __name__ == "__main__":
    main()
