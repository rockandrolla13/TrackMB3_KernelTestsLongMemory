"""Exploratory stock-duration placebo using the existing E1 test implementations.

Pairs use different trading days at least 20 sessions apart. Every sampled
window remains inside one day. Reused days make pair outcomes dependent, so
raw rejection fractions are descriptive, not binomial false-positive estimates.
"""

from __future__ import annotations

import os

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/EXP1-20261016"))
import mbcore as mb  # noqa: E402

SYMBOLS = ("AAPL", "MSFT", "NVDA")
PAIRS = (("AAPL", "MSFT"), ("AAPL", "NVDA"), ("MSFT", "NVDA"), ("AAPL", "AAPL"))
MIN_GAP_DAYS = 20
SEED = 2026102301
DATA = {}
FIELDS = (
    "pair", "status", "error", "cpu_s", "day_a", "day_b", "start_a", "start_b",
    "sd_x", "sd_y", "hsic_stat", "dcov_stat", "p_hsic_wb", "p_hsic_shift",
    "p_hsic_sub", "rej_hsic_sub", "p_dcov_wb", "p_dcov_sub", "rej_dcov_sub",
    "shift_check", "pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw",
    "ewc_nu", "t_ewc", "rej_ewc", "nw_lrv_nonpos",
)


def load_data() -> None:
    cache = ROOT / "data/raw/stock_duration_cache_2024"
    for symbol in SYMBOLS:
        with np.load(cache / f"{symbol}.npz") as z:
            x = z["dur"].astype(float)
            day = z["dur_day"].astype(int)
            days = z["days"].astype(str)
        first = np.searchsorted(day, np.arange(days.size + 1))
        if np.any(np.diff(first) < 2000):
            raise ValueError(f"A day has fewer than 2000 durations: {symbol}")
        DATA[symbol] = {"x": x, "days": days, "first": first}
    if not all(np.array_equal(DATA["AAPL"]["days"], DATA[s]["days"]) for s in SYMBOLS[1:]):
        raise ValueError("Stock trading calendars differ")


def draw(symbol: str, day: int, n: int, rng: np.random.Generator):
    record = DATA[symbol]
    lo, hi = record["first"][day : day + 2]
    start = int(rng.integers(lo, hi - n + 1))
    return record["x"][start : start + n], start


def run_chunk(task):
    cell, indexes, b_wild = task
    a, b, n, cidx = cell
    rows = []
    ndays = DATA[a]["days"].size
    for k in indexes:
        t0 = time.process_time()
        rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(cidx, k)))
        status, error, result = "ok", "", {}
        try:
            for _ in range(1000):
                i, j = (int(v) for v in rng.integers(0, ndays, size=2))
                if abs(i - j) >= MIN_GAP_DAYS:
                    break
            else:
                raise RuntimeError("Could not sample separated days")
            x, sa = draw(a, i, n, rng)
            y, sb = draw(b, j, n, rng)
            result = {
                "day_a": str(DATA[a]["days"][i]), "day_b": str(DATA[b]["days"][j]),
                "start_a": sa, "start_b": sb, "sd_x": float(x.std()), "sd_y": float(y.std()),
            }
            if result["sd_x"] == 0 or result["sd_y"] == 0:
                raise ValueError("Constant duration window")
            result.update(mb.kind_tests(x, y, rng, B=b_wild))
            result.update(mb.lind_tests(x, y))
        except Exception as exc:
            status, error = "error", repr(exc)[:300]
        rows.append({"pair": k, "status": status, "error": error,
                     "cpu_s": round(time.process_time() - t0, 4), **result})
    return cell, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=100)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--b-wild", type=int, default=500)
    parser.add_argument("--n", type=int, nargs="+", default=[500, 1000, 2000])
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    if not 1 <= args.workers <= 6 or args.pairs < 1 or args.b_wild < 10:
        parser.error("Expected 1–6 workers, positive pair count, and at least 10 wild draws")
    if any(n not in (500, 1000, 2000) for n in args.n):
        parser.error("n must be 500, 1000, or 2000")
    load_data()
    args.out.mkdir(parents=True, exist_ok=True)
    cells = [(a, b, n, idx) for idx, (a, b, n) in enumerate(
        (a, b, n) for a, b in PAIRS for n in sorted(set(args.n))
    )]
    manifest = {
        "source": "experiments/STOCK-20261023/duration_extraction_summary.json",
        "method": "different-day exploratory duration placebo; within-day windows",
        "seed": SEED, "min_gap_trading_days": MIN_GAP_DAYS,
        "pairs_per_cell": args.pairs, "wild_draws": args.b_wild,
        "cells": [{"a": a, "b": b, "n": n, "index": idx} for a, b, n, idx in cells],
        "note": "Reused days make sampled pairs dependent; rejection fractions are descriptive.",
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    tasks = []
    for cell in cells:
        a, b, n, _ = cell
        path = args.out / f"{a}-{b}_n{n}.csv"
        done = set()
        if path.exists():
            with path.open(newline="") as stream:
                done = {int(row["pair"]) for row in csv.DictReader(stream)}
        pending = [k for k in range(args.pairs) if k not in done]
        chunk = 5 if n == 2000 else 10
        tasks.extend((cell, pending[i : i + chunk], args.b_wild)
                     for i in range(0, len(pending), chunk))
    print(f"{len(tasks)} chunks, {sum(len(t[1]) for t in tasks)} pairs", flush=True)
    start = time.time()
    with mp.get_context("fork").Pool(args.workers, maxtasksperchild=50) as pool:
        for count, (cell, rows) in enumerate(pool.imap_unordered(run_chunk, tasks), 1):
            a, b, n, _ = cell
            path = args.out / f"{a}-{b}_n{n}.csv"
            fresh = not path.exists()
            with path.open("a", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction="ignore")
                if fresh:
                    writer.writeheader()
                writer.writerows(rows)
            if count % 10 == 0 or count == len(tasks):
                print(f"{count}/{len(tasks)} chunks, {time.time() - start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
