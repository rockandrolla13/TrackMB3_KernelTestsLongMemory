"""Stock-duration placebo under D0, D1, D2, copied from STOCK-20261023/run_stock_placebo.py.

Each pair re-seeds the original SeedSequence(SEED, spawn_key=(cell, pair)) once per adjustment.
Window draws depend only on day boundaries, so every adjustment sees the same windows and the
same wild multipliers. D0 must reproduce the original CSVs exactly.
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
ADJ = ("D0", "D1", "D2")
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


def load_data(adjs) -> None:
    for symbol in SYMBOLS:
        with np.load(ROOT / f"data/raw/stock_duration_cache_2024/{symbol}.npz") as z:
            series = {"D0": z["dur"].astype(float)}
            day = z["dur_day"].astype(int)
            days = z["days"].astype(str)
        with np.load(ROOT / f"data/raw/stock_diurnal_cache_2024/{symbol}.npz") as z:
            if not (np.array_equal(z["dur_day"], day) and np.array_equal(z["days"].astype(str), days)):
                raise ValueError(f"{symbol}: diurnal cache misaligned")
            for adj in adjs:
                if adj != "D0":
                    series[adj] = z[adj].astype(float)
        first = np.searchsorted(day, np.arange(days.size + 1))
        if np.any(np.diff(first) < 2000):
            raise ValueError(f"A day has fewer than 2000 durations: {symbol}")
        DATA[symbol] = {"days": days, "first": first, **series}
    if not all(np.array_equal(DATA["AAPL"]["days"], DATA[s]["days"]) for s in SYMBOLS[1:]):
        raise ValueError("Stock trading calendars differ")


def draw(symbol: str, day: int, n: int, rng: np.random.Generator):
    lo, hi = DATA[symbol]["first"][day : day + 2]
    return int(rng.integers(lo, hi - n + 1))


def one_pair(a, b, n, cidx, k, adj, b_wild):
    t0 = time.process_time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(cidx, k)))
    status, error, result = "ok", "", {}
    ndays = DATA[a]["days"].size
    try:
        for _ in range(1000):
            i, j = (int(v) for v in rng.integers(0, ndays, size=2))
            if abs(i - j) >= MIN_GAP_DAYS:
                break
        else:
            raise RuntimeError("Could not sample separated days")
        sa = draw(a, i, n, rng)
        sb = draw(b, j, n, rng)
        x = DATA[a][adj][sa : sa + n]
        y = DATA[b][adj][sb : sb + n]
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
    return {"pair": k, "status": status, "error": error,
            "cpu_s": round(time.process_time() - t0, 4), **result}


def run_chunk(task):
    cell, adj, indexes, b_wild = task
    a, b, n, cidx = cell
    return cell, adj, [one_pair(a, b, n, cidx, k, adj, b_wild) for k in indexes]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=500)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--b-wild", type=int, default=500)
    parser.add_argument("--n", type=int, nargs="+", default=[500, 1000, 2000])
    parser.add_argument("--adj", nargs="+", default=list(ADJ))
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error("Expected 1-4 workers")
    if any(a not in ADJ for a in args.adj):
        parser.error(f"adj must be in {ADJ}")
    load_data(args.adj)
    # Cell indexes must match the original run, which enumerated all three sizes per pairing.
    cells = [(a, b, n, idx) for idx, (a, b, n) in enumerate(
        (a, b, n) for a, b in PAIRS for n in (500, 1000, 2000)) if n in args.n]
    manifest = {
        "source": "experiments/STOCK-20261023 cache plus data/raw/stock_diurnal_cache_2024",
        "seed": SEED, "min_gap_trading_days": MIN_GAP_DAYS,
        "pairs_per_cell": args.pairs, "wild_draws": args.b_wild, "adjustments": args.adj,
        "cells": [{"a": a, "b": b, "n": n, "index": idx} for a, b, n, idx in cells],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"manifest_{'-'.join(args.adj)}_n{'-'.join(map(str, sorted(args.n)))}.json").write_text(json.dumps(manifest, indent=2) + "\n")
    tasks = []
    # Largest windows first, so the Table 7 cells finish early.
    for cell in sorted(cells, key=lambda c: -c[2]):
        a, b, n, _ = cell
        for adj in args.adj:
            path = args.out / adj / f"{a}-{b}_n{n}.csv"
            done = set()
            if path.exists():
                with path.open(newline="") as stream:
                    done = {int(row["pair"]) for row in csv.DictReader(stream)}
            pending = [k for k in range(args.pairs) if k not in done]
            chunk = 5 if n == 2000 else 10
            tasks.extend((cell, adj, pending[i : i + chunk], args.b_wild)
                         for i in range(0, len(pending), chunk))
    print(f"{len(tasks)} chunks, {sum(len(t[2]) for t in tasks)} pair-tests", flush=True)
    start = time.time()
    with mp.get_context("fork").Pool(args.workers, maxtasksperchild=50) as pool:
        for count, (cell, adj, rows) in enumerate(pool.imap(run_chunk, tasks), 1):
            a, b, n, _ = cell
            path = args.out / adj / f"{a}-{b}_n{n}.csv"
            path.parent.mkdir(parents=True, exist_ok=True)
            fresh = not path.exists()
            with path.open("a", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction="ignore")
                if fresh:
                    writer.writeheader()
                writer.writerows(rows)
            if count % 20 == 0 or count == len(tasks):
                print(f"{count}/{len(tasks)} chunks, {time.time() - start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
