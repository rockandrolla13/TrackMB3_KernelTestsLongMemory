"""Paired H=0.5 null check of centred and uncentred wild multipliers.

Only the multiplier centering changes. Each pair and its Markov multiplier draws
are shared between the two versions. This is a finite-sample control, not a proof
that either version works under long memory or on the trade-window placebos.
"""

from __future__ import annotations

import os

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import argparse
import csv
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXP1 = REPO / "experiments" / "EXP1-20261016"
sys.path.insert(0, str(EXP1))
import mbcore as mb  # noqa: E402

SEED = 2026102301
B = 500
FIELDS = ("n", "rep", "p_hsic_centred", "p_hsic_raw", "p_dcov_centred", "p_dcov_raw")


def multiplier_pair(n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    length = math.floor(math.sqrt(n))
    flip_probability = 0.5 * (1 - math.exp(-1 / length))
    first = np.where(rng.random(B) < 0.5, -1.0, 1.0)
    flips = (rng.random((n - 1, B)) < flip_probability).astype(np.int8)
    parity = np.concatenate((np.zeros((1, B), dtype=np.int64),
                             np.cumsum(flips, axis=0, dtype=np.int64)), axis=0)
    raw = first[None, :] * (1 - 2 * (parity & 1))
    return raw - raw.mean(axis=0, keepdims=True), raw


def pvalue(observed: float, samples: np.ndarray) -> float:
    return float((1 + np.count_nonzero(samples >= observed)) / (B + 1))


def one(task: tuple[int, int]) -> dict[str, float | int]:
    n, rep = task
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(n, rep)))
    x, y = mb.fgn(0.5, n, rng), mb.fgn(0.5, n, rng)
    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    kc = mb.double_center(mb.gauss_gram(dx, mb.median_bandwidth(dx)))
    lc = mb.double_center(mb.gauss_gram(dy, mb.median_bandwidth(dy)))
    ac, bc = mb.double_center(dx), mb.double_center(dy)
    centred, raw = multiplier_pair(n, rng)
    out: dict[str, float | int] = {"n": n, "rep": rep}
    for label, product in (("hsic", kc * lc), ("dcov", ac * bc)):
        observed = float(product.sum() / n)
        out[f"p_{label}_centred"] = pvalue(observed, mb.wild_boot_stats(product, centred))
        out[f"p_{label}_raw"] = pvalue(observed, mb.wild_boot_stats(product, raw))
    return out


def wilson(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return mid - half, mid + half


def summarize(rows: list[dict[str, str]]) -> dict:
    result: dict[str, object] = {"seed": SEED, "B": B, "null": "independent H=0.5 Gaussian fGn"}
    for n in sorted({int(row["n"]) for row in rows}):
        subset = [row for row in rows if int(row["n"]) == n]
        cells = {}
        for column in FIELDS[2:]:
            count = sum(float(row[column]) <= 0.05 for row in subset)
            lo, hi = wilson(count, len(subset))
            cells[column] = {"rejections": count, "replications": len(subset),
                             "rate": count / len(subset), "wilson95": [lo, hi]}
        result[str(n)] = cells
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", default="500,1000,2000")
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    ns = [int(part) for part in args.n.split(",")]
    path = HERE / "control_pairs.csv"
    existing: dict[tuple[int, int], dict[str, str]] = {}
    if path.exists():
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                existing[(int(row["n"]), int(row["rep"]))] = row
    tasks = [(n, rep) for n in ns for rep in range(args.reps) if (n, rep) not in existing]
    if tasks:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for row in pool.map(one, tasks, chunksize=4):
                existing[(int(row["n"]), int(row["rep"]))] = {key: str(row[key]) for key in FIELDS}
    rows = [existing[key] for key in sorted(existing)]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    summary = summarize(rows)
    summary["inputs"] = {
        str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__), EXP1 / "mbcore.py")
    }
    summary["output_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (HERE / "control_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"new_replications": len(tasks), "summary": summary}, indent=2), flush=True)


if __name__ == "__main__":
    main()
