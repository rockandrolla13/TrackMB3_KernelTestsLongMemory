"""EXP7 critical-value table for the self-normalised tests of the removed Hermite layers (DESIGN.md A2).

Usage: python -I sn_tables.py --workers 6  -> sn_table.json

For each n and each memory h on K.SN_GRID, the quantiles at K.SN_LEVELS of three statistics, on normal scores:
  split1: G_n(n//2)^2 of He_1(z), z from fGn(h)
  split2: G_n(n//2)^2 of He_2(z), z from fGn(H1), H1 = (h + 1)/2 above 1/2, so He_2 has memory h. H1 = 1/2 at h = 1/2.
  prod:   T of He_1(z_x) He_1(z_y), independent fGn(H1) with the same H1, so the product has memory h.
Seeds: SeedSequence(entropy=2026110174, spawn_key=(n, grid index, chunk)). Reads no experiment output.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import numpy as np  # noqa: E402

import exp7_core as K  # noqa: E402
from kmem.legacy import mbcore as mb  # noqa: E402

MASTER = 2026110174
NS = (500, 2000, 4000)
CHUNK = 500


def h1_of(h: float) -> float:
    return 0.5 if h == 0.5 else (h + 1.0) / 2.0


def chunk(args):
    n, ih, c, reps = args
    h = K.SN_GRID[ih]
    rng = np.random.default_rng(np.random.SeedSequence(entropy=MASTER, spawn_key=(n, ih, c)))
    out = np.empty((reps, 3))
    for i in range(reps):
        out[i, 0] = K.sn_split_stat(K.hermite_layers(mb.fgn(h, n, rng), 1)[:, 0])
        a = K.hermite_layers(mb.fgn(h1_of(h), n, rng))
        b = K.hermite_layers(mb.fgn(h1_of(h), n, rng), 1)
        out[i, 1] = K.sn_split_stat(a[:, 1])
        out[i, 2] = K.sn_mean_stat(a[:, 0] * b[:, 0])
    return n, ih, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--reps", type=int, default=40000)
    a = ap.parse_args()
    tasks = [(n, ih, c, CHUNK) for n in NS for ih in range(len(K.SN_GRID)) for c in range(a.reps // CHUNK)]
    store: dict = {}
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool:
        for k, (n, ih, out) in enumerate(pool.imap_unordered(chunk, tasks)):
            store.setdefault((n, ih), []).append(out)
            if (k + 1) % 200 == 0:
                print(f"{k + 1}/{len(tasks)} chunks, {time.time() - t0:.0f}s", flush=True)
    q = {}
    for j, kind in enumerate(("split1", "split2", "prod")):
        for n in NS:
            for lev in K.SN_LEVELS:
                q[f"{kind}|{n}|{lev}"] = [float(np.quantile(np.vstack(store[(n, ih)])[:, j], lev)) for ih in range(len(K.SN_GRID))]
    with open(os.path.join(HERE, "sn_table.json"), "w") as f:
        json.dump({"master": MASTER, "reps_per_point": a.reps, "grid": K.SN_GRID, "levels": K.SN_LEVELS, "quantiles": q}, f, indent=0)
    print(f"DONE {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
