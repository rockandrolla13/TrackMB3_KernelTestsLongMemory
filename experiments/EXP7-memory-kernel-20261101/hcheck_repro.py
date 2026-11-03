"""EXP7 check of paper/long-memory-kernel-tests/study-notes/hermite-centring-check/ (DESIGN.md A2).

Statistic A: hcheck.py as written. n * MMD^2 on normal scores, K on the scores, bandwidth^2 = median(d^2)/2, Q removes He_0..He_{m-1}.
Statistic B: the EXP7 herm form. K on raw x with the median bandwidth from raw x, Q from the normal scores (m = r + 1).
Data: fGn; SV as in hcheck_sv.py, exp(h/2) e with h = fGn(H); and NL-SV from kmem.legacy.gens (the EXP2 generator).
Inflation = mean statistic at H divided by mean statistic at H = 0.5. n = 800, 2,000 replications per cell, seed 2026110175.
Usage: python -I hcheck_repro.py --workers 4 -> hcheck_repro.csv
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import numpy as np  # noqa: E402
from numpy.polynomial.hermite_e import hermevander  # noqa: E402
from scipy.stats import norm, rankdata  # noqa: E402

from kmem.legacy import gens, mbcore as mb  # noqa: E402

N, MS, HS, DATA = 800, (1, 2, 3, 4), (0.5, 0.7, 0.8), ("fGn", "SV-hcheck", "NL-SV")
MASTER = 2026110175


def series(data, H, rng):
    if data == "fGn":
        return mb.fgn(H, N, rng)
    if data == "SV-hcheck":
        return np.exp(0.5 * mb.fgn(H, N, rng)) * rng.standard_normal(N)
    return gens.generate("NL-SV", H, N, rng)


def stats_one(x):
    n = x.size
    z = norm.ppf(rankdata(x) / (n + 1))
    a = np.r_[np.ones(n // 2), -np.ones(n - n // 2)] * (2 / n)
    d2 = (z[:, None] - z[None, :]) ** 2
    s2 = np.median(d2[np.triu_indices(n, 1)])
    ka = np.exp(-d2 / (2 * s2 / 2))
    dx = np.abs(x[:, None] - x[None, :])
    kb = mb.gauss_gram(dx, mb.median_bandwidth(dx))
    zb = norm.ppf((rankdata(x) - 0.5) / n)
    out = {}
    for m in MS:
        qa, _ = np.linalg.qr(hermevander(z, m - 1))
        pa = a - qa @ (qa.T @ a)
        out[("A", m)] = n * pa @ ka @ pa
        qb, _ = np.linalg.qr(hermevander(zb, m - 1))
        pb = a - qb @ (qb.T @ a)
        out[("B", m)] = n * pb @ kb @ pb
    return out


def run(args):
    data, H, c, reps = args
    rng = np.random.default_rng(np.random.SeedSequence(entropy=MASTER, spawn_key=(DATA.index(data), int(H * 100), c)))
    acc = {}
    for _ in range(reps):
        for k, v in stats_one(series(data, H, rng)).items():
            acc.setdefault(k, []).append(v)
    return data, H, acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--reps", type=int, default=2000)
    a = ap.parse_args()
    tasks = [(d, H, c, 100) for d in DATA for H in HS for c in range(a.reps // 100)]
    store = {}
    with mp.get_context("fork").Pool(a.workers) as pool:
        for data, H, acc in pool.imap_unordered(run, tasks):
            for k, v in acc.items():
                store.setdefault((data, H, *k), []).extend(v)
    rows = []
    for data in DATA:
        for st in ("A", "B"):
            for m in MS:
                base = np.array(store[(data, 0.5, st, m)])
                row = {"data": data, "statistic": st, "m": m, "mean_H0.5": base.mean()}
                for H in (0.7, 0.8):
                    v = np.array(store[(data, H, st, m)])
                    rng = np.random.default_rng(42)
                    b = v[rng.integers(0, v.size, (2000, v.size))].mean(axis=1) / base[rng.integers(0, base.size, (2000, base.size))].mean(axis=1)
                    row[f"inflation_H{H}"] = v.mean() / base.mean()
                    row[f"lo_H{H}"], row[f"hi_H{H}"] = float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))
                rows.append(row)
    with open(os.path.join(HERE, "hcheck_repro.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['data']:10s} {r['statistic']} m={r['m']}  H=0.7 x{r['inflation_H0.7']:.2f} [{r['lo_H0.7']:.2f}, {r['hi_H0.7']:.2f}]  "
              f"H=0.8 x{r['inflation_H0.8']:.2f} [{r['lo_H0.8']:.2f}, {r['hi_H0.8']:.2f}]")


if __name__ == "__main__":
    main()
