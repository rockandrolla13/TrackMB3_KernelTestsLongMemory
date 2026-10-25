"""Stage 2: unit checks of EXP2 generators (marginal, memory, periodicity, half balance).

Usage: python -I unit_checks.py > unit_checks.txt
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import math  # noqa: E402
import sys  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import gens  # noqa: E402
import mmdcal as M  # noqa: E402

SEED = 2026102301
N = 2000
R = 200


def hill(x: np.ndarray, frac: float = 0.05) -> float:
    a = np.sort(np.abs(x))[::-1]
    k = int(frac * a.size)
    return 1.0 / float(np.mean(np.log(a[:k] / a[k])))


def acf(x: np.ndarray, k: int) -> float:
    x = x - x.mean()
    return float(np.dot(x[k:], x[:-k]) / np.dot(x, x))


def season_ratio(x: np.ndarray) -> float:
    """Periodogram at frequency n/P divided by the mean of 10 neighbours on each side."""
    f = np.abs(np.fft.rfft(x - x.mean())) ** 2
    j = x.size // gens.PERIOD
    nb = np.r_[f[j - 11:j - 1], f[j + 2:j + 12]]
    return float(f[j] / nb.mean())


def main():
    m = math.floor(N ** 0.65)
    print(f"EXP2 unit checks: n = {N}, {R} series per row, seed {SEED}, local Whittle m = {m}")
    print("Means over series. hill = Hill tail index on top 5% of |X| (pooled per series, then averaged).")
    print("dmean, dsd = mean |half-1 minus half-2| of mean and sd. seas = periodogram ratio at 1/P.")
    hdr = ("gen", "H", "mean", "sd", "skew", "kurt", "hill", "Hlw_X", "Hlw_X2", "acf1", "acf10", "acf100",
           "acf1|X|", "acf100|X|", "seas", "dmean", "dsd")
    print(" ".join(f"{h:>9}" for h in hdr))
    gl = list(gens.MAIN) + ["SEAS-M-shared"]
    for gi, g in enumerate(gl):
        for hi, H in enumerate((0.5, 0.7, 0.8)):
            rows = []
            for r in range(R):
                rng = np.random.default_rng(np.random.SeedSequence(entropy=SEED, spawn_key=(gi, hi, r)))
                x = gens.generate(g, H, N, rng)
                z = (x - x.mean()) / x.std()
                h1 = N // 2
                rows.append([x.mean(), x.std(), np.mean(z ** 3), np.mean(z ** 4), hill(x),
                             M.local_whittle_d(x, m) + 0.5, M.local_whittle_d(x * x, m) + 0.5,
                             acf(x, 1), acf(x, 10), acf(x, 100), acf(np.abs(x), 1), acf(np.abs(x), 100),
                             season_ratio(x), abs(x[:h1].mean() - x[h1:].mean()), abs(x[:h1].std() - x[h1:].std())])
            v = np.mean(np.array(rows), axis=0)
            print(f"{g:>9} {H:>9}" + "".join(f" {u:9.3f}" for u in v), flush=True)
    # identity checks
    rng = np.random.default_rng(SEED)
    xi = rng.standard_normal(10 ** 6)
    t = gens.t3_from_gauss(xi)
    print("\nHT-S transform check (1e6 iid draws): var =", round(float(t.var()), 3),
          "(t3/sqrt3 has var 1; heavy tail makes this noisy)",
          "| P(|X|>3/sqrt3) =", round(float(np.mean(np.abs(t) > 3 / math.sqrt(3))), 5),
          "vs exact", round(float(2 * __import__('scipy').stats.t.sf(3, 3)), 5))
    print("HT-S oddness: max |G(x)+G(-x)| =", float(np.max(np.abs(gens.t3_from_gauss(xi[:1000]) + gens.t3_from_gauss(-xi[:1000])))))
    for H in (0.7, 0.8):
        print(f"NL-SV s_H at H={H}: {gens.sv_sd(H):.4f} (H=0.5: {gens.sv_sd(0.5):.4f}, should be 1)")
        psi = gens.arfima_psi(H - 0.5, 20 * N)
        print(f"HT-L psi at H={H}: sum psi^2 = {float(psi @ psi):.4f}, psi_0^2 share = {psi[0] ** 2:.3f}")
    print("Whole periods per half: n/2/P =", 500 // 2 / gens.PERIOD, "and", N // 2 / gens.PERIOD)


if __name__ == "__main__":
    main()
