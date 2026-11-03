"""Seeded path check for C2: does estimating F spoil the Hermite deflation?

This file IS Monte Carlo, unlike theory_checks.py. Nothing here enters theory_exact.csv.
Estimated normal scores have no exact expectation, so this is the only available check.

For independent fGn pairs it compares three things.
1. The Dehling-Taqqu reduction: z_t against (xi_t - xibar)/s.
2. The rank-one coefficient left in Q phi against the true driver, on the scale of the raw coefficient.
3. Mean T_n = tr(Qx Kx Qx Qy Ky Qy)/n against the exact E T_n with T known and sample centring.

Run:  OMP_NUM_THREADS=2 python -I theory_path_check.py > theory_path_check.out
      OMP_NUM_THREADS=2 python -I theory_path_check.py --note > theory_note_repro.out
"""
from __future__ import annotations

import math
import pathlib
import sys

import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import theory_checks as tc  # noqa: E402

SEED = 42
REPS = 60
SIG = tc.SIG


def fgn_path(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Davies-Harte circulant embedding. The embedding is nonnegative for fGn with H >= 1/2."""
    g = tc.rho_all(H, n)
    lam = np.fft.rfft(np.r_[g, g[-2:0:-1]]).real
    lam = np.maximum(lam, 0.0)
    m = 2 * n
    w = rng.standard_normal(m // 2 + 1) + 1j * rng.standard_normal(m // 2 + 1)
    w[0] = w[0].real * math.sqrt(2)
    w[-1] = w[-1].real * math.sqrt(2)
    return np.fft.irfft(np.sqrt(lam / 2) * w * math.sqrt(m), m)[:n] * 1.0


def scores(x: np.ndarray) -> np.ndarray:
    return ndtri((rankdata(x) - 0.5) / x.size)


def qkq(x: np.ndarray, z: np.ndarray, r: int) -> tuple[np.ndarray, np.ndarray]:
    K = np.exp(-((x[:, None] - x[None, :]) ** 2) / (2 * SIG ** 2))
    B = tc.he_norm(z, r).T
    Bq, _ = np.linalg.qr(B)
    KB = K @ Bq
    out = K - KB @ Bq.T - Bq @ KB.T + Bq @ (Bq.T @ KB) @ Bq.T
    return out, K


def main() -> None:
    rng = np.random.default_rng(SEED)
    lay = tc.gauss_layers()
    print("H     n      r  reduction_max  resid_rank1/raw  mean_T(est)  se      E_T(known T)  ratio")
    for H in (0.7, 0.8, 0.9):
        for n in (1000, 2000, 4000):
            rho = tc.rho_all(H, n + 60)
            for r in (1, 2):
                kr = lay.kappa(rho, r)
                exact = tc.et_centred(kr, kr, [n])[n]
                vals, red, res = [], [], []
                for _ in range(REPS):
                    xi, eta = fgn_path(H, n, rng), fgn_path(H, n, rng)
                    zx, zy = scores(xi), scores(eta)
                    red.append(float(np.abs(zx - (xi - xi.mean()) / xi.std()).max()))
                    Qx, Kx = qkq(xi, zx, r)
                    Qy, _ = qkq(eta, zy, r)
                    vals.append(float(np.sum(Qx * Qy) / n))
                    u = (xi - xi.mean()) / n
                    res.append(float(math.sqrt(max(u @ Qx @ u, 0.0)) / math.sqrt(max(u @ Kx @ u, 1e-300))))
                v = np.array(vals)
                print(f"{H:<4}  {n:<5d}  {r}  {np.median(red):13.3f}  {np.median(res):15.4f}  {v.mean():11.4f}  {v.std(ddof=1) / math.sqrt(REPS):6.4f}  {exact:12.4f}  {v.mean() / exact:5.3f}", flush=True)


def note_stat(x: np.ndarray, m: int) -> float:
    """Statistic of study-notes/hermite-centring-check/hcheck.py, unchanged in substance."""
    n = x.size
    z = ndtri(rankdata(x) / (n + 1))
    d2 = (z[:, None] - z[None, :]) ** 2
    s2 = np.median(d2[np.triu_indices(n, 1)])
    K = np.exp(-d2 / s2)
    B = np.polynomial.hermite_e.hermevander(z, m - 1)
    Q, _ = np.linalg.qr(B)
    a = np.r_[np.ones(n // 2), -np.ones(n - n // 2)] * (2 / n)
    pa = a - Q @ (Q.T @ a)
    return float(n * pa @ K @ pa)


def reproduce_note(reps: int = 600) -> None:
    """Means and standard errors for the note's table at n=800, with four times its replications."""
    rng = np.random.default_rng(20261024)
    n = 800
    for kind in ("fGn", "SV"):
        for m in (1, 2, 3, 4):
            means, ses = {}, {}
            for H in (0.5, 0.7, 0.8):
                vals = []
                for _ in range(reps):
                    xi = fgn_path(H, n, rng)
                    x = xi if kind == "fGn" else np.exp(0.5 * xi) * rng.standard_normal(n)
                    vals.append(note_stat(x, m))
                v = np.array(vals)
                means[H], ses[H] = v.mean(), v.std(ddof=1) / math.sqrt(reps)
            infl = {H: means[H] / means[0.5] for H in (0.7, 0.8)}
            err = {H: infl[H] * math.sqrt((ses[H] / means[H]) ** 2 + (ses[0.5] / means[0.5]) ** 2) for H in (0.7, 0.8)}
            print(f"{kind:4s} m={m}: inflation H=0.7 x{infl[0.7]:.2f} (se {err[0.7]:.2f}); H=0.8 x{infl[0.8]:.2f} (se {err[0.8]:.2f})", flush=True)


if __name__ == "__main__":
    if "--note" in sys.argv:
        reproduce_note()
    else:
        main()
