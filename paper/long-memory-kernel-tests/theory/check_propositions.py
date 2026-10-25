"""Numeric checks for theory/propositions.tex (Propositions A and B).

Run from anywhere:  OMP_NUM_THREADS=1 python -I check_propositions.py > check_output.txt
Single process. Uses the paper's fGn generator, rules, and HSIC helpers from kmem.legacy.mbcore
(read-only import). Seeds are fixed, so the output is reproducible.
"""
from __future__ import annotations

import math
import pathlib
import sys
import time

import numpy as np
from scipy.linalg import solve_toeplitz
from scipy.stats import norm

REPO = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from kmem.legacy import mbcore as mb  # noqa: E402

P, A = 50, 1.0
T0 = time.time()


def season(n: int, phase: int) -> np.ndarray:
    return A * np.sin(2.0 * np.pi * (np.arange(n) + phase) / P)


def c_limit(d: int) -> float:
    """Proposition A limit of r_n for phase gap d = U - V."""
    return (A * A / 2) * math.cos(2 * math.pi * d / P) / (1 + A * A / 2)


def hsic_v(x: np.ndarray, y: np.ndarray) -> float:
    """Biased HSIC V-statistic with the paper's median bandwidths (T_n / n)."""
    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    kc = mb.double_center(mb.gauss_gram(dx, mb.median_bandwidth(dx)))
    lc = mb.double_center(mb.gauss_gram(dy, mb.median_bandwidth(dy)))
    return float((kc * lc).sum() / x.size ** 2)


print("=" * 72)
print("Proposition A (i): sample Pearson r_n versus the limit c(U-V)")
print(f"P={P}, a={A}; u=0 fixed, v=-d so that U-V=d; 20 noise paths per cell")
print("=" * 72)
rng = np.random.default_rng(20261025)
print(f"{'H':>4} {'n':>6} {'d':>3} {'c(d)':>8} {'mean r_n':>9} {'sd r_n':>8} {'max|r-c|':>9}")
for H in (0.5, 0.8):
    for n in (2000, 20000):
        for d in (0, 5, 12, 20, 25):
            r = []
            for _ in range(20):
                x = mb.fgn(H, n, rng) + season(n, 0)
                y = mb.fgn(H, n, rng) + season(n, -d)
                r.append(np.corrcoef(x, y)[0, 1])
            r = np.array(r)
            c = c_limit(d)
            print(f"{H:4.1f} {n:6d} {d:3d} {c:8.4f} {r.mean():9.4f} {r.std(ddof=1):8.4f} "
                  f"{np.abs(r - c).max():9.4f}")
print(f"[{time.time() - T0:.0f}s]")

print()
print("P(c(U-V) = 0) for P=50: ", sum(abs(math.cos(2 * math.pi * d / P)) < 1e-12 for d in range(P)) / P)
cs = np.array([abs(c_limit(d)) for d in range(P)])
print("smallest |c(d)| over d:", f"{cs.min():.4f}", "at d =", [d for d in range(P) if abs(cs[d] - cs.min()) < 1e-12])

print()
print("=" * 72)
print("Proposition A (ii): Pearson NW and EWC rejection by phase gap (paper's lind_tests)")
print("H=0.5; n=2000 with 10 paths per gap d=0..49; n=20000 with 2 paths per gap")
print("=" * 72)
rng = np.random.default_rng(7)
for n, reps in ((2000, 10), (20000, 2)):
    rej_nw = np.zeros(P)
    rej_ewc = np.zeros(P)
    for d in range(P):
        for _ in range(reps):
            x = mb.fgn(0.5, n, rng) + season(n, 0)
            y = mb.fgn(0.5, n, rng) + season(n, -d)
            out = mb.lind_tests(x, y)
            rej_nw[d] += out["rej_nw"] == 1
            rej_ewc[d] += out["rej_ewc"] == 1
    rej_nw /= reps
    rej_ewc /= reps
    miss = [d for d in range(P) if rej_nw[d] < 1]
    print(f"n={n}: NW rejection averaged over uniform gap = {rej_nw.mean():.3f}; "
          f"EWC = {rej_ewc.mean():.3f}; gaps with NW rate < 1: {miss}")
print(f"[{time.time() - T0:.0f}s]")

print()
print("=" * 72)
print("Proposition A (iii): HSIC V-statistic T_n/n under SEAS-M with fixed gap")
print("H=0.5 and 0.8, gap d=12 (smallest |c|) and d=0; 3 paths per cell")
print("Wild and subsampling critical values are O(sqrt n); T_n = n V_n grows like n")
print("=" * 72)
rng = np.random.default_rng(11)
for H in (0.5, 0.8):
    for d in (12, 0):
        line = []
        for n in (1000, 2000, 4000):
            v = [hsic_v(mb.fgn(H, n, rng) + season(n, 0), mb.fgn(H, n, rng) + season(n, -d))
                 for _ in range(3)]
            line.append(f"n={n}: V_n={np.mean(v):.5f} (T_n={n * np.mean(v):7.2f}, 4*floor(sqrt n)={4 * math.floor(math.sqrt(n))})")
        print(f"H={H} d={d:2d}  " + " | ".join(line))
# Null reference without seasonality, for scale
for H in (0.5, 0.8):
    v = [hsic_v(mb.fgn(H, 2000, rng), mb.fgn(H, 2000, rng)) for _ in range(3)]
    print(f"no season, H={H}, n=2000: V_n={np.mean(v):.5f}")
print(f"[{time.time() - T0:.0f}s]")


def g_kern(mu: np.ndarray, s2: float, sig: float) -> np.ndarray:
    """E exp(-Z^2 / (2 sig^2)) for Z ~ N(mu, s2)."""
    v = sig * sig + s2
    return sig / np.sqrt(v) * np.exp(-mu * mu / (2 * v))


def mixture_median_dist(cent: np.ndarray) -> float:
    """Median of |X - X'| for X, X' iid from the mixture of N(cent_j, 1), equal weights."""
    from scipy.optimize import brentq
    mu = (cent[:, None] - cent[None, :]).ravel()
    sd = math.sqrt(2.0)
    cdf = lambda x: float(np.mean(norm.cdf((x - mu) / sd) - norm.cdf((-x - mu) / sd))) - 0.5
    return brentq(cdf, 1e-6, 20.0)


def hsic_limit(d: int) -> float:
    """Population HSIC of the limit measure mu_{u,v} (Gaussian mixture), median bandwidths."""
    j = np.arange(P)
    s = A * np.sin(2 * np.pi * j / P)
    q = A * np.sin(2 * np.pi * (j - d) / P)
    sx, sy = mixture_median_dist(s), mixture_median_dist(q)
    gk = g_kern(s[:, None] - s[None, :], 2.0, sx)
    gl = g_kern(q[:, None] - q[None, :], 2.0, sy)
    t1 = np.mean(gk * gl)
    t2 = np.mean(gk.mean(axis=1) * gl.mean(axis=1))
    t3 = gk.mean() * gl.mean()
    return float(t1 - 2 * t2 + t3)


print()
print("Limit HSIC h(d) of Proposition A (iii), closed form for the Gaussian-mixture limit")
hs = np.array([hsic_limit(d) for d in range(P)])
for d in (0, 5, 10, 11, 12, 13, 20, 25):
    print(f"  d={d:2d}: h(d) = {hs[d]:.6f}   n*h(d) at n=2000: {2000 * hs[d]:7.2f}   |c(d)| = {abs(c_limit(d)):.4f}")
print(f"  min_d h(d) = {hs.min():.6f} at d = {int(hs.argmin())}; all positive: {bool((hs > 0).all())}")
print(f"[{time.time() - T0:.0f}s]")

print()
print("=" * 72)
print("Proposition B:Var(mean of first half - mean of second half), n = 2m")
print("=" * 72)
rng = np.random.default_rng(3)
n, m, delta, R = 2000, 1000, 0.1, 10000
z95, z975 = norm.ppf(0.95), norm.ppf(0.975)
for H in (0.5, 0.7, 0.8):
    v_formula = (4 * m ** (2 * H) - (2 * m) ** (2 * H)) / m ** 2
    diffs = np.empty(R)
    for i in range(R):
        x = mb.fgn(H, n, rng)
        diffs[i] = x[:m].mean() - x[m:].mean()
    v_mc = diffs.var(ddof=1)
    se_v = v_mc * math.sqrt(2 / (R - 1))
    ratio = delta / math.sqrt(v_formula)
    asym = delta * m ** (1 - H) / math.sqrt(4 - 2 ** (2 * H))
    # GLS / Neyman-Pearson bound for location-invariant tests
    g = mb.fgn_acov(H, n)
    one = np.ones(n)
    c = np.r_[np.zeros(m), np.ones(m)]
    gi1 = solve_toeplitz(g, one)
    gic = solve_toeplitz(g, c)
    info = c @ gic - (one @ gic) ** 2 / (one @ gi1)
    r_gls = delta * math.sqrt(info)
    pow_lin2 = norm.cdf(ratio - z975) + norm.cdf(-ratio - z975)
    pow_np1 = norm.cdf(r_gls - z95)
    print(f"H={H}: Var formula={v_formula:.6f}  MC={v_mc:.6f} (+-{se_v:.6f})  sd={math.sqrt(v_formula):.4f}")
    print(f"       delta/sd = {ratio:.3f} (closed form {asym:.3f});  GLS ratio delta*sqrt(I) = {r_gls:.3f}")
    print(f"       two-sided half-mean z-test power = {pow_lin2:.3f};  NP bound (one-sided, any "
          f"location-invariant level-5% test) = {pow_np1:.3f}")
print(f"[{time.time() - T0:.0f}s total]")
