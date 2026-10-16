"""Validation of the generators and helpers before any simulation.

Checks (written to validation.out; exit code 1 on any FAIL):
  1. Davies-Harte embedding eigenvalues nonnegative for every (H, n) used.
  2. fGn sample autocovariance at lags 0,1,2,5,10,50 (known mean, averaged over reps) within
     4 Monte Carlo s.e. of the exact fGn autocovariance, for H in grid + 0.1.
  3. Var(sum_{t<=m} X_t) / m^{2H} within 4 s.e. of 1 for m in {10, 100, n}.
  4. Marginal kurtosis within 4 s.e. of 3.
  5. Student-t 0.975 quantiles against tabulated values; Phi via math.erf against known values.
  6. G2 marginal uniform on (-1,1) (mean 0, var 1/3); G3 balanced signs.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import math  # noqa: E402

import numpy as np  # noqa: E402

import common as C  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "validation.out")
lines = []
fails = 0


def log(s):
    print(s)
    lines.append(s)


def check(name, ok, detail):
    global fails
    log(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        fails += 1


rng = np.random.default_rng(np.random.SeedSequence([C.MASTER_SEED, 999]))

# 1. embedding eigenvalues
for h in (0.1,) + C.HURST_GRID[1:]:
    for n in (500, 1000, 2000, 4000, 8000, 9000):
        s = C.dh_sqrt_eigs(h, n)
        check(f"DH eigen >= 0 H={h} n={n}", bool(np.all(np.isfinite(s))), f"min sqrt-eig {s.min():.3e}")

# 2-4. moments
R = 4000
n = 1000
lags = (0, 1, 2, 5, 10, 50)
for h in (0.1,) + C.HURST_GRID:
    g = C.fgn_autocov(h, n)
    acf = np.zeros((R, len(lags)))
    sums = {10: np.zeros(R), 100: np.zeros(R), n: np.zeros(R)}
    kurt = np.zeros(R)
    for r in range(R):
        x = C.fgn(h, n, rng)
        for i, k in enumerate(lags):
            acf[r, i] = np.dot(x[k:], x[: n - k]) / (n - k)
        for m in sums:
            sums[m][r] = x[:m].sum()
        kurt[r] = np.mean(x ** 4)
    for i, k in enumerate(lags):
        est = acf[:, i].mean()
        se = acf[:, i].std(ddof=1) / math.sqrt(R)
        check(f"fGn acf H={h} lag={k}", abs(est - g[k]) < 4 * se + 1e-12,
              f"est {est:.5f} exact {g[k]:.5f} se {se:.5f}")
    for m, v in sums.items():
        ratio = v.var(ddof=1) / m ** (2 * h)
        se = ratio * math.sqrt(2.0 / (R - 1))
        check(f"fGn Var(S_m)/m^2H H={h} m={m}", abs(ratio - 1) < 4 * se, f"ratio {ratio:.4f} se {se:.4f}")
    est = kurt.mean()
    se = kurt.std(ddof=1) / math.sqrt(R)
    check(f"fGn E X^4 H={h}", abs(est - 3) < 4 * se, f"est {est:.4f} se {se:.4f}")

# 5. quantiles and Phi
for nu, tab in ((5, 2.570582), (10, 2.228139), (25, 2.059539), (30, 2.042272), (63, 1.998341),
                 (160, 1.974902)):
    q = C.t_quantile_975(nu)
    check(f"t_0.975 nu={nu}", abs(q - tab) < 2e-5, f"got {q:.6f} table {tab:.6f}")
for x, v in ((0.0, 0.5), (1.959963984540054, 0.975), (-1.0, 0.15865525393145707)):
    got = float(C.norm_cdf(np.array([x]))[0])
    check(f"Phi({x})", abs(got - v) < 1e-12, f"got {got:.15f}")

# 6. subordinated marginals
x2 = np.concatenate([C.generate("G2", 0.8, 2000, rng) for _ in range(200)])
check("G2 mean", abs(x2.mean()) < 0.03, f"{x2.mean():.4f}")
check("G2 var=1/3", abs(x2.var() - 1 / 3) < 0.01, f"{x2.var():.4f}")
check("G2 range", bool(x2.min() > -1 and x2.max() < 1), f"[{x2.min():.4f},{x2.max():.4f}]")
x3 = np.concatenate([C.generate("G3", 0.8, 2000, rng) for _ in range(200)])
check("G3 values", set(np.unique(x3)) == {-1.0, 1.0}, str(np.unique(x3)))
check("G3 balance", abs(x3.mean()) < 0.05, f"{x3.mean():.4f}")
for lam in (0.01, 0.1):
    y = np.array([C.generate("R1", lam, 2000, rng) for _ in range(100)])
    v = y.var(axis=1).mean()
    log(f"INFO  R1 lambda={lam}: mean within-path var {v:.3f}, lag-1 acf "
        f"{np.mean([np.corrcoef(p[1:], p[:-1])[0, 1] for p in y]):.3f}")

# 7. learners and multiplier (DML / K-hom helpers)
import dml  # noqa: E402
import khom  # noqa: E402

xs = rng.standard_normal(200000)
ds = (rng.random(xs.size) < dml.expit(0.5 * xs)).astype(float)
beta, ok = dml.logistic_irls(xs, ds)
check("IRLS recovers (0, 0.5)", ok and abs(beta[0]) < 0.02 and abs(beta[1] - 0.5) < 0.02, str(beta))
xg = rng.standard_normal(4000)
yg = xg + xg ** 2 + rng.standard_normal(4000)
xt = rng.standard_normal(20000)
pred = dml.gb_fit_predict(xg, yg, xt)
rmse = float(np.sqrt(np.mean((pred - (xt + xt ** 2)) ** 2)))
check("GB fits x + x^2 (rmse under N(0,1), n=4000)", rmse < 0.35, f"rmse {rmse:.3f}")
xg2 = rng.standard_normal(400)
pred2 = dml.gb_fit_predict(xg2, xg2 + xg2 ** 2 + rng.standard_normal(400), xt)
rmse2 = float(np.sqrt(np.mean((pred2 - (xt + xt ** 2)) ** 2)))
check("GB error falls with n (400 -> 4000)", rmse < rmse2, f"{rmse2:.3f} -> {rmse:.3f}")
w = khom.markov_rademacher(20000, 4, 20, rng)
for lag in (1, 20):
    ac = float(np.mean(w[lag:] * w[:-lag]))
    check(f"multiplier autocorr lag {lag}", abs(ac - math.exp(-lag / 20)) < 0.04,
          f"{ac:.3f} vs {math.exp(-lag / 20):.3f}")
check("multiplier values", set(np.unique(w)) == {-1.0, 1.0}, str(np.unique(w)))
# K-hom wild bootstrap statistic equals a^T K a with a = s for unit multipliers
xk = rng.standard_normal(200)
dk = np.abs(xk[:, None] - xk[None, :])
sig = khom.median_heuristic(dk)
iu = np.abs(xk[:, None] - xk[None, :])[np.triu_indices(200, 1)]
check("median heuristic = median over pairs", abs(sig - np.median(iu)) < 1e-12, f"{sig:.6f}")
kk = np.exp(-dk ** 2 / (2 * sig * sig))
direct = kk[:100, :100].mean() + kk[100:, 100:].mean() - 2 * kk[:100, 100:].mean()
sv = np.r_[np.full(100, 0.01), np.full(100, -0.01)]
check("MMD2 V-stat identity", abs(direct - sv @ kk @ sv) < 1e-12, f"{direct:.6e}")

log(f"\nTOTAL FAILS: {fails}")
with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
sys.exit(1 if fails else 0)
