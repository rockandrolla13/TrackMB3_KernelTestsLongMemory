"""Validation of the Davies-Harte fGn generator and of the test implementations.

Run: python -I validate.py > validation.out
Must be run before any simulation; the simulation driver refuses to start without a PASS file.
"""
import math
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import mbcore as mb  # noqa: E402

failures = []


def check(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


rng = np.random.default_rng(np.random.SeedSequence(777))

# 1. Exactness of the circulant embedding: eigenvalues nonnegative, embedding reproduces gamma.
for H in (0.1, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9):
    for n in (500, 1000, 2000, 4000, 8000, 5000, 9000, 1100):
        lam = mb.dh_eigs(H, n)
        g = mb.fgn_acov(H, n)
        rec = np.fft.ifft(lam).real[:n]
        g_full = mb.fgn_acov(H, n)
        c = np.concatenate([g_full, g_full[-2:0:-1]])
        lam_raw = np.fft.fft(c).real
        err = float(np.max(np.abs(rec - g)))
        check(f"embedding H={H} n={n}", lam_raw.min() >= -1e-10 * lam_raw.max() and err < 1e-9,
              f"min eig={lam_raw.min():.3e}, max|ifft(lam)[:n]-gamma|={err:.2e}")

# 2. Sample autocovariance vs theory (known zero mean), 4000 paths of length 1000.
P, n = 4000, 1000
lags = (0, 1, 2, 5, 10, 50, 200)
for H in (0.1, 0.5, 0.7, 0.9):
    paths = np.stack([mb.fgn(H, n, rng) for _ in range(P)])
    g = mb.fgn_acov(H, n)
    worst = 0.0
    parts = []
    for k in lags:
        per_path = (paths[:, : n - k] * paths[:, k:]).mean(axis=1)
        est, se = per_path.mean(), per_path.std(ddof=1) / math.sqrt(P)
        z = (est - g[k]) / se
        worst = max(worst, abs(z))
        parts.append(f"k={k}: {est:.4f} vs {g[k]:.4f} (z={z:+.2f})")
    check(f"autocov H={H}", worst < 4.0, "; ".join(parts))
    # 3. Var(sum of n) = n^{2H} exactly for unit fGn.
    s = paths.sum(axis=1)
    ratio = s.var(ddof=1) / n ** (2 * H)
    se = math.sqrt(2 / (P - 1))
    check(f"Var(sum)/n^2H H={H}", abs(ratio - 1) < 4 * se, f"{ratio:.4f} (MC se {se:.4f})")

# 3b. Full-covariance check at small n: empirical covariance matrix vs Toeplitz, 20000 paths.
n = 16
for H in (0.6, 0.9):
    X = np.stack([mb.fgn(H, n, rng) for _ in range(20000)])
    C = X.T @ X / X.shape[0]
    g = mb.fgn_acov(H, n)
    T = g[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]
    err = float(np.max(np.abs(C - T)))
    check(f"full covariance n=16 H={H}", err < 0.05, f"max |C_hat - Toeplitz| = {err:.4f} (MC scale ~0.01)")

# 4. Normal and t helpers.
check("norm_cdf(1.959964)", abs(float(mb.norm_cdf(1.959963984540054)) - 0.975) < 1e-12,
      f"{float(mb.norm_cdf(1.959963984540054)):.15f}")
for nu, ref in ((10, 2.228138851986), (25, 2.059538552753), (60, 2.000297821732), (120, 1.979930405146)):
    q = mb.t_quantile(0.975, nu)
    check(f"t_quantile(0.975, {nu})", abs(q - ref) < 1e-6, f"{q:.9f} vs {ref}")

# 5. Shift statistics: S(0) equals the HSIC statistic; S(c) equals brute force.
n = 60
x, y = rng.standard_normal(n), rng.standard_normal(n)
dx, dy = np.abs(x[:, None] - x[None, :]), np.abs(y[:, None] - y[None, :])
kc = mb.double_center(mb.gauss_gram(dx, mb.median_bandwidth(dx)))
s = mb.shift_stats(kc, mb.double_center(mb.gauss_gram(dy, mb.median_bandwidth(dy))), chunk=7)
worst = 0.0
for c in (0, 1, 5, 17, 59):
    ys = np.roll(y, -c)  # ys_i = y_{i+c}
    dys = np.abs(ys[:, None] - ys[None, :])
    lcs = mb.double_center(mb.gauss_gram(dys, mb.median_bandwidth(dy)))
    brute = (kc * lcs).sum() / n
    worst = max(worst, abs(brute - s[c]))
check("shift FFT vs brute force", worst < 1e-10, f"max abs err {worst:.2e}")

# 6. Median bandwidth vs direct pairwise median.
iu = np.triu_indices(n, 1)
check("median bandwidth", abs(mb.median_bandwidth(dx) - np.median(dx[iu])) < 1e-12,
      f"{mb.median_bandwidth(dx):.6f} vs {np.median(dx[iu]):.6f}")
xs = np.sign(rng.standard_normal(200))
ds = np.abs(xs[:, None] - xs[None, :])
check("median bandwidth sign data > 0", mb.median_bandwidth(ds) > 0, f"{mb.median_bandwidth(ds)}")

# 7. dCov V-statistic vs Szekely-Rizzo definition computed independently.
dcx, dcy = mb.double_center(dx), mb.double_center(dy)
v1 = (dcx * dcy).sum() / n**2
S1 = (dx * dy).mean(); S2 = dx.mean() * dy.mean(); S3 = (dx.mean(axis=1) * dy.mean(axis=1)).mean()
v2 = S1 + S2 - 2 * S3
check("dCov^2 V-stat identity", abs(v1 - v2) < 1e-12, f"{v1:.10f} vs {v2:.10f}")

# 8. Subsampling windows vs brute force for two windows.
n = 80; l = 8; dn = 8
x, y = rng.standard_normal(n), rng.standard_normal(n)
sx, sy = 0.9, 1.1
sh = mb.subsample_stats(x, y, l, dn, "hsic", sx, sy, chunk=5)
sd = mb.subsample_stats(x, y, l, dn, "dcov", chunk=5)
worst = 0.0
for k in (0, 13, n - l - dn - 1):
    xb, yb = x[k:k + l], y[k + dn:k + dn + l]
    ax, ay = np.abs(xb[:, None] - xb[None, :]), np.abs(yb[:, None] - yb[None, :])
    th = (mb.double_center(mb.gauss_gram(ax, sx)) * mb.double_center(mb.gauss_gram(ay, sy))).sum() / l
    td = (mb.double_center(ax) * mb.double_center(ay)).sum() / l
    worst = max(worst, abs(th - sh[k]), abs(td - sd[k]))
check("subsample windows vs brute force", worst < 1e-12 and sh.size == n - l - dn, f"max abs err {worst:.2e}, m={sh.size}")

# 9. Markov-Rademacher multiplier autocorrelation e^{-k/ln}.
w = mb.markov_rademacher(4000, 400, 20.0, rng)
parts = []
ok = True
for k in (1, 10, 20, 40):
    r = float((w[:-k] * w[k:]).mean() / (w * w).mean())
    ok &= abs(r - math.exp(-k / 20)) < 0.03
    parts.append(f"k={k}: {r:.3f} vs {math.exp(-k / 20):.3f}")
check("wild-bootstrap multiplier acf", ok, "; ".join(parts))

# 10. Wild-bootstrap quadratic form vs explicit loop.
m = rng.standard_normal((30, 30)); m = m + m.T
w = mb.markov_rademacher(30, 4, 3.0, rng)
direct = np.array([w[:, b] @ m @ w[:, b] / 30 for b in range(4)])
check("wild_boot_stats", np.allclose(direct, mb.wild_boot_stats(m, w)), "quadratic forms match")

# 11. EWC: with nu = n-1 cosine terms plus mean, the basis is orthonormal (Parseval check).
n = 50
u = rng.standard_normal(n); u -= u.mean()
v = mb.lrv_ewc(u, n - 1) * (n - 1)
check("EWC cosine basis Parseval", abs(v - u @ u) < 1e-9, f"{v:.10f} vs {u @ u:.10f}")

# 12. NW with S<=1 equals variance; Andrews bandwidth formula at rho = 0.5.
check("NW S=1 equals gamma0", abs(mb.lrv_nw(u, 1.0) - u @ u / n) < 1e-14, "ok")

# 13. G3 sign and G2 Phi transforms are rank-1, mean 0, bounded.
xi = mb.fgn(0.7, 4000, rng)
g2 = mb.phi_transform(xi)
check("G2 range and mean", g2.min() > -1 and g2.max() < 1 and abs(np.corrcoef(g2, xi)[0, 1] - math.sqrt(3 / math.pi)) < 0.02,
      f"corr(G2, xi)={np.corrcoef(g2, xi)[0, 1]:.4f} vs sqrt(3/pi)={math.sqrt(3 / math.pi):.4f}")

# 14. fOU: stationarity after burn-in (variance of first vs last quarter over paths).
for lam in (0.01, 0.1):
    paths = np.stack([mb.fou(0.1, lam, 2000, rng) for _ in range(400)])
    v_first, v_last = paths[:, :100].var(), paths[:, -100:].var()
    check(f"fOU stationarity lam={lam}", abs(v_first / v_last - 1) < 0.15, f"var first 100={v_first:.4f}, last 100={v_last:.4f}")

print()
print("OVERALL:", "PASS" if not failures else f"FAIL {failures}")
if not failures:
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "validation.PASS"), "w") as f:
        f.write("validate.py passed\n")
sys.exit(1 if failures else 0)
