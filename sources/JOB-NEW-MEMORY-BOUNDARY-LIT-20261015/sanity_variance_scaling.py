"""Formula sanity check (no simulation, no data).

For fractional Gaussian noise with Hurst H, rho(k) = 0.5(|k+1|^2H - 2|k|^2H + |k-1|^2H).
Computes n * Var(sample mean) exactly from autocovariances for three statistics:
  rank1:  mean of X_t                      (Hermite rank 1; e.g. a DML/AIPW state term tau(X_t)-theta)
  rank2:  mean of X_t^2 - 1               (Hermite rank 2; e.g. sample variance/covariance)
  cross:  mean of X_t * Y_t, X indep of Y (Hermite rank (1,1); e.g. HSIC/dCov cross term)
If n*Var stays bounded as n grows, sqrt(n) normal scaling holds; if it grows like n^a, it fails.
"""
import numpy as np

def rho(k, H):
    k = np.abs(k).astype(float)
    return 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))

def nvar(acf_fn, n):
    k = np.arange(1, n)
    return acf_fn(np.array([0]))[0] + 2 * np.sum((1 - k / n) * acf_fn(k))

ns = [10**3, 10**4, 10**5, 10**6]
print("H_X  H_Y | stat  | n*Var at n=1e3,1e4,1e5,1e6 | log10 growth per decade (last)")
for HX, HY in [(0.6, 0.6), (0.7, 0.7), (0.75, 0.75), (0.8, 0.8), (0.9, 0.9), (0.7, 0.9), (0.6, 0.9)]:
    rows = {
        "rank1": lambda k: rho(k, HX),
        "rank2": lambda k: 2 * rho(k, HX) ** 2,
        "cross": lambda k: rho(k, HX) * rho(k, HY),
    }
    for name, f in rows.items():
        v = [nvar(f, n) for n in ns]
        slope = np.log10(v[-1] / v[-2])
        print(f"{HX:.2f} {HY:.2f} | {name:5s} | " + " ".join(f"{x:9.3f}" for x in v) + f" | {slope:+.3f}")
print("Theory: rank1 slope 2H-1; rank2 slope max(0,4H-3); cross slope max(0,2HX+2HY-3) (log terms at boundary).")
