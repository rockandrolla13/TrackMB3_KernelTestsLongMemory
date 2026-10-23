"""Checks for mmdcal.py against brute-force definitions. Writes validate_mmdcal.PASS on success.

Usage: python -I validate_mmdcal.py
"""
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import math  # noqa: E402
import sys  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import mmdcal as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(1)
fails = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name, info)
    if not ok:
        fails.append(name)


# 1. SN statistic vs brute force in feature space (linear kernel -> feature = x itself; Gaussian via Gram)
n = 37
x = rng.standard_normal(n)
d = np.abs(x[:, None] - x[None, :])
k = np.exp(-d * d / 2)
p = M.prefix2(k)
m = n // 2
r = n - m


def nrm(a):
    return float(a @ k @ a)


def ind(a, b):
    v = np.zeros(n)
    v[a:b] = 1.0
    return v


s = np.concatenate([np.full(m, 1 / m), np.full(r, -1 / r)])
mmd = nrm(s)
v = 0.0
for kk in range(1, m + 1):
    v += nrm(ind(0, kk) - kk / m * ind(0, m))
for kk in range(m + 1, n + 1):  # 1-based k
    v += nrm(ind(kk - 1, n) - (n - kk + 1) / r * ind(m, n))
t_bf = n * mmd / (v / n ** 2)
t_fast = M.sn_stat(p, np.array([0]), n)[0]
check("sn_stat full sample", abs(t_bf - t_fast) < 1e-8 * abs(t_bf), f"{t_bf:.10g} {t_fast:.10g}")

# windowed SN equals SN on the sub-series with the same Gram entries
l = 15
st = np.arange(n - l + 1)
tw = M.sn_stat(p, st, l)
ok = True
for s0 in st:
    ksub = k[s0:s0 + l, s0:s0 + l]
    t1 = M.sn_stat(M.prefix2(ksub), np.array([0]), l)[0]
    ok &= abs(t1 - tw[s0]) < 1e-8 * abs(t1)
check("sn_stat windows", ok)

# 2. window_mmd vs direct
wm = M.window_mmd(p, n, l)
ok = True
for s0 in range(n - l + 1):
    m1 = l // 2
    sv = np.zeros(n)
    sv[s0:s0 + m1] = 1 / m1
    sv[s0 + m1:s0 + l] = -1 / (l - m1)
    ok &= abs(l * nrm(sv) - wm[s0]) < 1e-9
check("window_mmd", ok)

# 3. multipliers: csg chains independent across halves, centred; unc not centred; corr ~ exp(-1/l)
w = M.multipliers("csg", 400, 200, 10, rng, b=4000)
check("csg centred", np.allclose(w[:200].mean(0), 0) and np.allclose(w[200:].mean(0), 0))
c_cross = np.mean(w[199] * w[200])
check("csg halves independent at the boundary", abs(c_cross) < 0.06, f"corr={c_cross:.3f}")
w = M.multipliers("unc", 400, 200, 10, rng, b=4000)
c1 = np.mean(w[50] * w[51])
check("unc lag-1 corr ~ exp(-1/10)", abs(c1 - math.exp(-0.1)) < 0.03, f"{c1:.3f} vs {math.exp(-0.1):.3f}")
w = M.multipliers("old", 400, 200, 10, rng, b=4000)
c_cross = np.mean(w[199] * w[200])
check("old chain crosses the boundary", c_cross > 0.7, f"corr={c_cross:.3f}")

# 4. local Whittle / GPH on fGn recover H roughly (n = 4096, 40 reps)
for h in (0.5, 0.7, 0.9):
    est_lw, est_g = [], []
    for _ in range(40):
        z = M.generate("G1", h, 4096, rng)
        mm = math.floor(4096 ** 0.65)
        est_lw.append(M.local_whittle_d(z, mm) + 0.5)
        est_g.append(M.gph_d(z, mm) + 0.5)
    check(f"LW mean H at {h}", abs(np.mean(est_lw) - h) < 0.04, f"{np.mean(est_lw):.3f} sd {np.std(est_lw):.3f}")
    check(f"GPH mean H at {h}", abs(np.mean(est_g) - h) < 0.06, f"{np.mean(est_g):.3f} sd {np.std(est_g):.3f}")

# 5. processes: G5, G6 have Hermite rank 2 -> lag-1 autocorr of X ~ rho^2 (G5) for fGn
h = 0.8
rho1 = 0.5 * (2 ** (2 * h) - 2)
acs = []
for _ in range(200):
    z = M.generate("G5", h, 2000, rng)
    z = z - z.mean()
    acs.append(np.dot(z[1:], z[:-1]) / np.dot(z, z))
check("G5 lag-1 autocorr ~ rho1^2", abs(np.mean(acs) - rho1 ** 2) < 0.03, f"{np.mean(acs):.3f} vs {rho1 ** 2:.3f}")

# 6. rep_e2 / rep_h5 run and are finite
o = M.rep_e2(M.generate("G1", 0.7, 300, rng), rng)
check("rep_e2 finite", all(np.isfinite(v) for v in o.values()), str(sorted(o)))
o = M.rep_h5(M.generate("G1", 0.5, 300, rng), rng)
check("rep_h5 finite", all(np.isfinite(v) for v in o.values()))

if fails:
    print("FAILED:", fails)
    sys.exit(1)
open(os.path.join(HERE, "validate_mmdcal.PASS"), "w").write("all checks passed\n")
print("ALL PASS")
