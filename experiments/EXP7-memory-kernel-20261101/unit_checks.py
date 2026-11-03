"""EXP7 unit checks. Usage: python -I unit_checks.py  -> prints, exit 1 on failure.

1. Generators equal EXP6 (G0, lin, sq), gens.pair (process arms) and EXP4 (SEAS-V-slow) on the same seed.
2. raw_exact equals EXP6 kind_bw at kappa = 1 field by field on the same seed, difference 0.
3. raw_exact equals saved EXP6 main rows, same master and cell index, difference 0.
4. dCov and Pearson columns equal EXP6 kind_bw.
5. Fast path equals the full path on every kernel and several generators. The shift function at zero
   shift equals T on every transformed Gram matrix.
6. Split MMD: raw exact equals EXP3 mmd_ext (p_unc_05, r_unc_05); fast path equals explicit Gram matrices.
7. The increment kernel telescopes in split MMD, so it carries no information there.
8. growth_bundle equals the full T at n = 2000 and n = 8000.
9. Stress generators have the stated moments.
10. Self-normalised statistics: equal mmdcal.sn_stat on a linear kernel, reproduce Betken (2016) Table 1 and Lobato's 45.4.
11. Harris (2020): V_s on the fast path equal explicit segment Gram matrices, and the ratio matches the counting rule.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import csv  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
EXP3 = os.path.join(ROOT, "experiments", "EXP3-power-20261025")
EXP4 = os.path.join(ROOT, "experiments", "EXP4-deseason-20261025")
EXP6 = os.path.join(ROOT, "experiments", "EXP6-bandwidth-20261025")
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.append(EXP6)
sys.path.append(EXP3)
sys.path.append(EXP4)
import numpy as np  # noqa: E402

import bw_core as B6  # noqa: E402
import exp7_core as K  # noqa: E402
import gens_exp4 as G4  # noqa: E402  (EXP4, imports its own gens/mbcore copies)
import power_core as P3  # noqa: E402  (EXP3, imports its own gens/mbcore copies)
from kmem import memory_kernels as mk  # noqa: E402
from kmem.legacy import gens, mbcore as mb  # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {name} {detail}", flush=True)
    if not ok:
        fails.append(name)


def rng(*key):
    return np.random.default_rng(np.random.SeedSequence(entropy=20261101, spawn_key=tuple(key)))


def same(a, b):
    return float(np.max(np.abs(a - b)))


# ---------------------------------------------------------------- 1. generators
for alt, rho in (("null", 0.0), ("lin", 0.065), ("sq", 0.075)):
    a = K.ind_pair("G0", alt, 0.8, 0.8, 400, rho, rng(1))
    b = B6.ind_pair("G0", alt, 0.8, 400, rho, rng(1))
    check(f"ind_pair G0 {alt} equals EXP6", same(a[0], b[0]) == 0 and same(a[1], b[1]) == 0)
for g in ("HT-S", "HT-L", "NL-SV", "NL-TAR", "SEAS-M", "SEAS-V", "SEAS-M-shared"):
    a = K.ind_pair(g, "null", 0.8, 0.8, 400, 0.0, rng(2))
    b = gens.pair(g, 0.8, 400, rng(2))
    check(f"ind_pair {g} null equals gens.pair", same(a[0], b[0]) == 0 and same(a[1], b[1]) == 0)
a = K.ind_pair("SEAS-V-slow", "null", 0.8, 0.8, 400, 0.0, rng(3))
b = G4.pair("SEAS-V-slow", 0.8, 400, rng(3))
check("SEAS-V-slow equals EXP4 gens_exp4", same(a[0], b[0]) == 0 and same(a[1], b[1]) == 0)

# ---------------------------------------------------------------- 2. raw_exact and dCov vs EXP6 kind_bw
EXACT_MAP = {"hsic_stat": "t", "sx": "sx", "sy": "sy", "p_hsic_wb": "p_wb", "r_hsic_wb": "r_wb", "p_hsic_shift": "p_shift",
             "r_hsic_shift": "r_shift", "p_hsic_sub": "p_sub", "rej_hsic_sub": "rej_sub", "r_hsic_sub": "r_sub",
             "shift_check": "check"}
for n, gen, alt, rho, H in ((300, "G0", "lin", 0.1, 0.5), (300, "SEAS-V", "null", 0.0, 0.5), (2000, "G0", "null", 0.0, 0.8)):
    r6, r7 = rng(4), rng(4)
    x, y = B6.ind_pair(gen, alt, H, n, rho, r6)
    B6.ind_pair(gen, alt, H, n, rho, r7)
    want = B6.kind_bw(x, y, r6, kappas=(1.0,))
    got = K.ind_bundle(x, y, r7, kernels=("raw",), exact=True, hsplit=False)
    worst = max(abs(float(want[f"{a}_k1"]) - float(got[f"xraw_{b}"])) for a, b in EXACT_MAP.items())
    check(f"raw_exact equals EXP6 hsic_fields n={n} {gen} {alt}", worst == 0.0, f"maxdiff={worst}")
    for key in ("dcov_stat", "r_dcov_wb", "r_dcov_sub"):
        d = abs(float(want[key]) / float(got[key]) - 1.0)
        check(f"{key} equals EXP6 n={n} {gen} {alt}", d < 1e-9, f"rel diff={d:.1e}")
    for key in ("p_dcov_wb", "p_dcov_sub", "rej_dcov_sub", "t_nw", "t_ewc", "rej_nw", "rej_ewc", "r_nw", "r_ewc", "pearson_r"):
        check(f"{key} equals EXP6 n={n} {gen} {alt}", float(want[key]) == float(got[key]))

# ---------------------------------------------------------------- 3. saved EXP6 main rows
EXP6_MASTER = 2026102562
for idx, name, gen, H in ((0, "ind__G0__null0__H0.5__n2000", "G0", 0.5), (1, "ind__G0__null0__H0.8__n2000", "G0", 0.8),
                          (2, "ind__SEAS-V__null0__H0.5__n2000", "SEAS-V", 0.5)):
    saved = {int(r["rep"]): r for r in csv.DictReader(open(os.path.join(EXP6, "main", name + ".csv")))}
    worst, nv = 0.0, 0
    for rep in range(3):
        r = np.random.default_rng(np.random.SeedSequence(entropy=EXP6_MASTER, spawn_key=(idx, rep)))
        x, y = K.ind_pair(gen, "null", H, H, 2000, 0.0, r)
        got = K.ind_bundle(x, y, r, kernels=("raw",), exact=True, hsplit=False)
        for a, b in EXACT_MAP.items():
            worst = max(worst, abs(float(saved[rep][f"{a}_k1"]) - float(got[f"xraw_{b}"])))
            nv += 1
    check(f"raw_exact equals saved EXP6 rows {name} (3 reps, {nv} values)", worst == 0.0, f"maxdiff={worst}")

# ---------------------------------------------------------------- 5. fast path vs full path, every kernel
FAST = ("raw", "inc1", "inc50", "herm1", "herm2", "sym", "datadiff", "hp10", "frac_oracle")
for gen, H, n in (("G0", 0.5, 1000), ("G0", 0.95, 1000), ("HT-S", 0.8, 1000), ("NL-SV", 0.8, 1000), ("SEAS-V", 0.5, 1000),
                  ("SEAS-M-shared", 0.8, 1000), ("NL-TAR", 0.8, 1000), ("G0", 0.8, 2000)):
    r = rng(5, int(H * 100), n)
    x, y = K.ind_pair(gen, "null", H, H, n, 0.0, r)
    mult = K.Multipliers(n, K.B_WILD, r)
    sx, sy = K.median_bw_fast(x), K.median_bw_fast(y)
    fx, fy = K.gauss_features(x, sx), K.gauss_features(y, sy)
    for k in FAST:
        if k in ("sym", "datadiff"):
            ux, uy = K.own_series(k, x), K.own_series(k, y)
            gx, gy = K.gauss_features(ux, K.median_bw_fast(ux)), K.gauss_features(uy, K.median_bw_fast(uy))
        else:
            gx, gy = K.transform_features(k, fx, x, H), K.transform_features(k, fy, y, H)
        kc, lc = K.gram_full(k, x, sx, H), K.gram_full(k, y, sy, H)
        w = mult.get(kc.shape[0])
        full = K.rules_full(kc, lc, w)
        fast = K.rules_fast(K.compress(gx), K.compress(gy), w)
        e = abs(fast["t"] / full["t"] - 1.0)
        er = max(abs(fast["r_wb"] / full["r_wb"] - 1.0), abs(fast["r_shift"] / full["r_shift"] - 1.0))
        ok = e < 5e-5 and er < 5e-5 and fast["p_wb"] == full["p_wb"] and fast["p_shift"] == full["p_shift"] \
            and abs(full["check"]) < 1e-8 and abs(fast["check"]) < 1e-8
        check(f"fast=full {k} {gen} H={H} n={n}", ok, f"T rel={e:.1e} R rel={er:.1e} q={fast['q']} "
              f"p_wb {fast['p_wb']:.3f}/{full['p_wb']:.3f} shift check={full['check']:.1e}")
    # dcov_inc1 shortcuts
    ka, la = mk.increment(-np.abs(x[:, None] - x[None, :]), 1), mk.increment(-np.abs(y[:, None] - y[None, :]), 1)
    w = mult.get(n - 1)
    a, b = K.rules_full(ka, la, w), K.rules_distance(ka, la, w)
    e = max(abs(a[z] / b[z] - 1.0) for z in ("t", "r_wb", "r_shift"))
    check(f"distance shortcuts dcov_inc1 {gen} H={H} n={n}", e < 1e-9 and a["p_wb"] == b["p_wb"] and a["p_shift"] == b["p_shift"]
          and abs(a["check"]) < 1e-8, f"rel={e:.1e}")
    del ka, la, kc, lc

# inc1 subsampling on the fast path against explicit windows (raw is in tests/test_memory_kernels.py)
n = 1000
x, y = K.ind_pair("G0", "null", 0.9, 0.9, n, 0.0, rng(6))
sx, sy = K.median_bw_fast(x), K.median_bw_fast(y)
dxk, dyk = mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sx), mb.gauss_gram(np.abs(y[:, None] - y[None, :]), sy)
fxg, fyg = K.compress(K.gauss_features(x, sx)[1:] - K.gauss_features(x, sx)[:-1]), K.compress(K.gauss_features(y, sy)[1:] - K.gauss_features(y, sy)[:-1])
i, j = K.select_pairs(fxg[1], fyg[1])
m = n - 1
l, dn = math.floor(m ** 0.5), math.floor(0.1 * m)
fast = K.sub_fast(fxg[0], fyg[0], i, j, l, dn)
brute = np.empty(m - l - dn)
for k0 in range(m - l - dn):
    a = mk.increment(dxk[k0:k0 + l + 1, k0:k0 + l + 1], 1)
    b = mk.increment(dyk[k0 + dn:k0 + dn + l + 1, k0 + dn:k0 + dn + l + 1], 1)
    brute[k0] = (a * b).sum() / l
e = np.abs(fast - brute).max() / brute.max()
check("inc1 subsampling windows equal explicit Gram windows, H=0.9", e < 1e-4 and fast.size == m - l - dn, f"rel={e:.1e}")
q = K.q95(fast), K.q95(brute)
check("inc1 subsampling q95 agrees", abs(q[0] / q[1] - 1.0) < 1e-4, f"{q[0]:.5f} {q[1]:.5f}")

# ---------------------------------------------------------------- 6. split MMD
for gen, alt, par, H, n in (("G0", "null", 0.0, 0.8, 500), ("G0", "shift", 0.1, 0.5, 2000)):
    r3, r7 = rng(7), rng(7)
    x3 = P3.mmd_series(gen, alt, H, n, par, r3)
    x7 = K.mmd_series(gen, "mean" if alt == "shift" else alt, H, n, par, r7)
    check(f"mmd_series {gen} {alt} equals EXP3", same(x3, x7) == 0)
    want = P3.mmd_ext(x3, r3)
    got = K.mmd_bundle(x7, r7, exact=True, hsplit=False)
    worst = max(abs(want["mmd2"] - got["xraw_mmd2"]), abs(want["p_unc_05"] - got["xraw_p_wb"]),
                abs(want["r_unc_05"] - got["xraw_r_wb"]), abs(want["sigma"] - got["xraw_sigma"]))
    check(f"mmd raw exact equals EXP3 mmd_ext {gen} {alt} n={n}", worst == 0.0, f"maxdiff={worst}")
    e = abs(got["mmd2_raw"] / got["xraw_mmd2"] - 1.0)
    check(f"mmd fast raw equals exact {gen} {alt} n={n}", e < 1e-6 and got["p_wb_raw"] == got["xraw_p_wb"], f"rel={e:.1e}")

for gen, H in (("G0", 0.8), ("HT-S", 0.8)):
    n = 1000
    r = rng(8)
    x = K.mmd_series(gen, "null", H, n, 0.0, r)
    w = K.M.multipliers("unc", n, n // 2, math.floor(n ** 0.5), r, 500)
    s = K.contrast(n, n // 2)
    a = w * s[:, None]
    sigma = K.median_bw_fast(x)
    k = mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sigma)
    f = K.gauss_features(x, sigma)
    for name, kf in (("herm1", mk.hermite(k, x, 1)), ("herm2", mk.hermite(k, x, 2))):
        u = mk.hermite_basis(x, int(name[-1]))
        st, pv, rv, _ = K.mmd_fast_rule(f - u @ (u.T @ f), s, a)
        boot = np.einsum("ij,ij->j", a, kf @ a)
        st_full = float(s @ kf @ s)
        check(f"mmd {name} fast equals w'QKQw {gen}", abs(st / st_full - 1.0) < 1e-6 and pv == K.pboot(st_full, boot),
              f"rel={abs(st / st_full - 1.0):.1e}")
    for name in ("sym", "datadiff"):
        v = K.own_series(name, x)
        kf = mb.gauss_gram(np.abs(v[:, None] - v[None, :]), K.median_bw_fast(v))
        kc = mk.double_center(kf)
        if name == "sym":
            s2, a2 = s, a
        else:
            s2 = K.contrast(n - 1, n // 2 - 1)
            a2 = np.delete(w, n // 2 - 1, axis=0) * s2[:, None]
        st, pv, rv, _ = K.mmd_fast_rule(K.gauss_features(v, K.median_bw_fast(v)), s2, a2)
        st_full = float(s2 @ kf @ s2)
        boot = np.einsum("ij,ij->j", a2, kc @ a2)
        check(f"mmd {name} fast equals Gaussian kernel on the series {gen}", abs(st / st_full - 1.0) < 1e-6
              and pv == K.pboot(st_full, boot), f"rel={abs(st / st_full - 1.0):.1e}")

# ---------------------------------------------------------------- 7. telescoping
print("Telescoping: the increment kernel in split MMD. s' D K D' s = v' K v with v on three points.")
for n, H in ((500, 0.5), (2000, 0.5), (2000, 0.9)):
    x = mb.fgn(H, n, rng(9, n))
    d = np.abs(x[:, None] - x[None, :])
    k = mb.gauss_gram(d, mb.median_bandwidth(d))
    m = n - 1
    h = m // 2
    s = K.contrast(m, h)
    v = mk.telescoped_contrast(n, h)
    st = float(s @ mk.increment_raw(k, 1) @ s)
    bound = float(np.abs(v).sum() ** 2)
    plain = float(K.contrast(n, n // 2) @ k @ K.contrast(n, n // 2))
    check(f"telescoping identity n={n} H={H}", abs(st - float(v @ k @ v)) < 1e-12 and st <= bound,
          f"inc MMD2={st:.2e} l1 bound={bound:.2e} raw MMD2={plain:.2e} n*inc MMD2={n * st:.2e}")

# ---------------------------------------------------------------- 8. growth bundle
for n in (2000, 8000):
    x, y = mb.fgn(0.9, n, rng(10, n)), mb.fgn(0.9, n, rng(11, n))
    g = K.growth_bundle(x, y)
    dx = np.abs(x[:, None] - x[None, :])
    sx = mb.median_bandwidth(dx)
    kx = mb.gauss_gram(dx, sx)
    del dx
    dy = np.abs(y[:, None] - y[None, :])
    sy = mb.median_bandwidth(dy)
    ky = mb.gauss_gram(dy, sy)
    del dy
    t_raw = float((mb.double_center(kx) * ky).sum() / n)
    ix = mk.increment(kx, 1)
    del kx
    t_inc = float((ix * mk.increment_raw(ky, 1)).sum() / (n - 1))
    del ix, ky
    e1, e2 = abs(g["t_raw"] / t_raw - 1.0), abs(g["t_inc1"] / t_inc - 1.0)
    check(f"growth_bundle equals full T n={n} H=0.9", e1 < 1e-7 and e2 < 1e-7, f"raw rel={e1:.1e} inc1 rel={e2:.1e}")

# ---------------------------------------------------------------- 9. stress generators
reps = 300
r = rng(12)
vy, cxy, cxy_null, vy_scale, vy_flat = [], [], [], [], []
for _ in range(reps):
    a, b = K.ind_pair("G0", "factor", 0.5, 0.5, 2000, 0.0, r)
    cxy.append(np.mean((a - a.mean()) * (b - b.mean())))
    vy.append(np.var(a))
    a, b = K.ind_pair("FACTOR", "null", 0.5, 0.5, 2000, 0.0, r)
    cxy_null.append(np.corrcoef(a, b)[0, 1])
    a, b = K.ind_pair("G0", "scale", 0.8, 0.8, 2000, 0.2, r)
    vy_scale.append(np.mean(b * b))
    a, b = K.ind_pair("G0", "flat", 0.8, 0.8, 2000, 0.3, r)
    vy_flat.append(np.mean(b * b))
vy, cxy = float(np.mean(vy)), float(np.mean(cxy))
check("factor: Var X near 10/7", abs(vy - 10.0 / 7.0) < 0.05, f"{vy:.3f}")
check("factor: corr(X, Y) near 0.3", abs(cxy / vy - 0.3) < 0.02, f"{cxy / vy:.3f}")
check("factor null: mean corr near 0", abs(float(np.mean(cxy_null))) < 0.03, f"{float(np.mean(cxy_null)):.3f}")
check("scale: Var Y near 1", abs(float(np.mean(vy_scale)) - 1.0) < 0.03, f"{float(np.mean(vy_scale)):.3f}")
check("flat: Var Y near 1", abs(float(np.mean(vy_flat)) - 1.0) < 0.03, f"{float(np.mean(vy_flat)):.3f}")
xi = r.standard_normal(1_000_000)
for a_ in (0.1, 0.4):
    z = math.sqrt(1 - a_ * a_) * xi + a_ * (xi * xi - 1) / math.sqrt(2)
    check(f"shape alternative a={a_}: mean 0, variance 1", abs(z.mean()) < 0.005 and abs(z.var() - 1.0) < 0.01,
          f"mean {z.mean():.4f} var {z.var():.4f}")
x, y = K.ind_pair("G0", "lag", 0.5, 0.5, 2000, 0.5, rng(13))
check("lag: Y depends on X_{t-1}", np.corrcoef(x[:-1], y[1:])[0, 1] > 0.35 and abs(np.corrcoef(x, y)[0, 1]) < 0.1)

# ---------------------------------------------------------------- 10. self-normalised statistics (hsplit)
from kmem.legacy import mmdcal as M6  # noqa: E402

v = rng(14).standard_normal(1000)
d = abs(M6.sn_stat(M6.prefix2(np.outer(v, v)), np.array([0]), 1000)[0] / 16.0 - K.sn_split_stat(v))
check("G_n(n/2)^2 equals mmdcal.sn_stat / 16 on the linear kernel", d < 1e-9, f"diff={d:.1e}")
l = 40
for name, a_, b_ in (("mean", K.sn_mean_windows(v, l), [K.sn_mean_stat(v[k:k + l]) for k in range(v.size - l + 1)]),
                     ("split", K.sn_split_windows(v, l), [K.sn_split_stat(v[k:k + l]) for k in range(v.size - l + 1)])):
    e = float(np.max(np.abs(a_ / np.array(b_) - 1.0)))
    check(f"SN {name} windows equal direct computation", e < 1e-6, f"rel={e:.1e}")
k0 = 437
num = float(v[:k0].sum()) - k0 / v.size * float(v.sum())
s1 = np.cumsum(v[:k0] - v[:k0].mean())
s2 = np.cumsum(v[k0:] - v[k0:].mean())
g_def = abs(num) / math.sqrt(((s1 * s1).sum() + (s2 * s2).sum()) / v.size)
check("sn_split_stat equals Betken (2016) G_n(k) as written", abs(math.sqrt(K.sn_split_stat(v, k0)) - g_def) < 1e-10, f"{g_def:.5f}")
# Betken (2016, arXiv:1403.0265v2) Table 1: sup_k G_n(k) on ranks, [0.15, 0.85], n = 1000
BETKEN = {0.6: (6.182835, 7.276568, 9.785915), 0.7: (6.847260, 8.190125, 11.380584), 0.8: (7.767277, 9.495194, 13.021080),
          0.9: (8.520039, 10.333602, 14.544094)}
from scipy.stats import rankdata  # noqa: E402

for H, tab in BETKEN.items():
    r = rng(15, int(H * 100))
    t = np.array([K.betken_sup(rankdata(mb.fgn(H, 1000, r))) for _ in range(6000)])
    q = np.quantile(t, (0.90, 0.95, 0.99))
    e = max(abs(q[i] / tab[i] - 1.0) for i in range(3))
    check(f"Betken Table 1 reproduced, H={H}, n=1000, ranks", e < 0.06,
          f"mine {q[0]:.2f} {q[1]:.2f} {q[2]:.2f} / hers {tab[0]:.2f} {tab[1]:.2f} {tab[2]:.2f}")
# Lobato (2001): 95% point of the zero-mean SN statistic under short memory is 45.4
r = rng(16)
t = np.array([K.sn_mean_stat(r.standard_normal(1000)) for _ in range(20000)])
q = float(np.quantile(t, 0.95))
check("zero-mean SN statistic: 95% point near Lobato's 45.4, iid", abs(q / 45.4 - 1.0) < 0.05, f"{q:.1f}")
# hsplit columns do not disturb any other column, and the hsplit-only pass matches the full bundle
if os.path.exists(K.SN_TABLE_PATH):
    r1, r2 = rng(17), rng(17)
    x, y = K.ind_pair("G0", "null", 0.8, 0.8, 500, 0.0, r1)
    K.ind_pair("G0", "null", 0.8, 0.8, 500, 0.0, r2)
    a = K.ind_bundle(x, y, r1, HX=0.8, HY=0.8)
    b = K.ind_bundle(x, y, r2, kernels=("herm1", "herm2"), HX=0.8, HY=0.8, full=False)
    same_cols = ("t_herm1", "p_wb_herm1", "p_shift_herm1", "r975_shift_herm2", "rl1", "rl2", "rb1", "rb2", "T_12")
    check("hsplit-only pass equals the full bundle", all(float(a[c]) == float(b[c]) for c in same_cols))
    r1, r2 = rng(18), rng(18)
    x = K.mmd_series("G0", "null", 0.8, 500, 0.0, r1)
    K.mmd_series("G0", "null", 0.8, 500, 0.0, r2)
    a = K.mmd_bundle(x, r1)
    b = K.mmd_bundle(x, r2, kernels=("herm1", "herm2"))
    check("hsplit-only MMD pass equals the full bundle", all(float(a[c]) == float(b[c]) for c in ("mmd2_herm2", "p_wb_herm1", "r975_wb_herm2", "rl1", "rl2", "rb2")))
    # layer memories: the second layer of a long-memory Gaussian has H = 2H - 1, and SV has its memory in the second layer
    r = rng(19)
    h1 = np.mean([K.layer_memory(K.hermite_layers(mb.fgn(0.9, 4000, r))[:, 0]) for _ in range(40)])
    h2 = np.mean([K.layer_memory(K.hermite_layers(mb.fgn(0.9, 4000, r))[:, 1]) for _ in range(40)])
    # Theory gives H = 2H - 1 = 0.8 for the second layer. Local Whittle at m = n^0.65 on normal scores reads about 0.71 (DESIGN.md A4).
    check("layer memory on fGn(0.9): first layer near 0.9, second layer between 0.6 and the first", abs(h1 - 0.9) < 0.04 and 0.6 < h2 < h1, f"{h1:.3f} {h2:.3f}")
    tab = K.sn_table()
    for kind in ("split1", "split2", "prod"):
        y1 = tab["quantiles"][f"{kind}|2000|0.975"]
        check(f"table {kind}|2000|0.975 grows from H=0.5 to H=0.95", y1[-1] > y1[0], " ".join(f"{v:.0f}" for v in y1))
else:
    print("SKIP hsplit bundle checks: sn_table.json not built yet")

# ---------------------------------------------------------------- 11. Harris (2020) rules on real Gram matrices
n, big = 800, 39
x, y = K.ind_pair("G0", "lin", 0.8, 0.8, n, 0.3, rng(20))
sx, sy = K.median_bw_fast(x), K.median_bw_fast(y)
kx, ky = mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sx), mb.gauss_gram(np.abs(y[:, None] - y[None, :]), sy)
for name, ax, ay, m in (("raw", kx, ky, n), ("inc1", mk.increment_raw(kx, 1), mk.increment_raw(ky, 1), n - 1)):
    a = mk.double_center(ax[big:m - big, big:m - big])
    brute = np.array([(a * mk.double_center(ay[s + big:s + m - big, s + big:s + m - big])).sum() / (m - 2 * big) for s in range(-big, big + 1)])
    fx = K.gauss_features(x, sx)
    fy = K.gauss_features(y, sy)
    if name == "inc1":
        fx, fy = fx[1:] - fx[:-1], fy[1:] - fy[:-1]
    fast = K.harris_shifts(K.compress(fx), K.compress(fy), big)
    e = float(np.max(np.abs(fast / brute - 1.0)))
    mm = int(np.sum(np.delete(brute, big) >= brute[big]))
    check(f"Harris V_s equal explicit segment Gram matrices, {name}", e < 1e-4 and (K.harris_ratio(fast, big, 2.0) > 1) == (mm < 2.0), f"rel={e:.1e}, m={mm}")

print("ALL PASS" if not fails else f"{len(fails)} FAIL: {fails}")
sys.exit(1 if fails else 0)
