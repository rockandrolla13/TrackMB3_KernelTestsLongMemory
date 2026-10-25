"""MMD homogeneity calibrations for H5 and E2 (EXP1-20261016, claude-mb-kernel-exp1-20261016).

Exploratory simulation; not evidence for any accepted claim. numpy only.

The statistic is the one in khom.py (copied from the mmd-dml arm): biased (V-statistic) MMD^2
between the first and second half of one series, Gaussian kernel, median-heuristic bandwidth over
nonzero pairwise distances. Calibrations computed on the same Gram matrix in every replication:

wild bootstrap (B = 500, Rademacher block-Markov multipliers, corr(W_t, W_t+k) = exp(-k/l)).
  H5 compares all five forms below at l = floor(n^0.5); E2 uses the WILD_E2 list:
  old   one chain across both halves, centred within each half (the mmd-dml K-hom scheme)
  csg   published two-sample form (Chwialkowski, Sejdinovic & Gretton 2014, eq. 4): an
        independent chain for each half, each centred within its half
  unc   independent chain per half, NOT centred, applied to the pooled double-centred Gram
        matrix (the degenerate-kernel V_b1 form of the same paper)
  iid   i.i.d. Rademacher per point, centred within each half (H5 variant i)
  ustat csg multipliers, U-statistic form: zero diagonal in both observed and bootstrap
        statistics (H5 variant iv)
subsampling (sampling windows, khom.py rule: window statistic l * MMD^2 of the window's halves,
  full-sample bandwidth, compare with n * MMD^2): we store the 95% window quantile so that any
  rate correction (n/l)^e can be applied afterwards; reject iff n MMD^2 > q95 * (n/l)^e.
  sub_sqrt  l = floor(n^0.5)          (reference, n-rate)
  sub_bd    l = floor(n^((1 + D)/2 - 0.05)), D = min(1, 2 - 2 H_lw), H_lw local Whittle (Betken-Dehling rate)
self-normalised MMD with block sampling (H-free; Bai-Taqqu-Zhang 2016 / Betken-Wendler 2018 route):
  T = n MMD^2 / V,  V = n^-2 [ sum_{k<=m} ||S_{1,k} - (k/m) S_{1,m}||^2
                              + sum_{k>m} ||S_{k,n} - ((n-k+1)/(n-m)) S_{m+1,n}||^2 ]
  in the RKHS (S_{a,b} = sum of feature maps), m = floor(n/2). The ratio is free of the scale
  n^{2H-1}. Critical value: 95% quantile of the same statistic on every window of length l
  (split at floor(l/2)), l = floor(n^0.5) and floor(n^0.7).
memory estimates on the series itself: local Whittle and GPH, m = floor(n^0.65).
"""
from __future__ import annotations

import math

import numpy as np

import common as C

B_WILD = 500
ALPHA = 0.05


# ---------------------------------------------------------------- processes

def generate(process: str, h: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """G1 fGn; G3 sign(xi); G5 xi^2 - 1 (Hermite rank 2); G6 |xi| (Hermite rank 2)."""
    xi = C.fgn(h, n, rng)
    if process == "G1":
        return xi
    if process == "G3":
        return np.where(xi >= 0.0, 1.0, -1.0)
    if process == "G5":
        return xi * xi - 1.0
    if process == "G6":
        return np.abs(xi)
    raise ValueError(process)


# ---------------------------------------------------------------- memory estimators

def periodogram(x: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    n = x.size
    f = np.fft.rfft(x - x.mean())
    lam = 2.0 * np.pi * np.arange(1, m + 1) / n
    i = (np.abs(f[1:m + 1]) ** 2) / (2.0 * np.pi * n)
    return lam, i


def local_whittle_d(x: np.ndarray, m: int) -> float:
    """Robinson (1995) local Whittle estimate of d on Fourier frequencies 1..m, d in [-0.49, 1.49]."""
    lam, i = periodogram(x, m)
    if not np.all(i > 0):
        i = np.where(i > 0, i, np.min(i[i > 0]) if np.any(i > 0) else 1e-300)
    loglam = np.log(lam)
    mlog = loglam.mean()

    def r(d):
        return math.log(np.mean(np.exp(2.0 * d * loglam) * i)) - 2.0 * d * mlog

    grid = np.linspace(-0.49, 1.49, 199)
    vals = np.array([r(d) for d in grid])
    j = int(np.argmin(vals))
    a, b = grid[max(j - 1, 0)], grid[min(j + 1, grid.size - 1)]
    g = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = r(c), r(d)
    for _ in range(40):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = r(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = r(d)
    return 0.5 * (a + b)


def gph_d(x: np.ndarray, m: int) -> float:
    """Geweke-Porter-Hudak log-periodogram estimate of d on frequencies 1..m."""
    lam, i = periodogram(x, m)
    ok = i > 0
    reg = -2.0 * np.log(2.0 * np.sin(lam[ok] / 2.0))
    y = np.log(i[ok])
    reg = reg - reg.mean()
    return float(np.dot(reg, y - y.mean()) / np.dot(reg, reg))


# ---------------------------------------------------------------- multipliers

def markov_rademacher(n: int, b: int, l: float, rng: np.random.Generator) -> np.ndarray:
    """khom.py chain: keep W_{t-1} w.p. exp(-1/l), else fresh Rademacher. corr = exp(-k/l)."""
    refresh = rng.random((n, b)) >= math.exp(-1.0 / l)
    refresh[0] = True
    idx = np.where(refresh, np.arange(n)[:, None], 0)
    idx = np.maximum.accumulate(idx, axis=0)
    r = rng.choice(np.array([-1.0, 1.0]), size=(n, b))
    return np.take_along_axis(r, idx, axis=0)


def multipliers(form: str, n: int, h1: int, l: float, rng: np.random.Generator, b: int = B_WILD) -> np.ndarray:
    if form == "old":
        w = markov_rademacher(n, b, l, rng)
    elif form in ("csg", "unc"):
        w = np.vstack([markov_rademacher(h1, b, l, rng), markov_rademacher(n - h1, b, l, rng)])
    elif form == "iid":
        w = rng.choice(np.array([-1.0, 1.0]), size=(n, b))
    else:
        raise ValueError(form)
    if form != "unc":
        w[:h1] -= w[:h1].mean(axis=0)
        w[h1:] -= w[h1:].mean(axis=0)
    return w


def wild_p(k: np.ndarray, s: np.ndarray, stat: float, w: np.ndarray) -> float:
    a = w * s[:, None]
    boot = np.einsum("ij,ij->j", a, k @ a)
    return float((1.0 + np.sum(boot >= stat)) / (boot.size + 1.0))


# ---------------------------------------------------------------- prefix-sum helpers

def prefix2(k: np.ndarray) -> np.ndarray:
    n = k.shape[0]
    p = np.zeros((n + 1, n + 1))
    p[1:, 1:] = k.cumsum(0).cumsum(1)
    return p


def rect(p, a0, a1, b0, b1):
    return p[a1, b1] - p[a0, b1] - p[a1, b0] + p[a0, b0]


def window_mmd(p: np.ndarray, n: int, l: int) -> np.ndarray:
    """l * MMD^2 (first floor(l/2) vs rest) on every window of length l (khom.py rule)."""
    ks = np.arange(n - l + 1)
    m1 = l // 2
    m2 = l - m1
    aa = rect(p, ks, ks + m1, ks, ks + m1)
    bb = rect(p, ks + m1, ks + l, ks + m1, ks + l)
    ab = rect(p, ks, ks + m1, ks + m1, ks + l)
    return l * (aa / m1 ** 2 + bb / m2 ** 2 - 2.0 * ab / (m1 * m2))


def sn_stat(p: np.ndarray, starts: np.ndarray, l: int) -> np.ndarray:
    """Self-normalised MMD statistic T = l MMD^2 / V on windows [s, s + l) for each s in starts."""
    s = starts[:, None]
    m = l // 2
    r = l - m
    # numerator
    aa = rect(p, starts, starts + m, starts, starts + m)
    bb = rect(p, starts + m, starts + l, starts + m, starts + l)
    ab = rect(p, starts, starts + m, starts + m, starts + l)
    mmd = aa / m ** 2 + bb / r ** 2 - 2.0 * ab / (m * r)
    # first half: ||S_{1,k} - (k/m) S_{1,m}||^2, k = 1..m
    k = np.arange(1, m + 1)[None, :]
    f = k / m
    v1 = rect(p, s, s + k, s, s + k) - 2.0 * f * rect(p, s, s + k, s, s + m) + (f ** 2) * aa[:, None]
    # second half, from the right end: ||S_{k,l} - ((l-k+1)/r) S_{m+1,l}||^2 with j = l - k + 1 = 1..r
    j = np.arange(1, r + 1)[None, :]
    g = j / r
    e = s + l
    v2 = rect(p, e - j, e, e - j, e) - 2.0 * g * rect(p, e - j, e, s + m, e) + (g ** 2) * bb[:, None]
    v = (v1.sum(axis=1) + v2.sum(axis=1)) / l ** 2
    return l * mmd / v


# ---------------------------------------------------------------- one replication

def median_heuristic(d: np.ndarray) -> float:
    v = d[d > 0]
    return float(np.median(v)) if v.size else 1.0


def rep_h5(x: np.ndarray, rng: np.random.Generator) -> dict:
    """H5: wild-bootstrap variants at l = floor(sqrt(n)) on one series."""
    n = x.size
    h1 = n // 2
    d = np.abs(x[:, None] - x[None, :])
    k = np.exp(-(d * d) / (2.0 * median_heuristic(d) ** 2))
    del d
    s = np.empty(n)
    s[:h1] = 1.0 / h1
    s[h1:] = -1.0 / (n - h1)
    stat = float(s @ k @ s)
    l = math.floor(math.sqrt(n))
    out = {"mmd2": stat}
    for form in ("old", "csg", "iid"):
        out[f"p_{form}"] = wild_p(k, s, stat, multipliers(form, n, h1, l, rng))
    kc = k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()
    out["p_unc"] = wild_p(kc, s, stat, multipliers("unc", n, h1, l, rng))
    k0 = k.copy()
    np.fill_diagonal(k0, 0.0)
    stat_u = float(s @ k0 @ s)
    out["p_ustat"] = wild_p(k0, s, stat_u, multipliers("csg", n, h1, l, rng))
    return out


# (form, exponent) pairs; amended after H5 (events.jsonl parameters_amended): the published
# centred two-sample form (csg) was oversized like the old scheme at H = 0.5, the uncentred
# degenerate-kernel form (unc) was on size, so the length variants use unc.
WILD_E2 = (("old", 0.5), ("unc", 0.5), ("unc", 0.3), ("unc", 0.4))


def rep_e2(x: np.ndarray, rng: np.random.Generator) -> dict:
    """E2: every calibration on one series. Rejection columns use alpha = 0.05."""
    n = x.size
    h1 = n // 2
    d = np.abs(x[:, None] - x[None, :])
    sigma = median_heuristic(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    s = np.empty(n)
    s[:h1] = 1.0 / h1
    s[h1:] = -1.0 / (n - h1)
    stat = float(s @ k @ s)
    t_n = n * stat
    out: dict = {"mmd2": stat, "sigma": sigma, "t_n": t_n}

    kc = k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()
    for form, e in WILD_E2:
        l = max(1, math.floor(n ** e))
        out[f"p_{form}_{int(e * 10):02d}"] = wild_p(kc if form == "unc" else k, s, stat, multipliers(form, n, h1, l, rng))
    del kc

    m_lw = math.floor(n ** 0.65)
    d_lw = local_whittle_d(x, m_lw)
    out["h_lw"] = d_lw + 0.5
    out["h_gph"] = gph_d(x, m_lw) + 0.5

    p = prefix2(k)
    del k
    def q95(v):
        v = np.sort(v)
        return float(v[int(math.ceil((1 - ALPHA) * v.size)) - 1])

    l_sq = math.floor(math.sqrt(n))
    out["l_sqrt"] = l_sq
    out["q95_sqrt"] = q95(window_mmd(p, n, l_sq))
    dd = min(1.0, max(0.0, 2.0 - 2.0 * out["h_lw"]))
    l_bd = max(4, min(n - 1, math.floor(n ** ((1.0 + dd) / 2.0 - 0.05))))
    out["l_bd"] = l_bd
    out["q95_bd"] = q95(window_mmd(p, n, l_bd))

    out["sn_T"] = float(sn_stat(p, np.array([0]), n)[0])
    for e in (0.5, 0.7):
        l = math.floor(n ** e)
        out[f"sn_q95_{int(e * 10):02d}"] = q95(sn_stat(p, np.arange(n - l + 1), l))
    return out
