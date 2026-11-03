"""EXP7 core: generators, memory-kernel rule bundles, and a low-rank fast path (see DESIGN.md).

Three paths compute the same rules. unit_checks.py compares them on every kernel.

Exact path: mbcore functions on n x n Gram matrices (raw_exact, mmd_raw_exact). It reproduces EXP6
hsic_fields and EXP3 mmd_ext bit for bit on the same seed. It runs on the first replications of each
cell (columns xraw_*), so the fast path is also compared with it on the real cells.

Fast path: a Gaussian Gram matrix on one series has fast-decaying eigenvalues. A pivoted Cholesky
factor K = F F' gets the trace error below 1e-12 relative. Each memory kernel is a linear map A on
features, so A K A' = (A F)(A F)'. All three rules then reduce to sums over a few products of feature
columns.

Distance path: the distance kernel has slowly decaying eigenvalues, so dCov and dcov_inc1 keep their
n x n matrices. Two exact shortcuts remove most of their cost. The block-Markov multipliers are
piecewise constant, so w' M w needs only 2-d tail sums of M at the sign flips. Symmetric Gram matrices
repeat each circulant diagonal, so the shift statistic needs only half of them.

Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

from kmem import memory_kernels as mk
from kmem.legacy import gens
from kmem.legacy import mbcore as mb
from kmem.legacy import mmdcal as M

ALPHA = 0.05
B_WILD = 500
TOL_CHOL = 1e-12  # residual trace of the pivoted Cholesky factor, relative to the trace
TOL_PAIR = 1e-7  # share of sum_ij lambda_i lambda_j dropped when pruning feature-column pairs
INC_LONG = gens.PERIOD
HP_FRAC = 0.1

IND_KERNELS = ("raw", "inc1", "inc50", "herm1", "herm2", "sym", "datadiff", "hp10", "dcov_inc1", "frac_oracle")
SUB_KERNELS = ("raw", "inc1")
MMD_KERNELS = ("raw", "herm1", "herm2", "sym", "datadiff")

# hsplit: herm_r on the remainder plus a self-normalised test of the removed layers, Bonferroni (DESIGN.md A2)
HSPLIT_R = (1, 2)
IND_PAIRS = {1: ((1, 1),), 2: ((1, 1), (1, 2), (2, 1), (2, 2))}
MMD_LAYERS = {1: (1,), 2: (1, 2)}
SN_GRID = (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
SN_LEVELS = (0.95, 0.975, 0.9875, 0.99375)
SN_TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sn_table.json")
_SN_TABLE: dict | None = None

# Harris (2020) non-circular shift test (DESIGN.md A3)
HARRIS_KERNELS = ("raw", "inc1")
HARRIS_N_MIN = 19
HARRIS_STRICT_REPS = {500: 40, 2000: 40, 4000: 10}

# Stress alternatives
FACTOR_PHI = 0.99
FACTOR_VAR = 3.0 / 7.0  # corr(X, Y) = v / (v + 1) = 0.3 with unit-variance noise


def q95(v: np.ndarray) -> float:
    return float(np.quantile(v, 1.0 - ALPHA))


def pboot(t: float, b: np.ndarray) -> float:
    return float((1.0 + np.sum(b >= t)) / (1.0 + b.size))


# ---------------------------------------------------------------- generators

def seas_scale(n: int, rng: np.random.Generator) -> np.ndarray:
    u = int(rng.integers(gens.PERIOD))
    return 1.0 + gens.SEAS_B * gens.season(n, u)


def seas_v_slow(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """EXP4 SEAS-V-slow, same draw order as gens_exp4.generate: period 4n, phase independent per series."""
    p = 4 * n
    u = int(rng.integers(p))
    return (1.0 + gens.SEAS_B * np.sin(2.0 * np.pi * (np.arange(n) + u) / p)) * mb.fgn(H, n, rng)


def factor(n: int, rng: np.random.Generator) -> np.ndarray:
    """Stationary AR(1), coefficient 0.99, variance FACTOR_VAR, Gaussian iid innovations."""
    e = rng.standard_normal(n)
    out = np.empty(n)
    prev = e[0]
    for t in range(n):
        if t:
            prev = FACTOR_PHI * prev + math.sqrt(1.0 - FACTOR_PHI ** 2) * e[t]
        out[t] = prev
    return math.sqrt(FACTOR_VAR) * out


def scale_norm(gamma: float) -> float:
    """sd of exp(gamma |xi|) Z with xi, Z independent standard normal: E exp(2 gamma |xi|) = 2 e^{2 gamma^2} Phi(2 gamma)."""
    return math.sqrt(2.0 * math.exp(2.0 * gamma * gamma) * float(mb.norm_cdf(np.array(2.0 * gamma))))


def flat_ma_norm(H: float) -> float:
    """sd of xi_{t-1} + xi_t + xi_{t+1} for unit fGn(H)."""
    g = mb.fgn_acov(H, 3)
    return math.sqrt(3.0 + 4.0 * g[1] + 2.0 * g[2])


def ind_pair(gen: str, alt: str, HX: float, HY: float, n: int, par: float, rng: np.random.Generator):
    """Independence-arm pair. Null generators and lin/sq follow EXP2 and EXP3 draw orders."""
    if alt == "null":
        if gen == "G0":
            x = mb.fgn(HX, n, rng)
            return x, mb.fgn(HY, n, rng)
        if gen == "SEAS-V-slow":
            x = seas_v_slow(HX, n, rng)
            return x, seas_v_slow(HY, n, rng)
        if gen == "FACTOR":  # matched null of the common slow factor: two independent factors
            x = factor(n, rng) + mb.fgn(HX, n, rng)
            return x, factor(n, rng) + mb.fgn(HY, n, rng)
        if gen == "SCALE":  # matched null of the scale coupling: Y's scale comes from an independent copy
            x = mb.fgn(HX, n, rng)
            xp = mb.fgn(HX, n, rng)
            return x, np.exp(par * np.abs(xp)) * mb.fgn(HY, n, rng) / scale_norm(par)
        if HX != HY:
            raise ValueError("process arms use equal H")
        return gens.pair(gen, HX, n, rng)
    if HX != HY:
        raise ValueError("alternatives use equal H")
    H = HX
    if alt == "lin":  # EXP3
        xi = mb.fgn(H, n, rng)
        w = mb.fgn(H, n, rng)
        return xi, par * xi + math.sqrt(1.0 - par * par) * w
    if alt == "sq":  # EXP3
        xi = mb.fgn(H, n, rng)
        w = mb.fgn(H, n, rng)
        return xi, par * (xi * xi - 1.0) / math.sqrt(2.0) + math.sqrt(1.0 - par * par) * w
    if alt == "lag":  # Y_t = rho X_{t-1} + noise
        xi = mb.fgn(H, n + 1, rng)
        w = mb.fgn(H, n, rng)
        return xi[1:], par * xi[:-1] + math.sqrt(1.0 - par * par) * w
    if alt == "flat":  # Y_t = rho (X_{t-1} + X_t + X_{t+1}) / sd + noise
        xi = mb.fgn(H, n + 2, rng)
        w = mb.fgn(H, n, rng)
        ma = (xi[:-2] + xi[1:-1] + xi[2:]) / flat_ma_norm(H)
        return xi[1:-1], par * ma + math.sqrt(1.0 - par * par) * w
    if alt == "factor":  # X = m + e, Y = m + u, corr 0.3
        m = factor(n, rng)
        return m + mb.fgn(H, n, rng), m + mb.fgn(H, n, rng)
    if alt == "scale":  # Y = exp(gamma |X|) Z / sd
        xi = mb.fgn(H, n, rng)
        return xi, np.exp(par * np.abs(xi)) * mb.fgn(H, n, rng) / scale_norm(par)
    raise ValueError(alt)


def mmd_series(gen: str, alt: str, H: float, n: int, par: float, rng: np.random.Generator) -> np.ndarray:
    """Split-MMD series. mean: add delta to the second half (EXP3 'shift'). scale: multiply the second
    half by par. shape: second half becomes sqrt(1 - a^2) xi + a (xi^2 - 1)/sqrt(2), same mean and variance."""
    x = gens.generate(gen, H, n, rng)
    if alt == "null":
        return x
    x = x.copy()
    h = n // 2
    if alt == "mean":
        x[h:] += par
    elif alt == "scale":
        x[h:] *= par
    elif alt == "shape":
        xi = x[h:]
        x[h:] = math.sqrt(1.0 - par * par) * xi + par * (xi * xi - 1.0) / math.sqrt(2.0)
    else:
        raise ValueError(alt)
    return x


# ---------------------------------------------------------------- fast path: features

def median_bw_fast(x: np.ndarray) -> float:
    """Median pairwise distance by bisection on the sorted series, about 1e-13 relative to the range.

    Equals mbcore.median_bandwidth up to that tolerance when the median is positive. O(n log n).
    """
    s = np.sort(x)
    n = s.size
    pairs = n * (n - 1) // 2
    base = n * (n + 1) // 2

    def kth(k: int) -> float:
        lo, hi = 0.0, float(s[-1] - s[0])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if int(np.searchsorted(s, s + mid, side="right").sum()) - base >= k:
                hi = mid
            else:
                lo = mid
        return hi

    if pairs % 2:
        med = kth(pairs // 2 + 1)
    else:
        med = 0.5 * (kth(pairs // 2) + kth(pairs // 2 + 1))
    if med <= 0.0:
        return mb.median_bandwidth(np.abs(x[:, None] - x[None, :]))
    return med


def gauss_features(x: np.ndarray, sigma: float, tol: float = TOL_CHOL) -> np.ndarray:
    """Pivoted Cholesky factor F of the Gaussian Gram matrix: K = F F' with residual trace <= tol * n."""
    n = x.size
    cap = 64
    f = np.empty((n, cap))
    diag = np.ones(n)
    g = 1.0 / (2.0 * sigma * sigma)
    r = 0
    while r < n and diag.sum() > tol * n:
        if r == cap:
            cap *= 2
            f = np.concatenate([f, np.empty((n, cap - f.shape[1]))], axis=1)
        p = int(np.argmax(diag))
        c = np.exp(-g * (x - x[p]) ** 2)
        if r:
            c -= f[:, :r] @ f[p, :r]
        c /= math.sqrt(diag[p])
        f[:, r] = c
        diag -= c * c
        diag[diag < 0.0] = 0.0
        diag[p] = 0.0
        r += 1
    return f[:, :r]


def compress(f: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centre the columns and rotate them to principal axes. Returns (features, squared singular values)."""
    u, s, _ = np.linalg.svd(f - f.mean(axis=0, keepdims=True), full_matrices=False)
    return u * s, s * s


def select_pairs(lx: np.ndarray, ly: np.ndarray, tol: float = TOL_PAIR) -> tuple[np.ndarray, np.ndarray]:
    """Column pairs (i, j) with the largest lx_i ly_j. Dropped pairs carry at most tol of the total."""
    w = np.outer(lx, ly).ravel()
    order = np.argsort(-w)
    cum = np.cumsum(w[order])
    keep = min(int(np.searchsorted(cum, (1.0 - tol) * cum[-1])) + 1, w.size)
    return np.unravel_index(order[:keep], (lx.size, ly.size))


def window_sum(v: np.ndarray, l: int) -> np.ndarray:
    """Sums of v[k:k+l] along axis 0 for k = 0..len(v) - l."""
    c = np.vstack([np.zeros((1, v.shape[1])), np.cumsum(v, axis=0)])
    return c[l:] - c[:-l]


def sub_fast(px: np.ndarray, py: np.ndarray, i: np.ndarray, j: np.ndarray, l: int, dn: int) -> np.ndarray:
    """Sampling-window statistics T_{l,k}, k = 0..m-l-dn-1, X block at k, Y block at k + dn (mbcore rule)."""
    m = px.shape[0]
    nw = m - l - dn
    a, b = px[:m - dn], py[dn:]
    ux, iu = np.unique(i, return_inverse=True)
    uy, ju = np.unique(j, return_inverse=True)
    sa, sb = window_sum(a[:, ux], l), window_sum(b[:, uy], l)
    cab = window_sum(a[:, i] * b[:, j], l) - sa[:, iu] * sb[:, ju] / l
    return ((cab * cab).sum(axis=1) / l)[:nw]


def rules_fast(fx, fy, w: np.ndarray, sub: bool = False) -> dict:
    """Bare statistic, wild bootstrap, circular shift, optional subsampling, from compressed features.

    fx, fy = compress(...). T = (1/m) sum_ij Kc_ij Lc_ij; w is m x B and centred.
    """
    px, lx = fx
    py, ly = fy
    m = px.shape[0]
    A = math.floor(0.1 * m)
    i, j = select_pairs(lx, ly)
    p = px[:, i] * py[:, j]
    t = float((p.sum(axis=0) ** 2).sum() / m)
    a = p.T @ w
    bw = (a * a).sum(axis=0) / m
    del a, p
    ux, iu = np.unique(i, return_inverse=True)
    uy, ju = np.unique(j, return_inverse=True)
    xc = np.fft.irfft(np.conj(np.fft.rfft(px[:, ux], axis=0)[:, iu]) * np.fft.rfft(py[:, uy], axis=0)[:, ju], n=m, axis=0)
    s = (xc * xc).sum(axis=1) / m
    del xc
    shifts = s[A:m - A + 1]
    q_w, q_s = q95(bw), q95(shifts)
    out = {"t": t, "m": m, "q": i.size, "p_wb": pboot(t, bw), "r_wb": t / q_w,
           "p_shift": pboot(s[0], shifts), "r_shift": s[0] / q_s, "check": float(s[0] - t),
           "r975_wb": t / float(np.quantile(bw, 0.975)), "r975_shift": s[0] / float(np.quantile(shifts, 0.975))}
    if sub:
        sb = sub_fast(px, py, i, j, math.floor(m ** 0.5), math.floor(0.1 * m))
        q_b = q95(sb)
        out.update({"p_sub": float(np.mean(sb >= t)), "rej_sub": int(t > q_b), "r_sub": t / q_b})
    return out


def rules_full(kc: np.ndarray, lc: np.ndarray, w: np.ndarray) -> dict:
    """Full path on transformed Gram matrices. Same statistic and rules as rules_fast (no subsampling)."""
    m = kc.shape[0]
    A = math.floor(0.1 * m)
    mh = kc * lc
    t = float(mh.sum() / m)
    bw = mb.wild_boot_stats(mh, w)
    del mh
    s = mb.shift_stats(kc, lc)
    shifts = s[A:m - A + 1]
    q_w, q_s = q95(bw), q95(shifts)
    return {"t": t, "m": m, "q": -1, "p_wb": pboot(t, bw), "r_wb": t / q_w,
            "p_shift": pboot(s[0], shifts), "r_shift": s[0] / q_s, "check": float(s[0] - t)}


def wild_prefix(mh: np.ndarray, w: np.ndarray) -> np.ndarray:
    """(1/m) w_b' M w_b per column, for piecewise-constant multipliers. Equals mbcore.wild_boot_stats.

    With jumps J_a at positions f_a (the first row counts as a jump), w' M w = sum_ab J_a J_b Q(f_a, f_b),
    where Q(i, j) is the sum of M over s >= i, t >= j. The chains flip about m / (2 ln) times.
    """
    m, b = w.shape
    q = np.cumsum(np.cumsum(mh[::-1, ::-1], axis=0), axis=1)[::-1, ::-1]
    jumps = np.empty_like(w)
    jumps[0] = w[0]
    jumps[1:] = w[1:] - w[:-1]
    chain, pos = np.nonzero(jumps.T)
    cuts = np.searchsorted(chain, np.arange(b + 1))
    out = np.empty(b)
    for c in range(b):
        f = pos[cuts[c]:cuts[c + 1]]
        j = jumps[f, c]
        out[c] = j @ q[np.ix_(f, f)] @ j
    return out / m


def shift_stats_sym(kc: np.ndarray, lc: np.ndarray, chunk: int = 500) -> np.ndarray:
    """mbcore.shift_stats for symmetric matrices. Circulant diagonals d and n - d give the same correlation."""
    n = kc.shape[0]
    i = np.arange(n)
    half = (n - 1) // 2
    acc = np.zeros(n // 2 + 1, dtype=complex)

    def add(d: np.ndarray, weight: float) -> np.ndarray:
        cols = (i[None, :] + d[:, None]) % n
        kd = kc[i[None, :], cols]
        ld = lc[i[None, :], cols]
        return weight * (np.conj(np.fft.rfft(kd, axis=1)) * np.fft.rfft(ld, axis=1)).sum(axis=0)

    acc += add(np.array([0]), 1.0)
    if n % 2 == 0:
        acc += add(np.array([n // 2]), 1.0)
    for d0 in range(1, half + 1, chunk):
        acc += add(np.arange(d0, min(half + 1, d0 + chunk)), 2.0)
    return np.fft.irfft(acc, n=n) / n


def rules_distance(kc: np.ndarray, lc: np.ndarray, w: np.ndarray) -> dict:
    """rules_full with the two shortcuts. Same values up to rounding (unit_checks.py)."""
    m = kc.shape[0]
    A = math.floor(0.1 * m)
    mh = kc * lc
    t = float(mh.sum() / m)
    bw = wild_prefix(mh, w)
    del mh
    s = shift_stats_sym(kc, lc)
    shifts = s[A:m - A + 1]
    q_w, q_s = q95(bw), q95(shifts)
    return {"t": t, "m": m, "q": -1, "p_wb": pboot(t, bw), "r_wb": t / q_w,
            "p_shift": pboot(s[0], shifts), "r_shift": s[0] / q_s, "check": float(s[0] - t)}


def transform_features(kernel: str, f: np.ndarray, x: np.ndarray, H: float | None) -> np.ndarray:
    """Apply the kernel's linear feature map to the Gaussian features f of raw x."""
    if kernel == "raw":
        return f
    if kernel == "inc1":
        return f[1:] - f[:-1]
    if kernel == "inc50":
        return f[INC_LONG:] - f[:-INC_LONG]
    if kernel in ("herm1", "herm2"):
        u = mk.hermite_basis(x, int(kernel[-1]))
        return f - u @ (u.T @ f)
    if kernel == "hp10":
        u = mk.dct_basis(x.size, mk.n_low(x.size, HP_FRAC))
        return f - u @ (u.T @ f)
    if kernel == "frac_oracle":
        if H is None:
            raise ValueError("frac_oracle needs the true H")
        d = H - 0.5
        return f - f.mean(axis=0, keepdims=True) if d == 0.0 else mk.fracdiff_features(f - f.mean(axis=0, keepdims=True), d)
    raise ValueError(kernel)


def own_series(kernel: str, x: np.ndarray) -> np.ndarray:
    """Series on which sym and datadiff build their own Gaussian kernel."""
    return mk.even_map(x) if kernel == "sym" else mk.first_difference(x)


# ---------------------------------------------------------------- full Gram matrices of each kernel (checks)

def gram_full(kernel: str, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
    """Full-path transformed Gram matrix for one series. Used by unit checks and tests."""
    if kernel in ("sym", "datadiff"):
        u = own_series(kernel, x)
        d = np.abs(u[:, None] - u[None, :])
        return mk.double_center(mb.gauss_gram(d, mb.median_bandwidth(d)))
    if kernel == "dcov_inc1":
        return mk.increment(mk.distance_kernel(x), 1)
    k = mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sigma)
    if kernel == "raw":
        return mk.double_center(k)
    if kernel == "inc1":
        return mk.increment(k, 1)
    if kernel == "inc50":
        return mk.increment(k, INC_LONG)
    if kernel in ("herm1", "herm2"):
        return mk.hermite(k, x, int(kernel[-1]))
    if kernel == "hp10":
        return mk.highpass(k, HP_FRAC)
    if kernel == "frac_oracle":
        return mk.fracdiff(k, H - 0.5)
    raise ValueError(kernel)


# ---------------------------------------------------------------- multipliers shared across kernels

class Multipliers:
    """One block-Markov chain per replication. Shorter kernels use its first m rows, re-centred.

    The main chain is drawn first, with ln = floor(n^0.5), as in EXP6, so raw reproduces EXP6.
    ln = floor(m^0.5) can differ only for inc50, and then one extra chain is drawn. Draw order does not
    depend on the kernel set, so any kernel subset sees the same multipliers.
    """

    def __init__(self, n: int, B: int, rng: np.random.Generator):
        self.n = n
        self.ln = math.floor(n ** 0.5)
        self.w = mb.markov_rademacher(n, B, self.ln, rng)
        self.ln_long = math.floor((n - INC_LONG) ** 0.5)
        self.w_long = self.w if self.ln_long == self.ln else mb.markov_rademacher(n, B, self.ln_long, rng)
        if math.floor((n - 1) ** 0.5) != self.ln:
            raise ValueError(f"n = {n}: floor(sqrt(n - 1)) differs from floor(sqrt(n))")

    def get(self, m: int) -> np.ndarray:
        if m == self.n:
            return self.w
        base = self.w_long if m == self.n - INC_LONG else self.w
        v = base[:m]
        return v - v.mean(axis=0, keepdims=True)


# ---------------------------------------------------------------- Harris (2020) shift test

def harris_big_n(n: int) -> int:
    """N near 0.05 n with N + 1 a multiple of 20, so alpha (N + 1) is a whole number at alpha = 5%. At least 19."""
    return max(HARRIS_N_MIN, 20 * round(0.05 * n / 20) - 1)


def harris_shifts(fx, fy, big_n: int) -> np.ndarray:
    """V_s for s = -N..N (index s + N). V_s = (1/D) sum_ij Kc_ij Lc^{(s)}_ij on the segments X[N:m-N], Y[s+N:s+m-N].

    fx, fy = compress(...) on all m points. Centring X within its segment makes centring Y unnecessary. Feature
    pairs are pruned as in rules_fast. Linear correlation by FFT, without wrap-around.
    """
    px, lx = fx
    py, ly = fy
    m = px.shape[0]
    d = m - 2 * big_n
    i, j = select_pairs(lx, ly)
    ux, iu = np.unique(i, return_inverse=True)
    uy, ju = np.unique(j, return_inverse=True)
    a = px[big_n:m - big_n][:, ux]
    a = a - a.mean(axis=0, keepdims=True)
    size = m + d
    c = np.fft.irfft(np.conj(np.fft.rfft(a, size, axis=0)[:, iu]) * np.fft.rfft(py[:, uy], size, axis=0)[:, ju], n=size, axis=0)
    c = c[:2 * big_n + 1]
    return (c * c).sum(axis=1) / d


def harris_ratio(v: np.ndarray, big_n: int, m_limit: float) -> float:
    """R = V_0 / (the M-th largest V_s, s != 0), M = ceil(m_limit). R > 1 iff m < m_limit, m = #{s != 0 : V_s >= V_0}."""
    k = math.ceil(m_limit - 1e-9)
    rest = np.delete(v, big_n)
    return float(v[big_n] / np.partition(rest, rest.size - k)[rest.size - k])


def harris_rules(fx, fy, n: int) -> dict:
    """rhc: conservative at N near 0.05 n (primary). rh19: conservative at N = 19. rha: approximate at N near 0.05 n."""
    big = harris_big_n(n)
    vb = harris_shifts(fx, fy, big)
    v19 = vb if big == HARRIS_N_MIN else harris_shifts(fx, fy, HARRIS_N_MIN)
    return {"rhc": harris_ratio(vb, big, ALPHA * (big + 1)), "rh19": harris_ratio(v19, HARRIS_N_MIN, ALPHA * (HARRIS_N_MIN + 1)),
            "rha": harris_ratio(vb, big, ALPHA * (2 * big + 1))}


def harris_strict19(x: np.ndarray, y: np.ndarray) -> dict:
    """The literal form at N = 19 for raw and inc1: every bandwidth is a function of the segment, so V_s is a stationary sequence.

    The main path takes the bandwidths from the whole series, which is much cheaper but not literally a function of the segment.
    """
    n, big = x.size, HARRIS_N_MIN
    xs = x[big:n - big]
    fxs = gauss_features(xs, median_bw_fast(xs))
    out = {}
    vals = {"raw": [], "inc1": []}
    ax = {"raw": fxs - fxs.mean(axis=0, keepdims=True), "inc1": (fxs[1:] - fxs[:-1]) - (fxs[1:] - fxs[:-1]).mean(axis=0, keepdims=True)}
    for s in range(-big, big + 1):
        ys = y[s + big:s + n - big]
        f = gauss_features(ys, median_bw_fast(ys))
        for k, by in (("raw", f), ("inc1", f[1:] - f[:-1])):
            c = ax[k].T @ by
            vals[k].append(float((c * c).sum() / by.shape[0]))
    for k in vals:
        out[f"xh19_{k}"] = harris_ratio(np.array(vals[k]), big, ALPHA * (big + 1))
    return out


# ---------------------------------------------------------------- hsplit: tests of the removed Hermite layers

def bridge_sum(v: np.ndarray) -> float:
    """A = sum_{t=1}^{m} (S_t - (t/m) S_m)^2, the self-normaliser of Lobato (2001) and Shao (2010) times m^2."""
    m = v.size
    s = np.cumsum(v)
    b = s - (np.arange(1, m + 1) / m) * s[-1]
    return float(b @ b)


def bridge_sums_ranges(v: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """bridge_sum of v[a_i:b_i] for arrays of ranges, in O(len(v) + len(a)) by prefix sums."""
    n = v.size
    c = np.concatenate([[0.0], np.cumsum(v)])
    j = np.arange(n + 1, dtype=float)
    p0 = np.concatenate([[0.0], np.cumsum(c[1:])])
    p1 = np.concatenate([[0.0], np.cumsum(j[1:] * c[1:])])
    p2 = np.concatenate([[0.0], np.cumsum(c[1:] * c[1:])])
    l = (b - a).astype(float)
    base = c[a]
    sum_c, sum_jc, sum_c2 = p0[b] - p0[a], p1[b] - p1[a], p2[b] - p2[a]
    ss = sum_c2 - 2.0 * base * sum_c + l * base * base  # sum S_t^2
    ts = (sum_jc - a * sum_c) - base * (l * (l + 1.0) / 2.0)  # sum t S_t
    sl = c[b] - base
    return ss - 2.0 * (sl / l) * ts + (sl / l) ** 2 * (l * (l + 1.0) * (2.0 * l + 1.0) / 6.0)


def bridge_sums_windows(v: np.ndarray, l: int) -> np.ndarray:
    """bridge_sum of v[k:k+l] for k = 0..len(v) - l."""
    k = np.arange(v.size - l + 1)
    return bridge_sums_ranges(v, k, k + l)


def sn_mean_stat(u: np.ndarray) -> float:
    """Self-normalised statistic for a zero mean: T = n ubar^2 / (n^-2 sum_t (S_t - (t/n) S_n)^2)."""
    s = float(u.sum())
    return u.size * s * s / bridge_sum(u)


def sn_mean_windows(u: np.ndarray, l: int) -> np.ndarray:
    """sn_mean_stat on every window of length l."""
    c = np.concatenate([[0.0], np.cumsum(u)])
    s = c[l:] - c[:-l]
    return l * s * s / bridge_sums_windows(u, l)


def sn_split_stat(v: np.ndarray, k: int | None = None) -> float:
    """Betken (2016) G_n(k)^2 at a fixed split, with v in place of the ranks. k defaults to n // 2.

    G_n(k) = |sum_{i<=k} v_i - (k/n) sum_{i<=n} v_i| / ((1/n) sum_{t<=k} S_t(1,k)^2 + (1/n) sum_{t>k} S_t(k+1,n)^2)^(1/2).
    G_n(k)^2 is the T_n(k)^2 / V_n(k) of the CUSUM form in Shao and Zhang (2010) and Shao (2011).
    """
    n = v.size
    k = n // 2 if k is None else k
    num = float(v[:k].sum()) - (k / n) * float(v.sum())
    return n * num * num / (bridge_sum(v[:k]) + bridge_sum(v[k:]))


def betken_sup(v: np.ndarray, t1: float = 0.15, t2: float = 0.85) -> float:
    """Betken (2016) T_n(t1, t2) = sup_k G_n(k), k in [floor(n t1), floor(n t2)]. For the check against her Table 1."""
    n = v.size
    k = np.arange(math.floor(n * t1), math.floor(n * t2) + 1)
    c = np.concatenate([[0.0], np.cumsum(v)])
    num = c[k] - (k / n) * c[n]
    den = bridge_sums_ranges(v, np.zeros_like(k), k) + bridge_sums_ranges(v, k, np.full_like(k, n))
    return float(np.sqrt(np.max(n * num * num / den)))


def sn_split_windows(v: np.ndarray, l: int) -> np.ndarray:
    """sn_split_stat on every window of length l, split at l // 2."""
    h = l // 2
    n = v.size
    c = np.concatenate([[0.0], np.cumsum(v)])
    k = np.arange(n - l + 1)
    num = (c[k + h] - c[k]) - (h / l) * (c[k + l] - c[k])
    return l * num * num / (bridge_sums_windows(v, h)[k] + bridge_sums_windows(v, l - h)[k + h])


def hermite_layers(x: np.ndarray, r: int = 2) -> np.ndarray:
    """Columns He_1(z), ..., He_r(z) of the normal scores of x, each centred. n x r."""
    h = mk.hermite_design(mk.normal_scores(x), r)[:, 1:]
    return h - h.mean(axis=0, keepdims=True)


def layer_memory(v: np.ndarray) -> float:
    """Local Whittle H of one layer series, m = floor(n^0.65) (mmdcal), clipped to the table range."""
    d = M.local_whittle_d(v, math.floor(v.size ** 0.65))
    return float(min(max(d + 0.5, SN_GRID[0]), SN_GRID[-1]))


def sn_table() -> dict:
    global _SN_TABLE
    if _SN_TABLE is None:
        with open(SN_TABLE_PATH) as f:
            _SN_TABLE = json.load(f)
    return _SN_TABLE


def sn_cut(kind: str, n: int, level: float, h: float) -> float:
    """Critical value at memory h: interpolated in h on the log scale from the simulated table."""
    y = sn_table()["quantiles"][f"{kind}|{n}|{level}"]
    return float(math.exp(np.interp(h, SN_GRID, np.log(y))))


def hsplit_ind(x: np.ndarray, y: np.ndarray) -> dict:
    """SN tests of the products He_a(z_x) He_b(z_y), a, b <= 2. rl_r and rb_r exceed 1 iff the low-layer part rejects.

    hsplit_r spends alpha/2 on the low layers, split equally over r^2 products. rl uses the plug-in table, keyed on
    H_u = max(1/2, H_a + H_b - 1) with H_a, H_b the local Whittle memories of the two layers. rb is the H-free
    variant: the cut-off is a quantile of the same statistic on windows of length floor(n^0.5).
    """
    n = x.size
    l = math.floor(n ** 0.5)
    lx, ly = hermite_layers(x), hermite_layers(y)
    hx = [layer_memory(lx[:, 0]), layer_memory(lx[:, 1])]
    hy = [layer_memory(ly[:, 0]), layer_memory(ly[:, 1])]
    out = {"hl_x1": hx[0], "hl_x2": hx[1], "hl_y1": hy[0], "hl_y2": hy[1]}
    rl2, rb2 = 0.0, 0.0
    for a, b in IND_PAIRS[2]:
        u = lx[:, a - 1] * ly[:, b - 1]
        t = sn_mean_stat(u)
        win = sn_mean_windows(u, l)
        hu = max(0.5, hx[a - 1] + hy[b - 1] - 1.0)
        out[f"T_{a}{b}"] = t
        rl2 = max(rl2, t / sn_cut("prod", n, 0.99375, hu))
        rb2 = max(rb2, t / float(np.quantile(win, 0.99375)))
        if (a, b) == (1, 1):
            out["rl1"] = t / sn_cut("prod", n, 0.975, hu)
            out["rb1"] = t / float(np.quantile(win, 0.975))
    out["rl2"], out["rb2"] = rl2, rb2
    return out


def hsplit_mmd(x: np.ndarray) -> dict:
    """SN tests of a change at n // 2 in He_1(z) and He_2(z). hsplit_r spends alpha/2 on the low layers, split over r."""
    n = x.size
    l = math.floor(n ** 0.5)
    v = hermite_layers(x)
    h1, h2 = layer_memory(v[:, 0]), layer_memory(v[:, 1])
    t1, t2 = sn_split_stat(v[:, 0]), sn_split_stat(v[:, 1])
    w1, w2 = sn_split_windows(v[:, 0], l), sn_split_windows(v[:, 1], l)
    return {"hl_1": h1, "hl_2": h2, "T_1": t1, "T_2": t2,
            "rl1": t1 / sn_cut("split1", n, 0.975, h1),
            "rl2": max(t1 / sn_cut("split1", n, 0.9875, h1), t2 / sn_cut("split2", n, 0.9875, h2)),
            "rb1": t1 / float(np.quantile(w1, 0.975)),
            "rb2": max(t1 / float(np.quantile(w1, 0.9875)), t2 / float(np.quantile(w2, 0.9875)))}


# ---------------------------------------------------------------- independence bundle

def tagged(kernel: str, rules: dict) -> dict:
    return {f"{key}_{kernel}": v for key, v in rules.items()}


def kernel_columns(kernels) -> list[str]:
    out = []
    for k in kernels:
        out += [f"{key}_{k}" for key in ("t", "m", "q", "p_wb", "r_wb", "p_shift", "r_shift", "check")]
        if k in SUB_KERNELS:
            out += [f"{key}_{k}" for key in ("p_sub", "rej_sub", "r_sub")]
        if k in ("herm1", "herm2"):
            out += [f"{key}_{k}" for key in ("r975_wb", "r975_shift")]
        if k in HARRIS_KERNELS:
            out += [f"{key}_{k}" for key in ("rhc", "rh19", "rha", "xh19")]
    return out


HSPLIT_IND_COLS = ["hl_x1", "hl_x2", "hl_y1", "hl_y2", "T_11", "T_12", "T_21", "T_22", "rl1", "rl2", "rb1", "rb2"]
HSPLIT_MMD_COLS = ["hl_1", "hl_2", "T_1", "T_2", "rl1", "rl2", "rb1", "rb2"]

BASE_COLS = ["sx", "sy", "rx_rank", "ry_rank", "dcov_stat", "p_dcov_wb", "r_dcov_wb", "p_dcov_sub", "rej_dcov_sub",
             "r_dcov_sub", "pearson_r", "nw_bw", "t_nw", "p_nw", "rej_nw", "r_nw", "ewc_nu", "t_ewc", "rej_ewc",
             "r_ewc", "nw_lrv_nonpos"]
EXACT_COLS = [f"xraw_{k}" for k in ("sx", "sy", "t", "p_wb", "r_wb", "p_shift", "r_shift", "check", "p_sub", "rej_sub",
                                    "r_sub")]


def ind_columns(kernels) -> list[str]:
    return BASE_COLS + EXACT_COLS + HSPLIT_IND_COLS + kernel_columns(kernels)


def raw_exact(x, y, dx, dy, w, n: int) -> dict:
    """EXP6 hsic_fields at kappa = 1: same functions, same order, same arithmetic."""
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    sx, sy = mb.median_bandwidth(dx), mb.median_bandwidth(dy)
    r = rules_full(mb.double_center(mb.gauss_gram(dx, sx)), mb.double_center(mb.gauss_gram(dy, sy)), w)
    sub_h = mb.subsample_stats(x, y, ln, dn, "hsic", sx, sy)
    q_hb = q95(sub_h)
    r.update({"sx": sx, "sy": sy, "p_sub": float(np.mean(sub_h >= r["t"])), "rej_sub": int(r["t"] > q_hb),
              "r_sub": r["t"] / q_hb})
    return r


def ind_bundle(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, kernels=IND_KERNELS, B: int = B_WILD,
               HX: float | None = None, HY: float | None = None, exact: bool = False, full: bool = True,
               hsplit: bool = True, dist: bool = True, strict: bool = False) -> dict:
    """All independence rules for one pair on one set of multipliers.

    exact: also compute the exact raw path (xraw_* columns). full = False skips the n x n matrices: no
    dCov, dcov_inc1, Pearson or exact columns. dist = False skips dCov and dcov_inc1 only. hsplit: also test
    the removed Hermite layers. strict: also compute the literal Harris form at N = 19. No random draw follows
    the multipliers, so these switches and kernel subsets never change any other column.
    """
    n = x.size
    ln = math.floor(n ** 0.5)
    dn = math.floor(0.1 * n)
    out: dict = {}
    mult = Multipliers(n, B, rng)
    w = mult.w
    sx, sy = median_bw_fast(x), median_bw_fast(y)
    if full:
        if dist or exact:
            dx = np.abs(x[:, None] - x[None, :])
            dy = np.abs(y[:, None] - y[None, :])
        if exact:
            out.update({f"xraw_{k}": v for k, v in raw_exact(x, y, dx, dy, w, n).items() if k not in ("m", "q", "r975_wb", "r975_shift")})
        if dist:
            ac = mb.double_center(dx)
            bc = mb.double_center(dy)
            md = ac * bc
            del ac, bc
            t_d = float(md.sum() / n)
            bd = wild_prefix(md, w)
            del md
            sub_d = mb.subsample_stats(x, y, ln, dn, "dcov")
            q_dw, q_db = q95(bd), q95(sub_d)
            out.update({"dcov_stat": t_d, "p_dcov_wb": pboot(t_d, bd), "r_dcov_wb": t_d / q_dw,
                        "p_dcov_sub": float(np.mean(sub_d >= t_d)), "rej_dcov_sub": int(t_d > q_db),
                        "r_dcov_sub": t_d / q_db})
            if "dcov_inc1" in kernels:
                ka = mk.increment(-dx, 1)
                la = mk.increment(-dy, 1)
                out.update(tagged("dcov_inc1", rules_distance(ka, la, mult.get(n - 1))))
                del ka, la
        if dist or exact:
            del dx, dy
        lt = mb.lind_tests(x, y)
        nu = lt["ewc_nu"]
        out.update(lt)
        out["r_nw"] = abs(lt["t_nw"]) / 1.959963984540054 if np.isfinite(lt["t_nw"]) else float("nan")
        out["r_ewc"] = abs(lt["t_ewc"]) / mb._TQ[nu]
    out["sx"], out["sy"] = sx, sy
    fx = gauss_features(x, sx)
    fy = gauss_features(y, sy)
    out["rx_rank"], out["ry_rank"] = fx.shape[1], fy.shape[1]
    for k in kernels:
        if k == "dcov_inc1":
            continue
        if k in ("sym", "datadiff"):
            ux, uy = own_series(k, x), own_series(k, y)
            gx = gauss_features(ux, median_bw_fast(ux))
            gy = gauss_features(uy, median_bw_fast(uy))
        else:
            gx = transform_features(k, fx, x, HX)
            gy = transform_features(k, fy, y, HY)
        cx, cy = compress(gx), compress(gy)
        out.update(tagged(k, rules_fast(cx, cy, mult.get(gx.shape[0]), sub=k in SUB_KERNELS)))
        if k in HARRIS_KERNELS:
            out.update(tagged(k, harris_rules(cx, cy, n)))
    if hsplit:
        out.update(hsplit_ind(x, y))
    if strict:
        out.update(harris_strict19(x, y))
    return out


# ---------------------------------------------------------------- split-MMD bundle

MMD_EXACT_COLS = ["xraw_sigma", "xraw_mmd2", "xraw_p_wb", "xraw_r_wb"]


def mmd_columns(kernels) -> list[str]:
    out = ["sigma", "rank"] + MMD_EXACT_COLS + HSPLIT_MMD_COLS
    for k in kernels:
        out += [f"{key}_{k}" for key in ("mmd2", "p_wb", "r_wb")]
        if k in ("herm1", "herm2"):
            out.append(f"r975_wb_{k}")
    return out


def contrast(m: int, h: int) -> np.ndarray:
    s = np.empty(m)
    s[:h] = 1.0 / h
    s[h:] = -1.0 / (m - h)
    return s


def mmd_raw_exact(x: np.ndarray, s: np.ndarray, a: np.ndarray) -> dict:
    """EXP3 mmd_ext for the uncentred wild rule: same arithmetic. a = w * s."""
    d = np.abs(x[:, None] - x[None, :])
    sigma = M.median_heuristic(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    stat = float(s @ k @ s)
    kc = k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()
    del k
    boot = np.einsum("ij,ij->j", a, kc @ a)
    return {"xraw_sigma": sigma, "xraw_mmd2": stat, "xraw_p_wb": pboot(stat, boot), "xraw_r_wb": stat / q95(boot)}


def mmd_fast_rule(g: np.ndarray, s: np.ndarray, a: np.ndarray) -> tuple[float, float, float, float]:
    """MMD^2_V = s' K s and its uncentred wild draws a' Kc a, from features g of K. Centring does not change s' K s."""
    psi = g - g.mean(axis=0, keepdims=True)
    v = psi.T @ s
    stat = float(v @ v)
    b = psi.T @ a
    boot = (b * b).sum(axis=0)
    return stat, pboot(stat, boot), stat / q95(boot), stat / float(np.quantile(boot, 0.975))


def mmd_bundle(x: np.ndarray, rng: np.random.Generator, kernels=MMD_KERNELS, B: int = M.B_WILD,
               exact: bool = False, hsplit: bool = True) -> dict:
    """Split MMD^2_V with the tested uncentred wild rule (mmdcal form 'unc', l = floor(n^0.5)).

    All kernels use the same two chains. datadiff has one point fewer: it drops the last row of the first
    chain and splits at h - 1. herm_r uses w' Q K Q w with the contrast w.
    """
    n = x.size
    h1 = n // 2
    s = contrast(n, h1)
    l = max(1, math.floor(n ** 0.5))
    if math.floor((n - 1) ** 0.5) != l:
        raise ValueError(f"n = {n}: floor(sqrt(n - 1)) differs from floor(sqrt(n))")
    w = M.multipliers("unc", n, h1, l, rng, B)
    a = w * s[:, None]
    out: dict = mmd_raw_exact(x, s, a) if exact else {}
    sigma = median_bw_fast(x)
    f = gauss_features(x, sigma)
    out["sigma"], out["rank"] = sigma, f.shape[1]
    for kn in kernels:
        if kn == "raw":
            st, pv, rv, r975 = mmd_fast_rule(f, s, a)
        elif kn in ("herm1", "herm2"):
            u = mk.hermite_basis(x, int(kn[-1]))
            st, pv, rv, r975 = mmd_fast_rule(f - u @ (u.T @ f), s, a)
            out[f"r975_wb_{kn}"] = r975
        elif kn == "sym":
            v = mk.even_map(x)
            st, pv, rv, r975 = mmd_fast_rule(gauss_features(v, median_bw_fast(v)), s, a)
        elif kn == "datadiff":
            v = mk.first_difference(x)
            sd = contrast(n - 1, h1 - 1)
            st, pv, rv, r975 = mmd_fast_rule(gauss_features(v, median_bw_fast(v)), sd, np.delete(w, h1 - 1, axis=0) * sd[:, None])
        else:
            raise ValueError(kn)
        out.update({f"mmd2_{kn}": st, f"p_wb_{kn}": pv, f"r_wb_{kn}": rv})
    if hsplit:
        out.update(hsplit_mmd(x))
    return out


# ---------------------------------------------------------------- bare-statistic growth

GROWTH_COLS = ["t_raw", "t_inc1", "sx", "sy", "rx_rank", "ry_rank"]


def growth_bundle(x: np.ndarray, y: np.ndarray) -> dict:
    """T = (1/m) sum_ij Kc_ij Lc_ij for raw and inc1 only. Fast path, no pruning of feature pairs."""
    sx, sy = median_bw_fast(x), median_bw_fast(y)
    fx, fy = gauss_features(x, sx), gauss_features(y, sy)
    out = {"sx": sx, "sy": sy, "rx_rank": fx.shape[1], "ry_rank": fy.shape[1]}
    for name, gx, gy in (("raw", fx, fy), ("inc1", fx[1:] - fx[:-1], fy[1:] - fy[:-1])):
        px = gx - gx.mean(axis=0, keepdims=True)
        py = gy - gy.mean(axis=0, keepdims=True)
        c = px.T @ py
        out[f"t_{name}"] = float((c * c).sum() / px.shape[0])
    return out
