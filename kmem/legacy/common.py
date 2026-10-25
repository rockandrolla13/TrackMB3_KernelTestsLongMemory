"""Shared pieces for the K-hom (MMD) and DML simulations.

Job JOB-NEW-MEMORY-BOUNDARY-SIM-20261015, worker claude-mb-mmddml-20261015.
numpy only. Exploratory simulation; not evidence for any accepted claim.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import math

import numpy as np

MASTER_SEED = 20261015
PILOT_SEED = 2026101599  # pilot only; pilot output is not evidence
HURST_GRID = (0.5, 0.6, 0.7, 0.75, 0.8, 0.9)
Z975 = 1.959963984540054


# ---------------------------------------------------------------- fGn (Davies-Harte)

def fgn_autocov(h, n):
    k = np.arange(n, dtype=float)
    return 0.5 * (np.abs(k + 1) ** (2 * h) - 2 * np.abs(k) ** (2 * h) + np.abs(k - 1) ** (2 * h))


_DH_CACHE = {}


def dh_sqrt_eigs(h, n):
    """sqrt(lambda_j / M) for the circulant embedding of size M = 2n (Davies & Harte 1987).

    First row c = (g0, g1, ..., g_{n-1}, g_n, g_{n-1}, ..., g1). Raises if any eigenvalue is
    negative beyond rounding: then the method is not exact and must not be used.
    """
    key = (float(h), int(n))
    if key not in _DH_CACHE:
        g = fgn_autocov(h, n + 1)
        c = np.concatenate([g[: n + 1], g[n - 1:0:-1]])
        lam = np.fft.fft(c).real
        if lam.min() < -1e-8 * lam.max():
            raise ValueError(f"Davies-Harte embedding not nonnegative for H={h}, n={n}: {lam.min()}")
        lam = np.clip(lam, 0.0, None)
        _DH_CACHE[key] = np.sqrt(lam / (2 * n))
    return _DH_CACHE[key]


def fgn(h, n, rng):
    """Exact unit-variance fractional Gaussian noise of length n, Hurst h (Davies-Harte)."""
    if h == 0.5:
        return rng.standard_normal(n)
    m = 2 * n
    s = dh_sqrt_eigs(h, n)
    z = np.empty(m, dtype=complex)
    z[0] = rng.standard_normal()
    z[n] = rng.standard_normal()
    u = rng.standard_normal(n - 1)
    v = rng.standard_normal(n - 1)
    z[1:n] = (u + 1j * v) / math.sqrt(2.0)
    z[n + 1:] = np.conj(z[n - 1:0:-1])
    x = np.fft.fft(s * z)
    return x.real[:n]


# ---------------------------------------------------------------- processes

_erf = np.frompyfunc(math.erf, 1, 1)


def norm_cdf(x):
    return 0.5 * (1.0 + _erf(np.asarray(x) / math.sqrt(2.0)).astype(float))


def r1_series(lam, n, rng):
    """R1: log sigma_t as Euler fOU, H = 0.1: Y_t = (1 - lam) Y_{t-1} + g_t, burn-in ceil(10/lam)."""
    burn = int(math.ceil(10.0 / lam))
    g = fgn(0.1, n + burn, rng)
    y = np.empty(n + burn)
    prev = 0.0
    a = 1.0 - lam
    for t in range(n + burn):
        prev = a * prev + g[t]
        y[t] = prev
    return y[burn:]


def generate(process, param, n, rng):
    """process in {G1, G2, G3, R1}; param is H for G*, lambda for R1."""
    if process == "R1":
        return r1_series(param, n, rng)
    xi = fgn(param, n, rng)
    if process == "G1":
        return xi
    if process == "G2":
        return 2.0 * norm_cdf(xi) - 1.0
    if process == "G3":
        return np.where(xi >= 0.0, 1.0, -1.0)
    raise ValueError(process)


# ---------------------------------------------------------------- long-run variance

def andrews_bandwidth(u):
    """Andrews (1991) AR(1) plug-in bandwidth for the Bartlett kernel."""
    t = len(u)
    rho = float(np.dot(u[1:], u[:-1]) / np.dot(u[:-1], u[:-1]))
    rho = min(max(rho, -0.97), 0.97)
    alpha = 4.0 * rho ** 2 / ((1.0 - rho) ** 2 * (1.0 + rho) ** 2)
    return 1.1447 * (alpha * t) ** (1.0 / 3.0)


def nw_lrv(u):
    """Newey-West (Bartlett) long-run variance of demeaned u, Andrews bandwidth.

    Weights 1 - j/(S+1) for j = 1..floor(S). Returns (lrv, S).
    """
    t = len(u)
    s = andrews_bandwidth(u)
    L = min(int(math.floor(s)), t - 1)
    om = float(np.dot(u, u)) / t
    for j in range(1, L + 1):
        om += 2.0 * (1.0 - j / (s + 1.0)) * float(np.dot(u[j:], u[:-j])) / t
    return om, s


_EWC_BASIS = {}


def ewc_nu(t):
    return int(math.floor(0.4 * t ** (2.0 / 3.0)))


def ewc_lrv(u):
    """Equal-weighted cosine LRV (Lazarus, Lewis, Stock & Watson 2018), nu = floor(0.4 T^{2/3})."""
    t = len(u)
    nu = ewc_nu(t)
    if t not in _EWC_BASIS:
        j = np.arange(1, nu + 1)[:, None]
        tt = np.arange(1, t + 1)[None, :]
        _EWC_BASIS[t] = math.sqrt(2.0 / t) * np.cos(np.pi * j * (tt - 0.5) / t)
    lam = _EWC_BASIS[t] @ u
    return float(np.dot(lam, lam)) / nu, nu


# ---------------------------------------------------------------- Student-t quantile (no scipy)

def _t_pdf(x, nu):
    c = math.lgamma((nu + 1) / 2) - math.lgamma(nu / 2) - 0.5 * math.log(nu * math.pi)
    return np.exp(c - (nu + 1) / 2 * np.log1p(x * x / nu))


def t_cdf(x, nu, m=20001):
    """P(T_nu <= x) for x >= 0 by composite Simpson on [0, x]."""
    g = np.linspace(0.0, x, m)
    f = _t_pdf(g, nu)
    h = x / (m - 1)
    return 0.5 + h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())


_TQ = {}


def t_quantile_975(nu):
    if nu not in _TQ:
        lo, hi = 1.9, 15.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if t_cdf(mid, nu) < 0.975:
                lo = mid
            else:
                hi = mid
        _TQ[nu] = 0.5 * (lo + hi)
    return _TQ[nu]


# ---------------------------------------------------------------- Monte Carlo intervals

def wilson(k, n, z=Z975):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    w = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - w, c + w)
