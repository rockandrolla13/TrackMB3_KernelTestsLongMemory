"""EXP2 generators: process features switchable independently of memory (see DESIGN.md).

Every generator is driven by unit-variance fGn from mbcore.fgn (the E3 Davies-Harte code).
Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import special, stats

from . import mbcore as mb

T_DF = 3
SV_PHI, SV_TAU = 0.9, 1.0
TAR_LO, TAR_HI = 0.8, -0.4
BURN = 200
PERIOD = 50
SEAS_A, SEAS_B = 1.0, 0.5

MAIN = ("G0", "HT-S", "HT-L", "NL-SV", "NL-TAR", "SEAS-M", "SEAS-V")
SUBORDINATED = {
    "G0": "yes (rank 1)",
    "HT-S": "yes (instantaneous, rank 1)",
    "HT-L": "no (linear filter, t3 innovations)",
    "NL-SV": "multivariate only (function of two Gaussian processes)",
    "NL-TAR": "no (nonlinear path functional)",
    "SEAS-M": "no (deterministic periodic mean, random phase)",
    "SEAS-V": "no (deterministic periodic scale, random phase)",
    "SEAS-M-shared": "no (fixed common phase; diagnostic)",
}


def t3_from_gauss(xi: np.ndarray) -> np.ndarray:
    """F_t3^{-1}(Phi(xi)) / sqrt(3), computed on the upper tail for accuracy, made odd."""
    a = np.abs(xi)
    q = stats.t.isf(special.ndtr(-a), T_DF)
    return np.sign(xi) * q / math.sqrt(T_DF / (T_DF - 2.0))


_PSI: dict[tuple[float, int], np.ndarray] = {}


def arfima_psi(d: float, m: int) -> np.ndarray:
    key = (d, m)
    if key not in _PSI:
        j = np.arange(1, m + 1)
        psi = np.concatenate([[1.0], np.cumprod((j - 1 + d) / j)])
        _PSI[key] = psi / math.sqrt(float(np.dot(psi, psi)))
    return _PSI[key]


def arfima_t(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """ARFIMA(0, H - 1/2, 0) with t3 innovations, MA truncated at M = 20n, unit variance."""
    d = H - 0.5
    eps = rng.standard_t(T_DF, size=n + 20 * n if d != 0.0 else n) / math.sqrt(T_DF / (T_DF - 2.0))
    if d == 0.0:
        return eps
    psi = arfima_psi(d, 20 * n)
    L = eps.size + psi.size - 1
    nf = 1 << (L - 1).bit_length()
    y = np.fft.irfft(np.fft.rfft(eps, nf) * np.fft.rfft(psi, nf), nf)
    return y[eps.size - n:eps.size]


_SV_SD: dict[float, float] = {}


def sv_sd(H: float) -> float:
    """Exact sd of h_t = phi h_{t-1} + sqrt(1 - phi^2) xi_t, xi unit fGn(H) (400-lag truncation)."""
    if H not in _SV_SD:
        k = 400
        g = mb.fgn_acov(H, k)
        w = SV_PHI ** np.arange(k)
        i = np.arange(k)
        gam = g[np.abs(i[:, None] - i[None, :])]
        _SV_SD[H] = math.sqrt((1 - SV_PHI ** 2) * float(w @ gam @ w))
    return _SV_SD[H]


def ar1_filter(e: np.ndarray, phi) -> np.ndarray:
    out = np.empty_like(e)
    prev = 0.0
    for t in range(e.size):
        prev = phi * prev + e[t]
        out[t] = prev
    return out


def sv(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    xi = mb.fgn(H, n + BURN, rng)
    h = ar1_filter(math.sqrt(1 - SV_PHI ** 2) * xi, SV_PHI)[BURN:] / sv_sd(H)
    return np.exp(0.5 * SV_TAU * h) * rng.standard_normal(n)


def tar(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    xi = mb.fgn(H, n + BURN, rng)
    out = np.empty(n + BURN)
    prev = 0.0
    for t in range(n + BURN):
        prev = (TAR_LO if prev <= 0.0 else TAR_HI) * prev + xi[t]
        out[t] = prev
    return out[BURN:]


def season(n: int, phase: int) -> np.ndarray:
    return np.sin(2.0 * np.pi * (np.arange(n) + phase) / PERIOD)


def generate(gen: str, H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """One series of length n from generator `gen` at memory H."""
    if gen == "G0":
        return mb.fgn(H, n, rng)
    if gen == "HT-S":
        return t3_from_gauss(mb.fgn(H, n, rng))
    if gen == "HT-L":
        return arfima_t(H, n, rng)
    if gen == "NL-SV":
        return sv(H, n, rng)
    if gen == "NL-TAR":
        return tar(H, n, rng)
    if gen == "SEAS-M":
        u = int(rng.integers(PERIOD))
        return mb.fgn(H, n, rng) + SEAS_A * season(n, u)
    if gen == "SEAS-V":
        u = int(rng.integers(PERIOD))
        return (1.0 + SEAS_B * season(n, u)) * mb.fgn(H, n, rng)
    if gen == "SEAS-M-shared":
        return mb.fgn(H, n, rng) + SEAS_A * season(n, 0)
    raise ValueError(gen)


def pair(gen: str, H: float, n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Two independent series (sequential draws from one rng; independent phases for SEAS-*)."""
    return generate(gen, H, n, rng), generate(gen, H, n, rng)
