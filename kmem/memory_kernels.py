"""Kernels built to remove long memory: pure maps from Gram matrices to transformed Gram matrices.

Every transform is a linear map A on the features of a positive semidefinite Gram matrix K, so the
result A K A' stays positive semidefinite. Inputs are n x n arrays and 1-d series. No I/O and no state
beyond a small cache of fixed DCT bases.

Used by experiments/EXP7-memory-kernel-20261101. Exploratory; not evidence for any accepted claim.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import special, stats


def double_center(k: np.ndarray) -> np.ndarray:
    """H K H with H = I - 11'/n. Same arithmetic as mbcore.double_center."""
    rm = k.mean(axis=1, keepdims=True)
    cm = k.mean(axis=0, keepdims=True)
    return k - rm - cm + k.mean()


def distance_kernel(x: np.ndarray) -> np.ndarray:
    """K = -|x_i - x_j|. Conditionally positive definite, so D K D' and H K H are positive semidefinite."""
    return -np.abs(x[:, None] - x[None, :])


# ---------------------------------------------------------------- increments

def increment_raw(k: np.ndarray, lag: int = 1) -> np.ndarray:
    """D K D' with (D f)_t = f_{t+lag} - f_t, not centred. Output has m = n - lag points."""
    if not 1 <= lag < k.shape[0]:
        raise ValueError(f"lag must be in [1, n), got {lag}")
    return k[lag:, lag:] - k[lag:, :-lag] - k[:-lag, lag:] + k[:-lag, :-lag]


def increment(k: np.ndarray, lag: int = 1) -> np.ndarray:
    """H D K D' H: the increment kernel, double centred."""
    return double_center(increment_raw(k, lag))


def telescoped_contrast(n: int, h: int) -> np.ndarray:
    """Weights v on the n original points such that s' D = v', for the split contrast s on n - 1 increments.

    s has weight 1/h on the first h increments and -1/(m - h) on the other m - h. The sum telescopes,
    so v is zero except at t = 0, t = h and t = n - 1. Then s' D K D' s = v' K v.
    """
    m = n - 1
    if not 1 <= h < m:
        raise ValueError(f"h must be in [1, n - 1), got {h}")
    v = np.zeros(n)
    v[0] = -1.0 / h
    v[h] = 1.0 / h + 1.0 / (m - h)
    v[n - 1] = -1.0 / (m - h)
    return v


# ---------------------------------------------------------------- Hermite projection

def normal_scores(x: np.ndarray) -> np.ndarray:
    """z = Phi^{-1}((rank - 0.5) / n), ties get the average rank."""
    return special.ndtri((stats.rankdata(x) - 0.5) / x.size)


def hermite_design(z: np.ndarray, r: int) -> np.ndarray:
    """n x (r + 1) matrix with columns 1, He_1(z), ..., He_r(z) (probabilists' Hermite polynomials)."""
    if r < 0:
        raise ValueError("r >= 0")
    out = np.empty((z.size, r + 1))
    out[:, 0] = 1.0
    if r >= 1:
        out[:, 1] = z
    for j in range(1, r):
        out[:, j + 1] = z * out[:, j] - j * out[:, j - 1]
    return out


def orthonormal(design: np.ndarray) -> np.ndarray:
    """Orthonormal basis of the column space of a full-rank design."""
    q, _ = np.linalg.qr(design)
    return q


def hermite_basis(x: np.ndarray, r: int) -> np.ndarray:
    """Orthonormal basis of span{1, He_1(z), ..., He_r(z)}, where z are the normal scores of x."""
    return orthonormal(hermite_design(normal_scores(x), r))


def project_out(k: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Q K Q with Q = I - U U', U orthonormal. Costs O(n^2 p), not O(n^3)."""
    ku = k @ u
    out = k - u @ ku.T - ku @ u.T + u @ (u.T @ ku) @ u.T
    return 0.5 * (out + out.T)


def hermite(k: np.ndarray, x: np.ndarray, r: int) -> np.ndarray:
    """Q_r K Q_r, with Q_r removing 1, He_1(z), ..., He_r(z) for the normal scores z of this series x."""
    return project_out(k, hermite_basis(x, r))


# ---------------------------------------------------------------- high-pass filter

_DCT: dict[tuple[int, int], np.ndarray] = {}


def dct_basis(n: int, p: int) -> np.ndarray:
    """Orthonormal DCT-II vectors for frequency indices 0..p-1 (index 0 is the constant), n x p."""
    key = (n, p)
    if key not in _DCT:
        j = np.arange(p)[None, :]
        t = np.arange(n)[:, None]
        b = math.sqrt(2.0 / n) * np.cos(np.pi * j * (t + 0.5) / n)
        b[:, 0] = 1.0 / math.sqrt(n)
        _DCT[key] = b
    return _DCT[key]


def n_low(n: int, frac: float) -> int:
    """Number of DCT frequencies removed: floor(frac * n), at least the constant."""
    return max(1, math.floor(frac * n))


def highpass(k: np.ndarray, frac: float = 0.1) -> np.ndarray:
    """F K F, with F removing the lowest frac of DCT-II frequencies, the constant included."""
    n = k.shape[0]
    return project_out(k, dct_basis(n, n_low(n, frac)))


# ---------------------------------------------------------------- fractional difference

def fracdiff_weights(d: float, n: int) -> np.ndarray:
    """Coefficients pi_0..pi_{n-1} of (1 - B)^d: pi_0 = 1, pi_j = pi_{j-1} (j - 1 - d) / j."""
    j = np.arange(1, n)
    return np.concatenate([[1.0], np.cumprod((j - 1.0 - d) / j)])


def fracdiff_features(f: np.ndarray, d: float) -> np.ndarray:
    """(F f)_t = sum_{j=0}^{t} pi_j f_{t-j} per column. The filter is truncated at the sample start."""
    n = f.shape[0]
    pi = fracdiff_weights(d, n)
    nfft = 2 * n
    spec = np.fft.rfft(f, nfft, axis=0) * np.fft.rfft(pi, nfft)[:, None]
    return np.fft.irfft(spec, nfft, axis=0)[:n]


def fracdiff(k: np.ndarray, d: float) -> np.ndarray:
    """H F H K H F' H, with F the truncated fractional difference of order d. d = 0 returns H K H."""
    kc = double_center(k)
    if d == 0.0:
        return kc
    a = fracdiff_features(kc, d)
    return double_center(fracdiff_features(a.T, d))


# ---------------------------------------------------------------- series maps with their own Gaussian kernel

def even_map(x: np.ndarray) -> np.ndarray:
    """|x - median(x)|: the input of the sym kernel."""
    return np.abs(x - np.median(x))


def first_difference(x: np.ndarray) -> np.ndarray:
    """x_{t+1} - x_t: the input of the datadiff kernel."""
    return x[1:] - x[:-1]
