"""Low-rank features of Gaussian Gram matrices (the EXP7 fast path).

A Gaussian Gram matrix of one series has fast-decaying eigenvalues. A pivoted Cholesky factor K = F F' gets the
residual trace below TOL_CHOL relative to the trace. Rules then reduce to sums over products of feature columns;
select_pairs keeps the column pairs that carry all but TOL_PAIR of the total weight.
"""
from __future__ import annotations

import math

import numpy as np

TOL_CHOL = 1e-12  # residual trace of the pivoted Cholesky factor, relative to the trace
TOL_PAIR = 1e-7  # share of sum_ij lambda_i lambda_j dropped when pruning feature-column pairs


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
