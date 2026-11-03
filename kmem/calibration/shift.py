"""Circular shift test: all circular-shift statistics of symmetric centred Gram matrices by FFT."""
from __future__ import annotations

import numpy as np


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
