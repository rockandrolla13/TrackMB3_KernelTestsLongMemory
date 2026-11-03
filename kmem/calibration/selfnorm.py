"""Self-normalised statistics (Lobato 2001, Shao 2010, Betken 2016) and their simulated critical values."""
from __future__ import annotations

import json
import math

import numpy as np

SN_GRID = (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
SN_LEVELS = (0.95, 0.975, 0.9875, 0.99375)


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


class SNTable:
    """Simulated critical values of the self-normalised statistics (EXP7 sn_tables.py), read lazily from JSON."""

    def __init__(self, path: str):
        self.path = path
        self._data: dict | None = None

    @property
    def data(self) -> dict:
        if self._data is None:
            with open(self.path) as f:
                self._data = json.load(f)
        return self._data

    def cut(self, kind: str, n: int, level: float, h: float) -> float:
        """Critical value at memory h: interpolated in h on the log scale from the simulated table."""
        y = self.data["quantiles"][f"{kind}|{n}|{level}"]
        return float(math.exp(np.interp(h, SN_GRID, np.log(y))))
