"""K-hom: MMD^2 between the first and second half of one series.

Null cells: G1-G3 x H grid, R1 x lambda in {0.01, 0.1}; n in {500, 1000, 2000, 4000}.
Power cells (design: one check per family at H = 0.7): G1, H = 0.7, second half shifted by
POWER_SHIFT, all n.

Statistic: biased (V-statistic) MMD^2 with Gaussian kernel exp(-d^2 / (2 sigma^2)), sigma = median
of |x_i - x_j| over i < j excluding zero distances (median heuristic; zeros excluded so the
binary G3 series gets a usable bandwidth).

Calibrations:
  wild bootstrap  Chwialkowski, Sejdinovic & Gretton (2014) two-sample wild bootstrap. One
                  Rademacher block-Markov multiplier chain over the whole series: W_1 Rademacher;
                  for t > 1 keep W_{t-1} with probability exp(-1/l), else draw a fresh Rademacher,
                  so corr(W_t, W_{t+k}) = exp(-k/l). l = floor(n^0.5). W centred within each half.
                  B = 500 draws. p = (1 + #{boot >= stat}) / (B + 1).
  subsampling     Betken-Dehling sampling-window method: windows of length l = floor(n^0.5),
                  k = 1..n-l+1; window statistic l * MMD^2(first floor(l/2) points, rest) with the
                  full-sample bandwidth; reject if n * MMD^2 exceeds the ceil(0.95 m)-th order
                  statistic of the window statistics. (d_n has no role in a one-series test.)
"""
import math

import numpy as np

import common as C

N_GRID = (500, 1000, 2000, 4000)
B_WILD = 500
POWER_SHIFT = 0.3
ALPHA = 0.05
COLUMNS = ("mmd2", "sigma", "p_wild", "rej_wild", "p_sub", "rej_sub", "sub_q95", "t_n")
SIZE_KEY = "n"


def cell_list():
    cells = []
    for proc in ("G1", "G2", "G3"):
        for h in C.HURST_GRID:
            for n in N_GRID:
                cells.append({"process": proc, "param": h, "n": n, "shift": 0.0})
    for lam in (0.01, 0.1):
        for n in N_GRID:
            cells.append({"process": "R1", "param": lam, "n": n, "shift": 0.0})
    for n in N_GRID:
        cells.append({"process": "G1", "param": 0.7, "n": n, "shift": POWER_SHIFT})
    for i, c in enumerate(cells):
        c["cell_idx"] = i
        c["cell_id"] = f"khom_{c['process']}_{c['param']}_n{c['n']}" + ("_power" if c["shift"] else "")
    return cells


def median_heuristic(d):
    """Median over i < j of the nonzero |x_i - x_j|, from the full symmetric distance matrix d.

    Each pair appears twice in d, which leaves the median unchanged; the zero diagonal is dropped.
    """
    v = d[d > 0]
    return float(np.median(v)) if v.size else 1.0


def markov_rademacher(n, b, l, rng):
    refresh = rng.random((n, b)) >= math.exp(-1.0 / l)
    refresh[0] = True
    idx = np.where(refresh, np.arange(n)[:, None], 0)
    idx = np.maximum.accumulate(idx, axis=0)
    r = rng.choice(np.array([-1.0, 1.0]), size=(n, b))
    return np.take_along_axis(r, idx, axis=0)


def run_rep(cell, rng):
    n = cell["n"]
    x = C.generate(cell["process"], cell["param"], n, rng)
    h1 = n // 2
    x = x.copy()
    x[h1:] += cell["shift"]
    d = np.abs(x[:, None] - x[None, :])
    sigma = median_heuristic(d)
    k = np.exp(-(d * d) / (2.0 * sigma * sigma))
    del d
    s = np.empty(n)
    s[:h1] = 1.0 / h1
    s[h1:] = -1.0 / (n - h1)
    stat = float(s @ k @ s)

    l = int(math.floor(math.sqrt(n)))
    w = markov_rademacher(n, B_WILD, l, rng)
    w[:h1] -= w[:h1].mean(axis=0)
    w[h1:] -= w[h1:].mean(axis=0)
    a = w * s[:, None]
    boot = np.einsum("ij,ij->j", a, k @ a)
    p_wild = (1.0 + np.sum(boot >= stat)) / (B_WILD + 1.0)

    p = np.zeros((n + 1, n + 1))
    p[1:, 1:] = k.cumsum(0).cumsum(1)
    ks = np.arange(n - l + 1)
    m1 = l // 2
    m2 = l - m1

    def rect(a0, a1, b0, b1):
        return p[a1, b1] - p[a0, b1] - p[a1, b0] + p[a0, b0]

    aa = rect(ks, ks + m1, ks, ks + m1)
    bb = rect(ks + m1, ks + l, ks + m1, ks + l)
    ab = rect(ks, ks + m1, ks + m1, ks + l)
    t_sub = l * (aa / m1 ** 2 + bb / m2 ** 2 - 2.0 * ab / (m1 * m2))
    t_n = n * stat
    m = len(t_sub)
    q = np.sort(t_sub)[int(math.ceil((1 - ALPHA) * m)) - 1]
    return {
        "mmd2": stat,
        "sigma": sigma,
        "p_wild": float(p_wild),
        "rej_wild": int(p_wild <= ALPHA),
        "p_sub": float(np.mean(t_sub >= t_n)),
        "rej_sub": int(t_n > q),
        "sub_q95": float(q),
        "t_n": t_n,
    }
