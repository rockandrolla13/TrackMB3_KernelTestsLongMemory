"""Hermite split (EXP7 DESIGN.md A2): tests of the low Hermite layers that herm_r removes, by self-normalised statistics.

hsplit_r tests the remainder with herm_r and the removed layers with a self-normalised test, Bonferroni.
rl uses the plug-in critical values of a simulated table (SNTable), rb the H-free window quantiles.
"""
from __future__ import annotations

import math

import numpy as np

from kmem.calibration.selfnorm import SNTable, sn_mean_stat, sn_mean_windows, sn_split_stat, sn_split_windows
from kmem.layers import hermite_layers, layer_memory

HSPLIT_R = (1, 2)
IND_PAIRS = {1: ((1, 1),), 2: ((1, 1), (1, 2), (2, 1), (2, 2))}
MMD_LAYERS = {1: (1,), 2: (1, 2)}


def hsplit_ind(x: np.ndarray, y: np.ndarray, table: SNTable) -> dict:
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
        rl2 = max(rl2, t / table.cut("prod", n, 0.99375, hu))
        rb2 = max(rb2, t / float(np.quantile(win, 0.99375)))
        if (a, b) == (1, 1):
            out["rl1"] = t / table.cut("prod", n, 0.975, hu)
            out["rb1"] = t / float(np.quantile(win, 0.975))
    out["rl2"], out["rb2"] = rl2, rb2
    return out


def hsplit_mmd(x: np.ndarray, table: SNTable) -> dict:
    """SN tests of a change at n // 2 in He_1(z) and He_2(z). hsplit_r spends alpha/2 on the low layers, split over r."""
    n = x.size
    l = math.floor(n ** 0.5)
    v = hermite_layers(x)
    h1, h2 = layer_memory(v[:, 0]), layer_memory(v[:, 1])
    t1, t2 = sn_split_stat(v[:, 0]), sn_split_stat(v[:, 1])
    w1, w2 = sn_split_windows(v[:, 0], l), sn_split_windows(v[:, 1], l)
    return {"hl_1": h1, "hl_2": h2, "T_1": t1, "T_2": t2,
            "rl1": t1 / table.cut("split1", n, 0.975, h1),
            "rl2": max(t1 / table.cut("split1", n, 0.9875, h1), t2 / table.cut("split2", n, 0.9875, h2)),
            "rb1": t1 / float(np.quantile(w1, 0.975)),
            "rb2": max(t1 / float(np.quantile(w1, 0.9875)), t2 / float(np.quantile(w2, 0.9875)))}
