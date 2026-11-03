"""EXP7 layer-memory checks (DESIGN.md A6). Separate from the main run, its own seeds, reads no main-run file except in --part apply.

--part bias   Local Whittle (m = floor(n^0.65)) on He_1(z) and He_2(z) of fGn(H), z the normal scores, against theory:
              layer 1 has H, layer 2 has max(1/2, 2H - 1). Also He_2 of the Gaussian input itself, to separate the effect of
              the normal scores from the estimator. H in {0.6, ..., 0.95}, n in {500, ..., 16000}, 400 replications. -> layer_bias.csv
--part fresh  G0 nulls at H in {0.8, 0.85, 0.9}, n = 2000, 1,000 replications, independence and split MMD: the low-layer statistics and
              the remainder ratios. -> layer_true_rows/*.csv
--part apply  The fixed comparison on fresh rows and on main-run G0 rows: plug-in cut-off from the estimated layer memory against the
              cut-off from the true layer memory. -> layer_true_vs_est.csv
Seeds: SeedSequence(entropy=2026110176, spawn_key=(part, H index, n or cell, chunk)).
Usage: nice -n 19 python -I layer_memory_check.py --part bias --workers 4
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import numpy as np  # noqa: E402

import exp7_core as K  # noqa: E402
from kmem.legacy import mbcore as mb  # noqa: E402
from kmem.legacy import mmdcal as M  # noqa: E402

MASTER = 2026110176
BIAS_H = (0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
BIAS_N = (500, 1000, 2000, 4000, 8000, 16000)
TRUE_H = (0.8, 0.85, 0.9)
ALPHA = 0.05


def lw(v: np.ndarray) -> float:
    return M.local_whittle_d(v, math.floor(v.size ** 0.65)) + 0.5


def layer_true(h: float, a: int) -> float:
    return h if a == 1 else max(0.5, 2.0 * h - 1.0)


# ---------------------------------------------------------------- part bias

def bias_chunk(args):
    ih, n, c, reps = args
    h = BIAS_H[ih]
    rng = np.random.default_rng(np.random.SeedSequence(entropy=MASTER, spawn_key=(0, ih, n, c)))
    out = np.empty((reps, 3))
    for i in range(reps):
        xi = mb.fgn(h, n, rng)
        v = K.hermite_layers(xi)
        g2 = xi * xi - 1.0
        out[i] = (lw(v[:, 0]), lw(v[:, 1]), lw(g2 - g2.mean()))
    return ih, n, out


def part_bias(workers: int, reps: int) -> None:
    tasks = [(ih, n, c, 50) for ih in range(len(BIAS_H)) for n in BIAS_N for c in range(reps // 50)]
    store: dict = {}
    with mp.get_context("fork").Pool(workers) as pool:
        for ih, n, out in pool.imap_unordered(bias_chunk, tasks):
            store.setdefault((ih, n), []).append(out)
    rows = []
    for ih, h in enumerate(BIAS_H):
        for n in BIAS_N:
            v = np.vstack(store[(ih, n)])
            row = {"H": h, "n": n, "reps": v.shape[0], "theory_layer1": h, "theory_layer2": layer_true(h, 2)}
            for j, name, th in ((0, "layer1", h), (1, "layer2", layer_true(h, 2)), (2, "layer2_gauss_input", layer_true(h, 2))):
                e = v[:, j] - th
                se = e.std(ddof=1) / math.sqrt(e.size)
                row.update({f"{name}_mean": float(v[:, j].mean()), f"{name}_bias": float(e.mean()), f"{name}_bias_lo": float(e.mean() - 1.96 * se),
                            f"{name}_bias_hi": float(e.mean() + 1.96 * se), f"{name}_sd": float(v[:, j].std(ddof=1))})
            rows.append(row)
    with open(os.path.join(HERE, "layer_bias.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote layer_bias.csv", flush=True)


# ---------------------------------------------------------------- part fresh

IND_COLS = ["rep", "hl_x1", "hl_x2", "hl_y1", "hl_y2", "T_11", "T_12", "T_21", "T_22", "rl1", "rl2",
            "r975_shift_herm1", "r975_wb_herm1", "r975_shift_herm2", "r975_wb_herm2", "p_shift_herm1", "p_wb_herm1", "p_shift_herm2", "p_wb_herm2"]
MMD_COLS = ["rep", "hl_1", "hl_2", "T_1", "T_2", "rl1", "rl2", "r975_wb_herm1", "r975_wb_herm2", "p_wb_herm1", "p_wb_herm2"]


def fresh_chunk(args):
    fam, ih, c, reps = args
    h = TRUE_H[ih]
    rows = []
    for r in range(c * reps, (c + 1) * reps):
        rng = np.random.default_rng(np.random.SeedSequence(entropy=MASTER, spawn_key=(1 if fam == "ind" else 2, ih, r)))
        if fam == "ind":
            x, y = K.ind_pair("G0", "null", h, h, 2000, 0.0, rng)
            out = K.ind_bundle(x, y, rng, kernels=("herm1", "herm2"), full=False, hsplit=True)
        else:
            x = K.mmd_series("G0", "null", h, 2000, 0.0, rng)
            out = K.mmd_bundle(x, rng, kernels=("herm1", "herm2"), hsplit=True)
        rows.append({"rep": r, **out})
    return fam, ih, rows


def part_fresh(workers: int, reps: int) -> None:
    os.makedirs(os.path.join(HERE, "layer_true_rows"), exist_ok=True)
    tasks = [(fam, ih, c, 25) for fam in ("ind", "mmd") for ih in range(len(TRUE_H)) for c in range(reps // 25)]
    store: dict = {}
    with mp.get_context("fork").Pool(workers) as pool:
        for fam, ih, rows in pool.imap_unordered(fresh_chunk, tasks):
            store.setdefault((fam, ih), []).extend(rows)
    for (fam, ih), rows in store.items():
        rows.sort(key=lambda r: r["rep"])
        with open(os.path.join(HERE, "layer_true_rows", f"{fam}__H{TRUE_H[ih]}__n2000.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=IND_COLS if fam == "ind" else MMD_COLS, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    print("wrote layer_true_rows/", flush=True)


# ---------------------------------------------------------------- part apply (the fixed comparison)

def wilson(k: int, m: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = k / m
    d = 1 + z * z / m
    c = (p + z * z / (2 * m)) / d
    hw = z * math.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / d
    return c - hw, c + hw


def mcnemar(a: np.ndarray, b: np.ndarray) -> float:
    n10, n01 = int(np.sum(a & ~b)), int(np.sum(~a & b))
    m = n10 + n01
    if m == 0:
        return 1.0
    return min(1.0, 2 * sum(math.comb(m, i) for i in range(min(n10, n01) + 1)) / 2 ** m)


def load_cols(path: str) -> dict:
    rows = [r for r in csv.DictReader(open(path)) if r.get("status", "ok") == "ok"]
    return {k: np.array([float(r[k]) if r[k] != "" else np.nan for r in rows]) for k in rows[0] if k not in ("status", "error")}


def true_cuts_ind(d: dict, h: float) -> dict:
    """Low-layer ratios with the true layer memory in place of the estimate. Same table, same levels, same Bonferroni split."""
    ht = {1: layer_true(h, 1), 2: layer_true(h, 2)}
    rl1 = d["T_11"] / K.sn_cut("prod", 2000, 0.975, max(0.5, ht[1] + ht[1] - 1.0))
    rl2 = np.max([d[f"T_{a}{b}"] / K.sn_cut("prod", 2000, 0.99375, max(0.5, ht[a] + ht[b] - 1.0)) for a, b in K.IND_PAIRS[2]], axis=0)
    est_c11 = np.array([K.sn_cut("prod", 2000, 0.975, max(0.5, a + b - 1.0)) for a, b in zip(d["hl_x1"], d["hl_y1"])])
    return {"rl1_true": rl1, "rl2_true": rl2, "c11_est": est_c11, "c11_true": K.sn_cut("prod", 2000, 0.975, max(0.5, 2 * ht[1] - 1.0))}


def true_cuts_mmd(d: dict, h: float) -> dict:
    h1, h2 = layer_true(h, 1), layer_true(h, 2)
    rl1 = d["T_1"] / K.sn_cut("split1", 2000, 0.975, h1)
    rl2 = np.maximum(d["T_1"] / K.sn_cut("split1", 2000, 0.9875, h1), d["T_2"] / K.sn_cut("split2", 2000, 0.9875, h2))
    est_c1 = np.array([K.sn_cut("split1", 2000, 0.975, v) for v in d["hl_1"]])
    return {"rl1_true": rl1, "rl2_true": rl2, "c11_est": est_c1, "c11_true": K.sn_cut("split1", 2000, 0.975, h1)}


def part_apply() -> None:
    out = []
    sources = []
    for h in TRUE_H:
        sources.append(("fresh", "ind", h, os.path.join(HERE, "layer_true_rows", f"ind__H{h}__n2000.csv")))
        sources.append(("fresh", "mmd", h, os.path.join(HERE, "layer_true_rows", f"mmd__H{h}__n2000.csv")))
        sources.append(("main", "ind", h, os.path.join(HERE, "main", f"ind__G0__null0__H{h}__n2000.csv")))
        sources.append(("main", "mmd", h, os.path.join(HERE, "main", f"mmd__G0__null0__H{h}__n2000.csv")))
    for src, fam, h, path in sources:
        if not os.path.exists(path):
            continue
        d = load_cols(path)
        t = true_cuts_ind(d, h) if fam == "ind" else true_cuts_mmd(d, h)
        m = d["rl1"].size
        base = {"source": src, "family": fam, "H": h, "m": m,
                "layer1_est_mean": float(np.mean(np.r_[d["hl_x1"], d["hl_y1"]] if fam == "ind" else d["hl_1"])),
                "layer2_est_mean": float(np.mean(np.r_[d["hl_x2"], d["hl_y2"]] if fam == "ind" else d["hl_2"])),
                "layer1_true": layer_true(h, 1), "layer2_true": layer_true(h, 2),
                "cut_low1_est_median": float(np.median(t["c11_est"])), "cut_low1_true": float(t["c11_true"])}
        rules = [("low1", None), ("low2", None)] + ([(f"hsplit{r}|{ru}", (r, ru)) for r in (1, 2) for ru in ("shift", "wild")] if fam == "ind"
                                                    else [(f"hsplit{r}|wild", (r, "wild")) for r in (1, 2)])
        for name, spec in rules:
            if spec is None:
                r = int(name[-1])
                a, b = d[f"rl{r}"] > 1, t[f"rl{r}_true"] > 1
            else:
                r, ru = spec
                rem = d[f"p_{'shift' if ru == 'shift' else 'wb'}_herm{r}"] <= ALPHA / 2
                a, b = rem | (d[f"rl{r}"] > 1), rem | (t[f"rl{r}_true"] > 1)
            lo_a, hi_a = wilson(int(a.sum()), m)
            lo_b, hi_b = wilson(int(b.sum()), m)
            out.append({**base, "rule": name, "rate_est": float(a.mean()), "est_lo": lo_a, "est_hi": hi_a, "rate_true": float(b.mean()),
                        "true_lo": lo_b, "true_hi": hi_b, "only_est": int(np.sum(a & ~b)), "only_true": int(np.sum(~a & b)), "mcnemar_p": mcnemar(a, b)})
    with open(os.path.join(HERE, "layer_true_vs_est.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"wrote layer_true_vs_est.csv ({len(out)} rows)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=("bias", "fresh", "apply"), required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--reps", type=int, default=None)
    a = ap.parse_args()
    if a.workers > 4:
        sys.exit("at most 4 workers for this side job")
    if a.part == "bias":
        part_bias(a.workers, a.reps or 400)
    elif a.part == "fresh":
        part_fresh(a.workers, a.reps or 1000)
    else:
        part_apply()


if __name__ == "__main__":
    main()
