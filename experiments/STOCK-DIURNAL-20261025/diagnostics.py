"""Within-window variance-trend diagnostics for the n = 2000 placebo windows, and a D3 spot check.

Windows come from the original STOCK-20261023 pair CSVs (identical across adjustments).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/EXP1-20261016"))
sys.path.insert(0, str(HERE))
import mbcore as mb  # noqa: E402
import run_diurnal_placebo as rp  # noqa: E402

N = 2000
SEGMENTS = 8
SEED = 42


def segment_stats(w: np.ndarray) -> tuple[float, float]:
    """OLS slope of log segment variance on segment position (0..1), and max/min variance."""
    v = w.reshape(SEGMENTS, -1).var(axis=1)
    pos = np.linspace(0, 1, SEGMENTS)
    slope = np.polyfit(pos, np.log(v), 1)[0]
    return float(slope), float(v.max() / v.min())


def main() -> None:
    rp.load_data(("D1", "D2"))
    windows = []
    for a, b in rp.PAIRS:
        with (ROOT / f"experiments/STOCK-20261023/results/{a}-{b}_n{N}.csv").open() as s:
            for row in csv.DictReader(s):
                windows.append((a, int(row["start_a"])))
                windows.append((b, int(row["start_b"])))
    rng = np.random.default_rng(SEED)
    perms = [rng.permutation(N) for _ in windows]
    out = {"n": N, "segments": SEGMENTS, "windows": len(windows), "by_adjustment": {}}
    for adj in rp.ADJ:
        slopes, ratios, null_slopes = [], [], []
        for (sym, st), p in zip(windows, perms):
            w = rp.DATA[sym][adj][st:st + N]
            s, r = segment_stats(w)
            slopes.append(s)
            ratios.append(r)
            null_slopes.append(segment_stats(w[p])[0])
        slopes, null_slopes = np.abs(slopes), np.abs(null_slopes)
        cut = np.quantile(null_slopes, 0.95)
        out["by_adjustment"][adj] = {
            "median_abs_slope": float(np.median(slopes)),
            "median_abs_slope_permuted": float(np.median(null_slopes)),
            "share_exceeding_permuted_95pct": float(np.mean(slopes > cut)),
            "median_seg_var_max_over_min": float(np.median(ratios)),
        }
        print(adj, out["by_adjustment"][adj], flush=True)

    # D3: subtract each day's mean from D2. Tests are location invariant, so p-values must match.
    checks = []
    a, b, n, cidx = "AAPL", "MSFT", 500, 0
    for sym in {a, b}:
        first, x = rp.DATA[sym]["first"], rp.DATA[sym]["D2"]
        day = np.repeat(np.arange(first.size - 1), np.diff(first))
        means = np.bincount(day, weights=x) / np.diff(first)
        rp.DATA[sym]["D3"] = x - means[day]
    for k in range(20):
        r2 = rp.one_pair(a, b, n, cidx, k, "D2", 500)
        r3 = rp.one_pair(a, b, n, cidx, k, "D3", 500)
        keys = ("p_hsic_wb", "p_hsic_shift", "p_hsic_sub", "p_dcov_wb", "p_dcov_sub", "p_nw", "rej_ewc")
        checks.append(max(abs(float(r2[q]) - float(r3[q])) for q in keys))
    out["d3_spot_check"] = {"pairs": 20, "cell": "AAPL-MSFT n=500",
                            "max_abs_pvalue_difference_vs_D2": float(max(checks))}
    print(out["d3_spot_check"])
    (HERE / "diagnostics.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
