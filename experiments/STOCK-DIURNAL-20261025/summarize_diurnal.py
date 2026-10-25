"""Audit D0 against the original run and tabulate rejection percentages by adjustment."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import diagnostics as dg  # noqa: E402
import run_diurnal_placebo as rp  # noqa: E402

ORIG = ROOT / "experiments/STOCK-20261023/results"
RES = HERE / "results"
TESTS = {
    "HSIC wild": lambda r: float(r["p_hsic_wb"]) <= 0.05,
    "HSIC shift": lambda r: float(r["p_hsic_shift"]) <= 0.05,
    "HSIC sub": lambda r: int(r["rej_hsic_sub"]) == 1,
    "dCov wild": lambda r: float(r["p_dcov_wb"]) <= 0.05,
    "dCov sub": lambda r: int(r["rej_dcov_sub"]) == 1,
    "Pearson NW": lambda r: int(r["rej_nw"]) == 1,
    "Pearson EWC": lambda r: int(r["rej_ewc"]) == 1,
}


def read(path: Path) -> dict:
    with path.open(newline="") as s:
        return {int(r["pair"]): r for r in csv.DictReader(s)}


def main() -> None:
    rows, mismatches, cpu = [], {}, {}
    for adj in rp.ADJ:
        cpu[adj] = 0.0
        for a, b in rp.PAIRS:
            for n in (500, 1000, 2000):
                name = f"{a}-{b}_n{n}.csv"
                # Amendment 1: D0 reran only at n = 2000. Smaller D0 rows are the original run.
                rerun = adj != "D0" or n == 2000
                old = read(ORIG / name)
                new = read(RES / adj / name) if rerun else old
                if sorted(new) != list(range(500)) or any(r["status"] != "ok" for r in new.values()):
                    raise ValueError(f"{adj}/{name}: incomplete or failed rows")
                if any(new[k][f] != old[k][f] for k in new for f in ("day_a", "day_b", "start_a", "start_b")):
                    raise ValueError(f"{adj}/{name}: windows differ from the original run")
                if adj == "D0" and rerun:
                    bad = sum(any(new[k][f] != old[k][f] for f in rp.FIELDS if f != "cpu_s") for k in new)
                    mismatches[name] = bad
                if rerun:
                    cpu[adj] += sum(float(r["cpu_s"]) for r in new.values())
                good = list(new.values())
                rows.append({"adjustment": adj, "source": "rerun" if rerun else "original", "pairing": f"{a}-{b}", "n": n, "pairs": len(good),
                             **{t: round(100 * np.mean([f(r) for r in good]), 1) for t, f in TESTS.items()}})
    with (HERE / "summary.csv").open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # Do pairs with a within-window variance trend in both windows reject more? (n = 2000, D2)
    rp.load_data(("D2",))
    rng = np.random.default_rng(dg.SEED)
    null = []
    for _ in range(4000):
        sym = rp.SYMBOLS[int(rng.integers(3))]
        st = int(rng.integers(0, rp.DATA[sym]["D2"].size - dg.N))
        null.append(abs(dg.segment_stats(rng.permutation(rp.DATA[sym]["D2"][st:st + dg.N]))[0]))
    cut = float(np.quantile(null, 0.95))
    groups = {}
    for a, b in rp.PAIRS:
        for r in read(RES / "D2" / f"{a}-{b}_n2000.csv").values():
            ta = abs(dg.segment_stats(rp.DATA[a]["D2"][int(r["start_a"]):int(r["start_a"]) + dg.N])[0]) > cut
            tb = abs(dg.segment_stats(rp.DATA[b]["D2"][int(r["start_b"]):int(r["start_b"]) + dg.N])[0]) > cut
            groups.setdefault(int(ta) + int(tb), []).append(r)
    trend = {f"{k}_windows_trending": {"pairs": len(v), **{t: round(100 * np.mean([TESTS[t](r) for r in v]), 1)
                                                         for t in ("HSIC wild", "HSIC shift", "HSIC sub", "dCov sub", "Pearson NW")}}
             for k, v in sorted(groups.items())}

    diag = json.loads((HERE / "diagnostics.json").read_text())
    diag.update({"d0_regression_mismatched_rows": mismatches, "cpu_seconds_by_adjustment": cpu,
                 "rejection_by_trend_count_D2_n2000": trend, "trend_cut_abs_slope": cut})
    (HERE / "diagnostics.json").write_text(json.dumps(diag, indent=2) + "\n")

    fit = json.loads((HERE / "diurnal_fit.json").read_text())
    cols = list(TESTS)
    lines = ["# STOCK-DIURNAL tables", "",
             "Rejection percentages at a nominal 5% level. 500 pairs per row. Same windows and seeds in every adjustment.", ""]
    for n in (2000, 1000, 500):
        lines += [f"## n = {n}", "", "| Pairing | Adj | " + " | ".join(cols) + " |",
                  "| --- | --- | " + " | ".join("---:" for _ in cols) + " |"]
        for a, b in rp.PAIRS:
            for row in rows:
                if row["pairing"] == f"{a}-{b}" and row["n"] == n:
                    lines.append(f"| {a}-{b} | {row['adjustment']}{'' if row['source'] == 'rerun' else ' (orig)'} | " + " | ".join(f"{row[c]:.1f}" for c in cols) + " |")
        lines.append("")
    lines += ["## Diurnal scale range (max/min sd across the session)", "",
              "| Stock | Adj | Smooth sd max/min | 15-min bin sd max/min |", "| --- | --- | ---: | ---: |"]
    for sym, v in fit["symbols"].items():
        for adj in rp.ADJ:
            sc = v[f"scale_{adj}"]
            lines.append(f"| {sym} | {adj} | {sc['smooth_sd_max_over_min']:.3f} | {sc['bin_sd_max_over_min']:.3f} |")
    lines += ["", "## Within-window variance trend (n = 2000 windows, 8 segments)", "",
              "| Adj | Median abs slope | Same, permuted | Share above permuted 95% | Median seg-var max/min |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for adj, d in diag["by_adjustment"].items():
        lines.append(f"| {adj} | {d['median_abs_slope']:.3f} | {d['median_abs_slope_permuted']:.3f} | "
                     f"{100 * d['share_exceeding_permuted_95pct']:.1f}% | {d['median_seg_var_max_over_min']:.2f} |")
    lines += ["", "## D2, n = 2000: rejection by number of windows with a variance trend", "",
              "| Trending windows | Pairs | HSIC wild | HSIC shift | HSIC sub | dCov sub | Pearson NW |",
              "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for k, d in trend.items():
        lines.append(f"| {k.split('_')[0]} | {d['pairs']} | {d['HSIC wild']:.1f} | {d['HSIC shift']:.1f} | "
                     f"{d['HSIC sub']:.1f} | {d['dCov sub']:.1f} | {d['Pearson NW']:.1f} |")
    lines += ["", f"D0 regression: mismatched rows per cell = {sum(mismatches.values())} in total.",
              f"D3 spot check: max p-value difference from D2 = {diag['d3_spot_check']['max_abs_pvalue_difference_vs_D2']:.1e}.",
              "CPU seconds: " + ", ".join(f"{k} {v:,.0f}" for k, v in cpu.items()), ""]
    (HERE / "tables.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
