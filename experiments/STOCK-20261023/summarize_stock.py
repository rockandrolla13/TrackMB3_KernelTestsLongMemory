"""Audit saved stock-placebo pairs and summarize descriptive rejection rates."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
TESTS = {
    "HSIC wild": lambda row: float(row["p_hsic_wb"]) <= 0.05,
    "HSIC shift": lambda row: float(row["p_hsic_shift"]) <= 0.05,
    "HSIC subsampling": lambda row: int(row["rej_hsic_sub"]) == 1,
    "dCov wild": lambda row: float(row["p_dcov_wb"]) <= 0.05,
    "dCov subsampling": lambda row: int(row["rej_dcov_sub"]) == 1,
    "Pearson NW": lambda row: int(row["rej_nw"]) == 1,
    "Pearson EWC": lambda row: int(row["rej_ewc"]) == 1,
}


def main() -> None:
    manifest = json.loads((RESULTS / "manifest.json").read_text())
    cache = {}
    for symbol in ("AAPL", "MSFT", "NVDA"):
        with np.load(ROOT / f"data/raw/stock_duration_cache_2024/{symbol}.npz") as z:
            days = z["days"].astype(str)
            day_idx = z["dur_day"].astype(int)
        cache[symbol] = {"days": days, "first": np.searchsorted(day_idx, np.arange(len(days) + 1))}
    calendar = cache["AAPL"]["days"]
    if not all(np.array_equal(calendar, cache[s]["days"]) for s in cache):
        raise ValueError("Stock calendars differ")
    day_position = {day: i for i, day in enumerate(calendar)}

    summaries = []
    errors = []
    for cell in manifest["cells"]:
        a, b, n = cell["a"], cell["b"], cell["n"]
        path = RESULTS / f"{a}-{b}_n{n}.csv"
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        pair_ids = [int(row["pair"]) for row in rows]
        if sorted(pair_ids) != list(range(manifest["pairs_per_cell"])):
            errors.append(f"{path.name}: missing, duplicate, or extra pair IDs")
        good = [row for row in rows if row["status"] == "ok"]
        if len(good) != len(rows):
            errors.append(f"{path.name}: {len(rows) - len(good)} failed rows")
        left_days, right_days = Counter(), Counter()
        for row in good:
            da, db = row["day_a"], row["day_b"]
            i, j = day_position[da], day_position[db]
            if abs(i - j) < manifest["min_gap_trading_days"]:
                errors.append(f"{path.name}: pair {row['pair']} has short date gap")
            for symbol, day, start in ((a, i, int(row["start_a"])),
                                       (b, j, int(row["start_b"]))):
                lo, hi = cache[symbol]["first"][day : day + 2]
                if not (lo <= start and start + n <= hi):
                    errors.append(f"{path.name}: pair {row['pair']} crosses a day")
            left_days[da] += 1
            right_days[db] += 1
        rates = {name: round(100 * sum(test(row) for row in good) / len(good), 2)
                 for name, test in TESTS.items()} if good else {}
        summaries.append({
            "pairing": f"{a}-{b}", "n": n, "rows": len(rows), "valid": len(good),
            "distinct_a_days": len(left_days), "distinct_b_days": len(right_days),
            "maximum_reuse_a": max(left_days.values(), default=0),
            "maximum_reuse_b": max(right_days.values(), default=0),
            "rejection_percent": rates,
        })
    if errors:
        raise ValueError("Audit failed: " + "; ".join(errors[:20]))
    report = {
        "data": "XNAS.ITCH AAPL, MSFT, NVDA regular-session trades, January–June 2024",
        "interpretation": "Exploratory placebo; process independence is unverified; days are reused.",
        "nominal_percent": 5,
        "cells": summaries,
    }
    (HERE / "stock_placebo_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Stock trade-duration placebo", "", report["interpretation"], "",
             "| Pairing | n | Valid pairs | HSIC wild | HSIC shift | HSIC sub | dCov wild | dCov sub |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for cell in summaries:
        rates = cell["rejection_percent"]
        lines.append("| " + " | ".join([
            cell["pairing"], str(cell["n"]), str(cell["valid"]),
            *[f"{rates[name]:.1f}%" for name in
              ("HSIC wild", "HSIC shift", "HSIC subsampling", "dCov wild", "dCov subsampling")],
        ]) + " |")
    lines.extend(["", "Each window stays within one session, and paired days are at least 20 trading sessions apart.",
                  "Reused days and unknown process independence prevent false-alarm interpretation.", ""])
    (HERE / "stock_placebo_results.md").write_text("\n".join(lines))
    print(f"Audited {len(summaries)} cells and {sum(s['valid'] for s in summaries):,} pairs")


if __name__ == "__main__":
    main()
