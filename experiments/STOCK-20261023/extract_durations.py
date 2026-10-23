"""Build within-session stock trade-duration series from Databento DBN files.

The primary adjustment subtracts the pooled mean *log* duration in each 15-minute
New York clock bin. Raw gaps, event times, and bin labels remain in the cache.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import databento as db
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SYMBOLS = ("AAPL", "MSFT", "NVDA")
OPEN_MINUTE = 9 * 60 + 30
CLOSE_MINUTE = 16 * 60
BIN_MINUTES = 15


def main() -> None:
    manifest = json.loads((HERE / "manifest_six_months.json").read_text())
    if manifest["symbols"] != list(SYMBOLS) or manifest["schema"] != "trades":
        raise ValueError("Unexpected download manifest")
    by_symbol = {s: defaultdict(list) for s in SYMBOLS}
    day_summary = []

    for item in manifest["files"]:
        path = ROOT / item["path"]
        frame = db.DBNStore.from_file(path).to_df()
        if not {"symbol", "ts_event", "side", "sequence"}.issubset(frame.columns):
            raise ValueError(f"Missing required trade fields in {path}")
        local = frame["ts_event"].dt.tz_convert("America/New_York")
        minute = local.dt.hour * 60 + local.dt.minute
        frame["local_day"] = local.dt.date.astype(str).to_numpy()
        frame["clock_minute"] = minute.to_numpy()
        frame = frame.loc[(minute >= OPEN_MINUTE) & (minute < CLOSE_MINUTE)].copy()

        for (symbol, day), group in frame.groupby(["symbol", "local_day"], sort=True):
            if symbol not in by_symbol:
                continue
            group = group.sort_values(["ts_event", "sequence"])
            events = group.drop_duplicates("ts_event", keep="first")
            timestamps = events["ts_event"].astype("int64").to_numpy()
            if timestamps.size < 2:
                continue
            gaps = np.diff(timestamps)
            if np.any(gaps <= 0):
                raise ValueError(f"Nonpositive distinct-time gap for {symbol} {day}")
            bins = ((events["clock_minute"].to_numpy()[1:] - OPEN_MINUTE) // BIN_MINUTES)
            if np.any((bins < 0) | (bins >= 26)):
                raise ValueError(f"Invalid clock bin for {symbol} {day}")
            cell = by_symbol[symbol]
            cell["days"].append(day)
            cell["dur_ns"].append(gaps)
            cell["event_ts_ns"].append(timestamps[1:])
            cell["bin"].append(bins.astype(np.uint8))
            cell["day_idx"].append(np.full(gaps.size, len(cell["days"]) - 1, dtype=np.int16))
            day_summary.append({
                "symbol": symbol,
                "day": day,
                "trades_rth": len(group),
                "distinct_event_times": timestamps.size,
                "durations": gaps.size,
                "unspecified_side_share": round(float((group["side"] == "N").mean()), 6),
                "duplicate_time_share": round(1 - timestamps.size / len(group), 6),
                "median_positive_gap_seconds": round(float(np.median(gaps)) / 1e9, 9),
            })
        print(f"Read {path.name}: {len(frame):,} regular-session trade records", flush=True)
        del frame

    outdir = ROOT / "data/raw/stock_duration_cache_2024"
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {"source_manifest": str((HERE / "manifest_six_months.json").relative_to(ROOT)), "symbols": {}}
    for symbol in SYMBOLS:
        cell = by_symbol[symbol]
        days = np.asarray(cell["days"])
        if days.size == 0 or np.any(days[1:] <= days[:-1]):
            raise ValueError(f"Missing or unordered days for {symbol}")
        gaps = np.concatenate(cell["dur_ns"]).astype(np.int64)
        event_ts = np.concatenate(cell["event_ts_ns"]).astype(np.int64)
        bins = np.concatenate(cell["bin"]).astype(np.uint8)
        day_idx = np.concatenate(cell["day_idx"]).astype(np.int16)
        log_gaps = np.log(gaps.astype(float) / 1e9)
        sums = np.bincount(bins, weights=log_gaps, minlength=26)
        counts = np.bincount(bins, minlength=26)
        mean_log = sums / np.maximum(counts, 1)
        adjusted = (log_gaps - mean_log[bins]).astype(np.float32)
        output = outdir / f"{symbol}.npz"
        np.savez_compressed(
            output,
            dur=adjusted,
            dur_day=day_idx,
            days=days,
            event_ts_ns=event_ts,
            raw_dur_ns=gaps,
            clock_bin=bins,
            mean_log_duration_by_bin=mean_log,
        )
        summary["symbols"][symbol] = {
            "days": int(days.size),
            "first_day": str(days[0]),
            "last_day": str(days[-1]),
            "durations": int(gaps.size),
            "cache": str(output.relative_to(ROOT)),
            "clock_bin_counts": counts.tolist(),
        }
        print(f"{symbol}: {days.size} days, {gaps.size:,} durations", flush=True)

    summary["day_quality"] = day_summary
    (HERE / "duration_extraction_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
