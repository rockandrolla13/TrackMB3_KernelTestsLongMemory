"""Download the approved XNAS stock-trade pilot or six-month sample.

Raw Databento records stay under data/raw/ (gitignored). This script performs a
fresh cost check before any paid timeseries request and records file hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import databento as db


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data/raw/xnas_stock_trades_2024"
DATASET = "XNAS.ITCH"
SYMBOLS = ["AAPL", "MSFT", "NVDA"]
KEY_FILE = Path(
    os.environ.get(
        "DATABENTO_KEY_FILE",
        "/media/ak/10E1026C4FA6006E/GitRepos/marketdb/.databento_key",
    )
)
PERIODS = {
    "pilot": [("2024-01-08", "2024-01-13")],
    "six_months": [
        ("2024-01-01", "2024-02-01"),
        ("2024-02-01", "2024-03-01"),
        ("2024-03-01", "2024-04-01"),
        ("2024-04-01", "2024-05-01"),
        ("2024-05-01", "2024-06-01"),
        ("2024-06-01", "2024-07-01"),
    ],
}
DEFAULT_MAX_COST = {"pilot": 1.0, "six_months": 18.0}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=PERIODS)
    parser.add_argument("--max-usd", type=float)
    args = parser.parse_args()
    if not KEY_FILE.is_file():
        raise FileNotFoundError(f"Databento key file missing: {KEY_FILE}")

    periods = PERIODS[args.phase]
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    client = db.Historical(KEY_FILE.read_text().strip())
    quote = client.metadata.get_cost(
        dataset=DATASET,
        symbols=SYMBOLS,
        schema="trades",
        stype_in="raw_symbol",
        start=periods[0][0],
        end=periods[-1][1],
    )
    cap = DEFAULT_MAX_COST[args.phase] if args.max_usd is None else args.max_usd
    print(f"{args.phase}: quoted ${quote:.4f}; cap ${cap:.2f}", flush=True)
    if quote > cap:
        raise RuntimeError("Cost quote exceeds cap; no data requested")

    files = []
    for start, end in periods:
        target = DATA_DIR / f"xnas_trades_{start}_{end}.dbn.zst"
        if not target.is_file():
            temporary = DATA_DIR / f"xnas_trades_{start}_{end}.part.dbn.zst"
            if temporary.exists():
                temporary.unlink()
            client.timeseries.get_range(
                dataset=DATASET,
                schema="trades",
                symbols=SYMBOLS,
                stype_in="raw_symbol",
                start=start,
                end=end,
                path=temporary,
            )
            if not temporary.is_file() or temporary.stat().st_size == 0:
                raise RuntimeError(f"Empty Databento response for {start} to {end}")
            temporary.replace(target)
        item = {
            "start": start,
            "end_exclusive": end,
            "path": str(target.relative_to(ROOT)),
            "bytes": target.stat().st_size,
            "sha256": sha256(target),
        }
        files.append(item)
        print(f"{target.name}: {item['bytes']:,} bytes", flush=True)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "phase": args.phase,
        "dataset": DATASET,
        "schema": "trades",
        "symbols": SYMBOLS,
        "stype_in": "raw_symbol",
        "cost_quote_usd": quote,
        "cost_cap_usd": cap,
        "files": files,
    }
    output = ROOT / "experiments/STOCK-20261023" / f"manifest_{args.phase}.json"
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Manifest: {output}")


if __name__ == "__main__":
    main()
