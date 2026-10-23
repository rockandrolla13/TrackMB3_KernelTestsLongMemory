# Nasdaq stock trade-duration boundary test

This applied check uses Databento `XNAS.ITCH` trade records for AAPL, MSFT, and NVDA from January through June 2024. It does not supply a known independence null.

## Data and cost

- Source schema: `trades`, raw-symbol input, six monthly files. The exact requests, byte counts, and SHA-256 hashes are in `manifest_six_months.json`.
- The API quoted **$15.0041** before the six-month request. The separate five-session pilot was quoted **$0.4311**. These are pre-request quotes, not a claim about the final invoice.
- Raw DBN files and duration caches live under `data/raw/`, which Git ignores. Do not include licensed raw files in an Overleaf package.
- The pilot had 1,607,092 trade records. About a quarter of regular-session prints had unspecified aggressor side, so this arm uses durations only.

## Extraction

`extract_durations.py` keeps trades from 09:30 to 16:00 New York time. It sorts by event timestamp and feed sequence, collapses prints sharing one event timestamp, and computes positive gaps within each trading day. The adjusted value is the log gap in seconds minus that stock's pooled mean log gap in the relevant 15-minute clock bin. Raw gaps, timestamps, bins, and dates remain in the local cache.

`duration_extraction_summary.json` records 124 trading days for each stock and 9,680,821 AAPL, 6,566,505 MSFT, and 14,212,862 NVDA durations. Every clock bin has observations. Median day-level shares of unspecified side are about 19.5%, 25.3%, and 36.0%, respectively.

## Rebuild

From the repository root, with a valid local Databento key:

```sh
python experiments/STOCK-20261023/download_xnas.py pilot
python experiments/STOCK-20261023/download_xnas.py six_months
python experiments/STOCK-20261023/extract_durations.py
python experiments/STOCK-20261023/run_stock_placebo.py --pairs 500 --workers 6 --b-wild 500 --out experiments/STOCK-20261023/results
python experiments/STOCK-20261023/summarize_stock.py
```

The downloader checks the current API cost before requesting data. Its default caps are $1 for the pilot and $18 for the six-month sample. Re-running it with existing files computes hashes and does not request those files again.

## Test interpretation

Each pair uses two windows of 500, 1,000, or 2,000 adjusted durations. Each window stays in one trading day. The days are at least 20 trading sessions apart. The pairings are AAPL–MSFT, AAPL–NVDA, MSFT–NVDA, and AAPL–AAPL, with 500 pairs per cell. The test runner reuses the existing HSIC, distance-covariance, and Pearson implementations. Seed, pairing, dates, starts, and per-pair outputs are saved.

The day gap prevents window overlap. It does not prove the stock processes are independent. Reused days also make pair outcomes dependent. Treat rejection fractions as descriptive market placebos, not measured false-alarm probabilities.
