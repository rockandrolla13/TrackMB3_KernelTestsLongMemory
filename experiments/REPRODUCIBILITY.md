# Reproducibility of the saved experiment outputs

Written 2026-11-03 for the `kmem` refactor (`reviews/2026_11_03_refactoring_plan.md`, steps 0.1 and 0.3).

## What is pinned

- `experiments/OUTPUTS.sha256` holds the sha256 of every saved output file in EXP2 to EXP7 (878 files),
  taken on 2026-11-03 before any refactor step. Code, shell scripts, `DESIGN.md` and `RESULTS.md` are left
  out because they may be edited. Check with `sha256sum -c experiments/OUTPUTS.sha256` from the repo root.
  Raw main-run CSVs are not committed (EXP7 alone is 127 MB); this list pins them instead (decision D2).
- `tests/legacy.sha256` pins `kmem/legacy`, which must never change.
- `tests/repro/fixtures/EXP*.json` hold a committed sample of 439 saved rows (305 fast tier, 134 slow tier).
  Every saved main-run cell of EXP2 to EXP7 appears once.
- `experiments/EXP7-memory-kernel-20261101/sn_table.json` is committed because EXP7 rows need it at run time.
  `sn_tables.py` (master 2026110174) produced it.

## The gate

`PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -I -m pytest tests -q` runs the fast tier. Add
`KMEM_REPRO_TIER=slow` (or `all`) for the rest. A row passes when every saved column equals the recomputed
value (`float(saved) == float(new)`, or both NaN, or equal text), after EXP7's 12-digit formatter. The gate
skips only machine measurements (`cpu_s`, `wall_s`, `maxrss_mb`), the `error` text, and EXP4's `reused`
flag, which depends on how a run was chunked.

## EXP1 is outside the gate (decision D5)

`EXP1-20261016/manifest.json` records `run_e1.py` with sha256 `ef5f7e5c9e1ca06cce8ffcf8eef464270ae3b6452392bcbcfa98dfc61028f6fa`.
The file on disk hashes to `08308ecbbda53f2745cc23d4759fbc3cedb5986657c4e0c23f89bb9d2ca4219f`; it was edited
on 2026-10-23, after E1 ran. E1 also reads licensed market data from `e1/data_cache/`. EXP1 stays frozen
as it is, including its copies of the legacy core.
