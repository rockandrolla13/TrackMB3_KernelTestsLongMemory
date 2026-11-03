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

## Rerunning a past experiment

From the repo root (kmem is not installed in the shared environment):

    python -m kmem reproduce EXP5                    # the committed sample, through the spec (minutes)
    python -m kmem reproduce EXP6 --full --dry-run   # list the rows a full rerun writes
    python -m kmem reproduce EXP6 --full --out <folder> --workers 5   # full rerun, resumable, then compared
    python -m kmem run EXP7 --out <folder> --workers 8                # just run, with the saved seeds

Each `experiments/EXPn/spec.json` holds the final replications per cell that ended up on disk, so a full
rerun needs none of the staged scripts the original runs used: EXP4's `run_main.sh` and `run_topup.sh`, and
EXP5's `run_main.sh`, `run_part1b.sh`, `run_part2.sh` and `run_part3.sh`. Those stay as records of how the
saved outputs were made. EXP5's `run_part3.sh` ran its n = 8000 stages only before 21:25 or 21:35 UTC, so
rerunning it at another hour writes no n = 8000 rows; `kmem reproduce EXP5 --full` has no such dependence.
EXP7 needs `sn_table.json` first; its sha256 is in the spec and is checked before every run.

Full reruns cost (CPU-hours): EXP2 26.4, EXP3 9.9, EXP4 8.8, EXP5 8.4, EXP6 4.3, EXP7 43.8.

## The analysis gate (G-analysis)

`python -I tools/analysis_gate.py --check` reruns all 18 saved analysis jobs of EXP2 to EXP7 in a scratch copy
(about 3 minutes on 8 cores) and compares the 50 files they write with `tests/analysis_baseline.json`.
At the baseline (2026-11-03), 46 of the 50 regenerate byte for byte equal to the saved files. The four that do not:

- EXP4 `tables_pilot.md` and `verdicts_pilot.json`, and EXP7 `growth_q95_pilot.png`: stable, but the saved
  copies came from an earlier code or input version. The gate checks them against the regenerated hash.
- EXP7 `bandwidth_check.csv` (`review_checks.py`): unstable. Its simulation stacks chunks in completion
  order, so means differ in the last bits. The gate only checks that it is written.

## EXP1 is outside the gate (decision D5)

`EXP1-20261016/manifest.json` records `run_e1.py` with sha256 `ef5f7e5c9e1ca06cce8ffcf8eef464270ae3b6452392bcbcfa98dfc61028f6fa`.
The file on disk hashes to `08308ecbbda53f2745cc23d4759fbc3cedb5986657c4e0c23f89bb9d2ca4219f`; it was edited
on 2026-10-23, after E1 ran. E1 also reads licensed market data from `e1/data_cache/`. EXP1 stays frozen
as it is, including its copies of the legacy core.
