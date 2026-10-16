# Job JOB-NEW-MEMORY-BOUNDARY-SIM-20261015 — memory breakdown map simulation

Status: closed — all three families finished; results copied to TrackMB1–3  
Track: NEW-MEMORY-BOUNDARY  
Owner: Claude Code coordinator; three workers — `claude-mb-kernel-20261015` (`kernel-linear/`: K-ind HSIC/dCov, L-ind Pearson, G1–G4, R1), `claude-mb-mmddml-20261015` (`mmd-dml/`: K-hom MMD, DML), `claude-mb-hawkes-20261015` (`hawkes/`: HG on HK1–HK3); coordinator merges  
Date opened: 2026-10-15

## Question

Run the approved design `runs/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/simulation-design.md` (SHA-256 `7bb1d6b9b2775a7f74dd57ad1c88d55e59fb57d7051093b5dc070176c611559f`, Decision `D-ac49ecfe79`) and evaluate predictions P1–P5 and kill criteria KC1–KC5.

## Inputs

Approved design above; `runs/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/` (closed).

## Required skill and version

Route `theory_method` (proof-methodology `232534b5…`). Simulation run without a locked protocol: exploratory; design fixed and approved before any run. Environment: ResearchHarness venv, numpy only (no scipy/sklearn); estimators implemented in numpy. Workers log `skill_loaded` with actual hashes.

## Dependencies

`check-assignment NEW-MEMORY-BOUNDARY` allowed 2026-10-15 immediately before assignment.

## Output

Each worker writes only its subdirectory: scripts, raw cell results (CSV), `results.md`, `handoff.md`, `events.jsonl`. Coordinator writes `synthesis.md`.

## Completion condition

Every fixed cell or item reported, including failures; `parameters_fixed` before computing; scripts reproduce with `python -I` and fixed seeds; `finished` event with output hashes. Exploratory; no claim, protocol or Gate change.

## Required handoff

Author, via coordinator, in plain English.

## Scope

CPU: kernel-linear ≤ 8 cores, mmd-dml ≤ 6, hawkes ≤ 4; cap 48 CPU-hours per family with the design's timing-pilot and fallback rule. No commits, ledger writes or Gate actions.

## Execution checkpoints

Worker-owned `events.jsonl` in each worker directory; `checkpoint.md` if interrupted.

## State, result and next action

Assigned 2026-10-15.

Author amendment 2026-10-16: kernel-linear pilot projected ~160 CPU-h (2,000 reps) and ~106 CPU-h (design fallback) against the 48 CPU-hour cap. Author chose to keep the cap and drop n = 4000 for K-ind cells; n = 2000 runs at 2,000 reps; P1, KC1 and KC2 are judged at n = 2000. Other families unchanged.

Coordinator check, hawkes family (2026-10-16): `finished` event under `claude-mb-hawkes-20261015`; 92 output hashes match; scripts import only from their own directory; `analyze.py` rerun on a scratch copy reproduces `summary_cells.json` and `table_cells.md` byte for byte. Result ([handoff](JOB-NEW-MEMORY-BOUNDARY-SIM-20261015/hawkes/handoff.md) `d1faa555…`): 2,000 reps in all 36 cells, 26.6 CPU-h. Exponential kernel, correctly specified (HK1): size 4.3–5.5% up to ρ = 0.99. Power-law kernel, exponential model (HK2): 10–14% at ρ = 0.95 and 14–28% at ρ = 0.99, growing with N. Shared slow baseline (HK3): 87–100% at ρ = 0.5, 12–87% at ρ = 0.95. KC5 not triggered; P5 supported (size depends on ρ, kernel tail, baseline). Logged deviations: power-law burn-in slightly short at ρ = 0.99 (counts ≤ 4.9% below nominal); power-cell burn-in 44.4 vs 40; gradient-flag tolerance artefact at N = 50,000 (no reps dropped; ≤ 0.26 pp effect); first immigrant-share check replaced by `validation_supp.json`. HK3 run at one amplitude and period only.

Coordinator check, mmd-dml family (2026-10-16): `finished` event under `claude-mb-mmddml-20261015`; 440 output hashes match; scripts import only from their own directory; `analyze.py main` rerun on a scratch copy from cleared outputs reproduces `tables.md`, `analysis.json` and both summaries byte for byte. Result ([handoff](JOB-NEW-MEMORY-BOUNDARY-SIM-20261015/mmd-dml/handoff.md) `8a1447a2…`): MMD half-vs-half size 12–19% at H = 0.6 and 44–71% at H = 0.9 (no 3/4 line; wild bootstrap oversized 2–4 pp at H = 0.5). DML population-ATE coverage at T = 8000: monotone effect 84% (H 0.6) → 12% (H 0.9); even effect fine to H = 0.7, 88% at 0.75, 26% at 0.9; constant effect fine. Nuisance learning and fold choice move coverage ≤ 2 pp. Path-conditional ATE coverage ≥ 93.45% in every cell. KC3, KC4 not triggered; P2 holds in substance (knots 0.50–0.56), P3 fails strictly at H = 0.75 only, P4 holds. 27.1 CPU-h; 0 errors. Logged deviations: BLAS-threaded first pilot discarded; K-hom fallback 1,000 reps at n = 4000; GB settings; MMD subsampling window choice; analyze.py vectorised after parameters_fixed; Newey–West convention differs from kernel-linear (clipped ρ, weights 1 − j/(S+1)).

Coordinator check, kernel-linear family (2026-10-16): `finished` event under `claude-mb-kernel-20261015`; 345 output hashes match; scripts import only from their own directory; `analyze.py --cells cells --out cells` rerun on a scratch copy from cleared outputs reproduces `cell_summary.csv` and `knot_analysis.json` byte for byte. Result ([handoff](JOB-NEW-MEMORY-BOUNDARY-SIM-20261015/kernel-linear/handoff.md) `d019215b…`): knots for all seven tests at H_X + H_Y ≈ 1.4–1.6 (every CI contains 1.5); growth with n beyond the line weak for kernel tests up to n = 2000 (P1 not supported for kernels; Pearson-EWC supports it); shift test best calibrated; subsampling worst; KC1, KC2 not triggered. 29.6 CPU-h, 0 errors. Worker missed completion notices twice (≈ 2 h and 3.5 h delays).

## Job close, 2026-10-16

Status: all three families finished and checked. Results copied to the split tracks: DML → TrackMB1, Hawkes → TrackMB2, MMD and kernel-linear → TrackMB3. Exploratory simulation; no claim.
