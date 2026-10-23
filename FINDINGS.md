# Findings — TrackMB3 (as of 2026-10-16)

All exploratory. Evidence class: simulation (S), literature (F/S/A), derivation (D).

## Simulation: MMD homogeneity (`experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_mmd-dml/`, K-hom part, handoff `8a1447a2…`)

MMD², first half against second half of one series, on fGn (G1), bounded subordinated (G2), sign (G3) and the rough-vol arm (R1). 2,000 reps (1,000 at n = 4000, design fallback). Coordinator reran the analysis byte for byte.

- **Breaks under any long memory:** rejection at nominal 5% is 12–19% at H = 0.6 and 44–71% at H = 0.9, for all processes and both calibrations. There is no 3/4 line, as predicted for a rank-1 statistic. Knots sit at H ≈ 0.50–0.56.
- The wild bootstrap is oversized by 2–4 pp even at H = 0.5: a calibration flaw independent of memory.

## Simulation: HSIC / dCov / Pearson independence (`experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear/`, handoff `d019215b…`)

HSIC (wild bootstrap, shift test, subsampling), distance covariance (wild bootstrap, subsampling), Pearson (Newey–West, EWC fixed-b) on fGn, subordinated and sign processes, unequal pairs and the rough-vol arm. 2,000 reps per null cell, 294,000 replications, 0 errors, 29.6 CPU-h. Kernel n = 4000 dropped by author decision D-a28e4cb5fa; kernel verdicts at n = 2000. Coordinator reran the analysis byte for byte (`cell_summary.csv`, `knot_analysis.json`).

- **The line sits where predicted.** Size worsens sharply once H_X + H_Y passes about 1.5. Every test's knot falls between 1.4 and 1.6, and every knot CI contains 1.5. The H grid is coarse near 1.5, so the location is only that precise.
- **But P1's growth-with-n part is not supported for kernel tests.** Beyond the line, size barely grows with n up to 2000; it grows clearly only at H = 0.9. Pearson-EWC follows P1 exactly. Pearson-NW is already oversized (7–9%) just below the line.
- **Size at n = 2000:** H ≤ 0.7 mostly 5–9%. H 0.75–0.8: wild bootstrap 6–10%, subsampling 10–25%, NW 10–15%. H = 0.9: wild bootstrap 13–18%, subsampling 34–52%, NW 21–34%.
- **Best calibration:** the shift test (3.9–7.6% everywhere), then EWC for the linear test below the line. Subsampling, the only calibration with a long-memory proof, is the worst. The wild bootstrap is up to 1.6 pp oversized even on i.i.d. data.
- **Rough-vol arm:** slow mean reversion (λ = 0.01) fakes long memory (subsampling 43–46%, NW 18%), but the distortion shrinks with n.
- KC1 and KC2 not triggered. Power check (H = 0.7, r = 0.25): 99–100% for every test.

## Theory and literature (`sources/`)

- Betken–Dehling 2107.03041 (F): distance covariance with subsampling under LRD; normal limit below the line, non-central above. Zhai–Zhong–Wu 2609.30565 (F): sample covariance boundary. For two series the condition is H_X + H_Y vs 3/2.
- Chwialkowski–Sejdinovic–Gretton (F): the wild bootstrap needs τ(r) = O(r⁻⁶) mixing.
- Extending the boundary to HSIC is a worker inference (D), not a stated theorem.
- Even below the line in practice, Betken–Dehling's own test rejects 10–17% at a 5% level at H = 0.7 with n ≤ 1000.

## Open

SSRN and Scholar not searched. No real-data test.

## EXP1 (2026-10-16/17) — first exploratory experiments (`experiments/EXP1-20261016/`, appended 2026-10-18)

Job `JOB-NEW-MB-KERNEL-EXP1-20261016`; executions `claude-mb-kernel-exp1-20261016` (died 2026-10-17 09:47Z) and `claude-mb-kernel-exp1r-20261017` (resumed, same seeds, completed rows kept). Exploratory: simulation (S) and real-data placebo (R). 32.1 of 56 CPU-h. Full tables in `results.md`; plain-English summary in `handoff.md`.

- **H5 (S):** the MMD wild-bootstrap oversize at H = 0.5 is a centring problem. The old scheme and the published centred form give 6.1–9.1%; i.i.d. and uncentred multipliers are on size (4.4–5.3%).
- **E1 (R), LQD/HYG/JNK, 2023-05-30..2025-05-29, 500 independent placebo pairs per cell, 4 of 12,000 failed (constant window):** trade signs carry strong memory (per-window Ĥ_X + Ĥ_Y about 1.5–1.7; long-sample Ĥ is bandwidth-sensitive, 0.66–1.02). At n = 2000 HSIC wild rejects 8.0–9.2%, the shift test 5.8–6.8%, subsampling 17–23%, NW 8–10%, EWC 7–8%. This matches the simulated map at sums 1.5–1.6 and does not grow with n. Adjusted durations (sum about 1.33–1.40) reject more than memory alone predicts (HSIC wild 8.6–13.8%, subsampling 15–25%), rising with n in 2 of 4 pairings and with within-cell Ĥ in all 4. Against the plan's reads the result is mixed: neither "near 5%" nor a clean memory-driven rise.
- **E2 (S), MMD:** the Hermite-rank prediction holds. With the on-size wild bootstrap, ξ² − 1 and |ξ| keep size within 2 pp up to H = 0.7 (ξ² − 1 also at 0.75) and break at H ≥ 0.8, while rank-1 series break from H = 0.55. No calibration meets the fix rule (5 ± 2 pp in every G1/G3 cell up to H = 0.8 at n = 2000). Self-normalised MMD with block sampling is nearly flat in H (G1 6.0–8.5%, G3 4.7–7.5%) but oversized, with lower power (50% against 76–87%). Plug-in rate-corrected subsampling is conservative and weak; the Betken–Dehling block rule breaks like √n subsampling.
- **E3 (S), fine grid at H_X + H_Y = 1.3–1.7 (trigger fired), 2,000 reps:** the KC2 read (excess > 2 pp at sum ≤ 1.4) is met by HSIC and dCov subsampling (9.2% and 10.6% at (0.7, 0.7), flat in n) and by Pearson-NW. It is met only marginally by the wild bootstrap, which is driven by its own H = 0.5 oversize, and not at all by the shift test or EWC. Unequal pairs track equal pairs with the same sum, which supports the sum rule over the per-series rule. Knot CIs at n = 2000 barely narrow: subsampling [1.43, 1.575] and [1.375, 1.585], wild [1.40, 1.635].
- **Open for the author:** whether the subsampling KC2 read changes the 3/4 headline; whether the duration excess is a real-data hook; whether the SN-block MMD is worth pursuing as a partial fix.

**Author decision 2026-10-18 (coordinator note):** keep the 3/4 headline (H_X + H_Y = 3/2 line) for the wild-bootstrap, shift and EWC calibrations. Subsampling is reported separately as a calibration that fails at any memory level (KC2 fired: 9–11% at a sum of 1.4, not shrinking with n). Still open: the duration excess as a real-data hook; SN-block MMD as a partial fix.

**Author decision 2026-10-19 (coordinator note, appended 2026-10-22; for subsampling this replaces the 2026-10-18 note above; the ledger does not say whether the 2026-10-18 decision still stands for the wild-bootstrap, shift, EWC and Pearson-NW calibrations, so that is an open author question):** ledger decision `D-d359e15a48` (ResearchHarness `registry/events.jsonl`, 2026-10-19T21:22:27Z, topic `scope:NEW-MB-KERNEL:KC2-headline`): drop the 3/4 headline at least for subsampling tests, and pursue trade durations (rejection beyond memory alone) as the real-data hook. Rationale as recorded: KC2 fired in EXP1 (HSIC/dCov subsampling 9–11% at H_X + H_Y = 1.4, n = 2000, not shrinking with n). Left open: how far the drop extends beyond subsampling (Pearson-NW +3.1 pp and HSIC wild +2.3 pp also exceeded 2 pp); SN-block MMD as a partial fix. Exploratory evidence only; no claim or protocol.

**Author decision 2026-10-22 (coordinator note):** ledger decision `D-786003ad54` (2026-10-22T00:39:40Z, topic `scope:NEW-MB-KERNEL:KC2-headline-other-calibrations`): keep the 3/4 headline for the wild-bootstrap, shift, EWC and Pearson-NW calibrations. So `D-849b9e34d7` stands for those calibrations; `D-d359e15a48` governs subsampling (headline dropped) and the durations hook. Still open: SN-block MMD as a partial fix. Exploratory evidence only; no claim or protocol.
