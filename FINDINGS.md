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
