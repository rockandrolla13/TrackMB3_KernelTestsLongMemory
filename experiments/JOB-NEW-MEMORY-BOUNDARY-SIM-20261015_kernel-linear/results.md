# Memory breakdown map — kernel-linear arm (K-ind, L-ind): results

Job JOB-NEW-MEMORY-BOUNDARY-SIM-20261015 · worker `claude-mb-kernel-20261015` · 2026-10-16.
Design `runs/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/simulation-design.md` (SHA-256 `7bb1d6b9…559f`, verified),
as amended by author decision **D-a28e4cb5fa**. **Evidence class: exploratory simulation.** No claim, protocol or Gate change.

## What was run

- **K-ind** (HSIC: wild bootstrap, shift test, subsampling; dCov: wild bootstrap, subsampling) at n ∈ {500, 1000, 2000}.
  **n = 4000 was dropped by author decision D-a28e4cb5fa** (option C) because neither 2,000 reps nor the design's
  1,000-rep fallback fit the 48 CPU-hour cap: pilot projections were 160 and 106 CPU-hours. For K-ind, P1, KC1, KC2 and
  the knot fit are judged at **n = 2000**, with the bands unchanged. The pilot's n = 4000 data (20 reps per cell, `pilot/`) is
  timing data only. It is not evidence.
- **L-ind** (Pearson cross-correlation, NW-HAC and EWC fixed-b) at T ∈ {500, 2000, 8000}, the design's grid for
  linear statistics.
- Processes: G1–G3 × H ∈ {0.5, 0.6, 0.7, 0.75, 0.8, 0.9}, G4 (4 unequal pairs), R1 (λ = 0.01, 0.1). That makes 24
  settings, with X ⊥ Y under the null.
- **2,000 replications per null cell** (MC s.e. at 5% = 0.49 pp; ±0.96 pp half-width), 1,000 per power cell.
  Every cell ran in full, with **0 errors** out of 294,000 replications.
- CPU: 28.0 CPU-hours for the evidence cells plus 1.6 for the pilot, **29.6 in total** (cap 48). Up to 8 processes,
  BLAS single-threaded.
- Seeds: `SeedSequence(entropy=20261015002, spawn_key=(cell_idx, rep))`. The pilot used entropy 20261015001.
  The per-cell manifest is `cells/manifest.json`. Re-running `python -I run_sim.py --mode main --reps 2000
  --power-reps 1000 --exclude "__n4000" --out cells` into an empty directory reproduces every CSV.

## Implementation (fixed in `parameters_fixed` before any run)

- **fGn**: Davies–Harte minimal circulant embedding. Validation (`validation.out`) confirms nonnegative eigenvalues for
  every H and n used, an exact reproduction of γ(k) (error ≤ 5e-16), sample autocovariances within 2.2 MC s.e. of
  theory, Var(Σx) = n^{2H} within 2 s.e., and the full covariance matrix at n = 16.
- **HSIC**: Gaussian kernel, median-heuristic bandwidth per series, V-statistic n·HSIC. **dCov**: V-statistic n·dCov².
- **Wild bootstrap**: Rademacher block-Markov multipliers (flip w.p. (1−e^{−1/l})/2, so corr = e^{−k/l}),
  l = ⌊n^{1/2}⌋, empirically centred, 500 draws. This matches the mmd-dml worker's chain.
- **Shift test**: all circular shifts c ∈ [⌊0.1n⌋, n − ⌊0.1n⌋], computed exactly by FFT.
- **Subsampling** (Betken–Dehling): l = ⌊n^{1/2}⌋, d_n = ⌊0.1n⌋, T_l = l·V-stat, reject if T_n > 95% quantile.
- **Newey–West (L-ind)**: Bartlett weights **1 − j/S for j < S** (Andrews' k(j/S_T) convention), Andrews (1991) AR(1)
  plug-in S = 1.1447(α̂T)^{1/3}, **ρ̂ not clipped**, S capped at T − 1, N(0,1) critical value.
  *This differs from the mmd-dml arm*, which clips ρ̂ to ±0.97 and uses weights 1 − j/(S+1). Compare NW numbers across
  arms with that in mind.
- **EWC**: ν = ⌊0.4T^{2/3}⌋ cosine projections, t_ν critical value.

## Verdicts

**P1 (cross rule): not supported for any kernel test at n ≤ 2000. Supported for Pearson-EWC.**

| test | sum < 1.5: excess ≤ 2 pp at n_max | sum > 1.5: excess grows with n | P1 | knot (95% CI) |
|---|---|---|---|---|
| HSIC wild bootstrap | no | no | **not supported** | 1.545 [1.400, 1.590] |
| HSIC shift | yes | no | **not supported** | 1.525 [1.100, 1.600] |
| HSIC subsampling | no | no | **not supported** | 1.530 [1.420, 1.575] |
| dCov wild bootstrap | no | no | **not supported** | 1.545 [1.375, 1.595] |
| dCov subsampling | no | no | **not supported** | 1.475 [1.365, 1.590] |
| Pearson NW-HAC | no | yes | **not supported** | 1.475 [1.440, 1.535] |
| Pearson EWC | yes | yes | **supported** | 1.535 [1.465, 1.555] |

What the pattern says:

- **The line is where theory puts it.** For all seven tests, the knot in size excess against H_X + H_Y has a CI that
  contains 1.5, and the slope change is positive with a CI that excludes 0. Size worsens sharply beyond H_X + H_Y ≈ 1.5
  (H ≈ 0.75 per series). The x-grid is coarse (1.0, 1.2, 1.4, 1.45, 1.5, 1.6, 1.8), so "contains 1.5" means the knot
  is located between about 1.4 and 1.6. It is not pinned to 1.5.
- **The kernel tests fail P1 on the "grows with n" half.** At n ≤ 2000, size beyond the line is high but mostly flat in
  n. Only the H = 0.9 cells grow significantly: subsampling by +7 to +10 pp, and wild bootstrap by +3 to +4 pp in G2
  only. The n = 4000 cells, which would show growth best, were dropped. P1 fails here as stated. That is not evidence
  that the growth mechanism is absent.
- **Wild bootstrap fails the "≤ 2 pp below the line" half even without memory.** At H = 0.5 (i.i.d.), size is 5.6–6.6%
  at n = 2000 and about 8% at n = 500, falling with n. This is a finite-sample oversize of the wild bootstrap with
  l = ⌊√n⌋, unrelated to memory.
- **Subsampling fails it near the line.** It gives 7.4–11.1% at H = 0.7 and 11.2–13.6% for G4 (0.6, 0.85, sum 1.45).
- **Pearson NW fails it at H = 0.7 and for G4 (0.6, 0.85).** It gives 7.1–9.2% at T = 8000. EWC stays within 1.4 pp of
  nominal for every sum < 1.5 cell, and its size grows with T in every sum > 1.5 cell. That is the cross rule exactly.

**KC1 (no bite): not triggered.** Only 12 of 36 wild-bootstrap G1–G3 cells at n = 2000 lie in [3.5%, 6.5%]. At H = 0.9,
HSIC-WB is 13.5–18.3% and dCov-WB 13.5–18.2%. Part of the out-of-band count comes from the i.i.d.-level oversize noted
above, but the H ≥ 0.8 cells are clearly outside the band for memory reasons. Memory bites for the wild bootstrap.

**KC2 (wrong line): not triggered** for any K-ind calibration (headline: HSIC wild bootstrap). Every knot CI contains
1.5. Mean G1–G3 excess at H = 0.6 is 0.2–1.9 pp, below 5 pp, and every slope-change CI excludes 0.

**Which calibration holds best.**
- Kernel: the **shift test**. Its size is 3.9–7.6% everywhere, and 3.9–7.2% at n = 2000, even at H = 0.9. Its knot is
  weak (slope change 6 pp/unit) and its size does not grow with n by 2000.
- Linear: **EWC** holds below the line, then breaks with T beyond it (26% at H = 0.9, T = 8000).
- Worst: **subsampling**, which is the only calibration with a proof under LRD (Betken–Dehling). It reaches 34–52% at
  H = 0.9 and 10–25% at H = 0.75–0.8 (n = 2000), and it grows with n at H = 0.9.

**R1 (rough-vol arm, H = 0.1 fOU).** λ = 0.1 behaves like a short-memory series: sizes 5–7% at n = 2000.
λ = 0.01 (a near-unit-root AR) mimics long memory at small n:
- subsampling 43–46% at n = 2000 (52–56% at 500);
- NW 17.5% at T = 8000;
- wild bootstrap 11.5–11.7%;
- EWC 8.2%;
- shift 6.4%.

All of these fall with n, which is the opposite of the LRD cells beyond the line. In this arm, shrinking with n
separates spurious long memory from the real thing.

**Power check** (G2, H = 0.7, r = 0.25): all seven tests reject 99–100% at every n. The tests are not trivially
conservative. At this r the check cannot separate them.

## Deviations and things not done

1. K-ind n = 4000 not run (author decision D-a28e4cb5fa). Verdicts are at n = 2000 for K-ind.
2. L-ind used the design's T grid {500, 2000, 8000}, not the {500…4000} in the coordinator's message.
3. Where the design was silent, I chose these and fixed them before running: the wild-bootstrap l_n = ⌊√n⌋; the
   shift range [0.1n, 0.9n]; full-sample bandwidths in the HSIC subsamples; a subsampling rejection rule of
   T_n > 95% quantile.
4. Process: two avoidable delays. About 2 h passed between the pilot and the main run, and about 3.5 h between run
   completion and analysis, because I did not act on completion notices. No effect on results.
5. Pilot timings were taken under heavy machine load (load ~30 on 22 cores). They overstate uncontended cost.
6. Not done: power at other H or smaller r; any n above 2000 for kernel tests.

## Files

`mbcore.py` (generators, tests), `validate.py` → `validation.out`, `run_sim.py` (driver), `analyze.py` (analysis),
`cells/` (one CSV per cell, `manifest.json`, `cell_summary.csv`, `knot_analysis.json`), `pilot/` (timing only),
`results_tables.md` (generated tables, appended below), `events.jsonl`.

---

## Cell tables: empirical size at 5%, % [Wilson 95% MC interval]

Null cells (r = 0). Every cell that was run is listed. 'err' = replications that raised an error. MC s.e. at 5% nominal: 0.49 pp at 2,000 reps, 0.69 pp at 1,000.

Dropped by author decision (not run): 25 cells: kind__G1_H0.5__n4000, kind__G1_H0.6__n4000, kind__G1_H0.7__n4000, kind__G1_H0.75__n4000, kind__G1_H0.8__n4000, kind__G1_H0.9__n4000, kind__G2_H0.5__n4000, kind__G2_H0.6__n4000, kind__G2_H0.7__n4000, kind__G2_H0.75__n4000, kind__G2_H0.8__n4000, kind__G2_H0.9__n4000, kind__G3_H0.5__n4000, kind__G3_H0.6__n4000, kind__G3_H0.7__n4000, kind__G3_H0.75__n4000, kind__G3_H0.8__n4000, kind__G3_H0.9__n4000, kind__G4_HX0.6_HY0.85__n4000, kind__G4_HX0.7_HY0.8__n4000, kind__G4_HX0.7_HY0.9__n4000, kind__G4_HX0.6_HY0.9__n4000, kind__R1_lam0.01__n4000, kind__R1_lam0.1__n4000, kind__G2_H0.7__n4000__r0.25


### hsic_wb

| setting | H_X+H_Y | n=500 | n=1000 | n=2000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 8.1 [7.0, 9.4] | 7.5 [6.5, 8.8] | 6.5 [5.5, 7.6] | 2000 | 0 |
| G1_H0.6 | 1.20 | 8.5 [7.3, 9.8] | 6.2 [5.2, 7.3] | 7.0 [5.9, 8.1] | 2000 | 0 |
| G1_H0.7 | 1.40 | 8.8 [7.6, 10.1] | 8.0 [6.9, 9.3] | 6.7 [5.7, 7.9] | 2000 | 0 |
| G1_H0.75 | 1.50 | 10.7 [9.4, 12.1] | 9.2 [8.0, 10.5] | 7.1 [6.1, 8.3] | 2000 | 0 |
| G1_H0.8 | 1.60 | 11.0 [9.7, 12.4] | 10.0 [8.7, 11.3] | 9.1 [7.9, 10.4] | 2000 | 0 |
| G1_H0.9 | 1.80 | 15.7 [14.1, 17.3] | 15.6 [14.1, 17.3] | 16.2 [14.6, 17.8] | 2000 | 0 |
| G2_H0.5 | 1.00 | 8.2 [7.1, 9.5] | 6.3 [5.3, 7.5] | 6.4 [5.4, 7.6] | 2000 | 0 |
| G2_H0.6 | 1.20 | 7.5 [6.4, 8.7] | 6.9 [5.9, 8.1] | 5.6 [4.7, 6.7] | 2000 | 0 |
| G2_H0.7 | 1.40 | 8.3 [7.2, 9.6] | 7.0 [6.0, 8.2] | 7.0 [6.0, 8.2] | 2000 | 0 |
| G2_H0.75 | 1.50 | 9.8 [8.5, 11.1] | 8.8 [7.7, 10.2] | 10.1 [8.9, 11.5] | 2000 | 0 |
| G2_H0.8 | 1.60 | 11.2 [9.8, 12.6] | 9.8 [8.6, 11.2] | 10.0 [8.7, 11.3] | 2000 | 0 |
| G2_H0.9 | 1.80 | 14.4 [12.9, 16.0] | 14.8 [13.3, 16.4] | 18.3 [16.7, 20.1] | 2000 | 0 |
| G3_H0.5 | 1.00 | 7.8 [6.7, 9.1] | 7.0 [5.9, 8.1] | 5.6 [4.7, 6.7] | 2000 | 0 |
| G3_H0.6 | 1.20 | 6.7 [5.6, 7.8] | 7.1 [6.1, 8.3] | 8.0 [6.9, 9.3] | 2000 | 0 |
| G3_H0.7 | 1.40 | 8.2 [7.1, 9.5] | 7.7 [6.6, 9.0] | 5.2 [4.4, 6.3] | 2000 | 0 |
| G3_H0.75 | 1.50 | 8.4 [7.3, 9.7] | 7.1 [6.1, 8.3] | 6.2 [5.2, 7.3] | 2000 | 0 |
| G3_H0.8 | 1.60 | 9.1 [7.9, 10.4] | 9.4 [8.2, 10.8] | 8.2 [7.0, 9.4] | 2000 | 0 |
| G3_H0.9 | 1.80 | 12.2 [10.9, 13.8] | 12.0 [10.6, 13.5] | 13.5 [12.1, 15.1] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 9.3 [8.2, 10.7] | 8.0 [6.9, 9.3] | 7.4 [6.4, 8.7] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 10.1 [8.9, 11.5] | 9.2 [8.0, 10.5] | 8.0 [6.9, 9.3] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 10.8 [9.6, 12.3] | 9.4 [8.2, 10.8] | 10.2 [8.9, 11.6] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 10.5 [9.3, 12.0] | 7.8 [6.7, 9.0] | 7.0 [6.0, 8.3] | 2000 | 0 |
| R1_lam0.01 | — | 16.1 [14.6, 17.8] | 13.2 [11.7, 14.7] | 11.5 [10.1, 12.9] | 2000 | 0 |
| R1_lam0.1 | — | 9.0 [7.8, 10.3] | 7.2 [6.2, 8.5] | 6.8 [5.7, 7.9] | 2000 | 0 |

### hsic_shift

| setting | H_X+H_Y | n=500 | n=1000 | n=2000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 4.8 [3.9, 5.8] | 5.4 [4.5, 6.5] | 5.0 [4.1, 6.0] | 2000 | 0 |
| G1_H0.6 | 1.20 | 5.7 [4.7, 6.7] | 4.0 [3.3, 5.0] | 5.1 [4.3, 6.2] | 2000 | 0 |
| G1_H0.7 | 1.40 | 5.9 [4.9, 7.0] | 5.5 [4.6, 6.6] | 5.0 [4.1, 6.0] | 2000 | 0 |
| G1_H0.75 | 1.50 | 6.3 [5.3, 7.5] | 6.3 [5.3, 7.5] | 4.7 [3.8, 5.7] | 2000 | 0 |
| G1_H0.8 | 1.60 | 4.9 [4.0, 5.9] | 6.1 [5.1, 7.2] | 5.8 [4.9, 6.9] | 2000 | 0 |
| G1_H0.9 | 1.80 | 7.6 [6.6, 8.9] | 6.8 [5.7, 7.9] | 7.1 [6.1, 8.3] | 2000 | 0 |
| G2_H0.5 | 1.00 | 5.5 [4.5, 6.5] | 4.7 [3.9, 5.7] | 4.9 [4.0, 5.9] | 2000 | 0 |
| G2_H0.6 | 1.20 | 4.5 [3.7, 5.5] | 4.7 [3.8, 5.7] | 4.4 [3.6, 5.4] | 2000 | 0 |
| G2_H0.7 | 1.40 | 5.3 [4.4, 6.4] | 5.1 [4.2, 6.1] | 5.2 [4.3, 6.3] | 2000 | 0 |
| G2_H0.75 | 1.50 | 5.4 [4.5, 6.5] | 5.5 [4.6, 6.6] | 6.5 [5.5, 7.6] | 2000 | 0 |
| G2_H0.8 | 1.60 | 6.0 [5.1, 7.2] | 5.1 [4.2, 6.1] | 5.9 [4.9, 7.0] | 2000 | 0 |
| G2_H0.9 | 1.80 | 6.0 [5.1, 7.2] | 6.8 [5.7, 7.9] | 7.2 [6.1, 8.4] | 2000 | 0 |
| G3_H0.5 | 1.00 | 5.4 [4.5, 6.5] | 5.2 [4.4, 6.3] | 4.0 [3.2, 5.0] | 2000 | 0 |
| G3_H0.6 | 1.20 | 5.1 [4.2, 6.1] | 4.7 [3.9, 5.7] | 6.2 [5.2, 7.3] | 2000 | 0 |
| G3_H0.7 | 1.40 | 5.5 [4.5, 6.5] | 5.7 [4.7, 6.7] | 3.9 [3.1, 4.8] | 2000 | 0 |
| G3_H0.75 | 1.50 | 4.9 [4.0, 5.9] | 5.5 [4.5, 6.5] | 4.5 [3.6, 5.4] | 2000 | 0 |
| G3_H0.8 | 1.60 | 5.3 [4.4, 6.4] | 5.5 [4.6, 6.6] | 4.8 [3.9, 5.8] | 2000 | 0 |
| G3_H0.9 | 1.80 | 5.8 [4.8, 6.9] | 6.3 [5.3, 7.5] | 6.3 [5.3, 7.5] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 5.7 [4.8, 6.8] | 5.7 [4.8, 6.8] | 5.3 [4.4, 6.4] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 5.5 [4.6, 6.6] | 5.5 [4.6, 6.6] | 5.1 [4.2, 6.2] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 6.5 [5.5, 7.7] | 5.3 [4.4, 6.4] | 5.7 [4.8, 6.8] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 5.7 [4.7, 6.7] | 5.0 [4.1, 6.0] | 4.7 [3.8, 5.7] | 2000 | 0 |
| R1_lam0.01 | — | 7.3 [6.3, 8.6] | 6.2 [5.2, 7.3] | 6.4 [5.4, 7.6] | 2000 | 0 |
| R1_lam0.1 | — | 5.2 [4.3, 6.3] | 4.9 [4.0, 5.9] | 5.1 [4.2, 6.1] | 2000 | 0 |

### hsic_sub

| setting | H_X+H_Y | n=500 | n=1000 | n=2000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 6.6 [5.5, 7.7] | 6.0 [5.0, 7.1] | 5.5 [4.6, 6.6] | 2000 | 0 |
| G1_H0.6 | 1.20 | 8.1 [7.0, 9.4] | 5.9 [4.9, 7.0] | 6.3 [5.4, 7.5] | 2000 | 0 |
| G1_H0.7 | 1.40 | 10.0 [8.8, 11.4] | 9.8 [8.6, 11.2] | 8.7 [7.5, 10.0] | 2000 | 0 |
| G1_H0.75 | 1.50 | 14.9 [13.5, 16.6] | 14.5 [13.0, 16.1] | 11.8 [10.4, 13.2] | 2000 | 0 |
| G1_H0.8 | 1.60 | 18.9 [17.2, 20.7] | 20.2 [18.5, 22.0] | 20.6 [18.9, 22.4] | 2000 | 0 |
| G1_H0.9 | 1.80 | 39.4 [37.2, 41.5] | 44.4 [42.2, 46.6] | 48.6 [46.4, 50.8] | 2000 | 0 |
| G2_H0.5 | 1.00 | 6.9 [5.8, 8.0] | 6.6 [5.6, 7.8] | 6.1 [5.1, 7.2] | 2000 | 0 |
| G2_H0.6 | 1.20 | 7.4 [6.3, 8.6] | 6.8 [5.7, 7.9] | 6.0 [5.0, 7.1] | 2000 | 0 |
| G2_H0.7 | 1.40 | 11.6 [10.2, 13.0] | 10.3 [9.1, 11.8] | 10.2 [8.9, 11.6] | 2000 | 0 |
| G2_H0.75 | 1.50 | 15.0 [13.5, 16.6] | 15.2 [13.7, 16.8] | 17.1 [15.5, 18.8] | 2000 | 0 |
| G2_H0.8 | 1.60 | 22.1 [20.3, 23.9] | 21.1 [19.4, 23.0] | 23.4 [21.6, 25.4] | 2000 | 0 |
| G2_H0.9 | 1.80 | 39.7 [37.6, 41.9] | 44.4 [42.2, 46.5] | 49.9 [47.7, 52.1] | 2000 | 0 |
| G3_H0.5 | 1.00 | 7.4 [6.4, 8.7] | 6.5 [5.5, 7.6] | 5.0 [4.1, 6.0] | 2000 | 0 |
| G3_H0.6 | 1.20 | 7.0 [6.0, 8.3] | 7.0 [6.0, 8.2] | 7.8 [6.7, 9.1] | 2000 | 0 |
| G3_H0.7 | 1.40 | 10.7 [9.4, 12.1] | 9.1 [7.9, 10.4] | 7.4 [6.4, 8.7] | 2000 | 0 |
| G3_H0.75 | 1.50 | 13.1 [11.7, 14.6] | 11.8 [10.4, 13.2] | 10.5 [9.3, 12.0] | 2000 | 0 |
| G3_H0.8 | 1.60 | 16.0 [14.4, 17.6] | 17.8 [16.1, 19.5] | 15.8 [14.2, 17.4] | 2000 | 0 |
| G3_H0.9 | 1.80 | 27.3 [25.4, 29.3] | 28.7 [26.8, 30.8] | 34.1 [32.1, 36.2] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 12.6 [11.2, 14.1] | 11.7 [10.3, 13.1] | 11.2 [9.9, 12.7] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 14.1 [12.7, 15.7] | 13.4 [11.9, 14.9] | 12.8 [11.4, 14.3] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 20.6 [18.9, 22.5] | 18.6 [17.0, 20.4] | 22.1 [20.3, 24.0] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 15.3 [13.8, 16.9] | 13.5 [12.1, 15.1] | 11.6 [10.3, 13.1] | 2000 | 0 |
| R1_lam0.01 | — | 51.8 [49.7, 54.0] | 49.2 [47.1, 51.4] | 42.8 [40.6, 45.0] | 2000 | 0 |
| R1_lam0.1 | — | 9.0 [7.8, 10.3] | 7.0 [6.0, 8.3] | 6.1 [5.1, 7.2] | 2000 | 0 |

### dcov_wb

| setting | H_X+H_Y | n=500 | n=1000 | n=2000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 8.1 [6.9, 9.3] | 8.2 [7.0, 9.4] | 6.6 [5.5, 7.7] | 2000 | 0 |
| G1_H0.6 | 1.20 | 8.7 [7.5, 10.0] | 6.7 [5.6, 7.8] | 6.6 [5.6, 7.8] | 2000 | 0 |
| G1_H0.7 | 1.40 | 9.7 [8.4, 11.0] | 8.0 [6.8, 9.2] | 6.4 [5.4, 7.6] | 2000 | 0 |
| G1_H0.75 | 1.50 | 11.0 [9.7, 12.4] | 9.3 [8.2, 10.7] | 8.2 [7.0, 9.4] | 2000 | 0 |
| G1_H0.8 | 1.60 | 11.3 [10.0, 12.8] | 10.9 [9.6, 12.3] | 9.7 [8.4, 11.0] | 2000 | 0 |
| G1_H0.9 | 1.80 | 17.0 [15.4, 18.7] | 15.4 [13.9, 17.1] | 17.2 [15.6, 18.9] | 2000 | 0 |
| G2_H0.5 | 1.00 | 8.1 [6.9, 9.3] | 6.9 [5.8, 8.0] | 6.5 [5.5, 7.7] | 2000 | 0 |
| G2_H0.6 | 1.20 | 7.8 [6.7, 9.1] | 6.9 [5.9, 8.1] | 6.0 [5.0, 7.1] | 2000 | 0 |
| G2_H0.7 | 1.40 | 9.0 [7.8, 10.3] | 7.2 [6.1, 8.4] | 7.1 [6.1, 8.3] | 2000 | 0 |
| G2_H0.75 | 1.50 | 10.4 [9.1, 11.8] | 9.0 [7.8, 10.3] | 10.2 [8.9, 11.6] | 2000 | 0 |
| G2_H0.8 | 1.60 | 11.8 [10.5, 13.3] | 10.2 [8.9, 11.6] | 10.4 [9.2, 11.9] | 2000 | 0 |
| G2_H0.9 | 1.80 | 14.8 [13.4, 16.5] | 15.2 [13.6, 16.8] | 18.2 [16.6, 20.0] | 2000 | 0 |
| G3_H0.5 | 1.00 | 7.8 [6.7, 9.1] | 7.0 [5.9, 8.1] | 5.6 [4.7, 6.7] | 2000 | 0 |
| G3_H0.6 | 1.20 | 6.7 [5.6, 7.8] | 7.1 [6.1, 8.3] | 8.0 [6.9, 9.3] | 2000 | 0 |
| G3_H0.7 | 1.40 | 8.2 [7.1, 9.5] | 7.7 [6.6, 9.0] | 5.2 [4.4, 6.3] | 2000 | 0 |
| G3_H0.75 | 1.50 | 8.4 [7.3, 9.7] | 7.1 [6.1, 8.3] | 6.2 [5.2, 7.3] | 2000 | 0 |
| G3_H0.8 | 1.60 | 9.1 [7.9, 10.4] | 9.4 [8.2, 10.8] | 8.2 [7.0, 9.4] | 2000 | 0 |
| G3_H0.9 | 1.80 | 12.2 [10.9, 13.8] | 12.0 [10.6, 13.5] | 13.5 [12.1, 15.1] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 9.8 [8.5, 11.1] | 9.2 [8.0, 10.5] | 8.2 [7.1, 9.5] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 10.1 [8.8, 11.4] | 8.9 [7.7, 10.2] | 8.9 [7.7, 10.2] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 11.5 [10.2, 13.0] | 10.3 [9.1, 11.8] | 10.4 [9.2, 11.9] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 9.8 [8.6, 11.2] | 8.3 [7.2, 9.6] | 6.9 [5.9, 8.1] | 2000 | 0 |
| R1_lam0.01 | — | 16.5 [14.9, 18.2] | 13.7 [12.2, 15.2] | 11.7 [10.4, 13.2] | 2000 | 0 |
| R1_lam0.1 | — | 8.8 [7.6, 10.1] | 7.6 [6.5, 8.8] | 6.7 [5.7, 7.9] | 2000 | 0 |

### dcov_sub

| setting | H_X+H_Y | n=500 | n=1000 | n=2000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 5.5 [4.6, 6.6] | 5.9 [4.9, 7.0] | 5.2 [4.4, 6.3] | 2000 | 0 |
| G1_H0.6 | 1.20 | 7.6 [6.5, 8.8] | 5.6 [4.7, 6.7] | 6.2 [5.3, 7.4] | 2000 | 0 |
| G1_H0.7 | 1.40 | 11.4 [10.1, 12.9] | 10.7 [9.4, 12.1] | 9.5 [8.3, 10.9] | 2000 | 0 |
| G1_H0.75 | 1.50 | 17.2 [15.7, 19.0] | 17.0 [15.4, 18.7] | 15.0 [13.5, 16.6] | 2000 | 0 |
| G1_H0.8 | 1.60 | 21.6 [19.9, 23.5] | 23.8 [22.0, 25.8] | 24.2 [22.4, 26.1] | 2000 | 0 |
| G1_H0.9 | 1.80 | 43.1 [41.0, 45.3] | 49.6 [47.5, 51.8] | 52.2 [50.1, 54.4] | 2000 | 0 |
| G2_H0.5 | 1.00 | 6.7 [5.7, 7.9] | 7.0 [6.0, 8.2] | 5.9 [5.0, 7.1] | 2000 | 0 |
| G2_H0.6 | 1.20 | 7.9 [6.8, 9.2] | 7.5 [6.4, 8.7] | 6.2 [5.2, 7.3] | 2000 | 0 |
| G2_H0.7 | 1.40 | 13.0 [11.5, 14.5] | 11.7 [10.3, 13.1] | 11.1 [9.7, 12.5] | 2000 | 0 |
| G2_H0.75 | 1.50 | 17.6 [16.0, 19.4] | 16.7 [15.1, 18.4] | 18.4 [16.8, 20.2] | 2000 | 0 |
| G2_H0.8 | 1.60 | 24.2 [22.4, 26.1] | 23.6 [21.8, 25.5] | 25.4 [23.6, 27.4] | 2000 | 0 |
| G2_H0.9 | 1.80 | 41.3 [39.2, 43.5] | 46.0 [43.8, 48.2] | 51.6 [49.5, 53.8] | 2000 | 0 |
| G3_H0.5 | 1.00 | 7.4 [6.4, 8.7] | 6.5 [5.5, 7.6] | 5.0 [4.1, 6.0] | 2000 | 0 |
| G3_H0.6 | 1.20 | 7.0 [6.0, 8.3] | 7.0 [6.0, 8.2] | 7.8 [6.7, 9.1] | 2000 | 0 |
| G3_H0.7 | 1.40 | 10.7 [9.4, 12.1] | 9.1 [7.9, 10.4] | 7.4 [6.4, 8.7] | 2000 | 0 |
| G3_H0.75 | 1.50 | 13.1 [11.7, 14.6] | 11.8 [10.4, 13.2] | 10.5 [9.3, 12.0] | 2000 | 0 |
| G3_H0.8 | 1.60 | 16.0 [14.4, 17.6] | 17.8 [16.1, 19.5] | 15.8 [14.2, 17.4] | 2000 | 0 |
| G3_H0.9 | 1.80 | 27.3 [25.4, 29.3] | 28.7 [26.8, 30.8] | 34.1 [32.1, 36.2] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 13.5 [12.0, 15.0] | 13.7 [12.2, 15.2] | 13.6 [12.2, 15.2] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 15.8 [14.3, 17.5] | 15.7 [14.2, 17.4] | 14.8 [13.3, 16.4] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 23.5 [21.7, 25.4] | 23.2 [21.4, 25.0] | 25.2 [23.3, 27.1] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 17.1 [15.5, 18.8] | 14.8 [13.3, 16.4] | 13.6 [12.2, 15.2] | 2000 | 0 |
| R1_lam0.01 | — | 56.1 [54.0, 58.3] | 53.9 [51.8, 56.1] | 45.9 [43.7, 48.1] | 2000 | 0 |
| R1_lam0.1 | — | 9.1 [7.9, 10.4] | 6.7 [5.6, 7.8] | 5.7 [4.7, 6.7] | 2000 | 0 |

### pearson_nw

| setting | H_X+H_Y | n=500 | n=2000 | n=8000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 4.5 [3.7, 5.5] | 5.8 [4.9, 6.9] | 4.2 [3.4, 5.2] | 2000 | 0 |
| G1_H0.6 | 1.20 | 5.3 [4.4, 6.4] | 5.7 [4.8, 6.8] | 5.6 [4.7, 6.7] | 2000 | 0 |
| G1_H0.7 | 1.40 | 8.1 [7.0, 9.4] | 8.7 [7.5, 10.0] | 7.1 [6.1, 8.3] | 2000 | 0 |
| G1_H0.75 | 1.50 | 9.4 [8.2, 10.8] | 13.0 [11.5, 14.5] | 11.6 [10.3, 13.1] | 2000 | 0 |
| G1_H0.8 | 1.60 | 12.4 [11.0, 13.9] | 14.6 [13.1, 16.2] | 18.4 [16.7, 20.1] | 2000 | 0 |
| G1_H0.9 | 1.80 | 17.7 [16.1, 19.4] | 26.6 [24.7, 28.6] | 34.2 [32.2, 36.3] | 2000 | 0 |
| G2_H0.5 | 1.00 | 5.0 [4.1, 6.0] | 5.2 [4.3, 6.3] | 4.9 [4.0, 5.9] | 2000 | 0 |
| G2_H0.6 | 1.20 | 6.0 [5.1, 7.2] | 5.2 [4.3, 6.3] | 5.2 [4.4, 6.3] | 2000 | 0 |
| G2_H0.7 | 1.40 | 8.6 [7.5, 10.0] | 7.5 [6.4, 8.7] | 7.9 [6.8, 9.2] | 2000 | 0 |
| G2_H0.75 | 1.50 | 10.5 [9.2, 11.9] | 11.1 [9.7, 12.5] | 11.7 [10.3, 13.1] | 2000 | 0 |
| G2_H0.8 | 1.60 | 14.4 [13.0, 16.1] | 14.8 [13.3, 16.4] | 18.2 [16.6, 20.0] | 2000 | 0 |
| G2_H0.9 | 1.80 | 20.5 [18.8, 22.4] | 26.0 [24.1, 28.0] | 34.2 [32.1, 36.3] | 2000 | 0 |
| G3_H0.5 | 1.00 | 5.3 [4.4, 6.4] | 5.3 [4.4, 6.4] | 5.1 [4.3, 6.2] | 2000 | 0 |
| G3_H0.6 | 1.20 | 6.2 [5.2, 7.3] | 5.0 [4.1, 6.0] | 5.4 [4.5, 6.5] | 2000 | 0 |
| G3_H0.7 | 1.40 | 7.2 [6.2, 8.5] | 6.3 [5.4, 7.5] | 7.3 [6.2, 8.5] | 2000 | 0 |
| G3_H0.75 | 1.50 | 8.3 [7.2, 9.6] | 9.7 [8.5, 11.1] | 8.3 [7.2, 9.6] | 2000 | 0 |
| G3_H0.8 | 1.60 | 10.0 [8.8, 11.4] | 12.8 [11.4, 14.3] | 15.4 [13.9, 17.1] | 2000 | 0 |
| G3_H0.9 | 1.80 | 16.4 [14.8, 18.0] | 21.4 [19.7, 23.3] | 33.1 [31.1, 35.2] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 8.5 [7.4, 9.8] | 8.2 [7.0, 9.4] | 9.2 [8.1, 10.6] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 10.1 [8.8, 11.4] | 10.8 [9.5, 12.2] | 11.8 [10.4, 13.2] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 13.0 [11.6, 14.5] | 14.2 [12.7, 15.8] | 18.0 [16.4, 19.7] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 9.2 [8.0, 10.5] | 9.0 [7.8, 10.3] | 9.2 [8.1, 10.6] | 2000 | 0 |
| R1_lam0.01 | — | 25.0 [23.2, 26.9] | 21.0 [19.3, 22.8] | 17.5 [15.9, 19.3] | 2000 | 0 |
| R1_lam0.1 | — | 8.0 [6.9, 9.3] | 6.3 [5.3, 7.5] | 6.4 [5.4, 7.6] | 2000 | 0 |

### pearson_ewc

| setting | H_X+H_Y | n=500 | n=2000 | n=8000 | reps (n_max) | err |
|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | 5.0 [4.1, 6.0] | 6.0 [5.1, 7.2] | 4.2 [3.5, 5.2] | 2000 | 0 |
| G1_H0.6 | 1.20 | 5.2 [4.4, 6.3] | 5.1 [4.2, 6.2] | 5.1 [4.3, 6.2] | 2000 | 0 |
| G1_H0.7 | 1.40 | 5.1 [4.2, 6.2] | 6.5 [5.5, 7.7] | 5.3 [4.4, 6.4] | 2000 | 0 |
| G1_H0.75 | 1.50 | 5.1 [4.2, 6.1] | 8.8 [7.6, 10.1] | 7.5 [6.5, 8.8] | 2000 | 0 |
| G1_H0.8 | 1.60 | 7.1 [6.1, 8.3] | 8.8 [7.6, 10.1] | 11.2 [9.9, 12.7] | 2000 | 0 |
| G1_H0.9 | 1.80 | 10.3 [9.1, 11.8] | 16.7 [15.1, 18.3] | 26.2 [24.3, 28.2] | 2000 | 0 |
| G2_H0.5 | 1.00 | 4.4 [3.6, 5.4] | 4.9 [4.0, 5.9] | 4.8 [3.9, 5.8] | 2000 | 0 |
| G2_H0.6 | 1.20 | 4.8 [3.9, 5.8] | 4.6 [3.8, 5.6] | 4.5 [3.6, 5.4] | 2000 | 0 |
| G2_H0.7 | 1.40 | 6.0 [5.0, 7.1] | 5.5 [4.6, 6.6] | 6.0 [5.0, 7.1] | 2000 | 0 |
| G2_H0.75 | 1.50 | 5.8 [4.9, 6.9] | 7.5 [6.4, 8.7] | 8.0 [6.8, 9.2] | 2000 | 0 |
| G2_H0.8 | 1.60 | 7.8 [6.8, 9.1] | 8.9 [7.8, 10.3] | 11.9 [10.6, 13.4] | 2000 | 0 |
| G2_H0.9 | 1.80 | 12.5 [11.1, 14.0] | 17.4 [15.8, 19.2] | 25.8 [23.9, 27.7] | 2000 | 0 |
| G3_H0.5 | 1.00 | 4.7 [3.9, 5.7] | 5.2 [4.3, 6.3] | 5.1 [4.2, 6.1] | 2000 | 0 |
| G3_H0.6 | 1.20 | 5.9 [4.9, 7.0] | 5.1 [4.3, 6.2] | 5.0 [4.1, 6.0] | 2000 | 0 |
| G3_H0.7 | 1.40 | 6.3 [5.3, 7.5] | 5.1 [4.2, 6.2] | 5.5 [4.6, 6.6] | 2000 | 0 |
| G3_H0.75 | 1.50 | 6.0 [5.0, 7.1] | 6.5 [5.5, 7.6] | 5.9 [5.0, 7.1] | 2000 | 0 |
| G3_H0.8 | 1.60 | 7.0 [6.0, 8.2] | 7.6 [6.5, 8.8] | 10.2 [9.0, 11.7] | 2000 | 0 |
| G3_H0.9 | 1.80 | 8.7 [7.5, 10.0] | 13.7 [12.3, 15.3] | 21.7 [19.9, 23.6] | 2000 | 0 |
| G4_HX0.6_HY0.85 | 1.45 | 6.5 [5.5, 7.6] | 5.7 [4.8, 6.8] | 6.4 [5.4, 7.6] | 2000 | 0 |
| G4_HX0.7_HY0.8 | 1.50 | 5.9 [4.9, 7.0] | 6.3 [5.4, 7.5] | 7.3 [6.2, 8.5] | 2000 | 0 |
| G4_HX0.7_HY0.9 | 1.60 | 7.0 [6.0, 8.2] | 9.3 [8.1, 10.7] | 11.0 [9.7, 12.4] | 2000 | 0 |
| G4_HX0.6_HY0.9 | 1.50 | 6.2 [5.3, 7.4] | 6.0 [5.0, 7.1] | 6.8 [5.8, 8.0] | 2000 | 0 |
| R1_lam0.01 | — | 12.3 [10.9, 13.8] | 9.9 [8.7, 11.3] | 8.2 [7.1, 9.5] | 2000 | 0 |
| R1_lam0.1 | — | 6.2 [5.2, 7.3] | 4.9 [4.0, 5.9] | 5.5 [4.5, 6.5] | 2000 | 0 |

## Power check: G2, H = 0.7, Betken–Dehling 'linear' r = 0.25 (rejection %, [95% MC])

| test | n=500 | n=1000 | n=2000 | n=8000 |
|---|---|---|---|---|
| hsic_wb | 99.3 [98.6, 99.7] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | — |
| hsic_shift | 99.1 [98.3, 99.5] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | — |
| hsic_sub | 99.5 [98.8, 99.8] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | — |
| dcov_wb | 99.5 [98.8, 99.8] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | — |
| dcov_sub | 99.6 [99.0, 99.8] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | — |
| pearson_nw | 99.8 [99.3, 99.9] | — | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] |
| pearson_ewc | 99.5 [98.8, 99.8] | — | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] |

## One-knot spline on excess at the largest n (x = H_X + H_Y, 22 G1–G4 cells)

| test | n_max | knot | knot 95% CI | slope change c | c 95% CI | CI contains 1.5 |
|---|---|---|---|---|---|---|
| hsic_wb | 2000 | 1.545 | [1.400, 1.590] | 30.8 pp/unit | [15.7, 44.5] | yes |
| hsic_shift | 2000 | 1.525 | [1.100, 1.600] | 6.0 pp/unit | [2.3, 13.2] | yes |
| hsic_sub | 2000 | 1.530 | [1.420, 1.575] | 104.4 pp/unit | [58.3, 143.2] | yes |
| dcov_wb | 2000 | 1.545 | [1.375, 1.595] | 30.1 pp/unit | [15.4, 43.8] | yes |
| dcov_sub | 2000 | 1.475 | [1.365, 1.590] | 92.6 pp/unit | [53.8, 146.4] | yes |
| pearson_nw | 8000 | 1.475 | [1.440, 1.535] | 69.8 pp/unit | [60.4, 78.1] | yes |
| pearson_ewc | 8000 | 1.535 | [1.465, 1.555] | 62.1 pp/unit | [38.0, 72.9] | yes |

## P1 rule per test

| test | all sum<1.5 cells: excess at n_max ≤ 2 pp | all sum>1.5 cells: excess grows with n | P1 supported |
|---|---|---|---|
| hsic_wb | False | False | **False** |
| hsic_shift | True | False | **False** |
| hsic_sub | False | False | **False** |
| dcov_wb | False | False | **False** |
| dcov_sub | False | False | **False** |
| pearson_nw | False | True | **False** |
| pearson_ewc | True | True | **True** |

Growth with n (size at n_max minus size at n_min, pp; * = significant by the fixed rule):

| setting | sum | hsic_wb | hsic_shift | hsic_sub | dcov_wb | dcov_sub | pearson_nw | pearson_ewc |
|---|---|---|---|---|---|---|---|---|
| G1_H0.5 | 1.00 | -1.7 | +0.2 | -1.0 | -1.5 | -0.3 | -0.3 | -0.7 |
| G1_H0.6 | 1.20 | -1.5 | -0.5 | -1.8 | -2.1 | -1.3 | +0.3 | -0.1 |
| G1_H0.7 | 1.40 | -2.0 | -0.9 | -1.3 | -3.2 | -1.9 | -1.0 | +0.2 |
| G1_H0.75 | 1.50 | -3.6 | -1.7 | -3.2 | -2.8 | -2.2 | +2.2* | +2.5* |
| G1_H0.8 | 1.60 | -1.9 | +1.0 | +1.7 | -1.7 | +2.6 | +5.9* | +4.1* |
| G1_H0.9 | 1.80 | +0.5 | -0.6 | +9.2* | +0.2 | +9.1* | +16.5* | +15.9* |
| G2_H0.5 | 1.00 | -1.8 | -0.5 | -0.8 | -1.6 | -0.8 | -0.1 | +0.4 |
| G2_H0.6 | 1.20 | -1.9 | -0.1 | -1.4 | -1.8 | -1.8 | -0.8 | -0.3 |
| G2_H0.7 | 1.40 | -1.3 | -0.1 | -1.4 | -1.9 | -1.9 | -0.7 | +0.0 |
| G2_H0.75 | 1.50 | +0.4 | +1.1 | +2.1 | -0.2 | +0.8 | +1.2 | +2.1* |
| G2_H0.8 | 1.60 | -1.2 | -0.2 | +1.4 | -1.3 | +1.3 | +3.8* | +4.1* |
| G2_H0.9 | 1.80 | +3.9* | +1.1 | +10.2* | +3.4* | +10.3* | +13.6* | +13.2* |
| G3_H0.5 | 1.00 | -2.2 | -1.4 | -2.5 | -2.2 | -2.5 | -0.2 | +0.4 |
| G3_H0.6 | 1.20 | +1.3 | +1.1 | +0.8 | +1.3 | +0.8 | -0.8 | -0.9 |
| G3_H0.7 | 1.40 | -3.0 | -1.6 | -3.2 | -3.0 | -3.2 | +0.1 | -0.8 |
| G3_H0.75 | 1.50 | -2.3 | -0.4 | -2.6 | -2.3 | -2.6 | +0.0 | -0.1 |
| G3_H0.8 | 1.60 | -0.9 | -0.5 | -0.2 | -0.9 | -0.2 | +5.4* | +3.2* |
| G3_H0.9 | 1.80 | +1.3 | +0.5 | +6.8* | +1.3 | +6.8* | +16.8* | +13.0* |
| G4_HX0.6_HY0.85 | 1.45 | -1.9 | -0.4 | -1.3 | -1.5 | +0.2 | +0.7 | -0.1 |
| G4_HX0.7_HY0.8 | 1.50 | -2.1 | -0.5 | -1.4 | -1.2 | -1.1 | +1.7 | +1.4 |
| G4_HX0.7_HY0.9 | 1.60 | -0.7 | -0.8 | +1.5 | -1.1 | +1.7 | +5.0* | +4.0* |
| G4_HX0.6_HY0.9 | 1.50 | -3.5 | -1.0 | -3.7 | -2.9 | -3.5 | +0.1 | +0.6 |

## Kill criteria

- KC1 (no bite): triggered = **False** (36 wild-bootstrap G1–G3 cells at n = 2000; in band [3.5, 6.5]%: 12).
- KC2 (hsic_wb): knot CI excludes 1.5 = False; mean G1–G3 excess at H = 0.6 = 1.8 pp; slope-change CI contains 0 = False; triggered = **False**
- KC2 (hsic_shift): knot CI excludes 1.5 = False; mean G1–G3 excess at H = 0.6 = 0.2 pp; slope-change CI contains 0 = False; triggered = **False**
- KC2 (hsic_sub): knot CI excludes 1.5 = False; mean G1–G3 excess at H = 0.6 = 1.7 pp; slope-change CI contains 0 = False; triggered = **False**
- KC2 (dcov_wb): knot CI excludes 1.5 = False; mean G1–G3 excess at H = 0.6 = 1.9 pp; slope-change CI contains 0 = False; triggered = **False**
- KC2 (dcov_sub): knot CI excludes 1.5 = False; mean G1–G3 excess at H = 0.6 = 1.7 pp; slope-change CI contains 0 = False; triggered = **False**

Total CPU in these cells: 28.00 CPU-hours. Replications with errors: 0.
