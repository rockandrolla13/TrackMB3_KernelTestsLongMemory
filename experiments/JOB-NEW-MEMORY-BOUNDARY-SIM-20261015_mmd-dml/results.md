# Results: K-hom (MMD) and DML families

Job JOB-NEW-MEMORY-BOUNDARY-SIM-20261015 · worker `claude-mb-mmddml-20261015` · run 2026-10-15/16.
Design `runs/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/simulation-design.md`, SHA-256 `7bb1d6b9…559f`, verified.
**Evidence class: exploratory simulation.** It cannot support an accepted claim.

Every cell is listed in `main/tables.md`, with its Monte Carlo 95% (Wilson) interval. The per-cell
summaries are `main/summary_khom.csv` and `main/summary_dml.csv`. The raw per-replication rows are
`main/khom/*.csv` and `main/dml/*.csv`. Verdicts are in `main/analysis.json`.
Failures: 0 errored replications out of 255,000. There were 0 IRLS non-convergences and 0 clipped
propensities.

## Run facts

| | K-hom | DML |
|---|---|---|
| cells | 80 null + 4 power | 54 |
| reps per cell | 2,000 (n ≤ 2000); 1,000 at n = 4000 (fallback) | 2,000 |
| CPU-hours, main run | 18.86 | 6.15 |
| pilot CPU-hours (not evidence) | 0.52 | 0.12 |

There was also a first, discarded pilot (1.45 CPU-h, in `pilot_v1_blas_multithread/`). It
mis-timed because BLAS ran multithreaded. It was fixed and re-run before any main run; see the
parameters_amended event. Total CPU used is about 27.1 CPU-hours.

The re-run pilot projected 52.2 CPU-h for K-hom at 2,000 reps everywhere. That is over the
48 CPU-h cap, so the design's fallback applied: 1,000 reps at n = 4000. The main run then cost
less than projected. The pilot timings had been inflated by contention. In hindsight 2,000 reps
everywhere would have fitted. The rule was applied as written, before the run.

## K-hom: MMD² of the first half against the second half

Size minus 5% at n = 4000, in percentage points:

| process | H=0.5 | 0.6 | 0.7 | 0.75 | 0.8 | 0.9 |
|---|---|---|---|---|---|---|
| G1 wild bootstrap | +2.2 | +9.8 | +20.8 | +25.0 | +30.6 | +38.7 |
| G1 subsampling | −0.2 | +14.2 | +33.3 | +41.0 | +49.5 | +63.3 |
| G2 wild | +2.8 | +9.6 | +20.9 | +25.2 | +34.1 | +40.0 |
| G2 subsampling | +1.5 | +14.1 | +31.9 | +41.8 | +55.5 | +65.8 |
| G3 wild | +1.2 | +6.7 | +17.8 | +22.6 | +28.8 | +39.5 |
| G3 subsampling | 0.0 | +11.3 | +30.2 | +38.0 | +48.2 | +61.5 |

At n = 4000 the MC s.e. is about 0.7–1.6 pp. Each cell's interval is in `main/tables.md`.

- **It breaks at the first H above 0.5.** At H = 0.6 size is already 11.7–19.2% for every process and
  both calibrations. There is no sign of a 3/4 line.
- **Subsampling excess grows with n** for every H > 0.5 and every process (z from 2.1 to 9.5,
  n = 500 → 4000).
- **Wild bootstrap excess grows with n only from H ≈ 0.7.** At H = 0.6 it stays flat at about
  +7 to +10 pp over n = 500–4000 (z between −0.9 and 0.1).
- **Calibration defect under no memory.** At H = 0.5 (i.i.d.), the wild bootstrap is already
  oversized by +2 to +4.5 pp for G1 and G2 at every n. Under the 2 pp rule, that counts as
  "breaking" at H = 0.5. It comes from the multiplier bootstrap with l = ⌊√n⌋ on these MMD
  statistics, not from memory. Subsampling is on size at H = 0.5 (−0.5 to +1.5 pp).
- **Breakpoint (hockey stick, n = 4000).**

  | series | wild κ̂ [95% CI] | subsampling κ̂ [95% CI] |
  |---|---|---|
  | G1 | 0.514 [0.500, 0.543] | 0.503 [0.500, 0.523] |
  | G2 | 0.528 [0.500, 0.557] | 0.524 [0.505, 0.541] |
  | G3 | 0.550 [0.522, 0.574] | 0.526 [0.506, 0.544] |

  The CI contains 0.5 in 3 of 6 fits. In the other 3 it misses 0.5 by 0.005–0.022. The H grid
  has no point between 0.5 and 0.6, and the excess is concave in H. So a knot just above 0.5 is
  what a break at 0.5 looks like on this grid. Strictly, the pre-registered "CI contains 0.5" test
  passes 3 of 6.
- **R1 (rough-volatility mimic).** With λ = 0.01 it mimics strong memory at small n: wild
  +26.9 pp at n = 500 and subsampling +51.1 pp. The excess then fades with n: wild falls to
  −2.3 pp and subsampling to +14.5 pp at n = 4000. With λ = 0.1 the test is conservative
  (size about 0–2%).
- **Power check (G1, H = 0.7, shift +0.3).** Rejection rises with n: wild 0.55 → 0.89,
  subsampling 0.60 → 0.94. The tests are not trivially conservative. At this H, though, part of
  that rejection rate is the size excess.

## DML: AIPW interval for the ATE

Coverage at T = 8000 with the HAC Newey–West interval and the oracle (known-propensity) arm, in %:

| τ | target | H=0.5 | 0.6 | 0.7 | 0.75 | 0.8 | 0.9 |
|---|---|---|---|---|---|---|---|
| const | population = path | 95.5 | 94.5 | 95.5 | 95.2 | 95.0 | 94.8 |
| rank-1 (1+X) | population ATE | 95.2 | 84.0 | 54.0 | 40.9 | 29.8 | **11.9** |
| rank-1 | path-conditional | 97.4 | 97.6 | 98.6 | 98.4 | 99.2 | 99.2 |
| rank-2 (X²) | population ATE | 94.2 | 95.7 | 94.0 | 88.3 | 76.4 | 26.4 |
| rank-2 | path-conditional | 98.1 | 98.8 | 98.9 | 98.7 | 99.2 | 99.9 |

MC s.e. is about 0.5 pp near 95%. Every T, both variance estimators and all five arms are in
`main/tables.md`.

- **Nuisance learning and fold scheme do not matter here.** The oracle, the known propensity with
  GB outcomes, and the logistic propensity with GB outcomes agree within about 1–2 pp in every
  cell. So do i.i.d. and blocked folds. The failure is entirely the state term τ(X_t) − θ, as the
  decomposition predicts.
- **EWC fixed-b is somewhat less bad than NW, but it fails the same way.** For rank-1 at H = 0.9
  and T = 8000, coverage is 18.4% (EWC) against 11.9% (NW).
- **Why it fails.** For rank-1 at T = 8000, the ratio of the mean NW SE to the true SD of θ̂
  around 1 falls with H: 0.99, 0.72, 0.39, 0.28, 0.19, 0.08.
- **The population interval over-covers the path-conditional ATE (97–100%).** That is correct
  coverage for the path target, but the interval is wider than that target needs. It also carries
  the state variance, which is irrelevant to the path target.
- **Power check (τ const, H = 0.7, θ₀ = 1.1).** Rejection is 0.19 / 0.55 / 0.99 at T =
  500 / 2000 / 8000.

Breakpoints for the population ATE (coverage deficit at T = 8000, hockey stick):

| τ | arm, variance | κ̂ [95% CI] | prediction |
|---|---|---|---|
| rank-1 | oracle NW | 0.538 [0.529, 0.546] | 0.5 |
| rank-1 | oracle EWC | 0.562 [0.554, 0.570] | 0.5 |
| rank-2 | oracle NW | 0.767 [0.741, 0.772] | 0.75 |
| rank-2 | oracle EWC | 0.777 [0.772, 0.781] | 0.75 |

The other arms give the same values to ±0.002 (see `main/analysis.json`).

## Predictions and kill criteria

- **P2 (rank-1 breaks for every H > 0.5, and the excess grows with n): substantively confirmed; the
  strict knot test is mixed.**
  - MMD (both calibrations, all of G1–G3) and the rank-1 population-ATE interval (all arms, NW
    and EWC) are beyond 2 pp at every H ≥ 0.6.
  - For DML and for MMD subsampling, the excess grows significantly with n or T at every H > 0.5.
    For the MMD wild bootstrap it grows only from H ≥ 0.7. At H = 0.6 it is flat at about +8 pp.
  - Knot CIs contain 0.5 for 3 of 6 MMD fits and for none of the DML fits. Those fits put the knot
    at 0.53–0.56, between the grid points 0.5 and 0.6. Nothing in the data suggests a later line:
    DML rank-1 coverage is already 84% at H = 0.6 and falling with T.
- **P3 (rank-2 breaks only for H > 0.75): holds away from the boundary; at H = 0.75 there is
  moderate under-coverage.**
  - H ≤ 0.7: coverage is within about 1 pp of nominal (deficit −1.4 to +1.3 pp).
  - H = 0.8 and 0.9: coverage breaks and worsens with T (z 5–19).
  - H = 0.75 itself: deficit is 6.7 pp (NW) and 3.3–4.3 pp (EWC), with weak growth in T
    (z ≈ 1.2–2.6). Read literally, "only for H > 0.75" therefore fails at the boundary point.
    Theory puts a log-n variance inflation exactly at D = 1/2, which fits a slow drift rather than
    a clean pass.
  - The knot CI contains 0.75 only for the oracle NW fit. The others put κ at 0.767–0.777.
- **P4 (path-conditional coverage ≥ 93% at all H): confirmed.** The lowest cell across all 54
  cells, 5 arms and 2 variance estimators is 93.45% (τ const, H = 0.75, T = 2000; Wilson
  [92.3, 94.5]). That is an ordinary MC low near 95%. For heterogeneous τ, path coverage is
  96–100%.
- **KC3 (C1 killed if the rank-1 population-ATE NW coverage at H = 0.9, T = 8000 is ≥ 93%): not
  triggered.** Coverage is 11.95% [10.6, 13.4] (oracle). It is 12.0–12.1% in the other arms and
  18.4–18.6% with EWC. **C1 survives.**
- **KC4 (path-conditional coverage < 93% at H = 0.9 with known propensity): not triggered.** The
  minimum over τ and T is 94.4% (oracle NW), 94.85% (oracle EWC) and 94.45–94.75% for the
  known-propensity GB arms. No stop is required.

## Implementation notes the coordinator should know

- **Newey–West convention differs from kernel-linear.** This worker uses Andrews' AR(1)
  bandwidth S, ρ clipped to ±0.97, and weights 1 − j/(S+1) for j ≤ ⌊S⌋. The kernel-linear worker
  uses unclipped ρ, weights 1 − j/S for j < S, and S ≤ T − 1. EWC is the same in both.
- **Two small differences in the wild bootstrap, neither material.**
  - The multiplier chains match the kernel-linear worker's in distribution: corr e^{−k/l} with
    l = ⌊√n⌋.
  - The median heuristic differs only for G3. There, this worker drops zero distances and
    kernel-linear uses them unless the median is 0. Every calibration is scale-invariant for
    binary data, so this does not matter.
- **Subsampling for MMD** uses sampling windows of length ⌊√n⌋, with the window split into halves.
  The design's d_n (the lag between X and Y blocks) has no role in a one-series homogeneity test.
- **Interpretations of the design**, logged in parameters_fixed:
  - "Learned (logistic + gradient boosting)" was read as logistic for the propensity and gradient
    boosting for the outcome regressions.
  - "Known propensity" is reported both with true outcome regressions (oracle) and with GB
    outcomes.
  - One interval is scored against both estimands.
- **Gradient boosting** is a local numpy histogram GB: 255 bins, depth-2 trees, learning rate
  0.1, 100 rounds, minimum leaf 20. It replaces sklearn, which is not installed. Validation showed
  its error falls with sample size.
- **Validation:** `validation.out`, 0 failures. Davies–Harte autocovariances, Var(S_m)/m^{2H},
  kurtosis, t-quantiles, Φ, the learners and the multiplier were checked before use.

## Not done

- No path-conditional-specific (shorter) standard error was built. The design did not ask for one.
- No fGn H-estimation and no real data, as in the design.
