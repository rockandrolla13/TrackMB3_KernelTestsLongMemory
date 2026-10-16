# Proposed simulation — memory breakdown map (NOT RUN; for author approval)

Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 · plan item (5) · 2026-10-15.
Status: **proposal only.** Nothing here has been run. No data. Exploratory evidence class
(CLAUDE.md invariant 5). Approval would be a separate author decision and a separate job.

## Question it answers

Where does nominal 5% size (or 95% coverage) fail as memory grows, for three method families? Does
the failure point follow the Hermite-rank prediction in boundary-theorems.md §3? That prediction is:
H = 1/2 for rank-1 statistics, H = 3/4 for rank-2, H_X + H_Y = 3/2 for cross statistics, and no
H-line for Hawkes.

## Pre-registered predictions (what would be confirmed or refuted)

- **P1 (cross rule).** For HSIC/dCov independence tests and Pearson cross-correlation, size excess
  at the largest n stays ≤ 2 pp when H_X + H_Y < 1.5. It grows with n when H_X + H_Y > 1.5.
- **P2 (rank-1).** MMD homogeneity (first half vs second half of one series) and the population-ATE
  AIPW interval under monotone heterogeneity break for every H > 0.5. Excess grows with n.
- **P3 (rank-2).** The population-ATE interval under even heterogeneity (τ(x) = x² − 1) breaks
  only for H > 0.75.
- **P4 (C1).** Path-conditional ATE coverage stays nominal (≥ 93%) at all H.
- **P5 (Hawkes).** Granger-test size depends on branching ratio ρ, kernel tail and shared baseline.
  There is no threshold in an H-like index.

## Processes (fixed)

| ID | Process | Parameters |
|---|---|---|
| G1 | Gaussian fGn via Davies–Harte circulant embedding (exact) | H ∈ {0.5, 0.6, 0.7, 0.75, 0.8, 0.9} |
| G2 | Subordinated: X = 2Φ(ξ) − 1 (rank 1, bounded; as Betken–Dehling) | same H grid |
| G3 | Subordinated sign: X = sign(ξ) (rank 1; order-sign-like) | same H grid |
| G4 | Unequal pairs for the cross rule | (H_X, H_Y) ∈ {(0.6, 0.85), (0.7, 0.8), (0.7, 0.9), (0.6, 0.9)}; sums 1.45, 1.5, 1.6, 1.5 |
| R1 | Rough-volatility arm: log σ_t = fOU with H = 0.1, mean reversion λ ∈ {0.01, 0.1} per step, series = log σ_t | tests "spurious long memory" mimicry |
| HK1 | Bivariate Hawkes, exponential kernels, zero cross-excitation | ρ ∈ {0.5, 0.8, 0.95, 0.99} |
| HK2 | Bivariate Hawkes, power-law kernels φ(t) = ρ·α c^α/(t + c)^{1+α}, α = 0.5, zero cross-excitation | ρ as HK1 |
| HK3 | HK1 and HK2 with a shared slowly varying baseline μ(t) = μ₀(1 + 0.5 sin(2πt/T_b)), T_b = T/4 | ρ ∈ {0.5, 0.95} |

Independence nulls: X and Y drawn from independent generators. Seeds: fixed master seed, one stream per
cell (numpy SeedSequence spawn), recorded in the manifest.

## Statistics and calibrations (fixed)

| Family | Statistic | Calibrations compared |
|---|---|---|
| K-ind | HSIC, Gaussian kernels, median-heuristic bandwidth, lag 0 | (a) wild bootstrap (Chwialkowski et al. 2014, Rademacher block-Markov multiplier, 500 draws); (b) shift test (Chwialkowski–Gretton 2014); (c) subsampling (Betken–Dehling, l = ⌊n^{0.5}⌋, d_n = 0.1n) |
| K-ind | Distance covariance, lag 0 | (a), (c) |
| L-ind | Pearson cross-correlation | HAC Newey–West (Andrews bandwidth) and EWC fixed-b (Lazarus et al. 2018 rule) |
| K-hom | MMD², first half vs second half of one series (G1–G3, R1) | wild bootstrap; subsampling |
| DML | AIPW for ATE; state X_t from G1; D_t ~ Bernoulli(logit⁻¹(0.5 X_t)); Y_t = τ(X_t)D_t + X_t + ε_t, ε i.i.d. N(0,1) | effects: (i) τ = 1; (ii) τ = 1 + X_t (rank 1); (iii) τ = 1 + (X_t² − 1) (rank 2). Estimands: population ATE = 1 vs path-conditional (1/T)Σ τ(X_t). Nuisances: known propensity; and learned (logistic + gradient boosting) with 5-fold cross-fitting, as-i.i.d. folds vs blocked folds with buffer ⌊T^{1/2}⌋. Variance: HAC Newey–West and EWC fixed-b |
| HG | Likelihood-ratio test of zero cross-excitation, exponential-kernel model fitted by MLE (correct for HK1, misspecified for HK2) | χ²-boundary mixture (½χ²₀ + ½χ²₁ per cross term) |

## Sample sizes and replications

- Kernel families: n ∈ {500, 1000, 2000, 4000}. Linear and DML: T ∈ {500, 2000, 8000}.
  Hawkes: horizon chosen so expected events per stream ∈ {2,000, 10,000, 50,000}.
- Size and coverage: 2,000 replications per cell (Monte Carlo s.e. ≈ 0.49 pp at 5%, ≈ 0.49 pp at
  95% coverage). Fallback if the timing pilot exceeds the cap: 1,000 at the largest n (s.e. 0.69 pp).
- Power: not in scope for the first run (sizes first). One power check per family at H = 0.7
  (e.g. Betken–Dehling "linear" r = 0.25) to confirm the tests are not trivially conservative.

## Compute (estimate, to be checked by a timing pilot)

Rough single-core figures: HSIC with wild bootstrap at n = 4000 ≈ 5–10 s per replication. Kernel
arm ≈ 22 process settings × 2,000 reps at the largest n ≈ 60–120 CPU-hours. Smaller n and other
families are cheaper. Hawkes MLE at 50k events ≈ seconds per fit. **Proposed cap: 48 CPU-hours per
family.** A 1% timing pilot (20 reps per cell) is part of the approved run, and the run reduces
replications per the fallback rule if needed. The pilot is not evidence.

## Analysis (fixed before running)

- Primary outcome: empirical size (coverage) minus nominal, with Monte Carlo 95% intervals.
- Breakpoint test for P1–P3: per family, regress size excess at the largest n on H (or H_X + H_Y)
  with a one-knot linear spline. Estimate the knot by profile least squares with a bootstrap CI
  over cells. "Follows the predicted line" = knot CI contains the predicted value (0.5, 0.75 or 1.5),
  and excess grows with n only beyond it.
- Report every cell, including failures and non-convergent fits (invariant 4).

## Kill criteria (any one ends or reframes the direction)

- **KC1 — no bite.** For K-ind with wild bootstrap, size at n = 4000 lies within [3.5%, 6.5%] for
  every G1–G3 H ≤ 0.9. Then memory is a scope note, not a breakdown. Kill the "breakdown map"
  framing for kernel tests; Track E's G5 ranking stands.
- **KC2 — wrong line.** The K-ind knot CI excludes the cross-rule value, or excess is already
  > 5 pp at H = 0.6 with no slope change. Drop the H = 3/4 headline. Reframe as "calibration fails
  at any memory", if the data say so.
- **KC3 — C1 dead.** Population-ATE HAC coverage under (ii) at H = 0.9, T = 8000 is ≥ 93%. C1
  is killed.
- **KC4 — C1 mechanism wrong.** Path-conditional coverage < 93% at H = 0.9 with known propensity.
  The decomposition is false or incomplete; stop and diagnose before any write-up.
- **KC5 — C2 dead.** Hawkes LR size within ±1 pp of nominal across all ρ and baseline variants
  (as in the prior C2 proposal).
- **KC6 — pre-empted.** A paper mapping size across H for any two of the three families appears
  before the run. Stop and re-scope.

## What it does not do

No real data; no estimation of H on market series; no KCI/conditional or signature-kernel tests
(cost; add later if K-ind survives); no theory. Results are exploratory and cannot support an
accepted claim.

## Ownership

Overlaps Track E (kernel arm, G5), K and F-CAUSAL (DML arm), B2-COMP and the Track A ρ item
(Hawkes arm). Proposal: run as one shared simulator job under NEW-MEMORY-BOUNDARY and hand each
arm's result to its owning track. That ownership choice is the author's.
