# EXP2 (2026-10-23) — process feature × memory factorial under known nulls

Written before any generator check, pilot or main run. Simulation only. No market data.
Exploratory: rates are false-alarm rates of the tested implementations under the stated
generators, not validity results for any rule.

## Question

Where do the kernel tests used in the paper break because of

1. memory alone (Gaussian fGn, H > 1/2),
2. a process feature alone (heavy tails, nonlinear dynamics, or seasonality, with H = 1/2), or
3. a feature combined with memory?

Each feature is switched on one at a time and crossed with memory H ∈ {0.5, 0.7, 0.8}.
The Gaussian row (feature off) is re-run in the same code path.

## Hypotheses (per feature f, test-rule r, sample size n)

- H0_alone: rate(f, H=0.5) = rate(Gauss, H=0.5).
  H1_alone: rate(f, H=0.5) > rate(Gauss, H=0.5).
- H0_int: [rate(f, h) − rate(f, 0.5)] = [rate(Gauss, h) − rate(Gauss, 0.5)], h ∈ {0.7, 0.8}.
  H1_int: the left side differs (feature amplifies or dampens the memory effect).

## Common notation

- ξ_t: unit-variance fGn with Hurst H, exact Davies–Harte (`mbcore.fgn`, the E3 generator).
  Used for every generator in both test families, so the Gaussian row shares the code path.
  At H = 0.5, ξ is i.i.d. N(0, 1).
- All tests used here are invariant to a positive rescaling of a series (median-heuristic
  Gaussian kernel; distance covariance compared with its own resampled reference). So
  marginal scale does not matter. Generators are still standardized where it is cheap.
- "Memory" means the H of the driving Gaussian input. The unit check reports the
  local-Whittle Ĥ of the observed series and of its square.

## Generators (one per row)

| id | feature | equation | Gaussian-subordinated? |
|---|---|---|---|
| G0 | none (baseline) | X_t = ξ_t | yes, G = identity, Hermite rank 1 |
| HT-S | heavy tails | X_t = F_{t3}^{-1}(Φ(ξ_t)) / √3 | **yes**, instantaneous G(ξ_t); G is odd and increasing, so Hermite rank 1. The paper's rank theory applies. |
| HT-L | heavy tails | X_t = Σ_{j=0}^{M} ψ_j ε_{t−j} / (Σψ_j²)^{1/2}, ψ_j = Γ(j+d)/(Γ(d)Γ(j+1)), d = H − 1/2, ε_t i.i.d. t_3/√3, M = 20n | **no**. A linear long-memory (ARFIMA(0,d,0)) filter of non-Gaussian innovations. Outside the Hermite theory. Tail index stays 3 (linear filter of regularly varying noise). |
| NL-SV | nonlinear dynamics (volatility clustering) | X_t = exp(τ h_t / 2) ε_t, h_t = (φ h_{t−1} + √(1−φ²) ξ_t)/s_H, φ = 0.9, τ = 1, ε_t i.i.d. N(0,1) independent of ξ, burn-in 200 | **multivariate only**. X is a function of the two Gaussian processes (h, ε). Covered by multivariate Hermite expansions (Arcones 1994), **not** by the paper's univariate G(ξ_t) calculation. X is serially uncorrelated; memory sits in its volatility. Heavy-ish tails (kurtosis 3e^{τ²} ≈ 8.2). |
| NL-TAR | nonlinear dynamics (threshold AR) | X_t = a(X_{t−1}) X_{t−1} + ξ_t, a(x) = 0.8 if x ≤ 0, −0.4 if x > 0, X_{−200} = 0, burn-in 200 | **no** (in the paper's sense). X_t is a nonlinear functional of the whole Gaussian path, not G(ξ_t). The map is a contraction (Lipschitz 0.8), so a unique stationary solution exists; burn-in error ≤ 0.8^200. |
| SEAS-M | seasonality in mean | X_t = ξ_t + a sin(2π(t + U)/P), a = 1, P = 50, U ~ Uniform{0, …, P−1} independent of ξ | **no**. Additive deterministic periodic term with random phase. |
| SEAS-V | seasonality in variance | X_t = (1 + b sin(2π(t + U)/P)) ξ_t, b = 0.5, P = 50, U as above | **no**. Periodic scale. |

s_H in NL-SV is the exact standard deviation of the AR(1)-filtered fGn, computed from the
fGn autocovariance (truncated at 400 lags, φ^400 < 1e-18), so h_t has unit variance at every H.

Why these settings: t_3 matches the tail index often reported for returns and is finite in
mean (dCov needs a finite first moment). τ = 1 is a strong but common log-volatility scale.
φ = 0.9 gives clear short-memory clustering at H = 0.5, so "nonlinear dynamics, no long
memory" is a real cell. The TAR coefficients give an asymmetric, skewed, stationary process.
a = 1 gives seasonal variance 0.5 against noise variance 1. These are not calibrated to any
market series.

Diagnostic row (independence family only, excluded from the verdicts):

| id | equation | null status |
|---|---|---|
| SEAS-M-shared | both series X_t = ξ_t + sin(2πt/P), Y_t = η_t + sin(2πt/P), fixed phase 0 | X and Y are independent **processes** (the seasonal term is deterministic), but neither is stationary and the pooled pairs (X_t, Y_t) co-move. This is the situation a common clock pattern creates in market data. The known-null property holds only in the process-independence sense, so it is reported separately. |

## Null argument per family

**Split MMD (one series, first half vs second half, h1 = n/2).**
- G0, HT-S, HT-L, NL-SV, NL-TAR are strictly stationary after the stated burn-in (HT-L is
  stationary up to the MA truncation at M = 20n lags). The two halves have the same
  marginal law, so the population MMD is zero and any rejection is a false alarm.
- SEAS-M and SEAS-V use **both** a random phase and a whole number of periods per half.
  P = 50 divides n/2 = 250 and 1000. Random phase makes the process strictly stationary
  (shifting t by one maps U to U+1 mod P, which has the same law). Whole periods make each
  half contain every seasonal value equally often, so the half-averaged marginals, and
  hence the population mean embeddings that MMD compares, are equal for **every** phase,
  not only on average. So the null holds both unconditionally and conditionally on U.

**HSIC and dCov (two series).** X and Y are generated from the same generator and H with
independent random draws (independent innovations, independent volatility processes and
independent phases U_X, U_Y). Hence the processes X and Y are independent. Every pair
row is therefore a true null.

## Grid

- Main factorial: 7 generators × H ∈ {0.5, 0.7, 0.8} × n ∈ {500, 2000} = 42 cells per family.
  Two families (MMD; independence) = 84 cells. Independence uses H_X = H_Y.
- Diagnostic: SEAS-M-shared × 3 H × 2 n = 6 independence cells.
- Total: 90 cells. 1,000 replications per cell in the main run; 100 in the pilot.

At H = 0.5, HT-S and HT-L both reduce to i.i.d. t_3/√3. Their H = 0.5 cells are two
independent replications of the same law (different seeds and code paths), which gives a
free consistency check.

## Tests and rules (reused unchanged from EXP1)

Code copied from `../EXP1-20261016/` with hashes in `copied_from_EXP1.sha256` (they match
the EXP1 manifest). No test code is rewritten.

**Independence: `mbcore.kind_tests(x, y, rng, B=500)` and `mbcore.lind_tests(x, y)`, as in E3.**
- HSIC wild (`p_hsic_wb ≤ 0.05`): centred block-Markov Rademacher multipliers, ℓ = ⌊√n⌋, B = 500.
- HSIC shift (`p_hsic_shift ≤ 0.05`): circular shifts in [⌊0.1n⌋, n − ⌊0.1n⌋].
- HSIC sub (`rej_hsic_sub`): windows of length ⌊√n⌋, cross offset ⌊0.1n⌋.
- dCov wild (`p_dcov_wb ≤ 0.05`), dCov sub (`rej_dcov_sub`).
- Pearson NW and EWC are stored as a linear reference; they are not part of the verdicts.

**Split MMD: `mmdcal.rep_e2(x, rng)`, exactly the E2 replication.**
- Primary rule: uncentred wild multipliers, ℓ = ⌊n^0.5⌋ (`p_unc_05 ≤ 0.05`), the E2 rule.
- Secondary rules from the same call: √n subsampling (`t_n > q95_sqrt`) and self-normalized
  block MMD with ℓ = ⌊n^0.5⌋ (`sn_T > sn_q95_05`). The other stored E2 columns are kept but
  not read.

## Seeds

Replication r of cell c uses `np.random.SeedSequence(entropy=MASTER, spawn_key=(c, r))`,
one `default_rng` per replication (the EXP1 rule). The same rng then feeds the generator
and the test's multipliers. Masters:

| stage | family | master |
|---|---|---|
| unit check | all | 2026102301 |
| pilot | MMD | 2026102321 |
| pilot | independence | 2026102331 |
| main | MMD | 2026102322 |
| main | independence | 2026102332 |

Cell indices are fixed by the order of the cell list in `run_exp2.py` and are written to
each run's `manifest.json`.

## Pre-declared reading rule (fixed before the pilot)

Notation for one family, rule and n: B_h = Gaussian (G0) rate at H = h; F_h = rate of a
feature generator at H = h. Rates are rejection fractions at nominal 5% over 1,000 reps
(errors excluded and counted). Wilson 95% intervals as in EXP1 (`common.wilson`).

1. **Breaks alone (A):** F_0.5 − B_0.5 > 2 pp **and** their Wilson intervals do not overlap.
2. **Memory interaction at h ∈ {0.7, 0.8}:** I_h = (F_h − F_0.5) − (B_h − B_0.5).
   SE_h = √(Σ p(1−p)/R) over the four cells. **Amplifies (I+)** if I_h > 2 pp and
   I_h / SE_h > 1.96 at h = 0.7 or h = 0.8. **Dampens (I−)** if I_h < −2 pp and
   I_h / SE_h < −1.96 at either h (and no I+).
3. **Verdict per generator × rule, read at n = 2000** (n = 500 reported alongside):
   - "breaks alone" = A and not I+;
   - "breaks alone and amplifies memory" = A and I+;
   - "only with memory" = not A and I+;
   - "dampens memory" = I− and not I+ (A reported in the same line);
   - "no effect" = none of these.
4. **Per feature:** heavy tails, nonlinear dynamics and seasonality each have two generators.
   The feature verdict lists both generators' verdicts; it does not pool them. Where the two
   disagree, the subordinated/non-subordinated split is the reading to report.
5. The Gaussian memory row (B_h − B_0.5) is reported for context. It is the "memory alone"
   effect and is already known from E2/E3; EXP2 re-measures it, not tests it.
6. **Multiplicity:** about 6 generators × 8 rules × 2 h interaction tests per n. At a
   1.96 threshold, around two false "I" flags are expected by chance across the grid. The
   2 pp floor guards partly against this. Isolated flags are reported as such, not as findings.
7. **Suspicion check:** if a feature cell rejects *less* than nominal by more than 2 pp, or
   the two H = 0.5 heavy-tail cells (same law) differ by more than their intervals allow,
   it is flagged in results.md before any verdict is read.

## Unit checks (stage 2, saved to `unit_checks.txt`)

On 200 series of length 2000 per generator × H (seed 2026102301):
- marginal: mean, sd, skewness, sample kurtosis, Hill tail-index estimate (top 5%);
- memory: mean local-Whittle Ĥ of X and of X² (m = ⌊n^0.65⌋, the E2 estimator);
- lag-1, 10, 100 ACF of X and |X|;
- seasonality: mean periodogram at the seasonal frequency 1/P against neighbouring
  frequencies, and the per-half mean of X (whole-period check);
- for MMD-null sanity: the half-1 vs half-2 difference in mean and sd, averaged over series.

Expected: HT-S/HT-L Hill index ≈ 3; NL-SV kurtosis ≈ 8 and Ĥ(X) ≈ 0.5 but Ĥ(X²) > 0.5
when H > 0.5; SEAS rows show a periodogram spike at 1/P.

## Compute estimate (from EXP1 per-row CPU)

| family | per rep n = 500 | per rep n = 2000 | settings | CPU-h at 1,000 reps |
|---|---|---|---|---|
| independence (`kind_tests` + `lind_tests`) | 0.06 s | 1.0 s | 21 + 3 diagnostic | ≈ 7.1 |
| MMD (`rep_e2`) | 0.07 s | 1.53 s | 21 | ≈ 9.3 |
| **total** | | | | **≈ 16.5** |

At 16 worker processes, perfect scaling gives about 1.1 h wall. The machine has other load,
so 1.5–3 h is the working estimate. The pilot replaces this estimate with a measured one.

## Budget gate

At most 16 worker processes. If the pilot extrapolates the main run to ≤ 6 h wall, run it in
full under nohup, resumable by cell (rows already on disk are skipped). If not, cut in this
order and record the cut here: (1) drop the SEAS-M-shared diagnostic; (2) drop the secondary
MMD rules by replacing `rep_e2` with its unc-only part (would need a new function, so this is
the last resort); (3) run n = 2000 only.

## Amendments

(1) 2026-10-25, after the pilot, before the main run. NL-TAR memory estimates sit above the
nominal H. Mean local-Whittle estimate over 100 pilot reps (G0 for reference):

| H | G0 n=500 | G0 n=2000 | NL-TAR n=500 | NL-TAR n=2000 |
|---|---|---|---|---|
| 0.5 | 0.494 | 0.494 | 0.671 | 0.589 |
| 0.7 | 0.701 | 0.698 | 0.930 | 0.819 |
| 0.8 | 0.813 | 0.800 | 1.031 | 0.931 |

The upward bias appears even at H = 0.5 and shrinks from n = 500 to n = 2000. This pattern fits
short-range persistence from the TAR recursion (coefficient 0.8 below zero) leaking into the
low-frequency band that the estimator uses. It does not show extra long memory. Consequence: an
NL-TAR excess over G0 at fixed nominal H can come from that short-range persistence. Read NL-TAR
interactions as "nonlinear dynamics plus short-range persistence", not as a pure nonlinearity
effect. No change to the design or code. The main run proceeds as specified.

(2) 2026-10-25, pilot check. NL-SV and NL-TAR gave identical mmd_unc_wild rejection counts at
n = 2000 (11, 18, 22 of 100). The raw rows, seeds, p-values and rejecting replications all differ.
The match is a chance coincidence, not a bug.

