# EXP6 (2026-10-25) — bandwidth sensitivity of the Gaussian-kernel HSIC rules

Written before the pilot and the main run. Simulation only. No market data.
Exploratory: rates describe the tested implementations under the stated generators.

## Question

A reviewer notes that every HSIC result uses the median-heuristic Gaussian kernel. Kernel tests
are bandwidth-sensitive. The paper (§4.4, Table 5) reports size-adjusted HSIC power of about 14%
under the seasonal scale (SEAS-V, linear alternative) while dCov and Pearson keep more. The paper
says it does not know why. EXP6 asks two questions.

1. Is the size picture (which HSIC rules over-reject in which cell) stable when the bandwidth is
   scaled by κ ∈ {0.5, 1, 2}?
2. Does a different bandwidth recover size-adjusted HSIC power under SEAS-V? In other words, is
   the weakness an artefact of the median heuristic?

## Hypotheses

Size (per null cell, per HSIC rule):

- H0_size: the raw rejection rate at κ = 0.5 and κ = 2 equals the rate at κ = 1.
- H1_size: it differs. Read with a paired McNemar test on the same replications.

Reading rule, fixed now: the size picture is **stable** if every rule keeps its EXP2 verdict
class across κ. The classes are: valid (Wilson 95% interval contains 5%, or upper bound ≤ 7.5%),
mildly oversized (rate in (5%, 10%]), oversized (rate > 10%).

Power (per alternative cell, per HSIC rule and the HSIC oracle):

- H0_pow: size-adjusted power at κ = κ' equals size-adjusted power at κ = 1.
- H1_pow: it differs.

Reading rule, fixed now: a bandwidth **recovers** HSIC power under SEAS-V if some κ ≠ 1 raises
size-adjusted power of the HSIC oracle, or of the best HSIC rule, by more than 10 pp, and the
paired bootstrap 95% interval of the gain excludes 0. A gain to within 5 pp of dCov wild counts
as **full** recovery. Otherwise the weakness is **not a bandwidth artefact** over this κ range.

## Mechanism prediction (written before any EXP6 data)

Under SEAS-V, X_t = s_X(t) ξ_t and Y_t = s_Y(t) η_t, with s(t) = 1 + 0.5 sin(2π(t+U)/50) and
independent phases U_X, U_Y. Given the two phases, the deterministic scales s_X and s_Y are
correlated, with correlation cos(2π(U_X − U_Y)/50). So, within one sample, |X| and |Y| co-move
even under independence. That co-movement averages to zero only over phases. A Gaussian kernel
sees it, because the kernel reacts to second moments. Pearson does not see it, because the sign
of ξ is symmetric. So the null distribution of HSIC is widened, c* rises, and power falls.

For a large bandwidth, k(x, x') ≈ 1 − (x − x')²/(2σ²) + (x − x')⁴/(8σ⁴). Then σ⁴·HSIC tends to
a multiple of Cov(X, Y)², and the second-moment (scale) term shrinks like 1/σ² relative to it.
Prediction: κ = 2 raises size-adjusted HSIC power under SEAS-V-lin (towards Pearson). κ = 0.5
lowers it. Under G0-sq (Corr = 0, quadratic dependence) the order reverses: κ = 2 loses power.
To test the mechanism, each replication also records corr(X², Y²) and the median bandwidths.
If the null HSIC statistic tracks corr(X², Y²)² across SEAS-V replications, that supports the
mechanism. A pure "bandwidth inflated by the scale mixture" story predicts only a small change,
because Var(X) = 1 + 0.5²/2 = 1.125 under SEAS-V, about 6% on the median distance.

## Bandwidth

Median heuristic as in EXP2/EXP3: σ_X = median of off-diagonal |X_i − X_j| (mbcore.median_bandwidth).
EXP6 uses σ_X(κ) = κ·σ_X and σ_Y(κ) = κ·σ_Y. The same scaled σ is used in the full-sample Gram
matrices (wild, shift, oracle) and inside the subsampling windows (subsample_stats gets the
scaled sx, sy). dCov and Pearson do not depend on κ. They are computed once per replication as
references.

κ ∈ {0.5, 1, 2} is the primary grid. κ = 4 is an exploratory probe of the large-bandwidth limit.
The pilot timing decides whether κ = 4 enters the main run (amendment A1).

**Pairing.** All κ use the same (X, Y) and the same 500 block-Markov wild multipliers in one
replication. κ comparisons are therefore paired. At κ = 1 every field equals EXP3 `kind_ext` on
the same seed. `unit_checks.py` verifies this to 0 difference.

## Rules and size adjustment

Exactly EXP3 (experiments/EXP3-power-20261025/DESIGN.md). Per rule the ratio R = T / c_ref, with
c_ref the rule's own 95% reference quantile. Raw rule: p ≤ 0.05 (wild, shift) or T > q95 (sub).
Size-adjusted cut-off c* = the k-th order statistic of the null ratios with
#(null > c*) = ⌊0.05 m⌋ (EXP3 `cstar`). Oracle: the bare statistic n·HSIC_V(κ), n·dCov_V, |r|.
Uncertainty: Wilson 95% for counting error; bootstrap 95% (2,000 resamples, seed 42) resampling
null and alternative replications independently, recomputing c*. κ-versus-κ gaps use the same
resamples, so the pairing is kept.

## Cells (n = 2000 primary)

| id | arm | H | role |
|---|---|---|---|
| G0-null | X, Y independent fGn | 0.5, 0.8 | size; c* for G0 alternatives at H = 0.5 |
| SEAS-V-null | independent seasonal scale | 0.5 | size; c* for SEAS-V-lin |
| G0-lin | Y = ρξ + √(1−ρ²)w, ρ = 0.065 | 0.5 | power |
| G0-sq | Y = ρ(ξ²−1)/√2 + √(1−ρ²)w, ρ = 0.075 | 0.5 | power |
| SEAS-V-lin | as G0-lin, each series times its own seasonal scale, ρ = 0.065 | 0.5 | power |

Generators are EXP3 `power_core.ind_pair` (copied into `bw_core.py`, checked to be identical
on the same seed). Effect sizes are the EXP3 values, not re-tuned.

n = 500 (ρ doubled: 0.13, 0.15, as in EXP3) is a secondary block. It runs only if the n = 2000
block finishes inside the budget.

## Seeds and runs

- Seeds: SeedSequence(entropy=MASTER, spawn_key=(cell_idx, rep)); one default_rng per replication.
  MASTER = 2026102561 (pilot), 2026102562 (main). Not used by any other experiment.
- Regression mode: master 2026102542 with the EXP3 cell indices (0, 2, 6, 18, 24, 30). The first
  5 replications of each of the six cells are recomputed and compared with the EXP3 main rows.
  The κ = 1 fields must match exactly. This reuses the EXP3 master for a check only. No result
  in EXP6 tables comes from it.
- B = 500 wild draws. `python -I`. BLAS threads pinned to 1.
- Workers: at most 2 while EXP4 (`run_exp4.py`) runs, at most 5 after. Wall budget about 2.5 h.
- Pilot: 40 reps per n = 2000 cell (κ ∈ {0.5, 1, 2, 4}), for timing and a sanity reading.
- Main: 1,000 reps per cell if the projected wall time fits the budget. Otherwise 500 reps per
  cell, fixed by amendment A1 before the main run. The rep count is the same for all cells.

## Bias and leakage checklist (simulation version)

- Look-ahead: not applicable. Every statistic uses the full simulated path by design.
- Survivorship, universe selection, transaction costs: not applicable (no market data).
- Cherry-picking: κ grid, cells, effect sizes and reading rules are fixed above. Effect sizes
  come from EXP3. No κ is selected on the data. Picking the best κ per cell after the fact
  would be an oracle choice, and EXP6 reports it as such, not as a usable procedure.
- Threshold reuse: c* comes from null replications independent of the alternative replications.
- Skill registry paths named in the agent protocol (/.claude/skills/ml_stats, /scripts/...) do
  not exist on this machine (checked in EXP3). No registry skill is used.

## Code

`bw_core.py` (κ-variant of `kind_ext`, generators), `run_exp6.py` (driver, resumable),
`analyze_exp6.py` (size, size-adjusted power, intervals, tables), `unit_checks.py`.
Core imported from `kmem.legacy` (identical to the EXP2/EXP3 `mbcore.py`).

## Amendments

(A1 is added after the pilot, before the main run.)

### A1 (2026-10-25 19:00 UTC, after the pilot, before the main run): κ grid, reps, stop rule

**Pilot.** Master 2026102561, κ ∈ {0.5, 1, 2, 4}, n = 2000. Stopped after 48 replications
(8 per cell) instead of 40 per cell, because the machine is overcommitted (load ≈ 49 on 22 cores,
EXP4 still running with 16+ workers). Cost: 7.2 CPU-s per replication (median 7.1), 0 errors.
Pilot sanity reading (8 reps, not evidence): under the SEAS-V null, the HSIC statistic at
κ ≤ 2 ranks almost exactly with corr(X², Y²)² (Spearman ≈ 0.9). At κ = 4 it ranks with r²
instead, and raw HSIC wild and sub rejections fall. This matches the mechanism written above.

**Changes.**
1. κ = 4 enters the main run with all three rules, because the pilot suggests the qualitative
   change sits between κ = 2 and κ = 4. κ ∈ {0.5, 1, 2} stays the primary grid for the reading
   rules. κ = 4 is reported as a secondary column with the same reading rules.
2. Oracle-only κ ∈ {0.25, 8}: the bare statistic n·HSIC_V only, no rules (about 0.3 CPU-s each).
   It maps the two bandwidth limits for the oracle row. Unit check: the cheap path equals the
   full path to 7e-14 relative.
3. Replications: 1,000 per cell does not fit. Projected cost is about 8 CPU-s per replication,
   and only 2 workers are allowed while EXP4 runs. **Target 500 per cell.** Chunks are
   interleaved by replication index across the six cells, so a stopped run has balanced cells.
4. Stop rule, fixed now: the main run is killed at 20:45 UTC if it has not finished. Each cell is
   then truncated to the replication indices 0..m−1 that are complete in **all** six cells, and
   only those enter the analysis. If EXP4 finishes earlier, the run is restarted with 5 workers
   (resume skips completed replications, so seeds and rows are unchanged).
5. n = 500 block: dropped (no budget).
6. Size reading with m < 1,000: Wilson intervals widen to about ±2 pp at 5%. The size class
   rule is unchanged.
