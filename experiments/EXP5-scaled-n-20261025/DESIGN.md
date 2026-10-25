# EXP5 (2026-10-25) — null-SD-scaled MMD shift, and larger n for the memory-sum cells

Written before the pilot and the main run. Simulation only. No market data.
Exploratory: rates and power describe the tested implementations under the stated generators.

## Reviewer concerns

1. EXP3 fixed the split-MMD shift at δ = 0.1 (n = 2000). Under the null, the half-mean
   difference has SD 0.045 at H = 0.5 but 0.247 at H = 0.8. So δ = 0.1 is about 2.2 null SDs
   without memory and 0.4 null SDs at H = 0.8. Part of the power collapse may be less
   information in long-memory data, not a weakness of the kernel test.
2. The memory-sum cells (sums 1.4 and 1.6) stop at n = 2000. Does the excess at sum 1.4 stay
   flat, and does the sum-1.6 rate keep falling, at larger n?

## Part 1 — split MMD with a shift in null-SD units

**Null SD.** For unit-variance fGn(H) and n = 2m, the half-mean difference
D = mean(x_{m+1..n}) − mean(x_{1..m}) has SD sd_H = sqrt(4 m^{2H} − (2m)^{2H}) / m.
Checked against the exact fGn autocovariance (`exp5_core.sd_half_exact`) to 5 digits:

| H | n = 500 | n = 2000 |
|---|---|---|
| 0.5 | 0.0894 | 0.0447 |
| 0.7 | 0.2226 | 0.1469 |
| 0.8 | 0.3262 | 0.2472 |

**Shift.** x_t = ξ_t + δ_H 1{t ≥ n/2}, δ_H = c · sd_H(n). Two values of c:

- c = √5 ≈ 2.236. This gives δ = 0.100 at H = 0.5, n = 2000, and δ = 0.200 at n = 500, which are
  exactly EXP3's shifts. The H = 0.5 cells therefore replicate EXP3 on new seeds.
  δ = 0.328 (H = 0.7) and 0.553 (H = 0.8) at n = 2000.
- c = 1.5 (robustness). δ = 0.067, 0.220, 0.371 at n = 2000.

**Rules.** EXP3 split-MMD bundle, unchanged formulas: `mmd_unc_wild` (uncentred multiplier,
l = ⌊√n⌋, B = 500), `mmd_sub_sqrt` (window quantile, l = ⌊√n⌋), `mmd_sn_block`
(self-normalised, window quantile, l = ⌊√n⌋), and the oracle (bare MMD², size-adjusted).
Size adjustment as in EXP3: ratio R = T / c_ref, c* = empirical 95% quantile of R over the
1,000 matched null replications (same H, n), reject iff R > c*.

**Linear benchmark.** The half-mean difference itself, size-adjusted on the same null cells:
`lin_oracle` uses |D|, `lin_t` uses |D / SE_iid| (the naive two-sample Welch t, whose iid
standard error is wrong under memory, then size-adjusted). Because D is exactly Gaussian with
SD sd_H under the null, the exact size-adjusted power of |D| at 5% is
P(|Z + c| > 1.96): 60.9% at c = √5 and 32.3% at c = 1.5, the same at every H. This is what a
mean-based test can do at the scaled shift. It is printed next to the empirical value.

**Null cells.** Fresh nulls are run in the EXP5 code path, because the linear benchmark needs
D, which the EXP3 null CSVs do not store. The EXP3 null cells (same rules, ratio columns
present) serve as a regression check: fresh raw null rates and the c* values must agree with
EXP3 within Monte Carlo error.

**Hypotheses (per rule r, per c, at n = 2000; n = 500 reported the same way).**

- H0 (information story): adjusted power of r at H = 0.8 equals that at H = 0.5.
- H1 (rule or kernel weakness): adjusted power falls with H even at constant scaled shift.

**Reading rule, fixed now.** Let Δ_r = adj. power(H = 0.5) − adj. power(H = 0.8), with a
bootstrap 95% interval that resamples null and alternative replications (2,000 resamples,
seed 42, threshold error included).

- **flat**: |Δ_r| ≤ 10 pp.
- **collapse**: Δ_r > 10 pp, the interval excludes 0, and power at H = 0.8 is below half of the
  power at H = 0.5.
- **partial loss**: Δ_r > 10 pp and the interval excludes 0, but not a collapse.
- **gain**: Δ_r < −10 pp and the interval excludes 0.
- otherwise **unclear**.

The key question is answered from the oracle row (statistic) and the three rule rows
(critical-value rules). If the oracle is flat but a rule collapses, the loss is in the rule.
If the oracle also falls, the loss is in the MMD statistic. Gap to the linear benchmark is
reported as context. A rule that beats the linear oracle by a wide margin would be suspicious
for a pure mean shift, and is flagged.

## Part 2 — independence null at larger n

**Generator and rules.** The EXP1 E3 code path unchanged: independent fGn pairs, proc G4,
`mbcore.gen_series`, `mbcore.kind_tests` (B = 500), `mbcore.lind_tests`. Rules reported:
HSIC wild, HSIC shift, HSIC sub, dCov sub, Pearson NW, Pearson EWC (dCov wild is computed and
reported as an extra column).

**Cells.** H_X = H_Y ∈ {0.5, 0.7, 0.8}, i.e. memory sums 1.0, 1.4, 1.6.
n = 4000 with 1,000 replications. n = 8000 if the pilot timing allows, with fewer replications
(set in amendment A1). The comparison rows are EXP1 E3 at n = 500, 1000, 2000 (2,000 reps,
`experiments/EXP1-20261016/e3_summary.csv`) and, for the H = 0.5 baseline, the KL reference
row in `e3_tables.md` (the paper's Table 3 baseline column).

**Hypotheses (per rule).**

- Sum 1.4. H0: the rejection rate at n = 4000 equals the n = 2000 rate (the excess is flat).
  H1: it differs.
- Sum 1.6. H0: the rate at n = 4000 is not below the n = 2000 rate. H1: it keeps falling.

**Reading rule, fixed now.** Two-proportion z between n = 4000 (EXP5) and n = 2000 (E3).
|z| < 2 reads as "no change detected at this precision". For sum 1.6 a fall is read only if
z < −2. The trend over n = 500…4000 is also summarised by a weighted least-squares slope of the
rate on log2(n) (pp per doubling), weights 1/variance. At sum 1.4 an excess is "flat" if the
n = 4000 rate is within 2 SE of the n = 2000 rate and remains above 5% by more than 2 SE.
The MC standard error at 1,000 reps is about 0.9 pp near 9% and 1.3 pp near 20%, so changes
smaller than about 2.5 pp cannot be resolved. That limit is reported, not hidden.

## Seeds, resources, budget

- SeedSequence(entropy=2026102551, spawn_key=(cell_idx, rep)) for the pilot,
  2026102552 for the main run. Never reused from EXP1–EXP4.
- `python -I`, BLAS threads pinned to 1 before numpy import.
- Workers: at most 3 while EXP4 (`run_exp4.py`) is running, at most 6 after it finishes.
- Memory: an n = 4000 replication holds about 6 n×n float64 matrices (~0.8 GB). n = 8000 is
  about 3–4 GB per worker. `free -g` is checked before n = 8000 is launched.
- Wall budget about 3 h for pilot and main together. Rough prior estimate: Part 1 about
  2.5 CPU-h, Part 2 at n = 4000 about 8 CPU-h (≈10 s per replication). The pilot fixes the
  numbers. If Part 2 does not fit, the n = 4000 reps are cut for the H = 0.5 baseline first,
  then n = 8000 is dropped, before the sum 1.4 and 1.6 cells are cut. Any cut is an amendment.

## Pilot

- Part 1: 50 reps for the null and the c = √5 alternative at H = 0.5 and 0.8, n = 2000 only
  (4 cells). Checks: δ, timing, no errors, and raw power at H = 0.5 near EXP3's value.
- Part 2: 4 reps at n = 4000 (H = 0.5, 0.8) and 2 reps at n = 8000 (H = 0.8) for timing and
  peak memory.
- Pilot rows are never pooled with main rows.

## Code

`exp5_core.py` (generators and rule bundles; imports `kmem.legacy` without edits),
`run_exp5.py` (driver), `analyze_exp5.py` (size adjustment, intervals, tables, reading rules),
`unit_checks.py` (EXP5 `mmd_ext` equals EXP3 `power_core.mmd_ext` on the same seed for every
shared field; `ind_rep` equals the E3 code path on the same seed).

## Bias and leakage checklist (simulation version)

- Look-ahead: not applicable. Every statistic uses the full simulated path by design.
- Survivorship, universe selection, transaction costs: not applicable (no market data).
- Cherry-picking: c = √5 is set by EXP3's δ, not by looking at EXP5 output. c = 1.5 is the
  value suggested in the brief. Grid and reading rules fixed above.
- Threshold reuse: c* comes from null replications with seeds disjoint from the alternatives.
- Skill registry paths in the agent protocol (/.claude/skills/ml_stats, /scripts/...) do not
  exist on this machine. No registry skill was used.

## Amendments

(Added after the pilot, before the main run.)

### A1 (2026-10-25 18:58 UTC, after the pilot, before the main run): timing, reps, order

Pilot: seed master 2026102551, `pilot/`, `pilot.log`. 0 errors. Unit checks pass
(`unit_checks.txt`: EXP5 `mmd_ext` and `ind_rep` reproduce EXP3 and E3 to 0 difference).

Timing under shared load (EXP4 on 16 workers, load ≈ 40 on 22 cores):

| cell | CPU-s / rep | wall-s / rep | peak RSS / worker |
|---|---|---|---|
| MMD, n = 2000 | 0.96–1.02 | ≈ 1.1 | small |
| independence, n = 4000 | 9.0–12.0 | 10.4–13.9 | 1.05 GB |
| independence, n = 8000 | 37–47 | 47–58 | 4.0 GB |

MMD pilot sanity (50 reps, not evidence): raw power at H = 0.5, c = √5 is 50–54%, close to
EXP3's 38–46% on 1,000 reps. Null half-mean SD is 0.044 (H = 0.5) and 0.228 (H = 0.8),
against 0.045 and 0.247 from the formula.

Projected main run: Part 1 ≈ 2.8 CPU-h; Part 2 sums 1.4 and 1.6 at n = 4000, 1,000 reps each,
≈ 7 worker-h. That is ≈ 3.3 h wall on 3 workers and ≈ 1.7 h on 6. Decisions:

1. Order: Part 1 (all 18 MMD cells, 1,000 reps), then sums 1.4 and 1.6 at n = 4000 in steps of
   250 reps up to 1,000, then the H = 0.5 baseline at n = 4000 with **500 reps** (cut first,
   as planned). Workers are re-checked before each stage (3 while EXP4 runs, 6 after).
   Script: `run_main.sh`.
2. **n = 8000 is dropped from the main plan.** At ≈ 50 s per replication it would need
   ≈ 4 worker-h for 300 reps at one memory sum, which does not fit. If the n = 4000 stages
   finish with more than 45 min left, a small n = 8000 sum-1.6 cell (≤ 150 reps) may be added
   as amendment A2. It would be reported as low-precision.
3. If the wall budget runs out, the reps reached in each cell are reported as they stand.

### A2 (2026-10-25 19:18 UTC, during the main run, no results looked at): parallel stages

EXP4 finished at about 19:10 UTC. `run_main.sh` checked EXP4 with `pgrep -f run_exp4.py`, which
also matches other agents' shell watchers that contain that string, so it would have stayed at
3 workers. `run_main.sh` was stopped after its first stage had started (that stage, Part 1 reps
0–499, kept running). It is replaced by two scripts that run side by side, 6 workers in total:
`run_part1b.sh` (Part 1 reps 500–999, 3 workers, after the first stage ends) and
`run_part2.sh` (sums 1.4 and 1.6 at n = 4000 in steps of 125 reps to 1,000, then the
H = 0.5 baseline to 500 reps; 3 workers while a Part 1 pool runs, 6 after). Cells, seeds and
rep targets are unchanged. Resume logic skips reps already on disk, so no rep is run twice.
The machine is also loaded by other users' jobs (load ≈ 60 on 22 cores), so wall times are slower
than the pilot.

### A3 (2026-10-25 20:28 UTC, during the main run, Part 2 rates not yet looked at): n = 8000

Part 2 now runs at about 1 rep per 1.6 s on 6 workers (machine load fell), so the n = 4000 stages
should end near 21:00 UTC, about 50 min before the budget end. As allowed by A1 item 2,
`run_part3.sh` then runs n = 8000 at sum 1.6 (cell 105), and also at sum 1.4 (cell 104), each in
stages of 48 reps up to 144, 6 workers (≈ 4 GB each, 71 GB available). A stage starts only before
21:25 UTC (sum 1.6) or 21:35 UTC (sum 1.4). The n = 8000 rows are low precision (MC SE ≈ 3 pp
near 20%) and are reported as such. They do not enter the reading-rule z test, which uses
n = 4000 against n = 2000. They do enter the slope fit as extra points.
