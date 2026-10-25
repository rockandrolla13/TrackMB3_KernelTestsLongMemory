# EXP3 (2026-10-25) — size-adjusted power of the independence and split-MMD rules

Written before the pilot and the main run. Simulation only. No market data.
Exploratory: power figures describe the tested implementations under the stated generators.

## Question

The paper finds circular-shift HSIC closest to the nominal 5% size in every arm. A reviewer
objects that a rule which rarely rejects also looks good on size. EXP3 asks:

**After every rule is tuned to a true 5% false-alarm rate on the matched null cell, does
circular-shift HSIC keep power comparable to the other rules?**

## Hypotheses

For each alternative cell a (generator, dependence type, H, n):

- H0: size-adjusted power(HSIC shift) = size-adjusted power(best other HSIC rule) in cell a.
- H1: size-adjusted power(HSIC shift) < size-adjusted power(best other HSIC rule).

Reading rule, fixed now: shift is **competitive** in a cell if its size-adjusted power is at most
5 pp below the best other HSIC rule (wild, subsampling). It **loses** if the gap exceeds 5 pp and
the paired bootstrap 95% interval of the gap excludes 0. Otherwise the cell is **unclear**.
The same comparison is reported against dCov and Pearson rules, as context only.

## Size adjustment (oracle)

Every rule is reduced to one continuous number per replication, its **ratio**
R = T / c_ref, where T is the rule's statistic and c_ref is the rule's own 95% reference
quantile, computed on the same reference draws that give its p-value:

| rule | T | c_ref | raw rule (as EXP2) |
|---|---|---|---|
| HSIC wild | n·HSIC_V | 95% quantile of 500 block-Markov wild draws | p ≤ 0.05 |
| HSIC shift | zero-shift HSIC | 95% quantile of circular shifts ⌊0.1n⌋..n−⌊0.1n⌋ | p ≤ 0.05 |
| HSIC sub | n·HSIC_V | 95% quantile of window statistics, l = ⌊√n⌋ | T > c_ref |
| dCov wild, dCov sub | n·dCov_V | as for HSIC | as for HSIC |
| Pearson NW | \|t_NW\| | 1.96 | \|t\| > 1.96 |
| Pearson EWC | \|t_EWC\| | t_{ν,0.975} | \|t\| > t_ν |
| MMD uncentred wild | MMD²_V | 95% quantile of 500 uncentred multiplier draws, l = ⌊√n⌋ | p ≤ 0.05 |
| MMD sqrt subsampling | n·MMD² | 95% window quantile, l = ⌊√n⌋ | T > c_ref |
| MMD self-normalised block | T_SN | 95% window quantile of T_SN, l = ⌊√n⌋ | T > c_ref |

The raw rule rejects (approximately, exactly for the window rules) iff R > 1.
The size-adjusted rule rejects iff R > c*, where c* is the empirical 95% quantile of R over the
1,000 null replications of the matched null cell (same generator, H, n). The null rejection rate
is then 5% by construction.

Why the ratio and not the p-value: p-values from 500 wild draws or from windows tie at their
floor when the raw rule is badly oversized (SEAS-V: 73% raw size). A 5% quantile of tied
p-values is not identified. The ratio is continuous. As a robustness column the analysis also
reports the p-value version, with a randomized boundary so the null rate is exactly 5%.

**Oracle reference.** The analysis also size-adjusts the bare statistic (n·HSIC_V, n·dCov_V,
\|Pearson r\|, MMD²). This is the best a rule built on that statistic can do in that cell.
The gap between a rule and the oracle measures what the critical-value rule costs.

**Statement for the paper.** Size adjustment uses the true null of a simulated cell, which an
analyst never has. It answers "is the statistic-plus-rule informative at all?". It is not a
procedure that can be applied to data.

## Uncertainty

- Wilson 95% interval on the alternative rejection count, given c* (counting error only).
- Bootstrap 95% interval that also includes threshold error: resample null and alternative
  replications independently (2,000 resamples, seed 42), recompute c* and power.
- Rule-versus-rule gaps use the same resamples, so the pairing of rules on one replication is
  kept.

## Alternatives

ξ, w are independent unit-variance fGn(H) (Davies–Harte, `mbcore.fgn`).

| id | arm | equation | expectation |
|---|---|---|---|
| lin | independence, G0 | X = ξ, Y = ρ ξ + √(1−ρ²) w | Pearson most powerful; kernels close |
| sq | independence, G0 | X = ξ, Y = ρ (ξ² − 1)/√2 + √(1−ρ²) w | Corr(X,Y) = 0: Pearson ≈ 5%, kernels have power |
| lin-SEAS-V | independence, SEAS-V | X = s_X(t) ξ, Y = s_Y(t)(ρ ξ + √(1−ρ²) w), s(t) = 1 + 0.5 sin(2π(t+U)/50), U_X ⊥ U_Y | tests shift where it was the only valid rule |
| shift | split MMD, G0 | x_t = ξ_t + δ·1{t ≥ n/2} | mean change in the second half |

Size adjustment matches each alternative to its null cell: lin and sq → G0 independence null;
lin-SEAS-V → SEAS-V independence null; shift → G0 split-MMD null.

Effect sizes ρ_lin, ρ_sq, δ are chosen at n = 2000 so that size-adjusted HSIC wild (MMD
uncentred wild for δ) power at G0, H = 0.5 is near 50%. The pilot picks them (amendment A1).
They are then **fixed across H and generator**, so power at H = 0.7, 0.8 shows the cost of memory.
At n = 500 the effect is doubled (local alternatives, effect ∝ n^{−1/2}), so n = 500 power is
on a similar scale rather than near 5%.

## Grid and runs

- Null cells, fresh, in the EXP3 code path: independence G0 and SEAS-V; split MMD G0.
- Alternative cells: lin, sq, lin-SEAS-V, shift.
- H ∈ {0.5, 0.7, 0.8}, n ∈ {2000, 500}, 1,000 replications per cell. 42 cells.
- Why fresh nulls rather than the EXP2 rows: EXP2 rows store p-values but not reference
  quantiles, so the ratio is unavailable. The EXP2 rows serve as a regression check: fresh raw
  null rates must match EXP2 rates within Monte Carlo error. The pilot uses the EXP2 G0 null
  p-values only to pick ρ, δ.
- Seeds: SeedSequence(entropy=2026102542, spawn_key=(cell_idx, rep)) main; 2026102541 pilot.
- B = 500 wild draws. At most 10 worker processes (shared machine). Compute cap 3 h wall.
- Pilot: 100 reps per pilot cell, G0, H = 0.5, n = 2000, three candidate values per effect,
  plus one SEAS-V and one n = 500 cell for timing.

## Code

`power_core.py` (rule bundles with ratios and EXP3 generators), `run_exp3.py` (driver),
`analyze.py` (size adjustment, intervals, tables), `unit_checks.py`. `gens.py`, `mbcore.py`,
`mmdcal.py`, `common.py` are unmodified copies of the EXP2 files (hashes in
`copied_from_EXP2.sha256`). Unit check: `kind_ext` reproduces every `kind_tests` and
`lind_tests` field on the same seed to 0 difference.

## Bias and leakage checklist (simulation version)

- Look-ahead: not applicable. Every statistic uses the full simulated path by design.
- Survivorship, universe selection, transaction costs: not applicable (no market data).
- Cherry-picking: the generators, rules and reading rule are fixed above. Effect sizes are
  picked on G0, H = 0.5 only, before any other cell is seen.
- Threshold reuse: c* is estimated on null replications that are independent of the
  alternative replications (different seeds), so the size-adjusted power is not optimistic
  by reuse.
- Skill registry paths named in the agent protocol (/.claude/skills/ml_stats, /scripts/...) do
  not exist on this machine. No registry skill was used.

## Amendments

(A1 is added after the pilot, before the main run.)

### A1 (2026-10-25, after the pilot, before the main run): effect sizes and timing

Pilot: 100 reps per cell, G0, H = 0.5, n = 2000, seed master 2026102541 (`pilot/`, `pilot.log`).
Size adjustment in the pilot used the EXP2 main G0 null p-values (1,000 reps), 5% quantile.
Size-adjusted power of HSIC wild (MMD uncentred wild for δ):

| effect | value | HSIC wild adj. power | note |
|---|---|---|---|
| ρ_lin | 0.04 / 0.06 / 0.08 | 22% / 42% / 65% | Pearson NW 33% / 75% / 96% |
| ρ_sq | 0.06 / 0.09 / 0.12 | 23% / 69% / 90% | Pearson NW 2–6% (no power, as designed) |
| δ | 0.08 / 0.12 / 0.16 | 29% / 65% / 91% | MMD uncentred wild |

Linear interpolation to 50% gives ρ_lin ≈ 0.067, ρ_sq ≈ 0.077, δ ≈ 0.10.
**Fixed: ρ_lin = 0.065, ρ_sq = 0.075, δ = 0.10 at n = 2000; 0.13, 0.15, 0.20 at n = 500.**

Timing under the shared load: 2.2–2.6 CPU-s per independence rep at n = 2000 (EXP2 logged
1.3 s on an idle machine), 0.11 s at n = 500, 0.8 s per MMD rep at n = 2000. Projected main
run: about 12.5 CPU-h, about 75 min wall with 10 workers. Inside the 3 h cap, so the full grid
including n = 500 runs.

Pilot sanity reading (not evidence): under SEAS-V with ρ = 0.06, raw HSIC shift rejects 13%
while the oversized rules reject about 80%. Whether shift keeps power after size adjustment
is the main-run question.
