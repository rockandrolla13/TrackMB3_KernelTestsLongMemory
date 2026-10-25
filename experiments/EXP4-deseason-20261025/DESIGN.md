# EXP4 (2026-10-25) — does seasonal adjustment remove the EXP2 seasonal false alarms?

Written before the pilot and the main run. Simulation only. No market data.
Exploratory: rates are false-alarm rates of the tested implementations under the stated
generators and adjustments, not validity results for any rule.

## Question

EXP2 found that two independent series sharing a seasonal pattern mislead single-path tests.
SEAS-M (seasonal mean) gave about 91% rejection for HSIC, dCov and Pearson at H = 0.5, n = 2000.
SEAS-V (seasonal scale) gave HSIC wild 73.3%, HSIC sub 71.5%, dCov sub 45.6%, Pearson-NW 4.7%.
A reviewer asks whether this survives the standard first step of deseasonalising each series.

Sub-questions, fixed now:

- (a) Does mean adjustment alone (A1) fix SEAS-M?
- (b) Does SEAS-V need scale adjustment (A2)? This is the key question.
- (c) Does removing 50 bin means distort size under long memory (G0, H = 0.8)?
- (d) What does this imply for the paper's stock arm, which removed a pooled 15-minute
  clock-bin mean of log durations but not a scale?

## Hypotheses (per generator g, rule r, n, H)

- H0_fix: rate(g, A_k) = rate(G0, A_k) at the same H and n. Adjustment restores the G0 rate.
  H1_fix: rate(g, A_k) > rate(G0, A_k).
- H0_size: rate(G0, A_k) = rate(G0, A0). The adjustment itself does not change size.
  H1_size: the rates differ.

## Generators (copied unchanged from EXP2, hashes in `copied_from_EXP2.sha256`)

`gens.pair(gen, H, n, rng)`. Period P = 50, a = 1, b = 0.5. ξ, η independent unit fGn(H).

| id | series |
|---|---|
| G0 | X = ξ |
| SEAS-M | X = ξ + sin(2π(t + U)/50), U ~ Uniform{0..49}, independent per series |
| SEAS-V | X = (1 + 0.5 sin(2π(t + U)/50)) ξ, U independent per series |
| SEAS-M-shared | X = ξ + sin(2πt/50), fixed common phase (diagnostic, as in EXP2) |

All four are true process-independence nulls. SEAS-M-shared is not stationary (see EXP2 DESIGN).

## Adjustments (`deseason.py`), each estimated from the one observed series only

Bin b(t) = t mod P, counted from the first observation of the series.

| id | step |
|---|---|
| A0 | none (the EXP2 rate, re-run on the same draws for the paired comparison) |
| A1 | X̃_t = X_t − mean{X_s : b(s) = b(t)}, P = 50 known |
| A2 | A1, then divide by the bin sd of the A1 residuals (ddof = 1), P = 50 known |
| A3 | estimate P̂, then A2 with P̂ |

A1 is the per-path analogue of the clock-bin mean used in the paper's stock arm.
Bins have n/P = 10 (n = 500) or 40 (n = 2000) observations.

**A3 period rule.** Periodograms of x − x̄ and of (x − x̄)² − mean at Fourier indices
j = 10 … ⌊n/2⌋ (periods n/j from n/10 down to 2). Each is divided by its own sum over that range.
Take the largest normalised ordinate across the two, P̂ = round(n / j*). The squared series is
included because a seasonal scale leaves no peak in the periodogram of x itself. Recovery rate
(P̂ = 50) is recorded per series, with which periodogram won.

**Shortcut, exact by construction.** When P̂_X = P̂_Y = 50, A3 produces the same series as A2.
The test then reuses the A2 output (flag `reused = 1`). This is exact because every adjustment of
one replication uses the same data and the same test multipliers (seeds below).

## Tests (unchanged from EXP2: `mbcore.kind_tests(x, y, rng, B=500)` and `mbcore.lind_tests`)

HSIC wild, HSIC shift, HSIC sub, dCov wild, dCov sub, Pearson NW, Pearson EWC. Rules as in
EXP2 `analyze.py`. All tests are invariant to a positive rescaling of a whole series, so A2 acts
only through the bin-to-bin change of scale.

## Grid

4 generators × H ∈ {0.5, 0.7, 0.8} × n ∈ {500, 2000}, 1,000 replications per cell.
Each replication runs every adjustment on the same pair.
Planned adjustments: A0–A3 at n = 500. At n = 2000, A1–A3 only. The unadjusted n = 2000 rates
come from EXP2 `summary_main.csv` (independent seeds, same code). This saves about 9,000
test calls of about 1.2 s each.

## Seeds

Masters: pilot 2026102541, main 2026102542. Replication r of cell c draws the pair from
`SeedSequence(entropy=master, spawn_key=(c, r, 0))`. Each adjustment's tests use a fresh
`default_rng(SeedSequence(entropy=master, spawn_key=(c, r, 1)))`. So all adjustments of one
replication share data and multipliers (common random numbers), and paired differences are sharp.
Cell indices come from the order in `run_exp4.py` and are written to each manifest.

## Reading rule (fixed before the pilot)

Rates are rejection fractions at nominal 5% with Wilson 95% intervals (`common.wilson`).

1. **Fixed:** an adjustment fixes generator g for rule r at (H, n) if rate(g, A_k) − rate(G0, A_k)
   ≤ 2 pp, or their Wilson intervals overlap.
2. **Not fixed:** rate(g, A_k) − rate(G0, A_k) > 2 pp and the intervals do not overlap.
3. **Size distortion from adjustment (c):** compare rate(G0, A_k) with rate(G0, A0). At n = 500
   use the paired A0 on the same draws (McNemar-style discordant counts are reported). At n = 2000
   compare with EXP2 G0 (independent samples, two-sample z on proportions). A change > 2 pp with
   |z| > 1.96 is flagged.
4. Read at n = 2000. n = 500 alongside. The paper table uses H = 0.5 and 0.8.
5. SEAS-M-shared is a diagnostic. It is not used in rules 1–3.
6. Suspicion check: any adjusted G0 cell more than 2 pp below its A0 rate is flagged, as is any
   seasonal cell after adjustment more than 2 pp below the adjusted G0 rate.

## Compute estimate and budget gate

Measured 1.16 CPU-s per test call at n = 2000 and 0.10 s at n = 500 on an idle core
(EXP2 rows showed about 2.4 s at n = 2000 under load).
Test calls per replication at n = 2000: G0 3 (A3 never reuses), SEAS-M 2, SEAS-V 2,
SEAS-M-shared 2 (A3 reuses A2 when P̂ = 50 for both). Over 3 H × 1,000 reps: 27,000 calls,
about 9–18 CPU-h. At n = 500: about 4 calls × 12 cells × 1,000 = 48,000 calls, about 1.3 CPU-h.

Budget: 5 worker processes, 2 h wall, so about 10 CPU-h. The pilot (100 reps) measures the real
rate. Cut order if the main run would exceed about 1.6 h wall:
(1) drop SEAS-M-shared at n = 2000; (2) drop H = 0.7 at n = 2000; (3) reduce n = 2000 cells to
500 replications. Each cut is recorded below as an amendment before the main run.

## Amendments

(1) 2026-10-25, after the pilot, before the main run. **Exact equivalence.** The pilot showed, and
`unit_checks.py` confirms (`unit_checks.txt`, max difference about 4e-14), that the correct
known-period adjustment returns exactly the adjusted fGn core:
A1(ξ + s) = A1(ξ) and A2(ξ + s) = A2(ξ) for SEAS-M and SEAS-M-shared, and A2((1 + b s) ξ) = A2(ξ)
for SEAS-V. The reason is that s and 1 + b s are constant inside every bin of the true period.
So these adjusted rows equal the adjusted G0 rows **in law**. They need no separate simulation to
answer (a) and the A2 part of (b). Only three things carry new information: the adjusted G0 rows
(question c), SEAS-V under A1 (the wrong adjustment, question b), and A3 (whether P̂ = 50; pilot:
recovered in 100% of seasonal series at n = 2000).

(2) 2026-10-25, budget cut, before the main run. The pilot measured 1.9–2.6 CPU-s per test call at
n = 2000 under shared load, twice the idle estimate. The full n = 2000 grid would take about 3.4 h
wall at 5 workers. Cuts, in order of the gate plus the equivalence above:
- n = 2000: run G0 (A1, A2, A3) and SEAS-V (A1, A2, A3) at H ∈ {0.5, 0.8}, 1,000 reps.
  Drop H = 0.7 at n = 2000 (cut 2). Drop SEAS-M and SEAS-M-shared at n = 2000 (cut 1, and their
  adjusted rows equal adjusted G0 by (1)). Their pilot rows (100 reps) remain as a check.
- n = 500: full grid (4 generators × 3 H, A0–A3), 1,000 reps.

(3) 2026-10-25, added diagnostic, before the main run. **SEAS-V-slow** (`gens_exp4.py`): seasonal
scale with period 4n, so one window sees a quarter cycle. Independent random phase per series.
This mimics the stock arm, where a 2,000-duration window covers 7–15 minutes of a 390-minute day
and the windows start at independent random times. A per-window bin adjustment is impossible
there, so only A0 is run. Grid: H ∈ {0.5, 0.8}; n = 500 at 1,000 reps, n = 2000 at 300 reps.
Cell indices 24–29 are appended after the pilot cells, so pilot seeds are unchanged.
Read: if A0 on SEAS-V-slow rejects well above G0, a slow scale drift alone can mislead the
kernel tests, and the stock arm's missing scale adjustment matters more than the period-50 result
suggests.

(4) 2026-10-25, 14:10, during the main run. Machine load reached 47 on 22 cores (another agent's
run uses 10 workers). Each worker got about 60% of a core, and stage 1 projected to about 110 min
wall. Applied cut (3): the four n = 2000 cells run 500 replications instead of 1,000. The run was
stopped after 322 complete G0 H = 0.5 replications (no partial rows) and resumed with --reps 500.
Monte Carlo SE at n = 2000 is therefore about 1.0 pp at 5% and 2.2 pp at 50%.
