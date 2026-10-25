# EXP5 results (2026-10-25)

Simulation only. Exploratory. Design, hypotheses and amendments A1–A3 are in `DESIGN.md`.
Full tables are in `tables_main.md`. Per-cell numbers are in `summary_main.csv` (Part 1),
`summary_ind_main.csv` (Part 2) and `results_main.json` (reading-rule verdicts).

## Answer 1: the MMD power collapse is an information effect

The shift was scaled to a fixed number of null SDs of the half-mean difference.
At that scale, split-MMD power does not collapse with memory.
It stays level or rises.
No rule and no c meets the "collapse" reading rule in any cell.

At c = √5 and n = 2000, the oracle MMD detects 47.6% at H = 0.5 and 60.5% at H = 0.8.
The wild and subsampling rules move the same way (48% to 58–59%).
The self-normalised block rule is flat at 45%.
The linear benchmark is flat as designed: 66%, 56%, 59% (exact value 60.9% at every H).
At c = 1.5 the picture is the same at lower power: MMD rules 24–26% at H = 0.5 and 26–34% at H = 0.8.
n = 500 repeats the pattern.

So the EXP3 collapse came from the fixed shift δ = 0.1.
That shift is 2.2 null SDs without memory and 0.4 null SDs at H = 0.8.
It was not a weakness of the kernel or of the critical-value rules.

Two details matter for the paper.

1. Without memory, MMD loses ground to a mean test. At H = 0.5 it is about 18 points below the
   linear benchmark (48% against 66%). Under strong memory it matches the benchmark (about 59–61%).
   One reading is that long memory makes the rank-one (mean) part dominate the MMD statistic.
   This reading is not tested here.
2. The self-normalised rule costs power under memory. At H = 0.8 it is about 15 points below the
   oracle (45% against 60.5%). The wild and subsampling rules are within 2 points of the oracle.

Caveats. These are size-adjusted numbers. At H = 0.8 the raw rules are badly oversized
(wild 27%, subsampling 50% false alarms), so their raw power means little. At H = 0.8 the
scaled shift is large in marginal units (δ = 0.55 at c = √5). The H = 0.5 linear-benchmark
value 66.1% is above its exact value 60.9%. That gap is threshold noise from 1,000 null
replications, and the bootstrap interval [61.0, 72.5] covers the exact value.

Regression check. The fresh MMD nulls match the EXP3 nulls on raw size (all |z| < 1.5)
and on the size-adjustment cut-off.

## Answer 2: sum 1.4 stays flat, and sum 1.6 stops falling

**Sum 1.4.** The rate is flat from n = 500 to n = 4000 for every rule.
No change from n = 2000 to n = 4000 is detectable (all |z| ≤ 1.04).
The excess of the subsampling rules is still there at n = 4000.
HSIC subsampling is 8.2% and dCov subsampling is 9.4%.
Pearson NW is 7.4%. HSIC wild is 7.3%, against its own 5.8% at H = 0.5.

**Sum 1.6.** The fall from n = 500 to n = 2000 does not continue.
Every rule rejects more at n = 4000 than at n = 2000, by 0.6 to 4.0 points.
For dCov subsampling the rise is clear (22.9% to 26.9%, z = 2.34).
HSIC subsampling rises from 19.7% to 22.7% (z = 1.91).
Over n = 500 to 8000 the fitted trend is close to zero for every rule.
The paper's Δ column (n = 500 to 2000) should therefore not be read as a decline with n.
This is bad news for that reading, and good news for the theory: above the covariance balance,
the false-alarm rate does not shrink with n on this grid.

**n = 8000.** Only 144 replications per cell (amendment A3). The intervals are about ±6 points
near 20%, so these rows are low precision. At sum 1.6 the point estimates are higher again
(HSIC sub 25.0%, dCov sub 29.9%). At sum 1.4 several rates are lower (HSIC wild 2.1%), but the
intervals overlap the n = 4000 values.

**H = 0.5 baseline at n = 4000.** 500 replications (cut as planned). HSIC shift 4.8%,
HSIC sub 5.4%, HSIC wild 5.8%. The wild-rule oversize without memory keeps shrinking with n.

## Compute

8.5 CPU-h in total including the pilot (Part 1 1.7 h, Part 2 6.7 h). Wall time 18:52–21:22 UTC.
At most 3 workers while EXP4 ran, then 6. 0 errors. Unit checks pass (`unit_checks.txt`).
Peak memory 1.1 GB per worker at n = 4000 and 4.0 GB at n = 8000.

## Amendments in short

- A1: n = 8000 dropped from the main plan. H = 0.5 baseline at n = 4000 cut to 500 replications.
- A2: the run was split into two parallel scripts after EXP4 finished. The EXP4 check in the
  first script also matched other agents' shell watchers, so it would have stayed at 3 workers.
  Cells, seeds and targets did not change.
- A3: n = 8000 at sums 1.6 and 1.4, 144 replications each, added because time remained.

## Not applicable

The agent protocol asks for Sharpe, t-stat, drawdown and information ratio.
This is a size-and-power simulation with no returns, so those metrics do not exist here.
The skill registry paths named in the protocol do not exist on this machine.
