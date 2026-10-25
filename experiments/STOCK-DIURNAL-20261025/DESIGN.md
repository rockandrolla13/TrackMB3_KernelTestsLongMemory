# STOCK-DIURNAL: diurnal adjustment control for the stock-duration placebo

Written 2026-10-25, before any run.

## Question

A reviewer notes that the stock arm removes only a 15-minute clock-bin mean of log durations.
EXP4 shows that a shared seasonal scale misleads HSIC, and that a slow scale drift cannot be
removed per window. Do the Table 7 placebo fractions move toward 5% under a standard
multiplicative diurnal factor (Engle and Russell 1998) and a diurnal scale adjustment?

## Hypotheses

- H0: the rejection fractions under D1 and D2 equal those under D0, up to pair-level noise.
- H1: D1 or D2 lowers the HSIC wild, HSIC subsampling, and dCov subsampling fractions toward 5%.

## Adjustments

All smooth functions are fitted per stock on the pooled six-month sample, as the original bin
mean was. The time variable is seconds since 09:30 New York time at the event that ends the
duration. This matches the original bin assignment, which used the end event's clock minute.

- D0: the original cache value, log duration in seconds minus the pooled 15-minute bin mean.
  It is read unchanged from `data/raw/stock_duration_cache_2024/`. It is the regression check.
- D1: Engle-Russell multiplicative factor. Fit phi(t) = E[duration | t] by least squares on a
  cubic B-spline basis in t. Interior knots: 15, 30, 60, 120, 210, 300, 360 minutes after the open
  (7 knots, denser near the open where the pattern is steepest). Value = log(duration / phi(t)).
  We keep the paper's log scale.
- D2: D1 plus a diurnal scale. With r = log(duration / phi(t)), fit a smooth mean m(t) of r and a
  smooth variance v(t) of (r - m(t)) on the same basis. Value = (r - m(t)) / sqrt(v(t)).
- D3: D2 minus each day's mean. Every window stays inside one day. HSIC with the median bandwidth,
  distance covariance, and Pearson are all invariant to adding a constant to a whole window.
  D3 is therefore an exact no-op on the tests. We verify that on a 20-pair spot check and do not
  run it in full.

The fit uses per-second aggregates (sums and counts on the 23,400-second grid). The basis is then
evaluated at the exact event time. Aggregation only rounds the regressor to one second.

## Invariance note (stated before running)

The three tests are invariant to one affine map applied to a whole window. An adjustment
therefore changes a window only through the variation of phi, m, and v inside that window.
A 2,000-duration window spans a few minutes. We expect small changes, except for windows
near the open and close. A null result here would mean that pooled diurnal adjustment cannot
address within-window scale drift of the SEAS-V-slow kind.

## Placebo runs

- Same 12 cells: AAPL-MSFT, AAPL-NVDA, MSFT-NVDA, AAPL-AAPL x n in {500, 1000, 2000}, 500 pairs each.
- Same seed (2026102301), same SeedSequence spawn keys (cell index, pair index), same day and
  start draws, and same wild multipliers. Window draws depend only on day boundaries, so the
  windows are identical across adjustments. The rng is re-seeded for each adjustment.
- Same `mbcore.kind_tests` (B = 500 wild draws) and `mbcore.lind_tests` from experiments/EXP1-20261016.
- Regression check: D0 must reproduce `experiments/STOCK-20261023/results/*.csv` exactly in
  day, start, and every p-value and rejection flag.
- Workers: at most 2 while `run_exp4.py` is running, else at most 4.

## Diagnostics

1. Diurnal scale range: for each adjusted series, fit a smooth sd curve across the day and report
   max/min over the session. Also report max/min of the 26 clock-bin sds (model-free).
2. Within-window variance trend: for every n = 2000 window used in the placebo, split the window
   into 8 equal segments and regress log segment variance on segment position. Report the median
   |slope| and the share of windows whose |slope| exceeds the 95th percentile from the same window
   randomly permuted (5% if there is no trend). Also report max/min segment-variance ratio.

## Outputs

Code, this file, `results/<adj>/<pair>_n<n>.csv` (derived statistics only), `summary.csv`,
`diagnostics.json`, `tables.md`, `RESULTS.md`. The adjusted caches go to
`data/raw/stock_diurnal_cache_2024/`, which Git ignores. Nothing is committed.

## Caveats

Process independence is unknown and days recur across pairs. Fractions are descriptive placebo
rates, as in the paper. No transaction costs or trading metrics apply: this is a test-size study.

## Amendment 1 (2026-10-25, during the run, approved by the team lead)

The machine is shared and at most 2 workers are allowed, so the full run would take about 4-5 hours.
D0 now reruns only the four n = 2000 cells (the Table 7 cells). The D0 pilot (8 pairs at n = 500)
already matched the original CSVs in every field. D1 and D2 still run all 12 cells. The run was
stopped once and resumed. The resume logic skips completed pairs, and the rng is re-seeded per
pair, so the restart does not change any result. `run_all.sh` holds the two commands.
