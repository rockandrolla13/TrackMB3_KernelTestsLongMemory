# EXP8 (2026-11-08): memory kernels and a surrogate null on stock durations

Status: design only. Written before any EXP8 code or data run. No EXP8 rate has been computed.

## Questions

1. **Arm A (kernels on the stock placebo).** Do the EXP7 memory kernels change the stock placebo fractions? Do they remove the excess linked to drift inside a window? In the paper, HSIC subsampling rejects 11.0% of pairs with a variance trend in both windows and 7.5% of pairs with no trend.
2. **Arm B (surrogate null).** Each pair is a real window and a surrogate built with fresh randomness. In these pairs, how often does each rule reject? A surrogate is a randomised copy that keeps the source window's values and spectrum. This gives a false-alarm rate on realistic marginals and memory. The placebo cannot give one.

## Data

- **Source:** the local cache `data/raw/stock_duration_cache_2024/{AAPL,MSFT,NVDA}.npz`, covering AAPL, MSFT and NVDA from January to June 2024, with 124 days each.
- **Series:** `dur`, the log duration minus the stock's pooled clock-bin mean. This is the same series as the paper's main stock arm.
- **What stays local:** no download and no API call. Raw values never enter git, and only derived per-pair statistics are saved.

## Windows

Arm A and Arm B use the STOCK-20261023 sampler unchanged:
- one window per side, inside one day;
- paired days at least 20 trading sessions apart;
- pairings AAPL–MSFT, AAPL–NVDA, MSFT–NVDA and AAPL–AAPL.

Arm A uses the original seed `2026102301` and the original draw order. Its windows are therefore exactly the windows of the saved STOCK rows.

## Arm B surrogates

In each pair, X is the real window of stock a. Y is a surrogate of the real window of stock b.

- **S1, IAAFT** (iterated amplitude-adjusted Fourier transform, 100 iterations, then a final rank remap):
  - keeps Y's exact values and approximately its periodogram, and so its memory;
  - draws fresh random phases, so it cancels any time alignment with X.
- **S2, IAAFT plus drift:**
  - fit Y's variance trend first, using log segment variance on 8 segments, linearly in position;
  - divide Y by the fitted local scale, apply S1, then multiply the fitted scale back.
  - S2 keeps Y's slow scale drift and removes everything else that lines up with X. It tests whether drift alone, with independence, makes the rules reject.

**What "known null" means here.** Given the two source windows, the surrogate's phases are independent of X. Window-level features of Y, such as its marginal, spectrum and drift slope, still come from real data. A rejection is therefore a false alarm with respect to time alignment. That is the dependence HSIC is meant to detect.

## Cells

| Arm | Pairings | n | Pairs per cell | Cells |
|---|---|---|---|---|
| A, real–real placebo | 4 | 500, 2000 | 500 (as STOCK) | 8 |
| B, real–surrogate S1 | 4 | 500, 2000 | 1000 | 8 |
| B, real–surrogate S2 | 4 | 500, 2000 | 1000 | 8 |

Arm B seeds use master `2026110801`. Each pair gets `SeedSequence(master, spawn_key=(cell, pair))`. The cell list is append-only.

## Statistics and rules

- **Function:** `exp7_core.ind_bundle(x, y, rng, kernels=K, full=True, hsplit=True, strict=False)`.
- **Kernels:** `K` = raw, inc1, inc50, herm1, herm2, sym, datadiff, hp10 and dcov_inc1.
  - frac_oracle is dropped because it needs the true H.
- **Rules:** wild, circular shift and subsampling for raw and inc1, plus Harris for raw and inc1. hsplit1 and hsplit2 use both the plug-in rule and the block rule. dCov and Pearson are kept as baselines.
- **Per window, saved for the analysis:**
  - the local Whittle estimate of H, with bandwidth m = ⌊n^0.65⌋;
  - the variance-trend slope and its trend flag. The flag uses the rule in `experiments/STOCK-DIURNAL-20261025/diagnostics.py`:
    - fit the OLS slope of log variance on 8 segments;
    - permute each window once;
    - set the cut at the 95th percentile of the pooled permuted slopes in the cell.

## Predictions, written before any data

- **P0 (integrity).** Arm A reproduces the saved STOCK raw shift p-value and Pearson columns exactly. Those columns are deterministic given the window, so a mismatch means the windows differ. If P0 fails, Arm A stops.
- **P1 (S1, shift).** The raw and inc1 shift rules reject between 3.4% and 6.8% in every S1 cell at n = 2000. That is the Wilson band around 5% for 1000 pairs, widened by 0.4 points.
- **P2 (S1, wild and subsampling).** Raw wild and raw subsampling reject more than 5% at n = 2000, pooled over pairings. The size of the excess should rise with the mean estimated H of the pair.
- **P3 (S2 against S1).** For raw wild and raw subsampling, S2 rejects more than S1 at n = 2000, pooled over pairings, and the paired 95% interval of the difference excludes 0. The shift rules stay within the P1 band under S2.
- **P4 (Arm A, kernels).** At n = 2000, inc1 wild rejects less than raw wild, pooled over pairings. The paired McNemar test runs on the same windows. For inc1 subsampling, the gap between trend–trend and no-trend pairs is smaller than the raw gap of 3.5 points.
- **P5 (Harris).** Harris on raw and inc1 rejects at most 5% in every Arm B cell. The test is Wilson upper bound ≥ observed, read one-sided.

## Reading rules

- **Arm B fractions** are false-alarm rates with respect to time alignment, with the caveat above. Arm A fractions stay descriptive placebo fractions, as in the paper.
- **Pooling:** over pairings only. The n = 500 and n = 2000 cells are never pooled together.
- **Multiplicity:** there are five predictions and each gets a single pass or fail. No other comparison is called a finding. Anything else is labelled exploratory.
- **Wording:** a failed prediction is reported as failed, without rewording.

## Bias and leakage checklist

- The surrogates use only the window they replace, so no other window's data leaks in.
- The trend fit for S2 uses Y only.
- The windows, days and seeds are saved per pair, so any row can be rebuilt.
- The predictions and thresholds above are fixed before the run. Any change goes in a dated amendment below, stating what had been seen at that point.

## Compute

EXP7 took about 2.0 CPU-seconds per pair at n = 2000 and 0.27 at n = 500.
- There are 10,000 pairs in total.
- The estimate is about 3.5 CPU-hours, or about 15 minutes on 16 workers. This is rough, because surrogates and H estimates add time.
- The laptop must stay open.

## Code (to be written after approval)

- `experiments/EXP8-stock-surrogate-20261108/run_exp8.py`: the driver, with `run_chunk` and `job_for`.
- `spec.json`: run with `python -m kmem run EXP8`.
- `kmem/surrogates.py`: IAAFT and the drift variant, with unit tests:
  - the surrogate has Y's exact sorted values;
  - the periodogram correlation is above 0.95;
  - on simulated fGn, S1 at H = 0.8 gives shift-rule size within the P1 band, as a 1,000-rep check run before the real data.

## Amendments

### A1 (2026-11-08 02:05 UTC, after a 45-pair smoke run, before any rate was read)

- **P0, made precise.** Arm A must reproduce the saved STOCK days and window starts exactly. Its raw shift p-value must equal the saved value to 12 significant digits. Its Pearson r must agree to within 1e-12. Pearson differs at about 1e-14 because the two code paths do the floating-point sums in a different order. The smoke run passed this check on 6 of 6 pairs.
- **n = 500 trend fit.** 500 is not a multiple of 8, so the variance-trend fit drops the last 4 points. The S2 surrogate uses the same fit. At n = 2000 nothing changes. The trend-linked part of P4 is read at n = 2000 only, as in the paper.
- **Compute.** The smoke run measured 0.85 CPU-seconds per pair at n = 2000, so the run is shorter than estimated.

### A2 (2026-11-08 02:10 UTC, main run in progress, no main-run rate read): exact pass rules

The design left three pass rules loosely worded. Before any rate was read, they are fixed as follows.

- **Rejection definitions.** Each rule rejects exactly as in EXP7: `analyze_exp7.ind_rules` and `kmem.calibration.registry`. For example, wild rejects when p_wb ≤ 0.05, and subsampling rejects when rej_sub = 1.
- **P2.** "More than 5%" means the pooled Wilson 95% lower bound exceeds 5%, for raw wild and for raw subsampling. The H part is a logistic regression of rejection on the mean of H_x and H_y. It passes when the slope is positive with a one-sided p below 0.05. P2 passes only if all three parts pass.
- **P3.** S1 and S2 use different seeds, so their pairs are not paired. The interval is the Newcombe interval for a difference of two independent proportions. The shift part uses the P1 band in every S2 cell at n = 2000.
- **P4.** The McNemar test is paired on the same windows, using the exact binomial form from `kmem.analysis.mcnemar_p`. The test is one-sided in the predicted direction, at p < 0.05.
  - **Trend part, flags.** The trend flag is computed per cell: |slope| > the 95th percentile of the cell's pooled |slope_perm| over its x and y windows.
  - **Trend part, gap.** For each rule, the gap is the rate in trend–trend pairs minus the rate in no-trend pairs. The comparison uses this run's raw subsampling gap, pooled over pairings. The paper's 3.5 points is reported alongside it.
- **P5.** "At most 5%" means the cell's Wilson 95% lower bound is at or below 5%. That is, the rate is not significantly above 5%. The rule applies to raw Harris (rhc) and inc1 Harris in all 16 Arm B cells.
