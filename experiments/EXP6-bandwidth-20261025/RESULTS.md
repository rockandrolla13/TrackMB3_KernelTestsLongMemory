# EXP6 results — bandwidth sensitivity of the HSIC rules

Simulation only. n = 2000, 500 replications per cell, seed master 2026102562. Design and
amendment A1 are in DESIGN.md. Full tables are in tables_main.md. Rates are in percent.

## The answers

**1. The HSIC weakness under the seasonal scale is a bandwidth effect.** A wider kernel removes it.
Size-adjusted power of HSIC wild rises from 16% at the median bandwidth (κ = 1) to 43% at κ = 2
and 73% at κ = 4. dCov wild reaches 48% and Pearson NW 75% in the same cell. The oracle row moves
the same way (14%, 45%, 70%). So the weakness sits in the statistic, not in the critical-value rule.

**2. It is not a fault of the median heuristic as such.** The median bandwidth is not inflated by
the seasonal scale. Its median is 0.94 under SEAS-V and 0.95 under G0. The cause is what a
moderate bandwidth makes the kernel sensitive to (see Mechanism).

**3. No single bandwidth wins.** The bandwidth that rescues the seasonal cell destroys power
against the squared (non-linear) alternative. There HSIC wild falls from 37% at κ = 1 to 7% at
κ = 4. At κ = 8 the oracle has 4%, which is Pearson's level. A wide Gaussian kernel behaves like a
squared Pearson correlation and loses what makes it a kernel test.

**4. Circular-shift HSIC keeps its size at every bandwidth.** It is valid in all three null cells
at all four κ (4.0% to 5.8%). The other rules change with κ. The details follow.

## Size (raw null rejection, 500 reps, Wilson 95% about ±2 pp at 5%)

| null cell | rule | κ = 0.5 | κ = 1 | κ = 2 | κ = 4 |
|---|---|---|---|---|---|
| G0, H = 0.5 | wild | 7.2 | 7.8 | 5.0 | 5.4 |
| | shift | 5.4 | 5.4 | 4.4 | 5.0 |
| | sub | 6.0 | 6.6 | 5.8 | 5.4 |
| G0, H = 0.8 | wild | 6.4 | 8.0 | 8.6 | 8.6 |
| | shift | 4.0 | 4.8 | 5.0 | 5.2 |
| | sub | 12.2 | 19.2 | 24.2 | 26.4 |
| SEAS-V, H = 0.5 | wild | 72.0 | 74.4 | 56.8 | 6.6 |
| | shift | 5.0 | 5.6 | 5.8 | 4.8 |
| | sub | 71.8 | 73.6 | 52.0 | 6.8 |

- On the primary grid κ ∈ {0.5, 1, 2}, the picture is stable in kind. Shift is valid everywhere.
  Subsampling over-rejects under long memory. Wild and subsampling over-reject heavily under SEAS-V.
- Two changes are worth reporting. Subsampling at H = 0.8 gets worse as κ grows (12% to 26%,
  McNemar p < 0.001 for each κ ≠ 1). Under SEAS-V, wild and subsampling fall to 57% and 52% at κ = 2.
- At κ = 4 the SEAS-V over-rejection of wild and subsampling disappears (6.6% and 6.8%). This is
  the same mechanism as the power recovery. It is not a sign that these rules are fixed in general.
  At κ = 4 they lose the non-linear alternative.

## Size-adjusted power (c* from the matched null, bootstrap 95% in tables_main.md)

| alternative | rule | κ = 0.5 | κ = 1 | κ = 2 | κ = 4 | dCov wild | Pearson NW |
|---|---|---|---|---|---|---|---|
| G0 linear | wild | 23.6 | 53.8 | 80.4 | 81.8 | 77.4 | 81.4 |
| | shift | 22.8 | 55.0 | 79.0 | 81.8 | | |
| | oracle | 22.8 | 55.0 | 78.8 | 81.2 | | |
| G0 squared | wild | 22.0 | 37.2 | 41.0 | 7.2 | 35.4 | 2.8 |
| | shift | 20.0 | 39.8 | 36.4 | 6.8 | | |
| | oracle | 20.4 | 39.0 | 37.4 | 6.8 | | |
| SEAS-V linear | wild | 10.2 | 15.6 | 42.6 | 72.8 | 48.2 | 75.2 |
| | shift | 10.2 | 15.0 | 41.2 | 69.0 | | |
| | oracle | 10.2 | 13.8 | 44.6 | 70.0 | | |

- Every κ = 2 and κ = 4 gain in the SEAS-V cell is 26 to 60 pp, and every paired bootstrap
  interval excludes 0. By the pre-registered rule this is **full recovery** at κ = 4 for all HSIC
  rules, and for HSIC sub and the oracle already at κ = 2.
- κ = 0.5 loses power in every alternative cell.
- The median bandwidth is also far from best for the linear G0 alternative (55% against 80% at
  κ = 2). It is best, or tied best, only for the squared alternative.
- Picking κ per cell after seeing these numbers is an oracle choice. It is not a usable procedure.

## Mechanism

Under SEAS-V each series is X_t = s_X(t) ξ_t with a periodic scale s(t) and its own random phase.
Within one sample, the two deterministic scales are correlated, with correlation equal to the
cosine of the phase gap. So |X| and |Y| co-move even though X and Y are independent. That
co-movement averages to zero only across samples.

A Gaussian kernel at a moderate bandwidth reacts to second moments. So under the SEAS-V null its
statistic is driven by this scale co-movement. The data show it directly. Across SEAS-V null
replications, the rank correlation between n·HSIC and corr(X², Y²)² is 0.82 to 0.89 for κ ≤ 2.
Its rank correlation with r² is about 0. Under the G0 null the first correlation is about 0.05.
The scale co-movement widens the null distribution of HSIC, which pushes c* up and power down.

For a large bandwidth, the Gaussian kernel expands as 1 − (x − x')²/(2σ²) + O(σ⁻⁴). HSIC then tends
to a multiple of Cov(X, Y)², and the second-moment term shrinks like 1/σ². The data agree. At
κ = 8, n·HSIC under the SEAS-V null ranks with r² at 0.96, and with corr(X², Y²)² at 0.00. The
oracle power at κ = 8 equals the Pearson oracle in both the seasonal cell (69% against 70%) and the
squared cell (4% against 3%).

So the bandwidth sets a trade-off. A moderate bandwidth sees non-linear dependence, but it also
sees the spurious scale co-movement. A wide bandwidth ignores both. The median heuristic sits on
the sensitive side of this trade-off. The same scale effect explains why wild and subsampling
over-reject under SEAS-V at κ ≤ 2 and stop over-rejecting at κ = 4.

## Checks

- κ = 1 reproduces EXP3 exactly. Thirty EXP3 replications recomputed on the EXP3 seeds match all
  780 stored values with zero difference (unit_checks.txt).
- κ = 1 rates on the new seeds agree with EXP3 within Monte Carlo error. All raw-size |z| ≤ 1.1,
  and every EXP3 power value lies inside the EXP6 bootstrap interval (regression_vs_exp3.txt).
- 0 errors in 3,000 main replications. CPU: 4.34 h main, 0.10 h pilot, about 0.05 h regression.
- The n = 500 block was dropped for budget (A1). Effects at other H values were not tested.

## Limits

- One seasonal design (period 50, amplitude 0.5) and one effect size per alternative.
- κ is applied to the median heuristic. Other bandwidth rules (for example, power-maximising
  choices on held-out data) were not tested.
- 500 replications per cell, not 1,000. Size intervals are about ±2 pp at 5%.

## Files without a writing script (added 2026-11-03, refactor decision D11)

`results.json` and `regression_vs_exp3.txt` were written by hand from the analysis output; no script in this
folder writes them. They are records, not regenerable outputs, and `tools/analysis_gate.py` does not check them.
Their hashes are pinned in `experiments/OUTPUTS.sha256`.
