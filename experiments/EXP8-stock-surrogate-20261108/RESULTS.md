# EXP8 results (2026-11-08)

Run: `run_exp8.py --mode main --out main --workers 16`. It produced 20,000 pairs in 24 cells in 1,464 s wall time, and all rows have status ok. The simulated check (`--mode sim`) ran first, before any real data. Numbers below are rejection percentages at a nominal 5%, with Wilson 95% intervals. Pooled numbers pool the four pairings at n = 2000. Full tables are in `tables.md` and `summary.csv`. The verdicts are in `results.json`, produced by `analyze_exp8.py`.

## Verdicts on the pre-registered predictions

| | Prediction | Verdict | Evidence |
|---|---|---|---|
| P0 | Arm A windows are the saved STOCK windows | **Pass** | 0 mismatches in 4,000 pairs: days, starts, shift p-values to 12 digits, Pearson to within 1e-12 |
| P1 | On the surrogate null, the shift rules stay in 3.4–6.8% | **Pass** | raw 4.5–5.1%, inc1 4.1–4.9%, over all 4 pairings |
| P2 | On the surrogate null, raw wild and raw subsampling exceed 5%, rising with H | **Fail** (2 of 3 parts pass) | wild 6.3% [5.6, 7.1], sub 7.6% [6.8, 8.5]. The H slope is not significant (one-sided p = 0.43) |
| P3 | Keeping the drift (S2) raises raw wild and raw subsampling | **Fail** | wild +0.8 points [−0.3, 1.9], sub +1.0 [−0.2, 2.2]. The shift rules held under S2 (4.6–6.2%) |
| P4 | inc1 wild rejects less than raw wild, and inc1 shrinks the trend gap | **Fail** (trend part passes) | wild 6.6% to 5.7%, McNemar one-sided p = 0.11. The trend gap falls from 3.1 to 0.8 points |
| P5 | Harris rejects at most 5% on the surrogate null | **Pass** | 2.2–3.2% in all 16 cells |

## What the run shows

1. **Circular shifts hold their size on real stock durations, with a known null.** On the surrogate null, raw HSIC with shifts gives 4.9% [4.3, 5.6] pooled at n = 2000. Wild multipliers give 6.3% and subsampling 7.6%. These are true false-alarm rates with respect to time alignment. The paper had only placebo fractions for this arm.
2. **The stock windows have mild memory.** The mean local Whittle H is 0.61 (sd 0.05). That narrow range is why the H trend in P2 could not be detected.
3. **Drift alone adds little.** Keeping each surrogate's slow scale drift raised raw wild and raw subsampling by about one point. Both intervals include 0. The within-window drift that the paper links to the stock fractions is therefore at most a small part of them.
4. **The increment kernel fixes subsampling on real data.** In Arm A at n = 2000, inc1 subsampling rejects 4.3% [3.5, 5.3], against 8.8% for raw. On the surrogate null it rejects 4.2% [3.6, 4.9]. This was not a pre-registered prediction, so it is exploratory.
5. **Harris is conservative here.** It rejects 2.2–3.2% on the surrogate null. Its size holds. EXP8 does not measure its power.

## Exploratory comparisons (not pre-registered)

- **Placebo against known null.** Arm A (real against real) and S1 (real against surrogate) give similar raw rates at n = 2000: wild 6.6% and 6.3%, subsampling 8.8% and 7.6%. The placebo excess in the paper is therefore close to what the rules produce under a true null on these data. Most of it looks like calibration error rather than real cross-stock dependence. The Y windows differ between the arms, so this is a descriptive comparison.
- **At n = 500,** most wild and subsampling rates are higher, for example raw wild at 7.9–9.6%. The shift rules stay at 4.5–5.7%.
- **Hermite split, hsplit2.** With shifts it gives 4.8–5.3%. Its block rule rises to 5.9–8.2%.

## Caveats

- **What the surrogate keeps.** It keeps each window's values, spectrum and, in S2, its drift. It cannot keep higher-order time structure such as volatility clustering. A rejection on S1 is a false alarm with respect to time alignment only.
- **Reused days.** Days are reused across pairs, as in STOCK-20261023. Pair outcomes are therefore dependent. Arm B samples its days the same way, so it reuses days too.
