# STOCK-DIURNAL results (2026-10-25)

Full tables: `tables.md`. Summary: `summary.csv`. Diagnostics: `diagnostics.json`.

1. **Regression check passed.** D0 reproduces every saved row of `experiments/STOCK-20261023/results` at n = 2000 (0 mismatches).
2. **H0 not rejected.** D1 and D2 move no cell by more than 0.6 points at any n. The diurnal adjustment does not explain the placebo fractions.
3. **D2 removes the daily scale pattern.** Smooth sd max/min falls from 1.33-1.45 to 1.007-1.013.
4. **D2 does not remove the within-window variance trend.** Share of n = 2000 windows with a trend beyond the permutation 95th percentile: 45.6% (D0), 45.4% (D1), 43.7% (D2). The no-trend rate is 5%.
5. **Trend goes with rejection (descriptive).** D2, n = 2000, pooled over the four pairings: HSIC sub rejects 7.5% of pairs with no trending window, 8.8% with one, and 11.0% with two. Pearson NW: 4.6%, 6.5%, 6.4%.

D3 (day-mean removal) is an exact no-op: max p-value difference from D2 = 4.4e-16 on 20 pairs.
Process independence is unknown and days recur, so all fractions are descriptive placebo rates.
CPU: D0 3,794 s, D1 3,939 s, D2 3,844 s.
