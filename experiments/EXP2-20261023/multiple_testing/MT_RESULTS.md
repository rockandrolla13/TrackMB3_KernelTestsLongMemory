# EXP2 multiple-testing re-analysis

Re-analysis only. No new simulation. Script: `mt_analysis.py` (run with `python -I multiple_testing/mt_analysis.py` from the EXP2 folder). It imports the frozen `analyze.py` and checks that it rebuilds every verdict in `verdicts_main.json` exactly.

## Answer

- 28 of the 50 non-null verdicts survive Holm unchanged. At n = 2000 the count is 15 of 26.
- Four more verdicts keep part of their claim. For example, "breaks alone, dampens memory" keeps "breaks alone".
- No verdict that read "no effect" becomes non-null.
- The heavy-tail "8 points" claim survives only as a level gap at H = 0.8. As an interaction it fails Holm (adjusted p = 0.12).
- The heavy-tail "3 points" distance-covariance claim fails. Its raw p is 0.11 before any correction.

## The test family

Each verdict combines three tests. The tests run per rule, generator and n.

1. One "alone" test: rate(F, H=0.5) against rate(G0, H=0.5). Pooled two-proportion z. One-sided as in DESIGN.md H1. A two-sided version is also reported.
2. Two interaction tests, at h = 0.7 and h = 0.8. These use the same z and SE as `analyze.interaction`. They are two-sided because the reading rule flags both directions.

The grid has 8 verdict rules, 6 feature generators and 2 sample sizes. That gives 96 verdicts and 288 tests. The MMD arm holds 108 tests and the independence arm holds 180.

Several §4.3 claims compare levels at H = 0.8, not interactions. The verdict family does not contain those tests. So the script adds a second family of 96 level tests, F(0.8) against G0(0.8), two-sided. It also reports an "augmented" family of all 384 tests.

## Results across the 288 tests

- 88 tests carry the uncorrected flag (|effect| > 2 pp and |z| > 1.96).
- 92 have raw p < 0.05.
- 48 survive Holm. Bonferroni gives the same 48.
- 62 survive Benjamini–Yekutieli (BY). 76 survive Benjamini–Hochberg (BH).
- Holm needs a raw p below about 0.0002. BY needs a raw p below about 0.0016.

The tests that fail are almost all interaction flags with |z| between 2 and 3.5. The seasonal "alone" effects and the stochastic-volatility effects have p-values far below any threshold. They are not affected.

## Which family to use

I recommend the overall family of 288 tests. The reason is that §4.3 draws one conclusion from both arms ("the MMD arm agrees with this picture"). A family per arm treats the arms as separate papers. Per-arm Holm rescues only 2 more verdicts (30 of 50). Both rescued verdicts are self-normalized MMD flags that the paper does not cite. So the choice does not change any claim in the text.

Holm is the primary method. It controls the chance of any false finding under any dependence between tests. The tests here are dependent in two ways. First, the G0 cells enter every interaction test for their rule and n. Second, all rules in one cell use the same replications. BY controls the false discovery rate under the same arbitrary dependence. BH assumes positive dependence, which these tests do not guarantee. It is reported for reference only.

## Is the independent-binomial SE right for the interaction?

Yes. `run_exp2.py` gives every cell its own index (`cell_list` enumerates family × generator × H × n). Each replication draws from `SeedSequence(entropy=master, spawn_key=(cell_idx, rep))`. The four cells in one interaction contrast have four different indices. So their random streams are independent, and the variance of the contrast is the sum of the four binomial variances. The two arms also use different master seeds. One caveat applies. At rates near 0% (seasonal split MMD) the normal approximation is weak. No §4.3 claim depends on those interaction tests.

Alone tests: a Fisher exact test agrees with the z test on every alone comparison except one borderline cell (mmd_sub_sqrt, n = 2000, NL-TAR: z p = 0.045, Fisher p = 0.055). That cell fails Holm under both tests.

## §4.3 claims (n = 2000)

1. **Filtered t3 adds about 8 points to HSIC subsampling at H = 0.8.** As a level gap (25.9% against 17.5%) it survives. Holm p = 0.0003 in the level family and 0.0016 in the augmented family. As an interaction (+7.4 pp) it fails Holm (p = 0.12) and survives BY (p = 0.017).
2. **Filtered t3 adds about 3 points to distance-covariance subsampling at H = 0.8.** Fails. The level gap is +3.1 pp with raw p = 0.11. The interaction is +6.7 pp but fails Holm (p = 0.41) and BY (p = 0.051). The interaction is larger only because HT-L starts below G0 at H = 0.5.
3. **"No excess at H = 0.5" for filtered t3.** True for excess. But distance-covariance subsampling sits 3.6 pp *below* G0 (3.2% against 6.8%). The two-sided p is 0.0002. Holm in the two-sided family gives 0.050, right at the threshold. BY gives 0.006. The text should not imply HT-L matches G0 here.
4. **Threshold AR adds 11–12 points to both subsampling rules at H = 0.8.** Survives. Level gaps are +12.3 and +10.7 pp. Interactions are +10.9 and +11.1 pp. Every Holm p is below 0.001.
5. **Stochastic volatility keeps HSIC and dCov subsampling at or below 6.2% at H = 0.8.** Survives. Holm p < 1e-6 for both level gaps and both interactions.
6. **Stochastic volatility breaks MMD square-root subsampling alone (13.3%).** Survives. Holm p = 3e-8.
7. **Seasonal mean and seasonal scale break the kernel tests alone.** Survives for all four kernel rules. Holm p < 1e-80. Circular shifts show no effect (two-sided raw p 0.63 and 0.76).
8. **Seasonal series make split MMD conservative at H = 0.5.** Survives as a two-sided test. Rates sit 3.3–5.1 pp below G0, and every Holm p is below 5e-5. The one-sided "breaks alone" test cannot detect this direction by design.
9. **Uncentred wild MMD at 26–31% at H = 0.8 across non-seasonal generators.** Correct as a description (26.3–30.6%). No interaction with G0 survives Holm. The uncorrected HT-L "only with memory" flag is gone.

Verdicts that change under Holm and are cited in or near §4.3:

- hsic_sub, n = 2000, HT-L: "only with memory" becomes "no effect".
- dcov_sub, n = 2000, HT-L: "only with memory" becomes "no effect" (also under BY).
- hsic_wild and dcov_wild, n = 2000, NL-TAR: "only with memory" becomes "no effect".
- hsic_wild and dcov_sub, n = 2000, SEAS-V: "breaks alone; dampens memory" keeps only "breaks alone".

The full list of changed verdicts is in `mt_summary.json` under `survival`.

## Files

- `mt_results.csv`: one row per test (288 verdict tests plus 96 level tests). Raw p and Holm, Bonferroni, BH and BY p for the overall, per-arm, n = 2000-only, two-sided-alone and augmented families.
- `mt_claims.csv`: the §4.3 claims mapped to their tests and adjusted p-values.
- `mt_summary.json`: family counts and verdict survival under each method.
