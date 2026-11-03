# Adversarial review of EXP7 (memory kernels)

Reviewer: adversarial review agent, 2026-11-03. Read-only review. Scripts and outputs are in the session scratchpad:
`/tmp/claude-1000/-media-ak-10E1026C4FA6006E-GitRepos-TrackMB3-KernelTestsLongMemory/13e492c5-ed83-4c16-a7a0-b488e851d1ef/scratchpad/review/`
(`recompute.py`, `recompute2.py`, `mmd_size.py`, `trend.py`, `binom_boot.py`, `theory_lambda.py`, `theory_layers.py`, with `.out` logs).

## Verdict in one paragraph

Every number in RESULTS.md that I recomputed matches the raw rows exactly. I found no code error that changes a rate. One headline claim is wrong in scope. Point 1 ("under every memory kernel the false-alarm rate does not climb with H") is false for split MMD. RESULTS.md leaves out the pre-registered split-MMD size results, and there several memory kernels do climb with H. The strict "removes memory" verdicts are also close to a coin flip: a rule at exactly 5% fails "valid in all 24 cells" with probability 0.50. So the yes/no split between kernels under shift carries no information. The "falling inc1 null" comes from the median bandwidth, not from the kernel. THEORY's fixed-bandwidth limit does not describe the EXP7 statistic at high H.

## Critical

### C1. Split-MMD size is pre-registered, omitted from RESULTS.md, and contradicts headline point 1

- Where: `RESULTS.md:7` ("Under every memory kernel it does not"), `RESULTS.md:57` ("Every memory kernel is memory-flat under wild"). DESIGN.md:163 applies "removes memory" to split MMD on its G0 cells. DESIGN.md:144 lists 28 split-MMD null cells. RESULTS.md reports split-MMD power only (`RESULTS.md:99-108`). It never reports split-MMD size or verdicts. The numbers exist in `tables_main.md:756-842`.
- Evidence: my recomputation (`mmd_size.py`) agrees with `tables_main.md`. Uncentred wild, n = 2000, H = 0.5 → 0.9:

| kernel | G0 | HT-S | memory-flat (G0 / HT-S) |
|---|---|---|---|
| raw | 5.4 → 38.4 | 4.7 → 35.3 | no / no |
| datadiff | 4.3 → 4.3 | **3.0 → 17.4** | yes / **no** |
| sym | 5.9 → 8.9 | **4.2 → 11.6** | yes / **no** |
| herm1 | 5.4 → 7.9 | **3.6 → 9.4** | **no** / **no** |
| herm2 | 4.8 → 6.5 | 5.3 → 7.5 | **no** / yes |
| hsplit1 | 4.5 → 7.0 | **4.1 → 8.0** | yes / **no** |
| hsplit2 | 4.6 → 6.7 | 4.4 → 6.2 | yes / yes |

  On HT-S, datadiff, sym, herm1 and hsplit1 all exceed their own H = 0.5 baseline beyond the Wilson intervals. They also rise from n = 500 to n = 2000 at H = 0.9: datadiff 9.3 → 17.4, sym 7.9 → 11.6, herm1 5.2 → 9.4, hsplit1 4.4 → 8.0. That is the memory signature the study tests for. On G0, herm1 and herm2 fail memory-flat through a rise at H = 0.7 (3.7 → 7.7, 3.6 → 6.8).
- What it changes: headline point 1 holds for Gaussian independence testing only. For split MMD, only hsplit2 holds on both generators. datadiff holds on G0 but breaks on heavy tails, at 17.4% at H = 0.9. RESULTS.md must report the split-MMD size table and verdicts, and restrict point 1.

## Major

### M1. The strict "removes memory" verdict is a coin flip for a rule at exactly 5%

- Where: `RESULTS.md:9`, `RESULTS.md:38-55`, `RESULTS.md:125-126`, `RESULTS.md:162`. The rule is in `analyze_exp7.py:115-121` and `analyze_exp7.py:314`.
- Evidence (`binom_boot.py`, exact binomial): for a rule whose true rate is exactly 5%, one cell fails "valid" with probability 0.0284 at 1,000 replications and 0.0292 at 300. The chance that all 24 cells pass is **0.497**. RESULTS.md:162 says "about 1%" per cell, which is wrong by a factor of about 3. Pooled over the 24 cells, every memory kernel under shift sits at 4.95% to 5.34% (`trend.py`): sym 4.95, hp10 4.98, datadiff 5.10, dcov_inc1 5.17, hsplit2 5.17, herm2 5.20, inc1 5.25, herm1 5.32, hsplit1 5.34. Simulating 24 binomial cells at each pooled rate, P(24/24 valid) is 0.17 to 0.56.
- What it changes: point 3 says "under shift, herm1, sym and hsplit1 meet the strict rule". That does not separate them from inc1, herm2, datadiff, hp10 or dcov_inc1. All of these behave like 5% rules. The "Refuted as written" for P1 and P2 should rest on the wild half only (pooled 7.7% and 7.2%). The shift half (21/24 and 22/24) is consistent with exact size. The table at `RESULTS.md:38-55` should also give pooled rates with intervals.

### M2. "Memory-flat" is a weak test, and it hides one real trend

- Where: `RESULTS.md:7`, `RESULTS.md:57`. The rule is in `analyze_exp7.py:303-312`. A cell counts only if its Wilson interval does not overlap the H = 0.5 interval.
- Evidence: at n = 2000, a baseline of 6.4% needs a memory cell near 10% to register. At n = 4000 (300 replications), a baseline of 7.7% needs about 16%. I fitted a logistic model of rejection on H, with n fixed effects, on the 18 equal-H Gaussian cells (`trend.py`):
  - **Good news for the study:** no H trend for inc1, herm1, herm2, sym, datadiff or hsplit1, under wild or shift. Every p-value is above 0.2. So the reading "the wild excess is the rule's own error" is supported by a stronger test than the one used.
  - **hp10 wild climbs with H:** slope +0.62 per unit of H, p = 0.0016. That survives Bonferroni over 21 rules. Pooled, it goes from 7.6% (H ≤ 0.6) to 8.8% (H ≥ 0.8). dcov_inc1 wild has p = 0.021, which does not survive.
  - hsplit2 wild falls with H (p = 0.0006). It is conservative at high H.
- What it changes: "every memory kernel" in point 1 is false for hp10 under wild, by about 1 pp. Replace the overlap rule with the trend test in the write-up, or report both.

### M3. Headline point 1 promotes the secondary reading after the primary reading failed

- Where: `RESULTS.md:7` against `DESIGN.md:161-171`. DESIGN says the informative reading is under wild and subsampling. It also says memory-flat "is weaker than 'removes memory' and is reported as such".
- Evidence: under wild and subsampling, no memory kernel passes "removes memory" (24 Gaussian cells). Point 1 opens with "Changing the kernel removes the memory effect" and does not say that the pre-registered criterion failed. Point 3 says it later.
- What it changes: point 1 should begin "No memory kernel passes the pre-registered rule under wild or subsampling." The memory-flat result follows as the weaker reading. Together with C1, point 1 then reads "for Gaussian independence tests, under the weaker memory-flat reading".

### M4. The falling inc1 null is a median-bandwidth effect, and THEORY's limit does not describe EXP7's statistic

- Where: `RESULTS.md:8`, `RESULTS.md:75`. THEORY.md:14 (bottom line 1), `THEORY.md:136` (Theorem A, `|E T − Λ| ≤ C/n`), `THEORY.md:142` (Λ values), `THEORY.md:179` (median-bandwidth step, marked Sketch).
- Evidence (`theory_lambda.py`): THEORY's exact identity is right at a fixed bandwidth. I reproduce Λ = 1.166 … 0.0714 to four digits at σ = 0.954, the median heuristic for N(0, 1). But EXP7 uses the sample median bandwidth. At high H, that bandwidth is far below its limit and moves slowly. At H = 0.95, mean σ̂ is 0.650 at n = 500 and 0.732 at n = 8000, against 0.954 in the limit. Plugging the mean σ̂ into THEORY's fixed-σ formula reproduces the EXP7 mean of T_inc1 within 5% in every (H, n) cell:

| H, n | EXP7 mean T_inc1 | formula at mean σ̂ | THEORY Λ (fixed σ) |
|---|---|---|---|
| 0.9, 500 | 0.337 | 0.335 | 0.208 |
| 0.9, 8000 | 0.265 | 0.269 | 0.208 |
| 0.95, 500 | 0.250 | 0.239 | 0.071 |
| 0.95, 8000 | 0.171 | 0.168 | 0.071 |

- What it changes:
  1. The fall in inc1's 95th percentile (×0.65 at H = 0.95) is explained by the bandwidth. It is not a property of the increment kernel. It will continue past n = 8000, towards about 0.07. RESULTS.md should give this mechanism.
  2. "A falling null makes a rule conservative" (`RESULTS.md:75`) does not follow. Wild and shift recompute their reference on the same sample with the same bandwidth. A falling unconditional quantile says nothing about their size.
  3. THEORY's Theorems A, B and D, and C7 for raw, are proved at a fixed bandwidth. EXP7's σ̂ converges at about the rate of the sample variance, n^(2H−2), which is n^(−0.1) at H = 0.95. The step marked Sketch at `THEORY.md:179` is therefore not a technicality. At H ≥ 0.9 no practical n is near the fixed-bandwidth regime. THEORY should state this rate.
  4. The raw slope shortfall (0.538 against 0.8 at H = 0.95, `RESULTS.md:73-75`) may come partly from the same drift. RESULTS.md says only that growth "has not reached its limiting rate". I did not test this.

## Minor

1. **Harris literal-form check does not cover the primary rule.** `exp7_core.py:494-513` checks the segment-bandwidth form only at N = 19. The primary rule uses N = 99 and 199 at n = 2000 and 4000. "2 of 5,280 decisions differ" (`RESULTS.md:138`) is therefore evidence about harris19 and about n = 500. The largest ratio difference is 7.7% (`results_main.json`, `harris_max_rel`). "As its theorem promises" (`RESULTS.md:11`) overstates what is shown for the main columns.
2. **hsplit decision and ratio use different conventions.** The decision uses `p_shift_herm ≤ 0.025` with the +1 correction (`analyze_exp7.py:176`). The ratio R uses `t / np.quantile(·, 0.975)` (`analyze_exp7.py:67-70`, `exp7_core.py:290`). R > 1 is not exactly the rule. DESIGN.md:290 says "up to the p-value convention". It affects size-adjusted hsplit power only marginally.
3. **Matched nulls are approximate for sq, lag and flat.** All four of lin, sq, lag and flat use the G0 null (`analyze_exp7.py:196-207`). Under sq, Y carries a He_2(ξ) term, so its marginal law and memory differ from fGn. Under lag and flat, Y is not exactly fGn. DESIGN.md:133 says so for flat only. With ρ ≤ 0.13 the effect should be small. I did not measure it.
4. **Seeds and cut-offs: checked, no finding.** `spawn_key = (cell_idx, rep)` with unique cell names (`run_exp7.py:104-110, 132`). No seed is reused across cells. So c* comes from replications independent of the alternative, as DESIGN.md:214 states. Within a cell, all kernels share data and multipliers. The bootstrap keeps that pairing (`analyze_exp7.py:333-334`), which is correct. One caveat: lin, sq, lag and flat share one G0 null cell per H. So their c* errors are common, and comparisons across these alternatives are correlated.
5. **Percentile intervals are skewed, and several verdicts sit on the boundary.** Example: Harris factor at H = 0.7 has point 78.8%, bootstrap mean 76.7% and interval [64.2, 83.9]. Borderline verdicts include raw shift − Harris, lin H = 0.5: +5.1 [+0.6, +9.4]. Others are inc1 − raw, lin H = 0.9: +6.8 [+0.8, +12.6]; inc1 sq H = 0.7: +10.7 [+0.9, +18.4]; Harris − raw shift, scale H = 0.8: +6.0 [+0.2, +10.5]. There are several hundred gap verdicts and no multiplicity control. Read single "gains" or "loses" labels near the boundary as tentative.
6. **Seasonal omissions.** hsplit1 shift rejects 15.0% and 15.2% at SEAS-V-slow. It inherits herm1's remainder failure (19.2%, 19.7%). `RESULTS.md:117` attributes hsplit's seasonal failures to the low layer, but here the remainder is the cause. herm1 shift rejects 9.3% and 11.8% at SEAS-M-shared. Harris rejects 15.5% at SEAS-M-shared, H = 0.8. `RESULTS.md:118` gives only H = 0.5.
7. **Rank expectations are "refuted" on the wrong evidence.** `RESULTS.md:129` cites shift validity, but DESIGN.md:165 says shift passes by itself. The growth of T for herm1 and sym was never measured, because growth cells hold raw and inc1 only (`exp7_core.py:873-883`). The better evidence is the absent H trend under wild (herm1 p = 0.57, sym p = 0.95, `trend.py`). The refutation stands on that.
8. **THEORY power proxy against EXP7.** `THEORY.md:226` says raw is "slightly ahead up to H = 0.9" at n = 2000. On lin at H = 0.9, EXP7 at matched size has inc1 shift +6.8 pp over raw shift [+0.8, +12.6]. datadiff is +32.1 and dcov_inc1 +38.1. THEORY should note that EXP7 contradicts its proxy at H = 0.9.
9. **THEORY X = ξ² retention figures are underspecified.** `THEORY.md:77` says c_2 keeps 10.8%, 0.76% and 0.011% at r = 1, 2, 3. The plain L² residual share of He_2(ξ) after removing He_1..He_r(T) is 30.8%, 2.7% and 0.046% (`theory_layers.py`, 4M draws). THEORY's quantity is probably kernel-weighted, but THEORY.md does not define it. I could not reproduce it, and I do not claim it is wrong.
10. **Small wording slips.** `RESULTS.md:10` gives "4% to 10%" for factor and flat blindness, but datadiff factor at H = 0.9 is 10.3%. `RESULTS.md:114` calls 12.3% (herm2, δ = 0.5, H = 0.8) "blind".

## What I checked and found correct

- **Headline numbers.** I recomputed from `main/*.csv` with my own Wilson, size-class and c* code: section 1 (all 20 rows × 6 H), the verdict table (valid counts, removes-memory, memory-flat and worst cells for all 28 rules), growth ratios and slopes, the section 3 power table (48 entries), the split-MMD power table, process arms, FACTOR-null rates, P2 McNemar (p = 0.0067 and 0.125), P3 (9 of 11), P6 layer means (0.518 and 0.892) and the Harris range (1.6% to 3.8%). All match RESULTS.md to the printed digit.
- **Kernels.** These match DESIGN.md:
  - Q_r K Q_r uses Q from a QR of [1, He_1..He_r] of the same series' normal scores, with K on raw x.
  - D K D' uses lag 1, and lag 50 = `gens.PERIOD`. Centring follows the transform.
  - The bandwidth comes from raw x for inc, herm, hp10 and frac_oracle, and from the series itself for sym and datadiff.
  - The frac_oracle order is d = H − 1/2 per series.
- **Rules.**
  - Shift references run over c ∈ [⌊0.1m⌋, m − ⌊0.1m⌋], with p = (1 + #≥)/(1 + N).
  - Wild uses centred block-Markov multipliers, and shorter kernels take the first m rows, re-centred.
  - Subsampling uses no +1 and a strict "> q95", as in the paper.
  - Harris uses segments X[N:m−N] and Y[s+N:s+m−N], s = −N..N, N = 19/99/199. Its count excludes s = 0, and the M-th-largest ratio is equivalent to m < M.
  - The Bonferroni splits are 0.025/0.025, with 4 × 0.00625 (HSIC) and 2 × 0.0125 (MMD).
  - The self-normalised statistics match the Lobato zero-mean form and Betken's G_n(k)².
- **Audit.** 6,600 exact replications, largest relative T error 2.1e-6, 0 decision flips, 0 error rows. The stop rule did not trigger.
- **THEORY.** The Λ values reproduce exactly at a fixed σ. The SV table (correlation of He_p(ξ) with He_q(T)) reproduces to three decimals (0.488, −0.113, 0.052, −0.030 in row 1). The level-shift lower bound m^(2−2H)/(4 − 2^(2H)) checks algebraically, and gives 37.57 at H = 0.8, n = 16000.
