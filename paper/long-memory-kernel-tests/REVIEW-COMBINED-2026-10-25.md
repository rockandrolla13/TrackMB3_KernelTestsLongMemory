# Combined review, 25 Oct 2026: AISTATS + JFec

Manuscript: `paper/long-memory-kernel-tests/paper.tex` (version with EXP2 added), 15 pages
(main text to p. 12, references p. 14), 16 references. Read-only review. No code was changed.

## Bottom line

- **AISTATS:** borderline / weak reject. EXP2 makes the empirical story stronger, but the paper
  still has no theorem, no size-adjusted power, and does not meet the format rules (12 main pages
  against 8, not anonymised, no checklist).
- **JFec:** weighted score 5.95 / 10, so **major revision (probably reject)**. The economic
  interpretability gate scores 3, which is the minimum to pass. The ML special-issue deadline
  (1 Sep 2026) has passed. The regular issue is the realistic JFec target.
- **Top three fixes, in priority order:**
  1. Add size-adjusted power for every rule on the same cells. Without it, "circular shift is
     best" cannot be told apart from "circular shift has low power".
  2. Add a deseasonalised control to EXP2: estimate and remove the seasonal mean and scale, then
     re-test. This shows whether the seasonal-scale failure goes away under standard practice.
  3. Correct the EXP2 verdicts for multiple testing. About 190 interaction tests use an
     uncorrected z > 1.96 cut-off, so roughly 10 "effects" are expected by chance alone.

---

# Part A. AISTATS review

**1. Summary and Contributions:**

The paper measures how often three kernel tests (split-sample MMD, Gaussian-kernel HSIC, distance
covariance) reject true nulls when the data have long memory. It contributes (i) false-alarm
rates under known nulls across seven critical-value rules, (ii) a Hermite-rank covariance
calculation that predicts different memory rates for split MMD (rank one fails at any H > 1/2)
and for independence tests (the memory-sum line H_X + H_Y = 3/2), (iii) a new factorial
experiment (EXP2, §3.4 and §4.3) that crosses memory with heavy tails, nonlinear dynamics and
seasonality, and (iv) a descriptive placebo on Nasdaq trade durations. The main new finding is
that a seasonal scale shared by two independent series drives HSIC to 73.3% rejection with no
memory, while Pearson correlation stays at 4.7% and circular shifts at 5.0% (Table 4). I reviewed
the full manuscript and appendix. There is no supplement or reproducibility checklist.

**2. Strengths:**

- **The nulls are known.** Every simulation rejection is a real false alarm, and the paper keeps
  this apart from the stock placebo, which it labels descriptive throughout (§3, Table 2).
- **The assumption audit is unusually careful** (§2.3). It works out why each published guarantee
  excludes fGn. Examples are the τ-dependence lower bound τ(r) ≥ (2/π)^{1/2}|ρ(r)| against the
  Chwialkowski et al. rate conditions, and the Bradley bound ruling out strong mixing.
- **EXP2 isolates one feature at a time**, with a Gaussian row run through the same code path and
  a design record written before the pilot. Pearson acts as a linear control, so the paper can
  separate kernel-specific failures (seasonal scale) from failures shared by every test (seasonal
  mean).
- **The seasonal results come with a mechanism.** The phase-gap argument explains why the
  seasonal-mean failure does not shrink with n, and why circular shifts survive it.
- **Scope is stated honestly.** The TAR short-memory confound is reported with local-Whittle
  estimates (0.59 against 0.50), and the Hermite analogy for SEAS-V is labelled informal.

**3. Weaknesses:**

*Major*

1. **No formal result.** The Hermite calculation is a heuristic about scalar averages, and the
   paper says so (Table 1, §2.1). AISTATS reviewers will read the paper as purely empirical. That
   is acceptable only if the empirical evidence is decisive, and point 2 shows it is not yet.
2. **No size-adjusted power.** Circular shifts come out best in every arm, including both seasonal
   forms. A rule that rarely rejects also gets this result. The one power cell (§4.1) compares
   rules at different null sizes, and the paper says it is not size-adjusted. Without power, the
   paper cannot recommend a rule.
3. **Is the seasonal-scale result a test failure or a modelling choice?** Given the two phases,
   the paths do share a time-varying scale. Practitioners deseasonalise before testing. A reviewer
   will ask whether the 73% survives a standard estimated seasonal adjustment. If it does not, the
   finding narrows to "kernel tests need deseasonalised inputs", which is still useful but smaller.
   A deseasonalised control is missing.
4. **Multiple testing in EXP2.** The verdicts use an uncorrected interaction cut-off of |z| > 1.96
   across about 190 tests. The seasonal effects survive any correction, since their z values are
   above 8. Some secondary claims do not. The filtered t3 excess for HSIC subsampling at n = 2000
   has z = 3.47, which is below a Bonferroni cut-off of about 3.6. The "about 8 points" claim in
   §4.3 needs this caveat or a correction.
5. **The memory excess still falls with n.** Table 3 shows Δ < 0 at sum 1.6 for six of seven rules.
   The earlier review raised this, and it is unresolved. Readers will ask whether the excess is a
   finite-sample effect that disappears as n grows.

*Minor*

- Local-Whittle estimation is used in §4.3 without a citation (Robinson 1995).
- There is one bandwidth rule (the median heuristic). Kernel tests are known to be sensitive to it.
- EXP2 uses one parameter setting per generator and n ≤ 2000.
- With four contributions, the stock arm (descriptive only) takes space away from the core story.
- Only 16 references. The separate literature review lists 90 works, including Shao's
  self-normalisation papers and Boniece–Horváth–Trapani (arXiv 2608.26422), which are missing.

**4. Detailed Evaluation:**

- **Formatting & Structure:** fails as submitted. The main text runs about 12 pages against the
  usual 8. The draft is not in the AISTATS template, not anonymised (it gives a public GitHub link),
  and has no reproducibility checklist. The structure itself is clear: theory → design → results
  → discussion → limitations.
- **Theoretical Grounding:** the calculations are correct as stated, and the assumption audit is
  strong. There are no theorems, and the paper makes no claim to have any. The SEAS-V Hermite
  analogy is now labelled informal, which is right.
- **Methodological Novelty:** moderate. Mapping where the rules fail is new for MMD and
  Gaussian-kernel HSIC. The seasonal-scale finding is new as far as I can tell, but its
  seasonal-mean companion is Yule (1926), and the paper says so. The novelty of the seasonal-scale
  claim is assessed against the manuscript's own positioning only. A reviewer may know
  spurious-dependence results for squared series under shared volatility.
- **Empirical Validation:** strong on null size (1,000–2,000 reps, Wilson intervals, a baseline
  row). Weak on power, on sample-size range (n ≤ 2000) and on parameter range. There is no
  multiple-testing control in EXP2.
- **Clarity & Transparency:** good. Claims are scoped, limitations are explicit, and the audit
  trail is complete. Some paragraphs in §2.3 are dense.

**5. Constructive Feedback & Questions for Authors:**

Revisions, in priority order:
1. Add size-adjusted power for all rules on matched alternatives: a mean shift for MMD, and a
   weak linear and a weak squared dependence for HSIC.
2. Add a deseasonalised EXP2 arm that removes the estimated seasonal mean and scale.
3. Apply Holm or Bonferroni to the EXP2 verdicts and report which ones survive.
4. Extend the key cells (sum 1.4 and 1.6, SEAS-V) to n = 4000 or 8000.
5. Cut to 8 pages. Move the stock arm and §2.3 detail to the appendix.

Questions:
- Does the SEAS-V failure persist after an estimated seasonal-scale adjustment?
- At matched size, does circular-shift HSIC keep power against a small dependence?
- Does the HSIC subsampling excess at sum 1.4 stay flat beyond n = 2000?
- Would a rank-based or Hermite-split kernel (Discussion) also resist the seasonal-scale failure?

---

# Part B. JFec review

```
JFec Review — Special Issue on Machine Learning and AI in Economics
Title: When do kernel tests raise false alarms on long-memory data?
Manuscript classification (Phase 0): none of the four tags fits; closest is "structural/inferential
  methodology" (a size study of dependence tests). Flag: this is not an ML-in-economics paper.
Reviewer recommendation: major revision (probably reject) for the special issue
Weighted score: 5.95 / 10
```

## 1. Summary

The manuscript studies the null size of kernel two-sample and independence tests (split MMD, HSIC,
distance covariance) when observations have long memory. Using known-null simulations, it shows
that the size depends heavily on the critical-value rule. A Hermite-rank calculation organises the
failures by feature rank and memory sum. A new factorial experiment separates memory from heavy
tails, stochastic volatility, threshold dynamics and seasonality. It finds a seasonal-scale failure
that hits kernel tests but not Pearson correlation. A placebo on AAPL, MSFT and NVDA trade
durations reports descriptive rejection fractions.

## 2. Contribution and significance

The contribution is a methodological warning, supported by simulation, for financial
econometricians who use kernel dependence tests on persistent data: realised volatility,
durations, order flow. People who apply HSIC or dCov to market data, and authors of long-memory
HAC and self-normalised inference, would cite it. It is not ML used as a tool inside an economic
model, which is what the special issue asks for.

## 3. Strengths

- Dependence is the subject of the paper, not an afterthought. The mixing audit in §2.3 is better
  than most JFec submissions.
- Pearson with Newey–West and EWC long-run variances is the right econometric benchmark, and it is
  present in every independence table.
- Causal and predictive honesty is excellent. Nothing is claimed beyond size under the stated
  generators.
- The intraday-duration application is the right object for JFec, and the paper treats it
  cautiously.
- Reproducibility is strong: seeds, manifests, hashes, and a design record written before the
  pilot.

## 4. Major concerns

1. **Economic interpretability (special-issue gate, score 3).** No economic object is estimated.
   The duration placebo gives rejection fractions with no link to a pricing, liquidity or risk
   quantity. To pass the gate, connect the seasonal-scale result to the diurnal pattern in
   durations and volatility (Andersen–Bollerslev 1997; Engle–Russell 1998). Then show what a
   false "dependence" finding would do to a concrete decision, such as a cross-asset
   liquidity-spillover test or pairs selection.
2. **The intraday seasonality in the application is only partly handled.** §3.5 subtracts a
   15-minute clock-bin mean from log durations. EXP2 now shows that a shared seasonal *scale* is
   what misleads HSIC. The stock arm needs a clock-bin scale adjustment, or a standard diurnal
   factor as in ACD modelling, followed by a re-run. Otherwise a JFec referee will read Table 5
   through this confound.
3. **No power, so no recommendation.** JFec readers want to know which test to use. As in AISTATS
   point 2, size-adjusted power is required before circular shifts can be recommended.
4. **Asymptotics are missing where JFec expects them.** A size study is acceptable, but the house
   style expects at least one formal result. The cheapest candidate is a proposition stating that,
   under the random-phase SEAS-M null, the sample cross-correlation converges to a non-degenerate
   limit set by the phase gap. That would make the nonsense-correlation mechanism formal.
5. **Multiple testing.** About 190 uncorrected interaction tests are behind the EXP2 verdicts
   (see AISTATS point 4). JFec referees will ask for Holm or Romano–Wolf adjustment.

## 5. Minor concerns

1. Cite Robinson (1995) for local Whittle, and Taqqu (1975) / Dobrushin–Major (1979) for Hermite
   rank and non-central limits.
2. Cite the long-memory HAC literature, such as Abadir–Distaso–Giraitis (2009) on memory- and
   autocorrelation-consistent (MAC) estimation, beside NW and EWC.
3. Cite Shao (2010) for self-normalisation. The SN-MMD rule uses the idea.
4. State the regime coverage of the January–June 2024 stock sample, and note that it has no
   stress episode.
5. The SV parameters (φ = 0.9, τ = 1) should be compared with values estimated on equity returns.

## 6. Phase-1 audit table

| Claim | Estimand | Conditioning | Regime | Evaluation | Verdict |
|---|---|---|---|---|---|
| Split MMD rejects 27.1% at H = 0.8 (abstract) | null rejection probability | marginal over generator draws | finite-sample (n = 2000) | Monte Carlo, known null | valid |
| HSIC subsampling 9.2% at (0.7, 0.7) | same | same | finite-sample | Monte Carlo | valid |
| H_X + H_Y = 1.5 is a useful axis | covariance summability | — | asymptotic heuristic | indirect (grid plots) | plausible but unverified |
| Seasonal scale drives HSIC to 73.3%, Pearson 4.7% | null rejection probability | marginal over phases | finite-sample | Monte Carlo, known null | valid |
| Filtered t3 adds about 8 points (§4.3) | interaction in rejection | marginal | finite-sample | Monte Carlo, uncorrected z | plausible but unverified |
| TAR excess is not pure nonlinearity | attribution | — | finite-sample | local-Whittle diagnostic | valid as a caveat |
| Circular shifts closest to 5% | size | marginal | finite-sample | Monte Carlo | valid for size, silent on power |
| Stock rejections 3.4–10.2% | placebo fraction | unknown null | finite-sample | descriptive | valid as descriptive |

No claim is causal. There is no in-sample versus out-of-sample issue, because there is no
forecasting.

## 7. Score

| # | Criterion | Weight | Score | Weighted |
|---|---|---|---|---|
| 1 | Economic interpretability | 15% | 3 | 0.45 |
| 2 | Rigor under dependence | 15% | 7 | 1.05 |
| 3 | Causal/predictive honesty | 12% | 9 | 1.08 |
| 4 | Novelty at ML–econ interface | 10% | 5 | 0.50 |
| 5 | Benchmark adequacy | 10% | 5 | 0.50 |
| 6 | OOS discipline (here: design pre-registration) | 8% | 7 | 0.56 |
| 7 | Empirical significance | 8% | 4 | 0.32 |
| 8 | Multiple-testing awareness | 7% | 4 | 0.28 |
| 9 | Reproducibility | 8% | 9 | 0.72 |
| 10 | Presentation | 7% | 7 | 0.49 |
| | **Total** | | | **5.95** |

Criterion 1 sits exactly at the gate (an automatic reject below 3). Criterion 2 is safe.

## 8. Required next steps

```
MUST-FIX
  □ Size-adjusted power for all rules on matched alternatives (Major 3)
  □ Clock-bin scale adjustment of stock durations, then re-run Table 5 (Major 2)
  □ Link to an economic decision: spillover test or pairs selection on durations (Major 1)
  □ Holm/Romano–Wolf correction of EXP2 verdicts (Major 5)
SHOULD-FIX
  □ Proposition on the random-phase cross-correlation limit (Major 4)
  □ Deseasonalised EXP2 control arm
  □ Add citations: Robinson 1995, Shao 2010, Abadir–Distaso–Giraitis 2009,
    Andersen–Bollerslev 1997, Engle–Russell 1998
NICE-TO-HAVE
  □ n = 4000–8000 for sum 1.4/1.6 and SEAS-V cells
  □ A stress-period stock sample (for example March 2020)
```

## 9. Recommendation

As an inference-under-dependence paper, this is careful and honest work. For the ML/AI special
issue it misses the economic-interpretability bar, and the deadline has passed. For the JFec
regular issue, the must-fix list would make it competitive. The most valuable additions are the
stock re-run with a scale adjustment, which ties the new EXP2 result to a real microstructure
pattern, and size-adjusted power.
