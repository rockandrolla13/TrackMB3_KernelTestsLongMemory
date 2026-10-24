# Review — JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023

Reviewer `opp-memory-review-20261023`, 2026-10-23. Independent contrary review of `runs/JOB-OPPORTUNITIES-MEMORY-20261023/`.
Discovery output hashes matched its `finished` event before this review started.
Skill: `research-gap-finder` SKILL.md `e81cb1de…25919` and gap template `a5c625dd…f34448`, both as pinned.
This review makes no decision, no novelty claim and no experiment. A search that finds nothing is a search result.

Notation. Betken–Dehling write D for the covariance decay exponent, so D = 2 − 2H. Düker writes d for the fractional parameter, so d = H − 1/2.
The sum line H_X + H_Y = 3/2 is the same as D_X + D_Y = 1 and d_X + d_Y = 1/2.
"Mixed region" means one series has H above 3/4 and the pair is still below the sum line.

## Summary of verdicts

| ID | Discovery | Reviewer verdict |
|---|---|---|
| MB-G1-HSIC-MIXED-NULL | ranked 1, narrowed, moderate | **narrow** further. The linear case is already a theorem. |
| MB-G2-LOWFREQ-KERNEL | ranked 2, open, low, next_priority 1 | **narrow**. The obstruction above the line is now source-backed. Rank 2 and next_priority 1 look too high. |
| MB-UNIFORM-SIZE-LIMITS | ranked 3, narrowed, low | **provisional**. This is a framing check, not a literature gap. One source is now read and changes the wording. |
| MB-G3-LMSV-RANK1 | rejected, answered | **reject** (agree), with a corrected citation and narrower wording |
| MB-HAWKES-NO-H-LINE | rejected, narrowed | **reject** (agree), with one overstatement corrected |
| MB-CROSS-CLASS-MAP | provisional, uncertain | **provisional** (agree). The sweep is done for two anchors and still finds nothing. |

---

## MB-G1-HSIC-MIXED-NULL — narrow

**Closest work and section.**
1. Betken & Dehling, arXiv:2107.03041v3, Remark 2.1(i) and Theorem 2.2. Read in full (PDF).
   The quote is accurate: "Theorem 2.2 does not cover the case D_ξ+D_η>1, but requires either D_ξ<1/2 or D_η<1/2."
   It also says "in the case when D_ξ+D_η<1, our results only cover Gaussian processes".
2. **Düker, "Limit theorems in the context of multivariate long-range dependence", arXiv:1704.08609v4 (2020), Theorem 4.1 and Remark 4.2.** Read in full. Discovery did not cite it.
   Theorem 4.1 gives joint functional limits for all sample cross-covariances X_k,n X_l,n+ℓ of a multivariate linear process. Components may be short- or long-range dependent.
   Each pair is sorted by d_k + d_l. A pair with d_k + d_l < 1/2 gets a √N Brownian limit. A pair with d_k + d_l ∈ (1/2, 1) gets a matrix Rosenblatt limit.
   Remark 4.2: "The sum of two dependence parameters d_k+d_l determines if the corresponding component of the sample autocovariance process is long- or short-range dependent."
   This index set allows one component with d_k > 1/4 (H > 3/4) while the cross pair has d_k + d_l < 1/2. That is the mixed region.
3. Lévy-Leduc & Taqqu, "Hermite ranks and U-statistics", Metrika 77 (2014) 105–136, DOI 10.1007/s00184-013-0474-4. This is the HAL item hal-01197638.
   I identified it through the HAL API. Its abstract treats Σ h(X_i, X_j) for **one** stationary Gaussian sequence. Full text is blocked (HAL bot check, Springer login).
4. Dehling, Rooch & Wendler, arXiv:1404.0551, Section 1–2. Read in full. It treats two-sample U-processes built from segments of **one** LRD process. It cites Lévy-Leduc–Taqqu for "different notions of a two-dimensional Hermite rank".
5. Bai & Taqqu, arXiv:1211.0576, Theorem 3.6 and Conjecture 3.7 (partial text). The setting is functions of **one** common Gaussian series. Its mixed SRD/LRD case is proved only for ranks 1 and 2.

**Direct answer or partial result.** Partly answered. For the identity transform with finitely many lags, the mixed-region null limit of the cross-covariance is a theorem: Düker Theorem 4.1 gives a √N Gaussian limit.
This covers the linear benchmark test that Betken–Dehling compare against. It does not cover their distance-covariance statistic. That statistic integrates nonlinear characteristic-function features and needs a Hilbert-space-valued limit.
Neither Lévy-Leduc–Taqqu nor Dehling–Rooch–Wendler nor Bai–Taqqu treats two independent sequences with different memory. Discovery's fear that the HAL item pre-empts G1 looks unlikely from its abstract. That is abstract-level only.

**Reviewer inference, not from a source.** Under the null, the leading terms of centred HSIC are products H_k(ξ)H_l(η) with k, l ≥ 1. Their covariance decays like |i−j|^−(kD_ξ + lD_η).
In the mixed region D_ξ + D_η > 1, so every such term is summable. The expected limit is therefore the same √n Gaussian-chaos form as Betken–Dehling Theorem 2.2(i).
If this holds, the mixed case is a proof gap (tightness of an RKHS-valued process) with a predictable answer. It is not a case with an unknown limit. That lowers its importance as a stand-alone contribution.

**Contrary evidence.** Düker (above) answers the linear version. Forward citations: Semantic Scholar lists one citing paper for Betken–Dehling, Ernst–Viens–Yan 2510.23563, which treats non-stationary paths. OpenAlex lists zero. No newer paper closes the nonlinear case.

**Remaining question.** Prove the mixed-region √n limit for a fixed bounded characteristic kernel (nonlinear features, RKHS-valued). Separately, prove the above-line limit for non-Gaussian marginals. The second is the less predictable and probably more valuable part.

**Read depth.** Betken–Dehling full text. Düker full text. Dehling–Rooch–Wendler full text. Bai–Taqqu partial. Lévy-Leduc–Taqqu abstract only (blocked).

**Verdict: narrow.** Keep it as the strongest supported gap, but restate it as the nonlinear-kernel case only. Cite Düker for the linear case. The non-Gaussian above-line sub-case may deserve the higher priority.

---

## MB-G2-LOWFREQ-KERNEL — narrow

**Closest work and section.**
1. Müller & Watson, "Low-Frequency Analysis of Economic Time Series", Handbook of Econometrics vol. 7 draft (May 2024), §8.2, Theorem 2. Read in full.
   It states: "The inference results of this chapter crucially depend on the large sample Gaussianity of the suitably scaled cosine transform".
   Theorem 2 proves that CLT for **linear** processes v_t = Σ a_s ε_{t−s} with martingale-difference innovations.
2. Düker Theorem 4.1 (above). The sample cross-covariance of two components with d_k + d_l ∈ (1/2, 1) has a non-Gaussian matrix Rosenblatt limit.
3. Betken & Dehling Theorem 2.4 and Section 5. Their subsampling block length ⌊n^γ⌋ satisfies the theorem's condition for every D. But they choose the statistic's √n normalisation "Under a corresponding assumption" that D ≈ 0.6. So their calibration is not H-free, and it does not pre-empt G2.

**Direct answer or partial result.** Discovery's worry that the route "may fail above the line" was its own reasoning. It now has source support:
- The low-frequency method needs Gaussian cosine averages (Müller–Watson §8.2).
- Above the sum line the unweighted cross-covariance sums are non-Gaussian (Düker Theorem 4.1, linear case). Betken–Dehling Theorem 2.2(ii) shows the same for distance covariance.

No source shows that cosine-*weighted* averages of rank-1 × rank-1 cross features are non-Gaussian above the line. That last step is still inference. It is the standard expectation from non-central limit theory, though.

**Contrary evidence.** Of 86 Semantic Scholar records citing Müller–Watson 2018, none extends low-frequency inference to kernels, nonlinear features or independence tests. Fresh searches found no fixed-q, Hermite-rank-2 or Rosenblatt combination. No pre-emption exists, but the obstruction is now clearer.

**Remaining question.** Two narrower questions survive:
- (a) Below the sum line, including the mixed region, does Gaussian low-frequency calibration of HSIC cross features stay valid as the pair approaches the line? This is where HAC-type estimates are fragile (Müller 2007).
- (b) Above the line, is there a non-Gaussian low-frequency analogue whose limits can be bounded over a d range without estimating d?

(b) is hard and may not exist. The "one cheap derivation decides it" framing is too optimistic. The derivation will very likely show a non-Gaussian, d-dependent limit, which comes close to discovery's own stopping criterion.

**Read depth.** Müller–Watson handbook §8.2 full. Düker full. Betken–Dehling full. NBER w23186 not re-read.

**Verdict: narrow.** Restate it as "low-frequency calibration below the sum line, with a stated obstruction above it". I disagree with next_priority 1. A self-normalised block-sampling route (Bai–Taqqu–Zhang 1512.00820, already known to the programme) handles non-Gaussian limits for scalar statistics. It may be the better H-free comparator.

---

## MB-UNIFORM-SIZE-LIMITS — provisional

**Closest work and section.**
1. Preinerstorfer & Pötscher, arXiv:1304.1383v2 (Econometric Theory 2016), Section 3, Assumption 1, Theorem 3.3, Remark 3.5(i), and the discussion of Pötscher 2002 in §3.2.2. **Read in full this time.** Discovery read only the abstract.
   - The size-one / zero-power results "rest on Assumption 1, i.e., C ⊇ C_AR(1)". This means the covariance model contains AR(1) correlations with ρ → ±1.
   - Remark 3.5(i) weakens this to a model containing correlation matrices that converge to e₊e₊′ (all-ones) and e₋e₋′.
   - The setting is finite-sample F-type HAR tests in linear regression.
   - On Pötscher 2002 they write: it "only discusses ill-posedness of the problem of estimating f(0) … and does not make any statements regarding consequences of this ill-posedness for autocorrelation robust tests".
2. Pötscher 2002, Econometrica 70(3). Still abstract only. The working-paper URL returns HTML, and CORE redirects.
3. Müller 2007 Theorem 1. Not re-opened. Discovery's full-text quote that it is not "a statement about the impossibility of valid inference" fits the Preinerstorfer–Pötscher reading.

**Direct answer or partial result.** The impossibility mechanism is near-unit-root concentration, which is short memory. It is not fractional d.
*Reviewer inference:* ARFIMA(0,d,0) correlations tend to 1 at every lag as d → 1/2. So a memory-agnostic class that lets d approach 1/2 would also meet the Remark 3.5 condition. A class with d bounded below 1/2 would not.
The answer to the candidate's question is therefore conditional. It depends on the parameter range the claim allows, not on "memory-agnostic" as such.

**Contrary evidence.** No uniform-in-d inference result found in a fresh extended search.

**Remaining question.** Record each result's model class against the simulation's d grid. Pötscher 2002 §§2–3 is still unread.

**Read depth.** Preinerstorfer–Pötscher full text. Pötscher 2002 abstract only.

**Verdict: provisional.** This is a framing note and a reading task, not a literature gap. In the skill's terms it is a coverage shortfall ("文献读得还不够"), not a 真正的 gap. I recommend no rank.
I disagree with two pieces of discovery wording:
- The evidence matrix says Pötscher 2002 "suggests that any calibration needing an estimate of H … cannot be uniformly informative". Preinerstorfer–Pötscher say Pötscher 2002 concerns estimating f(0) itself and says nothing about tests.
- The gap analysis implies these results bound "H-free" claims. In fact they bind only if the claimed class reaches the nonstationary edge.

---

## MB-G3-LMSV-RANK1 — reject (agree), narrower wording

**Closest work and section.** Kulik & Soulier, arXiv:1109.5298 (Adv. Appl. Probab. 2012), partial text (first 100k of 150k characters).
- Section 2, eq. (4): Y_i = σ(X_i)Z_i. Z is heavy-tailed with E|Z₀|^α = ∞. In the LMSV case η and Z are independent.
- Section 4: "τ_p = τ(σ^p) be the Hermite rank of the function σ^p" — the rank governing partial sums of |Y|^p.
- Section 5.1: τ*_p(s) is "the Hermite rank of the bivariate function K*_{p,s}". This governs **sample covariances** of |Y|^p. For σ = exp, τ*_p(s) = 1. For σ(x) = x², τ*_p(s) = 2.
- Section 5, LMSV: "the Hermite rank of the function K depends only on the function σ".
- Example 5.2 (σ = x², leverage): rank 1 "except if E[Z₀η₀]=0, which holds in the LMSV case".

**Direct answer.** The conclusion survives. Discovery cites the wrong quantity, though. The §5.1 sentence is about the covariance kernel K*, not "every power |Y|^p".
For partial sums the relevant rank is τ_p = rank(σ^p). For σ = exp it is 1 by direct calculation, since E[e^{pX}X] = p·e^{p²/2} ≠ 0. That is reviewer arithmetic, not a quoted line.
Because rank depends "only on the function σ" in the LMSV case, the rank statement carries over to finite-variance Z whenever E|Z|^p < ∞. Limit rates and regimes still differ from the infinite-variance paper. I agree with discovery's derivation.

**Contrary evidence and narrowing.** Kulik–Soulier also show that rank 2 *can* hold. With a symmetric volatility function such as σ(x) = x² and no leverage, τ* = 2. Leverage then drops it to 1 (Example 5.2).
So the accurate statement is: rank-2 protection is not generic. It fails under exponential volatility, can hold under symmetric σ without leverage, and leverage removes it.
"Fails in the standard model" is fair only if exp-LMSV is named as that model. The real-data rank is empirical. A rank test exists (Beran–Möhrle–Ghosh, JCGS 2016, abstract only). The one applied study found tests log returns, not |r|.

**Transfer to the track's setting.** G3 asks about the break point of a single-series kernel test, so the H = 1/2 conclusion fits.
For a two-series cross statistic, rank-1 features give the sum line H_X + H_Y = 3/2, not H = 1/2. Any later wording should say which statistic it means.

**Read depth.** Partial full text (100k of 150k). Surgailis–Viano 2002 located but not opened.

**Verdict: reject as a stand-alone gap (agree).** Fix the citation to §4 τ_p for partial sums, or §5.1 τ* for sample covariances. Keep the σ-dependence and leverage caveats.

---

## MB-HAWKES-NO-H-LINE — reject (agree), one overstatement

**Closest work.** Jaisson & Rosenbaum, arXiv:1504.03100 (AAP 2016), abstract only. A nearly unstable Hawkes process with kernel tail x^−(1+α), α ∈ (1/2, 1), has a scaling limit "with associated Hurst parameter H=α−1/2".
Filimonov–Sornette 1308.6756 (abstract): kernel misspecification biases the branching ratio upward. Diouane–Schoenberg–Mohler (NSF PAR, search snippet): spurious cross-excitation near criticality.

**Direct answer.** Discovery's setting says Hawkes memory is "not a Gaussian Hurst index". That is too strong. Heavy-tailed near-critical Hawkes processes do carry a Hurst index through the kernel tail.
That index lies in (0, 1/2), which is rough, not long memory. It does not arise from Gaussian subordination, so no Hermite-rank 3/4 line follows. The rejection stands.
Discovery's candidate has empty `supporting_evidence` and abstract-only contrary evidence. The rejection rested mostly on the programme's own exploratory simulation. Jaisson–Rosenbaum gives it a theoretical anchor.

**Contrary search.** No size study of Hawkes Granger tests indexed by H was found.

**Verdict: reject (agree).** Reword it as "the Hawkes memory index is set by kernel tail and criticality, and it is rough near criticality. No Hermite mechanism applies."

---

## MB-CROSS-CLASS-MAP — provisional (agree)

**Sweep done this review.**
- Betken–Dehling: Semantic Scholar shows 1 citing paper, not a map. OpenAlex shows 0.
- Müller–Watson 2018: 86 citing papers, none on kernel, DML or Hawkes.
- Kulik–Soulier: 15 citing papers, none a map.
- Zhai–Cheng–Wu 2609.33868: 0 citations returned. The metadata call hit a rate limit, so I cannot tell whether the index holds the paper.
- A fresh search found only single-class size studies, such as change-point tests at H 0.6–0.9 (1506.07296) and LRD subsampling (1509.05720).

**Partial result relevant to the map.** The cross-term sum rule (H_X + H_Y = 3/2) is already a theorem for linear cross-covariances (Düker Remark 4.2). The map should cite it as an established result, not present it as a finding.

**Remaining.** Google Scholar, SSRN and a working Zhai–Cheng–Wu citation lookup are still missing.

**Verdict: provisional (agree).** No paper found is a search result, not proof of absence.

---

## Disagreements with discovery (explicit)

1. **G1.** Discovery says nothing closes any part of the mixed case. Düker Theorem 4.1 (1704.08609) closes the linear, identity-transform version, and discovery did not cite it. G1 should be restated for nonlinear bounded kernels only.
2. **G1.** Discovery treats "Hermite ranks and U-statistics" as possible pre-emption with unknown authors. It is Lévy-Leduc & Taqqu, Metrika 2014. Its abstract is single-sequence, so pre-emption of the two-sequence case looks unlikely (abstract-level only).
3. **G1.** In my reading, the mixed-region limit is predictable (√n, Gaussian chaos), so importance is lower than discovery states. This is reviewer inference.
4. **G2.** The "may fail above the line" point is no longer only worker inference. Müller–Watson §8.2 and Düker Theorem 4.1 support it. Given that, next_priority 1 ("cheapest decisive check") is too high. The outcome is largely predictable.
5. **UNIFORM.** Pötscher 2002 does not support the evidence-matrix wording about "any calibration needing an estimate of H". Preinerstorfer–Pötscher state that Pötscher 2002 concerns estimating f(0) only. Their own mechanism is AR(1) near-unit-root concentration. Recommend provisional, no rank.
6. **G3.** The quoted §5.1 line is about τ*, the sample-covariance kernel rank, not "every power |Y|^p". The conclusion holds via §4 τ_p. Kulik–Soulier also show rank 2 is possible for symmetric σ without leverage, which discovery omitted.
7. **Hawkes.** "Not a Hurst index" is overstated. Near-critical heavy-tailed Hawkes has H = α − 1/2 < 1/2. The rejection still holds.

## Coverage limits of this review

- Semantic Scholar rate-limited some calls (HTTP 429). Google Scholar and SSRN were not used.
- Lévy-Leduc–Taqqu 2014 and Pötscher 2002 are inaccessible, so I have abstracts only.
- Kulik–Soulier and Bai–Taqqu 2013 were read for the first 100k characters only.
- The reviewer inferences in G1, G2 and UNIFORM are labelled as such. They need a specialist check.
