# Deep literature review — kernel tests under long memory (TrackMB3, NEW-MB-KERNEL)

Execution `claude-mb3-lit-20261016` · Agent 2 · 2026-10-16.
Starts from `hypotheses/hypotheses-and-claims.md` (SHA-256 `b937424a…afb04`, verified).
Modes: deep-research `systematic-review` for the core question (`systematic-review.md`) and `lit-review`
for the surrounding questions (this file). Gap analysis: `gap-analysis.md`. Comparison of the closest
papers: `summary-3w.md`. Records: `evidence-matrix.md`, `bibliography.bib`, `search-log.md`.

**Evidence labels.** F = full-text section read by me. S = summarising fetch, or a verbatim sentence
read inside another full text. A = abstract or snippet. I = identity checked only. D = my own inference.
Programme simulation evidence is called "track S" to avoid a clash.
**No novelty claims.** "Not found" means not found in this bounded scan (SSRN unreachable, Scholar
partial — see `search-log.md`).

---

## Part A. Answers to the 18 hand-off questions

### For H1, H1b, H1c (MMD / homogeneity)

**Q1. Any paper with the null limit of MMD (or RKHS two-sample / kernel change-point) under
Gaussian-subordinated LRD, a rank-one χ² law, or the n^{2H−1} divergence?**
Not found. Searched: arXiv field queries (MMD, energy distance, two-sample U-statistic, kernel change
point, each with LRD terms), WebSearch (open and journal domains), Scholar, and forward citations of
Dehling–Rooch–Wendler (8 citers, all scalar change-point or U-statistic work, none on MMD).
The nearest is Dehling, Rooch & Wendler (Statistics 2017; arXiv:1404.0551, **F**): two-sample
U-statistic processes Σ_{i≤[nt]}Σ_{j>[nt]} h(X_i, X_j) under subordinated Gaussian LRD, with the
Hermite rank m of h deciding the limit; non-central when Dm < 1. It is written for change-point kernels
(Gauss test h = y − x, Wilcoxon), not for MMD. Energy-distance or MMD tests under LRD: not found.
*Effect on H1:* the rank-one χ² law and κ_H remain worker derivations (D). Pre-emption risk for that
law: unchanged (medium). For "MMD breaks under LRD" as a statement: still likely to be read as a
corollary of DRW plus Dehling–Taqqu (1989/1991) U-statistic theory (I, not read).

**Q2. Subsampling with an estimated rate for U/V-statistics under LRD; verify the R items.**
- The hand-off's "Hall, Horowitz & Jing (1998)" is **wrong**. The LRD sampling-window paper is **Hall,
  Jing & Lahiri (1998), Statistica Sinica 8(4):1189–1204** (I; from the Betken–Wendler reference list,
  F). It covers subordinated Gaussian sample means.
- **Jach, McElroy & Politis (2012), JTSA 33(1):96–111** (A, Crossref): a self-normalised mean for
  heavy-tailed LRD; its law is approximated "by subsampling without the knowledge of tail or/and memory
  parameters". So it avoids rate estimation by self-normalising, rather than estimating the rate.
  Corrigendum JTSA 37(5) 2016.
- **Bertail, Politis & Romano (1999), JASA 94:569–579** (A): estimate the unknown rate τ_n by
  subsampling, then use it. This is the general rate-correction recipe.
- **Betken & Wendler (2018), Statistica Sinica** (F): subsampling valid for general statistics of
  subordinated Gaussian LRD with l = O(n^{(1+D)/2−ε}); applied to a **self-normalised** statistic so
  that no nuisance parameter is estimated. They note: "for other statistics, one needs an unknown scaling
  … one has to estimate the scaling parameters before applying subsampling".
- **Zhang, Ho, Wendler & Wu (2013), SPA 123:2323–2339** (I; described in Zhai–Cheng–Wu, F, and ZKW, F):
  block sampling for LRD means, generalised to nonlinear transforms of linear processes.
- **Zhai, Cheng & Wu (arXiv:2609.33868, Sep 2026)** (F): block sampling **with an estimated
  normalisation** for realised conformal coverage, proved consistent across the Gaussian, critical
  (H = 3/4) and Rosenblatt regimes (Thm 7). This is the closest existing instance of H4a's
  "rate-corrected subsampling across a Hermite phase transition", for a scalar indicator statistic.
- No application to an MMD / HSIC / U-statistic test under LRD was found.
*Effect on H1b/H4a:* "rate-corrected subsampling" is known for scalars and is being used across the 3/4
transition by the Wu group. Kernel version: not found. Pre-emption risk for the *recipe*: high.
The (n/l)^{2H−1} collapse of size (H1b) as an explanatory index: not found.

**Q3. Does an even transform (|X|, X²) change the size of a kernel two-sample / change-point test
under LRD?**
Not found for kernel tests. Two related results:
- Zhai–Cheng–Wu (F): for |W| scores of Gaussian linear LRD residuals, the coverage indicator has
  covariance h^{2−4β}, i.e. the rank-2 line H = 3/4. This is H1c's mechanism, for a different statistic.
- **Contrary caution:** Bai & Taqqu, "How the instability of ranks under long memory affects
  large-sample inference", Statistical Science 33 (2018) (A), and "Sensitivity of the Hermite rank",
  SPA 129 (2019) (I): ranks other than 1 are statistically unstable; assuming a higher rank "usually
  results in underestimating the order of the fluctuations". For H1c this means a nearly-even G (any
  asymmetry, a non-zero mean shift, skewed marginals) brings back a small rank-1 term that dominates at
  large n. The clean |ξ| design will likely confirm H1c, but real |r_t| need not.
- Tewes (EJS 2017; A) treats change-point tests where the Hermite rank itself changes.

### For H2, H2b, H3 (HSIC / dCov cross rule)

**Q4. Betken–Dehling v3 in full.** (F; PDF SHA-256 `5bc4b26c…52b50df`, identical to the prior job's copy.)
(a) *Mixed case.* Remark 2.1(i): "Theorem 2.2 does not cover the case D_ξ + D_η > 1, but requires either
D_ξ < 1/2 or D_η < 1/2." No supplement section or later paper found that closes it. The only citer found
by Semantic Scholar is Ernst–Viens–Yan (2025, fBm paths; A).
(b) *Form of the Thm 2.2(ii) limit.* For G₁ = G₂ = id and D_ξ + D_η < 1:
n^{(D_ξ+D_η)/2 −1}L^{−1/2}… (φ̂_{XY} − φ̂_Xφ̂_Y)(s,t) → [∫_{[−π,π)²} K(x,y) Z_X(dx)Z_Y(dy)] × s t e^{−(s²+t²)/2},
with K the Fourier kernel of the partial-sum, Z_X, Z_Y independent random spectral measures. Corollary
2.1(ii): n^{D_ξ+D_η}L^{−1}… dcov_n → [∫ K dZ_X dZ_Y]² × ∫∫ s²t² e^{−(s²+t²)} w(s,t) ds dt.
The same double Wiener–Itô integral is the limit of the normalised empirical cross-covariance
(their Thm 4.1). **So for Gaussian series above the line, the dCov null is a squared scalar
cross-correlation variable times a constant — H3 points 1–2, for dCov with G = id, are in this paper.**
For non-linear G (e.g. signs) and for HSIC, not covered.
(c) *HSIC with general kernels.* Discussed only as related work (Sejdinovic et al. 2013 equivalence for
a specific kernel; HSIC under mixing). No HSIC result under LRD.
Betken–Dehling–Kroll (arXiv:2112.14091; A): α-mixing only. `betken2024bootstrap` in the prior bib is the
same paper. Neither covers LRD.

**Q5. HSIC / kernel independence tests under long memory, 2024–2026.**
Not found. Diz-Castro, Febrero-Bande & González-Manteiga (arXiv:2604.28104v2, Oct 2026; F scan): HSIC
under **L1-near-epoch dependence on a β-mixing process**; no LRD. AutoHSIC (Li, Xu & Zhou,
arXiv:2605.22025; A): serial independence, short memory. Ghoshal–Shao (EJS 2026): strong mixing.

**Q6. Bai & Taqqu on multivariate Hermite rank.**
Bai & Taqqu, JTSA 34(6) 2013 (A): limits of vectors of normalised sums of functions of LRD Gaussian
series can be Gaussian, Hermite, or a mix; mixed ranks 1 and 2 treated; the general mixed case is a
conjecture. Rank instability (Stat. Sci. 2018; A) as in Q3. Arcones (AoP 1994; I) gives limit theorems
for functionals of Gaussian **vector** sequences — the natural tool when X and Y are built from
*dependent* latents (H7's conditional null). Not read. For the unconditional H2 (ξ ⊥ η) the cross term's
expansion has only mixed terms H_a(ξ)H_b(η), a, b ≥ 1, so its covariance is Σ c_{ab}² a! b! ρ_ξ^a ρ_η^b
(D). Rank instability then matters only through the marginal centring terms (H2b's mechanism).

**Q7. Robinson & Hidalgo (1997) exact statement.**
Full text not reached (Project Euclid fetch returned an empty page). Abstract (A, via OpenAlex): "When
sufficiently strong spectral singularities in the error and a regressor coincide at the same frequency,
least squares need no longer be n^{1/2}-consistent … our class of estimates is n^{1/2}-consistent and
asymptotically normal." The **explicit sum condition** is in **Tsay & Chung (2000), J. Econometrics
96:155–182** (S): "as long as [the variables'] orders of integration sum up to a value greater than 1/2,
the t-ratios become divergent and spurious effects occur" (verbatim, read inside Ventosa-Santaulària et
al. 2022, F). VSVL (F) add: at d_x + d_y = 1/2 the t-ratio diverges at (log T)^{1/2}; below 1/2 it is
O_P(1) but not necessarily N(0,1) in finite samples.
*Effect on H2:* the linear sum rule (Pearson / OLS) is classical. H2's contribution can only be the
kernel version.

### For H4 (calibration)

**Q8. Shao self-normalisation.**
- Shao (2010, JRSSB 72(3):343–366; A): SN confidence intervals, weak dependence.
- Shao (2011, JTSA 32(6):598–606; A + S): SN CUSUM change-in-mean test "in the possible presence of LRD".
  Full text paywalled. Bai, Taqqu & Zhang (2016, SPA; F) state: "As observed by Shao (2011),
  self-normalization itself is not able to fully avoid the problem of estimating the nuisance
  parameters, as the asymptotic distribution at least depends on the unknown Hurst index H." **The R item
  is confirmed (S).** Whether Shao tabulated critical values in H: not verified (full text not reached).
- Betken (2016, JTSA 37(6):785–809; F): SN Wilcoxon change-point test; the limit is a functional of fBm
  B_H; **critical values tabulated in H** (Table 1, 5% level): 7.28 (H = 0.6), 8.19 (0.7), 9.50 (0.8),
  10.33 (0.9). With the true H, size is 4.8–5.5% for n ≥ 50 (Table 2). What is pivotal: the statistic is
  free of the scale d_n and of the Hermite coefficient, but **not of H**.
- Shao (2015, JASA review; A): lists "inference for long memory time series" among topics; full text not
  reached; no kernel or U-statistic SN under LRD visible in the abstract.
*Effect on H4b:* point 2's premise ("limit depends on H only; tabulate in H") is right. Its claim that
"the dependence on H is mild" is **contradicted** by Betken's table: the 5% critical value rises 42%
from H = 0.6 to 0.9. A plug-in Ĥ error of ±0.1 moves the critical value by about 10%.
The established fix is **SN + subsampling / block sampling** (Betken–Wendler 2018; Bai–Taqqu–Zhang 2016;
Jach et al. 2012), which needs no Ĥ.

**Q9. Ghoshal–Shao.** (F; arXiv v2 SHA-256 `b4607abe…6621`.)
- **Same work** as the EJS paper: Crossref lists "Resampling-free inference for time series via RKHS
  embedding", Ghoshal & Shao, Electronic Journal of Statistics 2026, doi 10.1214/26-EJS2488.
- Assumptions: Assumption 1 is a Hilbert-space FCLT for the RKHS-embedded series with a **Brownian
  motion** limit in H(K). Proposition 1 gives sufficient conditions: stationary **strong (α-) mixing**
  with Σ α(m)^{δ/(2+δ)} < ∞ and E|K(Y,Y)|^{1+δ/2} < ∞. (The v1 abstract said "absolute regularity";
  v2 says strong mixing. The hand-off's "β-mixing" is the v1 wording.)
- Long memory or fractional limits: **not mentioned anywhere** in the text (grep for long-range, long
  memory, fractional, Hurst: zero hits). Future work named: high-dimensional or non-stationary series.
*Effect on H4b:* pre-emption risk is not that Ghoshal–Shao already did LRD. It is that both halves are
public (GS for the kernel projection + SN; Bai–Taqqu–Zhang / Betken–Wendler for H-free SN under LRD), and
Shao's group has written on both. Risk: **high**.

**Q10. Block bootstrap under LRD.**
- Lahiri (1993, SPL 18:405–413; I + S): as summarised by Betken–Wendler (F), the moving block bootstrap
  "might fail under long range dependence".
- Kim & Nordman (2011, Sankhya A 73:79–109; I): identity verified; content not read.
- Tewes (2018, JTSA 39:28–53; A): the block-bootstrapped empirical process converges to a semi-degenerate
  limit whose random part is always Gaussian, "thus the bootstrap might fail when the original empirical
  process accomplishes a noncentral limit theorem". Supports H4c above the line.
- Zhang, Lahiri & Nordman (2022, AoS 50:3619–3646; A): block-resampling variance estimators under
  Gaussian-subordinated LRD depend on the non-linearity "beyond Hermite ranks"; block order n^{1/2} often
  reasonable; data-driven block rule.
- ZZW / ZKW (F): overlapping block bootstrap valid for covariance statistics for H < 3/4 only.
- Zhai–Cheng–Wu (F, Table 2): MBB with b = n^{2/3} is oversized at β = 0.75 (13–17%) and β = 0.6
  (32–44%) for the coverage statistic.
- No block bootstrap for **degenerate** U/V-statistics under LRD found.

**Q11. Wild / dependent wild bootstrap; null oversize at independence.**
- Shao (2010, JASA 105:218–235; A): dependent wild bootstrap for smooth functions of means, weak
  dependence.
- Leucht & Neumann (2013, JMVA 117:257–280; I): DWB for degenerate U/V-statistics; their Section 5
  finding is reported in CSG 2014 (below).
- **Chwialkowski, Sejdinovic & Gretton (NeurIPS 2014; F):** with the centred multipliers W̃_t = W_t − mean
  ("artificial degeneration", their B₂ / V_{b2} and the two-sample MMD_{k,b} in eq. 4), Type I error
  stays above 5%: 0.070, 0.105, 0.100 (MCMC) and 0.145, 0.120, 0.114 (audio, rising n), "although it can
  be observed to drop with increasing sample size". The uncentred V_{b1} gives 0.012–0.060. The authors
  say this is "consistent with findings in [Leucht–Neumann], Section 5". In their two-sample form each
  sample gets its **own independent** multiplier chain.
*Effect on H5:* candidate cause 1 (within-half centring) is directly supported by the method's authors
(F); cause 3 (one chain across both halves) departs from the published two-sample form. H5 is a hygiene
item, as stated; it is not new.

**Q12. Zhai–Zhong–Wu Section 6 and Zhai–Kolar–Wu.**
- ZZW §6 (F): "The theory relies on temporal dependence β > 3/4. For ultra-long-memory time series with
  1/2 < β ≤ 3/4, the asymptotic distribution will be Rosenblatt-type, and is beyond the scope of this
  paper." **No announcement of forthcoming work.** Future work named: entrywise precision inference.
  No kernel or non-linear extension.
- ZKW **verified**: IEEE Trans. Inf. Theory 72(6):4246–4296 (2026), doi 10.1109/TIT.2026.3685246;
  preprint arXiv:2604.16219 (F). Condition 1 needs β > 3/4; "Gaussian approximation holds only if
  β > 3/4. When 1/2 < β ≤ 3/4, the … distribution … transitions from Gaussian to the Rosenblatt
  distribution"; simulations at β = 0.55 "fail again". No kernel content.
- **But** the same group (Zhai, Cheng & Wu, two arXiv papers, 27 Sep 2026) already crosses the 3/4 line
  for conformal coverage with Rosenblatt limits and estimated-normalisation block sampling (F).
*Effect on H2/H7:* the scoop risk is real but specific: the group has the tools to cross 3/4, and
names it as a limitation, yet shows no kernel or HSIC interest in these three papers.

### For H7, H7b (kernel CI / Granger)

**Q13. CI or Granger tests with validity conditions stated in terms of residual vs conditioning-variable
memory.** Not found. GCM (Shah–Peters 2020; A), GKCM (PMLR 2026; A), GTCM (Shi–Shah 2026; A), dGCM
(2025; A): i.i.d. or weak dependence. Linear Granger tests robust to long memory: Bauer & Maynard (JoE
2012; R/A), Chen (2006) spectrum-based Granger allowing long memory in the transfer function (A).
ZZW (F) is Gaussian-linear CI with β > 3/4 on the whole process, not split into residual vs conditioning
memory.

**Q14. Residual autocorrelation of finite-order AR approximations to ARFIMA.**
Poskitt (2007, AISM 59:697–725; I) "Autoregressive approximation in nonstandard situations: the
fractionally integrated and non-invertible cases" and Poskitt (2008, JTSA, sieve bootstrap for FI; I)
exist. Content not read (no open abstract). Bhansali, Giraitis & Kokoszka (2006, JMVA; I) fit FAR(p, d)
models; Shao (2009, F) uses that FAR sieve for prewhitening. **Reading gap**, not a literature gap.

### For H8 (hand to F-KERNEL)

**Q15. Zhai–Cheng–Wu.** Two papers, both arXiv 27 Sep 2026:
- "Conformal Coverage of Time Series: Validity and Inference" (arXiv:2609.33868; F). Split conformal,
  adjacent calibration/test. FDM-based non-asymptotic coverage bounds; CLT for realised coverage under
  summable indicator FDM. Long-memory part: Gaussian linear residuals W_t = Σ a_j ε_{t−j}, a_j ~ c j^{−β},
  score |W|. The score and its coverage indicator have covariance ~ h^{2−4β}. Realised coverage is
  √n-Gaussian for 3/4 < β < 1 (H < 3/4), Gaussian with √(log n / n) at β = 3/4, and n^{1−2β} times a
  Rosenblatt functional for 1/2 < β < 3/4 (Thm 6). Block sampling with estimated normalisation is
  consistent in all three (Thm 7).
- "Valid and Efficient Split Conformal Regression for Time Series" (arXiv:2609.33866; A): CQR and median
  conformal; interval-length accuracy, sharper for Gaussian linear long memory.
*Effect on H8:* points 2–3 for **symmetric (|residual|) scores** are already proved there (rank 2 → the
3/4 line). Signed / one-sided scores (rank 1, any H > 1/2) are not treated in what I read. Kernel-weighted
scores (point 4): not found. Pre-emption for H8: **high**; hand to F-KERNEL with these two papers.

### For H9 (finance)

**Q16. Order-sign exponents after 2010; durations, spreads; spurious-memory tests.**
- Sato & Kanazawa (PRL 131:197401, 2023; A + S): Tokyo Stock Exchange 2012–2020, account-level data;
  confirm γ ≈ α − 1 (LMF). A search summary reports α ≈ 3/2, which gives γ ≈ 0.5, H = 1 − γ/2 ≈ 0.75
  (S; number not read in the paper). Order signs again sit **on** the line.
- Tóth, Palit, Lillo & Farmer (JEDC 51, 2015; A): LSE; persistence "out to … tens of thousands of
  orders"; the exponent could not be extracted from the PDF (pdftotext xref error).
- Crypto / futures order-sign exponents 2020–2026: not found by WebSearch (W13).
- Durations and spreads with CIs; spurious long-memory tests on these series: not reached in budget.
  Qu (2011) spurious-memory test exists (A, via an R package page).
- SSRN: not searchable (HTTP 403). Hwu & Tsay (SSRN 7546859, 2026; A, title only): aggregation-induced
  long memory and spurious inference in DiD — adjacent.

**Q17. SSRN and Scholar for Q1, 5, 8, 9, 12.** Scholar was reached and queried 10 times; 4 pages came
back empty of results, others noisy and recency-biased. SSRN direct search: HTTP 403; SSRN-restricted
WebSearch returned only arXiv records. **The SSRN gap is not closed.** Scholar is partly closed.

### For H10 (prewhitening)

**Q18. Fractional prewhitening before independence / cross-correlation tests.**
- **Shao (2009), Econometric Theory 25:195–210 (F):** a generalised portmanteau test of independence
  between two stationary series. Each series is prewhitened in the frequency domain with an estimated
  spectral density (parametric, or a FAR(p, d) sieve as in Bhansali et al. 2006); the kernel-weighted
  sum of squared residual cross-correlations (Hong 1996 form) is N(0,1) under independence, "each time
  series is allowed to possess short memory, long memory or anti-persistence". Builds on Chen & Deo
  (2004, ET) for goodness of fit.
- Tsay (Econometric Theory 16(3), 2000; A, search snippet only): in trend regressions with
  fractionally integrated errors, the t-ratio of the fractionally-differenced (FD) estimator "never
  diverges", so FD removes the spurious effect there. Different setting; noted as a filtering precedent.
*Effect on H10:* fractional prewhitening is **known for linear cross-correlation tests** with a
theorem. A kernel (HSIC) version on prewhitened series: not found. H10 stays a sensible baseline, and the
linear precedent is citable support. Referee "obvious" risk: high.

---

## Part B. Lit-review strands (surrounding literature)

### B1. Hermite-rank limit theory for U- and V-statistics
Classical and well identified, but **not re-read** here (I): Breuer–Major (1983), Dobrushin–Major (1979),
Taqqu (1979), Dehling–Taqqu (1989) empirical process with an application to U-statistics,
Dehling–Taqqu (1991) bivariate symmetric statistics of LRD observations, Hsing–Wu (2004) weighted
U-statistics for stationary processes (covers LRD linear processes), Lévy-Leduc et al. (2011)
U-processes under LRD, Beutner–Zähle (2014) continuous-mapping approach to U/V-statistics (limits can
combine several Hermite ranks — S, via Betken–Wendler), Arcones (1994) vector Gaussian functionals,
Bai–Taqqu (2013) multivariate limits. Read in full: Dehling–Rooch–Wendler (2017).
**Take-away.** H0 is the standard map. A referee will ask why kernel statistics need a new proof. The
honest answer is narrow: MMD/HSIC are degenerate V-statistics in the i.i.d. null, and the LRD expansion
of their leading term (rank 1 for MMD homogeneity, (1,1) for HSIC) has to be computed per kernel and
per G. That is a corollary-type contribution unless the above-line law and a calibration come with it.

### B2. Self-normalisation under LRD
Established for scalar statistics: CUSUM (Shao 2011), Wilcoxon (Betken 2016), LMSV change in volatility
(Betken–Kulik 2017, F), rank-based tests (Betken–Wendler 2022, A), self-normalised linear processes with
heavy tails (Peligrad–Sang; A), heavy-tailed SN means (McElroy–Politis 2007; Jach et al. 2012). In all
of them the SN limit is a functional of fBm (or a Hermite process) and **depends on H**. Two ways out:
tabulate in H and plug in Ĥ (Betken 2016), or subsample / block-sample the SN statistic (Betken–Wendler
2018; Bai–Taqqu–Zhang 2016, which covers Hermite rank ≥ 1 including Rosenblatt limits). Kernel SN under
LRD: not found. Ghoshal–Shao (2026) is the kernel SN method under mixing.

### B3. Subsampling under LRD
Validity for subordinated Gaussian processes: Hall–Jing–Lahiri (1998), Nordman–Lahiri (2005, ET) for
linear processes, Zhang–Ho–Wendler–Wu (2013), Betken–Wendler (2018) for general statistics with
l = O(n^{(1+D)/2−ε}), Bai–Taqqu (2017, AoS) with b = o(n) irrespective of memory. Betken–Dehling's
Thm 2.4 applies Betken–Wendler to dCov; it presumes a statistic that converges after normalisation, and
the normalisation differs on the two sides of the line (Rem 2.1(i)). So subsampling with the n-rate is
valid only where the n-rate is right — exactly H4a's premise (F). Rate estimation: Bertail–Politis–Romano
(1999). Estimated normalisation across regimes: Zhai–Cheng–Wu (2026).

### B4. Block and wild bootstrap under LRD
Block bootstrap: may fail for non-central limits (Lahiri 1993 per Betken–Wendler; Tewes 2018); valid for
rank-2 covariance statistics below 3/4 (ZZW/ZKW); block-size behaviour depends on non-linearity beyond
the Hermite rank (Zhang–Lahiri–Nordman 2022). Wild bootstrap: no LRD theory found; the τ-mixing
conditions of CSG 2014 exclude all H > 1/2 (prior job F); centred multipliers over-reject even under
mixing (CSG 2014 Table 1, F).

### B5. Finance memory estimates (update to `FME`)
New: Sato–Kanazawa (2023) puts TSE order-sign γ near 0.5 (H ≈ 0.75; S). Nothing new for durations,
spreads, volume or RV within budget. The earlier table (`sources/…/finance-memory-estimates.md`) stands;
its caveats (estimator sensitivity of 0.1–0.3 in d) stand too.

---

## Part C. Corrections to the hand-off file

| Item in hypotheses file | Correction | Evidence |
|---|---|---|
| H1b, H4a: "Hall, Horowitz & Jing (1998)" | Hall, **Jing & Lahiri** (1998), Statistica Sinica 8(4) | I (BW 2018 references, F) |
| H4b: Ghoshal–Shao "absolute regularity — β-mixing" | v2 and EJS: **strong (α-) mixing**, Σα^{δ/(2+δ)} < ∞; Brownian FCLT | F |
| H4b: "dependence on H is mild" | Not mild for SN-Wilcoxon: 5% CV +42% from H = 0.6 to 0.9 | F (Betken 2016 Table 1) |
| H4b: "Shao's LRD SN limit depends on H" (R) | Confirmed | S (BTZ 2016 quoting Shao 2011) |
| H7: Zhai–Kolar–Wu "not verified" | Verified, IEEE TIT 72(6):4246–4296, doi 10.1109/TIT.2026.3685246 | Crossref + F |
| H8: "Zhai–Cheng–Wu (R)" | Two arXiv papers 2609.33866 / 2609.33868; Rosenblatt regime for \|W\| coverage proved | F / A |
| H3: "Betken–Dehling Thm 2.2(ii) gives this for G = id" | Confirmed, explicitly as a squared scalar (Cor 2.1(ii)) | F |
| H2 "Robinson–Hidalgo product-memory rule (lit-A)" | Explicit sum rule is Tsay–Chung (2000) | S (verbatim) |
| H10 "I do not know a long-memory kernel version" | A long-memory **linear** version exists: Shao (2009) | F |
| H5 cause list | Cause 1 supported by CSG 2014's own Table 1 | F |
