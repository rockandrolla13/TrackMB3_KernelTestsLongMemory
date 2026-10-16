# Systematic review — calibration of kernel two-sample, independence and conditional-independence tests under long-range dependence

Execution `claude-mb3-lit-20261016` · Agent 2 · 2026-10-16. Mode: deep-research `systematic-review`
(protocol `references/systematic_review_protocol.md`, applied by one worker; the skill's agents could not
be spawned, so their roles — question, protocol, search, verification, appraisal, synthesis, devil's
advocate — are done inline and marked). Evidence labels F / S / A / I as in `search-log.md`.
**Bounded review. Not finding a paper is not evidence that it does not exist.** No novelty claim.

## 1. Protocol

### 1.1 Question (PICOS adapted to methods research)

| Element | This review |
|---|---|
| Population (data) | Stationary time series with long-range dependence (ρ(k) ~ k^{−D}, D ∈ (0,1); H ∈ (1/2,1)), or explicitly non-mixing series; Gaussian-subordinated, linear, or general |
| Intervention (method) | Kernel or distance-based tests: MMD / energy distance (two-sample, homogeneity, change point), HSIC / distance covariance / characteristic-function tests (independence), kernel or residual-covariance CI tests (KCI, GCM-type, kernel Granger) — and, as the linear comparator, cross-covariance / regression t-tests |
| Comparator | The same tests under i.i.d. or mixing data |
| Outcome | Null limit (rate, law) and validity of a calibration (wild bootstrap, block bootstrap, subsampling, self-normalisation, shift/permutation, prewhitening); empirical size |
| Study design | Theorems (refereed or preprint) and simulation studies; any year; any language found |

Primary question: *Where do these tests lose size under long memory, and is the loss governed by the
Hermite rank of the statistic's leading term (rank 1: any H > 1/2; rank 2: H = 3/4; cross (1,1):
H_X + H_Y = 3/2)?*

### 1.2 Eligibility

Include if (I1) the paper studies a test statistic of the kinds above, and (I2) it gives theory or
simulation under LRD / long memory / an explicitly non-mixing class.
Exclude: **E1** only i.i.d. or mixing / NED-on-mixing dependence; **E2** not a test of distributions
or dependence (e.g. comparing memory parameters); **E3** content could not be determined (no access,
abstract silent); **E4** off-question; **E5** scalar rank/CUSUM statistic (reviewed in the lit-review
strand in `review.md`, not in this systematic strand).

### 1.3 Sources and strategy

arXiv API (field search), WebSearch (open web; domain filters for Project Euclid / IMS / Bernoulli,
ScienceDirect, Wiley, T&F, Cambridge, OUP, NeurIPS / PMLR / JMLR / OpenReview, SSRN), Google Scholar
(parsed result pages), Crossref and OpenAlex (verification, abstracts), Semantic Scholar (forward
citations of Betken–Dehling, Dehling–Rooch–Wendler, Betken–Wendler, ZZW). Prior-track sources in
`sources/` were reused as seeds and re-checked where the hand-off asked. Exact strings and counts:
`search-log.md`. Date range: no lower limit; up to 2026-10-16.

Core Boolean concept (adapted per engine, since arXiv is AND-of-terms and WebSearch is free text):
`("long-range dependen*" OR "long memory" OR fractional OR Hurst OR "non-mixing")
AND (HSIC OR "Hilbert-Schmidt" OR "distance covariance" OR "maximum mean discrepancy" OR MMD OR
"energy distance" OR "characteristic function" OR "conditional independence" OR "Granger" OR
"cross-correlation") AND (test OR bootstrap OR subsampling OR "self-normali*")`.

### 1.4 Appraisal and synthesis plan

RoB 2 and ROBINS-I are for trials and observational studies. They do not fit theorems. I appraise each
included paper on four methods-specific domains (devised for this review, stated before scoring):
(a) **assumption realism** for market series (subordination, Gaussianity, linearity); (b) **proof
status** (refereed vs preprint); (c) **simulation grid near the boundary** (are H values close to the
predicted line simulated?); (d) **calibration feasibility** (does the paper give a usable critical value
under LRD?). Meta-analysis: **not feasible** — there is no common effect size; outcomes are limit
theorems and a handful of size tables on non-comparable designs. A structured narrative synthesis
(SWiM-style grouping by statistic type, then by calibration) is given instead. GRADE is replaced by a
certainty statement per finding (high / moderate / low), based on (a)–(d) and consistency.

Protocol registration: not registered (internal exploratory review); recorded here before screening
counts were finalised. Deviation: protocol written at the start of this session but refined after the
first five searches (E5 added when scalar change-point papers kept appearing).

## 2. PRISMA 2020 flow (text figure)

```
IDENTIFICATION
  Records identified, core question (first-source attribution):
    arXiv API ............................ 17
    WebSearch (open + 9 domain filters) .. 27
    Google Scholar (parsed) ............... 3
    Citation chasing (Semantic Scholar) ... 3
    Prior-track seeds (sources/) .......... 11
                                   total  61
  Duplicates removed ..................... 13
                                           │
SCREENING                                  ▼
  Records screened (title/abstract) ...... 48
  Excluded at title/abstract ............. 32
     E1 mixing / i.i.d. only ......... 24
     E2 not a distribution/dependence test 2
     E3 content not determinable ...... 3
     E4 off-question .................. 1
     E5 scalar statistic (→ review.md)  2
                                           │
ELIGIBILITY                                ▼
  Reports sought for full text ........... 16
  Reports not retrieved .................. 5  (Robinson–Hidalgo 1997: Project Euclid fetch empty;
                                               Chen 2006, Bauer–Maynard 2012, Ernst et al. 2025,
                                               Zhai–Cheng–Wu split conformal: not attempted, budget;
                                               all five kept at abstract level, flagged A)
  Full texts assessed .................... 11 (10 direct + Tsay–Chung 2000 via a verbatim
                                               quote inside VSVL 2022, ScienceDirect paywalled)
  Excluded after full text ............... 3  (E1: Ghoshal–Shao 2026; Diz-Castro et al. 2026;
                                               Chwialkowski et al. 2014 — kept as calibration
                                               comparators in the synthesis)
                                           │
INCLUDED                                   ▼
  Included in narrative synthesis ........ 13
     full text (F) ...................... 8  (Betken–Dehling; DRW; ZZW; ZKW; Shao 2009;
                                              Zhai–Cheng–Wu coverage; Ventosa-Santaulària et al.;
                                              Tsay–Chung via verbatim quote in VSVL = S*)
     abstract only, flagged (A) ......... 5  (Robinson–Hidalgo; Chen 2006; Bauer–Maynard;
                                              Ernst–Viens–Yan; Zhai–Cheng–Wu split conformal)
  Meta-analysis .......................... 0  (not feasible)
```
*Tsay–Chung 2000 is counted as included at S level: its key sentence was read verbatim inside a full
text (VSVL 2022). It is listed once among the 8.

Whole-review budget (both strands): screened 103 / 120, extracted 38 / 40, full text 15 / 20.

## 3. Screening table (core strand, 48 records)

| # | Record | Decision | Reason / note |
|---|---|---|---|
| C1 | Betken & Dehling, arXiv:2107.03041v3 | **Include (F)** | dCov under subordinated Gaussian LRD |
| C2 | Betken, Dehling & Kroll, arXiv:2112.14091 | E1 | α-mixing |
| C3 | Ghoshal & Shao, EJS 2026 / arXiv:2601.13468 | FT → E1 | strong mixing; kept as comparator |
| C4 | Diz-Castro et al., arXiv:2604.28104 | FT → E1 | NED on β-mixing |
| C5 | Dehling, Rooch & Wendler, Statistics 2017 | **Include (F)** | two-sample U-processes under LRD |
| C6 | Zhai, Zhong & Wu, arXiv:2609.30565 | **Include (F)** | GGM CI test, Gaussian linear LRD |
| C7 | Zhai, Kolar & Wu, IEEE TIT 2026 | **Include (F)** | covariance inference under LRD (block sampling) |
| C8 | Li, Xu & Zhou, AutoHSIC, arXiv:2605.22025 | E1 | serial independence, short memory |
| C9 | Wang, Li & Zhu, arXiv:1804.09866 | E1 | β-mixing |
| C10 | Chwialkowski & Gretton, ICML 2014 | E1 | φ/β-mixing |
| C11 | Chwialkowski, Sejdinovic & Gretton, NeurIPS 2014 | FT → E1 | τ-mixing; kept as comparator (H5) |
| C12 | Shen et al., TMLR 2024 (1908.06486) | E1 | block permutation; stationarity + mixing-type assumptions |
| C13 | Kernel tests for non-stationary random processes (2010.00271) | E1 | non-stationarity, not LRD |
| C14 | Ernst, Viens & Yan, arXiv:2510.23563 | **Include (A)** | independence of fBm / BM paths |
| C15 | Permutation-free kernel independence test (JMLR 2023) | E1 | i.i.d. |
| C16 | Shao, Econometric Theory 2009 | **Include (F)** | long memory allowed; prewhitened cross-correlation |
| C17 | Hong, Biometrika 1996 | E1 | short memory |
| C18 | Haugh, JASA 1976 | E1 | ARMA |
| C19 | Beering, J. Nonparametric Stat. 2025 | E3 | title only (HTTP 403, no abstract) |
| C20 | Shi & Shah, GTCM, arXiv:2609.20772 | E1 | weak dependence |
| C21 | dynamic GCM, arXiv:2504.21647 | E1 | non-stationary, not LRD |
| C22 | GKCM, PMLR 323 (2604.03721) | E1 | i.i.d. |
| C23 | Shah & Peters, GCM, AoS 2020 | E1 | i.i.d. |
| C24 | Manten et al., signature-kernel CI, ICLR 2025 | E1 | i.i.d. replicate paths (prior F) |
| C25 | MMD of ordinal patterns, arXiv:1210.4903 | E1 | no LRD theory |
| C26 | Kernel change-point, arXiv:2206.01853 | E1 | piecewise i.i.d. |
| C27 | Li & Xie, M-statistic, NeurIPS 2015 | E1 | i.i.d. |
| C28 | Jiang & Bodenham, arXiv:2609.22545 | E1 | online kernel CPD; no LRD in snippet |
| C29 | GOF via distance covariance, arXiv:1903.00708 | E1 | mixing |
| C30 | Wang & Hong (2018) CF-based CI test | E1 | weak dependence (secondary snippet) |
| C31 | Betken et al., arXiv:2408.10570 | E1 | two-sample i.i.d. |
| C32 | Betken, Micali & Ruiz Marín, arXiv:2601.14223 | E3 | "broad class of stationary processes"; LRD not stated |
| C33 | Zhai, Cheng & Wu, arXiv:2609.33868 | **Include (F)** | coverage test under Gaussian linear LRD |
| C34 | Tsay & Chung, JoE 2000 | **Include (S)** | spurious regression, sum rule |
| C35 | Ventosa-Santaulària et al., CSTM 2022 | **Include (F)** | multivariate spurious regression |
| C36 | Robinson & Hidalgo, AoS 1997 | **Include (A)** | LS under coinciding singularities |
| C37 | Chen (2006), Granger via spectrum estimator, long memory | **Include (A)** | prior seed; linear |
| C38 | Bauer & Maynard, JoE 2012 | **Include (A/R)** | linear Granger robust to long memory |
| C39 | Hwu & Tsay, SSRN 7546859 (2026) | E4 | DiD inference; adjacent |
| C40 | Two-sample test for long-memory parameters, arXiv:0907.1787 | E2 | compares d |
| C41 | ECF two-sample test, CSDA 2008 | E1 | i.i.d. |
| C42 | Fast ECF two-sample tests, 2016 | E1 | i.i.d. |
| C43 | Empirical correlation of two AR(1), arXiv:2310.08575 | E2 | short memory distribution theory |
| C44 | Greaves-Tunnell PhD thesis (UW, 2021) | E3 | title only |
| C45 | Barber & Pananjady, arXiv:2510.02471 | E1 | β-mixing (conformal; H8 context) |
| C46 | Zhai, Cheng & Wu, arXiv:2609.33866 | **Include (A)** | conformal length under LRD |
| C47 | Betken & Kulik, arXiv:1706.06351 | E5 | scalar SN change-point (read in review.md strand) |
| C48 | Dehling, Rooch & Taqqu, SJS 2013 | E5 | Wilcoxon change-point (review.md strand) |

E3 = C19, C32, C44. E4 = C39. E5 = C47, C48. E1 = the remaining 24 title/abstract exclusions.

## 4. Appraisal of included studies (methods-specific domains)

| Study | (a) Assumption realism | (b) Proof status | (c) Grid near line | (d) Usable calibration | Overall |
|---|---|---|---|---|---|
| Betken–Dehling | Medium: subordinated Gaussian (i); Gaussian only (ii) | Preprint, "submitted to Bernoulli" | Low: H ∈ {0.6, 0.7} only | Subsampling, needs correct normalisation; block cap depends on D | Some concerns |
| DRW 2017 | Medium | Refereed | none (theory) | none | Low concern (theory only) |
| ZZW 2026 | Low: Gaussian linear | Preprint | Low: β ∈ {2, 0.9, 0.55} | block bootstrap, H < 3/4 only | Some concerns |
| ZKW 2026 | Low: Gaussian linear | Refereed (IEEE TIT) | Low: β = 0.9, 0.55 | block sampling, H < 3/4 | Low–some concerns |
| Shao 2009 | Medium: linear, parametric or sieve spectrum | Refereed | moderate | yes (N(0,1)) | Low concern (linear only) |
| Zhai–Cheng–Wu 2026 (coverage) | Low: Gaussian linear, \|W\| score | Preprint | Good: β ∈ {0.9, 0.75, 0.6, 0.7} | yes, estimated normalisation | Some concerns (preprint) |
| VSVL 2022 / TC 2000 | Medium: independent I(d) | Refereed | moderate | none (shows failure) | Low concern |
| A-level includes | n/a | — | — | — | Not appraisable |

## 5. Structured narrative synthesis

**Grouping.** By statistic type (rank 1 two-sample/homogeneity; cross (1,1) independence; rank 2
quadratic; conditional), then by calibration.

### 5.1 Rank-1 two-sample / homogeneity statistics (MMD, energy distance)
- No paper found that derives the null law of MMD or energy distance under LRD.
- The machinery exists: two-sample U-statistic processes under LRD with Hermite rank m of the kernel and
  non-central limits when Dm < 1 (DRW 2017, F). It is written for change-point kernels (Gauss test,
  Wilcoxon), not for MMD's degenerate-in-i.i.d. structure. Extending it to MMD is worker inference (D).
- Certainty that "MMD loses √n-calibration at any H > 1/2 when the rank-1 coefficient is non-zero":
  **moderate** — theory by analogy (D) plus the track's own exploratory simulation, not literature.

### 5.2 Cross (1,1) independence statistics (dCov, HSIC, Pearson)
- **Linear statistics: the sum rule is classical.** Tsay & Chung (2000): OLS t-ratios between
  independent I(d) series diverge "as long as [the] orders of integration sum up to a value greater
  than 1/2" (S, verbatim). VSVL (2022, F) refine: at d_x + d_y = 1/2 divergence at (log T)^{1/2}; below,
  bounded but not N(0,1) in finite samples. Robinson–Hidalgo (1997, A): LS loses √n when singularities
  coincide. **Certainty: high.**
- **dCov:** Betken–Dehling (F) prove the √n Gaussian regime only when **both** D > 1/2 (both H < 3/4),
  and the non-central regime only for G = id. The mixed case is explicitly uncovered (Rem 2.1(i)).
  Their simulation never goes above H = 0.7 and shows 7–17.5% size at 5% for n ≤ 1000.
  **Certainty that the dCov boundary is H_X + H_Y = 3/2 in general: low–moderate** (proved for equal
  sides only as two separate one-sided regimes; mixed case open).
- **The above-line law is a squared scalar (G = id).** Corollary 2.1(ii) of Betken–Dehling (F): the
  normalised dCov converges to [∫ (Hermite kernel) Z_X(dx) Z_Y(dy)]² × ∫∫ s²t² e^{−(s²+t²)} w(s,t) ds dt.
  The same double Wiener–Itô integral is the limit of the normalised empirical cross-covariance
  (their Thm 4.1, F). So, for Gaussian series above the line, dCov and the squared cross-covariance share
  one random variable. This is H3's point 1–2 for dCov with G = id. **Certainty: high (for that case).**
- **HSIC with general kernels under LRD:** not found. Diz-Castro et al. (2026) and Ghoshal–Shao (2026)
  are mixing/NED only (F). **Gap in the bounded scan.**
- **Prewhitening (linear).** Shao (2009, F): estimate each series' spectral density (parametric, or
  FAR(p,d) sieve) and test residual cross-correlations; the statistic is N(0,1) under independence for
  short memory, long memory or anti-persistence. **Certainty: high (linear).** No kernel analogue found.

### 5.3 Rank-2 quadratic statistics (covariance, |W| indicators)
- ZZW (F) and ZKW (F, refereed) prove Gaussian approximation and block bootstrap / block sampling only for
  β > 3/4 (H < 3/4) and say the limit is Rosenblatt-type above; both fail at β = 0.55 in simulation.
- Zhai–Cheng–Wu (F): coverage indicators of |W| scores have covariance h^{2−4β}; Gaussian for H < 3/4,
  log-normalised at 3/4, Rosenblatt above. **Block sampling with an estimated normalisation** is
  consistent in all three regimes (Thm 7); MBB and HAC are oversized 13–44% at β ≤ 0.75 (Table 2);
  block sampling itself is 13.9% at n = 4,096 and 5.7% at n ≈ 10⁶ for β = 0.6 (contrary finite-sample
  evidence kept). **Certainty: moderate–high** (one preprint, clear theorems and simulation).

### 5.4 Conditional-independence and Granger tests
- Kernel / GCM-type CI tests: all found are i.i.d. or weak dependence (E1).
- Linear/Gaussian: ZZW GGM test (β > 3/4). Linear Granger tests robust to long memory exist
  (Bauer–Maynard, Chen 2006; A only).
- **No kernel CI test under LRD found. Certainty about the gap itself: low** (bounded scan, SSRN
  unsearched).

### 5.5 Calibration methods (cross-cutting; detail in `review.md`)
| Calibration | Evidence under LRD | Certainty |
|---|---|---|
| Wild / dependent-wild bootstrap | No LRD theory found. CSG 2014 (F): centred multipliers over-reject (0.07–0.145) even under mixing | moderate (for the finite-sample defect) |
| Moving / overlapping block bootstrap | valid below the rank-2 line for covariances (ZZW/ZKW); Lahiri 1993 and Tewes 2018: fails or "might fail" for non-central limits; ZCW Table 2: oversized at β ≤ 0.75 | moderate |
| Subsampling with known rate | Betken–Dehling Thm 2.4, Betken–Wendler, Bai–Taqqu 2017: valid if the statistic is correctly normalised | high (scalars, dCov below line) |
| Subsampling / block sampling with estimated rate or SN | Bertail–Politis–Romano 1999; Zhang–Ho–Wendler–Wu 2013; Jach et al. 2012; Bai–Taqqu–Zhang 2016; ZCW 2026 | high for scalars; none for kernel statistics |
| Self-normalisation alone | limit depends on H (Shao 2011 via BTZ, F; Betken 2016 Table 1, F) | high |
| Prewhitening | Shao 2009 (linear) | high (linear only) |

## 6. Devil's-advocate checks (inline, Checkpoint 2)

- *Cherry-picking?* Contrary items retained: Betken–Dehling's below-line oversize; VSVL's below-line
  non-normality; ZCW's 13.9% block-sampling size at n = 4,096; Bai–Taqqu's rank instability; Ghoshal–
  Shao's strong size accuracy under mixing (it is the competitor, not a failure).
- *Confirmation bias toward the Hermite map?* The map is confirmed for linear, covariance, dCov-Gaussian
  and coverage-indicator statistics. For HSIC and MMD it is still inference. The review does **not**
  upgrade those to "known".
- *Absence-as-evidence?* Every "not found" is tied to named failed or partial sources (SSRN 403; Scholar
  4/10 empty pages; S2 429). These are reading gaps, not gaps in the literature.
- *Strongest counter-argument:* a referee could say the whole kernel story is a corollary of
  Hsing–Wu (2004), Dehling–Taqqu (1991) and Arcones (1994) limit theory for U/V-statistics and vector
  Gaussian functionals. These were identified but **not read** (I). That is the largest unread risk.

## 7. Limitations

Single reviewer; no dual screening. SSRN not searched; Scholar partial; several key classics read only
second-hand (Shao 2011, Lahiri 1993, Leucht–Neumann 2013, Kim–Nordman 2011) or not at all (Hsing–Wu
2004, Dehling–Taqqu 1991, Arcones 1994, Poskitt 2007). Proofs were not checked. AI-assisted review:
all searches, extraction and synthesis were done by an AI agent; identifiers were machine-verified.
