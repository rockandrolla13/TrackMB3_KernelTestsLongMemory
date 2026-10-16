# Search log — workstream `memory-kernel`

Execution `claude-mkc-kernel-20261015` · 2026-10-15 (UTC evening).
Tools: WebSearch (standard; two "extended"), WebFetch (arXiv abstract pages, arXiv listing search,
CrossRef API), curl + pdftotext for full-text section reads (PDFs kept in the session scratchpad, not in
the repo). No data downloaded. No Semantic Scholar/OpenAlex calls.

Skills applied: deep-research lit-review protocol (search → verify → synthesise → devil's-advocate
checks done inline), paper-lit-finder (2025–26 recent-paper scan, programme-aware limitation notes),
research-gap-finder template (gap-map.md structure).

## Reuse read first (R, not counted in budget)

- research_assistant/docs/topics/track_e_kernels_causality/: E1_kernel_tests_dependence.md (full),
  GAP_ANALYSIS.md (full), LIT_REVIEW.md (taxonomy §3), E2 (summary table, grep for memory terms).
- research_assistant/docs/topics/conformal_causal_fin_ts/: S2 (long-memory entries #10–12, gaps),
  S8 (summary table, gaps).
- ResearchHarness programme/catalog.yaml: B1, B2-HTE, B2-COMP, E, F-CAUSAL, F-KERNEL, G, J, K.
- runs/JOB-A-HAWKES-CHECK/handoff.md; runs/JOB-A-P2-CALIB-20261015/handoff.md + calibration.md (grep).

## Queries (in order)

| # | Engine | Query | Useful hits |
|---|---|---|---|
| Q1 | WebSearch | kernel independence test long-range dependence HSIC distance covariance long memory arXiv 2025 | 2107.03041, 2112.14091 (known) |
| Q2 | WebSearch | MMD two-sample test long-range dependent time series non-central limit theorem | 1404.0551 |
| Q3 | WebSearch | signature kernel fractional Brownian motion rough paths Hurst estimation 2025 arXiv | 2506.16156, 2512.09588, 2512.04646, 2311.13071 |
| Q4 | WebSearch | nonparametric RKHS estimation Hawkes process excitation kernel power-law 2025 arXiv | 2411.00621, 1112.1838, 1401.0903 |
| Q5 | WebSearch | Hawkes process power-law kernel estimation near-critical nonparametric 2025 2026 | 1302.1405, 1403.5227 (old only) |
| Q6 | WebSearch | propagator model nonparametric estimation impact kernel kernel regression 2025 | 2301.05157 |
| Q7 | WebSearch | signature kernel two-sample test stochastic processes time series dependence 2025 2026 | 2506.01718, 1810.10971, 2305.04625, 2407.19848 |
| Q8 | WebSearch | conditional mean embedding dynamical systems Koopman operator kernel non-Markovian memory 2025 | 1712.01572, 2405.07432, 2002.03689 |
| Q9 | arXiv listing | Hawkes kernel RKHS | 2510.08916, 2411.00621 |
| Q10 | arXiv listing | long memory kernel test | 47 results, none on topic (off-field titles not screened); 2403.11773 noted |
| Q11 | arXiv listing | propagator price impact kernel estimation | 2510.06879, 2309.02994, 2301.05157 |
| Q12 | arXiv listing | long-range dependence independence test | 2107.03041, 2110.08089 |
| Q13 | arXiv listing | signature kernel rough | 2510.16728, 2509.07893, 2606.28869, 2605.08844, 2502.20392, 2501.06758, 2501.09683 |
| Q14 | arXiv listing | Hawkes nonparametric estimation kernel | 2510.08916, 2506.24065 |
| Q15 | arXiv listing | kernel instrumental variable time series | nothing on topic |
| Q16 | arXiv listing | signature kernel conditional independence | 2402.18477, 2502.08470 |
| Q17 | arXiv listing | Hawkes Granger causality nonparametric | 1605.06759, 2102.02741 (old) |
| Q18 | WebSearch (extended) | kernel conditional independence test time series "long memory" OR "long-range dependence" 2026 | **2609.30565**, 2505.11085 |
| Q19 | WebSearch | kernel proxy causal learning time series dependent data beta-mixing 2025 2026 | 2105.04544, 2602.13531 |
| Q20 | WebSearch | wild bootstrap kernel test fails long memory fractional ARIMA HSIC size distortion | nothing direct |
| Q21 | WebSearch | nearly unstable Hawkes process power-law kernel estimation rough volatility statistical inference 2025 | 2312.08784, Jaisson–Rosenbaum |
| Q22 | WebSearch | maximum mean discrepancy estimation long memory fractional process minimum distance estimator | nothing kernel-based |
| Q23 | WebSearch | kernel conformal prediction long-range dependence time series coverage 2026 | 2405.16828, 2510.02471 |
| Q24 | WebSearch | Hawkes process goodness-of-fit test kernel MMD or signature 2025 | 2410.05008, 2006.09439, 1709.08516 |
| Q25 | WebSearch | causal discovery multivariate Hawkes processes nonparametric kernel Granger 2025 2026 | 1607.06333 (old) |
| Q26 | WebSearch | signature kernel causal discovery stochastic processes 2025 2026 follow-up Manten Salvi Kilbertus | none beyond 2402.18477 |
| Q27 | WebSearch | fading memory RKHS universal approximation time series kernel 2025 | 2403.11945, 2212.14641, 2408.07386 |
| Q28 | WebSearch (extended) | "Hawkes" "power-law kernel" estimation consistency MLE heavy-tailed 2025 OR 2026 | 2503.14814, 2510.08085, 2504.16737, 2512.07563 |
| Q29 | WebSearch | kernel mean embedding causal effect time series dynamic treatment 2025 2026 Gretton | 2503.05024, 1809.08560 |
| Q30 | WebSearch | HSIC independence test fractional Gaussian noise Hermite rank non-central limit | 1902.06441, 2206.11607 (i.i.d.) |
| Q31 | WebSearch | learning memory kernel generalized Langevin equation RKHS nonparametric 2025 | 2402.11705 |
| Q32 | WebSearch | nonparametric IV RKHS dependent data mixing time series kernel IV rates 2025 | 2411.19653, 2511.21603, 2603.04681 |
| Q33 | WebSearch | signature roughness test Hurst parameter rough volatility estimation 2025 2026 | 2607.06340, 2401.02819, 2508.07151, 2509.20015 |
| Q34 | WebSearch | self-normalization long-range dependence robust test ... kernel or U-statistic | 1403.0265, 1509.05720, 2004.06574 |
| Q35 | WebSearch | Betken long-range dependence U-statistics kernel two-sample change-point subsampling 2024 2025 | same as Q34 |
| Q36 | WebSearch | Hawkes long memory variance of counts power-law kernel mixing rate polynomial beta-mixing 2025 | 2607.20838 |
| Q37 | WebSearch | Koopman operator regression RKHS learning rates beta-mixing Kostic Pontil | 2205.14027, 2302.02004 |
| Q38 | WebSearch | kernel ridge regression dependent data long memory learning rate RKHS | JMLR 17-306 (kernel PLS LRD), ICML 2024 KRR τ-mixing, 2510.15363 |
| Q39 | WebSearch | kernel Granger causality realized volatility long memory nonlinear 2025 2026 | 0711.2643, 2609.22007 (known) |
| Q40 | WebSearch | self-normalized test long memory distance covariance OR HSIC OR kernel independence | none new |
| Q41 | WebSearch | Hawkes near criticality branching ratio estimation inference kernel nonparametric 2026 | 1810.03730; unverified working paper (below) |
| V1–V12 | WebFetch arXiv abs / CrossRef | identifier checks: 2510.06879, 2510.08916, 2609.30565, 2607.06340, 2606.28869, 2511.21603, 2607.20838, 1401.0903, 1504.03100, cond-mat/0311053, 2205.14027, 2512.09588, 2212.14641, 2504.16737, 2603.04681; CrossRef 10.1111/1468-0262.00418; Bernoulli 10.3150/21-BEJ1338 (via search) | all resolved |

## Budget used

- **Screened: 70 records** (titles/abstracts assessed for relevance; off-field titles in arXiv listing
  noise — e.g. Q10's GPU/LLM papers — not counted). Cap 70: at cap.
- **Extracted: 25** (records with fields filled in field-map.md: the 12 full-text reads plus
  2609.30565-linked, 1401.0903, 1504.03100, 2312.08784, 2607.20838, 2301.05157, 1810.10971, 2006.14794,
  2205.14027, 2105.04544, 2511.21603, 2402.11705, 1912.05737/BEJ1338, cond-mat/0311053, 2607.06340,
  minus overlaps). Cap 25: at cap.
- **Full-text reads: 12** (sections: abstract/intro, assumptions, main theorems; proofs not checked):
  2107.03041, 1408.5404, 2510.06879, 2510.08916, 2411.00621, 2609.30565, 2402.18477, 2403.11945,
  2601.11233, 2405.16828, 1404.0551, 2606.28869. Cap 12: at cap.

## Inclusion / exclusion

Included: theorem-level or method papers on RKHS or memory-kernel methods with stated dependence
assumptions; empirical finance papers that estimate memory kernels. Excluded: pricing/hedging uses of
signatures (2501.06758, 2501.09683, 2508.07151 — noted only), deep-learning Hawkes without theory,
smoothing-kernel (non-RKHS) papers (2603.04681), quantum/reservoir ML (2602.13531).

## Leads seen but not verified (not cited)

- A report that Hawkes MLE is biased near criticality (search snippet attributing it to a 2018 J. Stat.
  Comput. Simul. paper; phbern.ch repository link). Not opened.
- A 2025 working paper on branching-ratio estimation for Hawkes processes with an Itô-semimartingale
  baseline, including a near-criticality test (keio.ac.jp PDF link, search snippet). Not opened.
  → pass both to the memory-causal worker / coordinator.

## Engine notes

- WebSearch returns mostly pre-2024 papers for long-memory queries; the 2025–26 update depended on arXiv
  listing search, which is AND-of-terms and misses synonyms.
- SSRN, JFEc/JoE/JBES and Bernoulli/SPA were not searched directly. This is the main reading-gap risk.
