# Evidence matrix — kernel tests under long memory

Execution `claude-mb3-lit-20261016` · 2026-10-16. Template: deep-research
`templates/literature_matrix_template.md` (extended matrix + convergence summary), adapted to
mathematical-statistics evidence. "Level" in that template is a medical design pyramid. For theorems it
does not fit, so I grade **Fitness** instead (the template's discipline-relative rule):
**T** = refereed theorem; **Tp** = theorem in a preprint (unrefereed); **Sim** = simulation evidence in the
paper; **Sec** = statement read second-hand.
Class: **F / S / A / I** as in `search-log.md`. Keys refer to `bibliography.bib`.
"+" supports, "−" contradicts, "~" partly / qualifies, "·" not addressed. Hypothesis IDs from
`hypotheses/hypotheses-and-claims.md`.

## Extended matrix (38 extracted records)

| # | Key | Class | Fitness | Dependence setting | Statistic | Calibration | What it shows (relevant part) | Hypotheses |
|---|---|---|---|---|---|---|---|---|
| 1 | betken2021ldr | F | Tp + Sim | two independent subordinated Gaussian LRD series | dCov (ECF process in L²(w)) | subsampling, l_n = O(n^{(1+min D)/2−ε}) | (i) both D>1/2: √n Gaussian process, n·dcov → ∫∫\|Z\|²w; (ii) G=id, D_ξ+D_η<1: limit = [double Wiener–Itô integral]² × ∫∫s²t²e^{−(s²+t²)}w — a squared scalar times a constant; mixed case not covered (Rem 2.1(i)); sims H∈{0.6,0.7} only, size 7–17.5% | H2 ~, H2b + (gap confirmed), **H3 + (G=id, dCov)**, H4a ~ |
| 2 | betken2024bootstrap | A | Tp | α-mixing | dCov | block bootstrap | no LRD | H4c · |
| 3 | dehling2017twosample | F | T | subordinated Gaussian LRD | two-sample U-statistic process | none (limit theory) | Hermite rank m of kernel h; when Dm<1 the process is dominated by the rank-m term and the normalisation uses d_n² ~ n^{2−Dm} (eq. 9) — non-central limit | H1 + (machinery), H0 + |
| 4 | zhai2026ggmlong | F | Tp + Sim | high-dim Gaussian linear, β>3/4 | covariance-based GGM statistic | overlapping block bootstrap l=n^{2/3} | valid for H<3/4; "beyond the scope" for 1/2<β≤3/4 (Rosenblatt); fails at β=0.55; **no forthcoming announcement** | H0 + (rank-2 line), H2 ~, H7 ~ |
| 5 | zhai2026simultaneous | F | T (IEEE TIT) + Sim | Gaussian linear, β>3/4 | sample covariance / precision | block sampling (finite-sample) | Gaussian approximation only if β>3/4; states Rosenblatt above; fails at β=0.55 | H0 +, H4a ~ |
| 6 | ghoshal2026rkhs | F | T (EJS) + Sim | strong mixing, Σα^{δ/(2+δ)}<∞ | RKHS split–project–SN (GOF, change point, two-series independence) | self-normalisation, Brownian FCLT | pivotal limits under mixing; **no mention of long memory**; future work: high-dim / non-stationary | H4b ~ (method exists, LRD absent) |
| 7 | dizcastro2026hsic | F | Tp | L1-NED on β-mixing | HSIC (topological spaces) | asymptotic | no LRD | H2 · |
| 8 | chwialkowski2014wild | F | T (NeurIPS) + Sim | τ-mixing (τ(r)=O(r^{−6−ε}) for HSIC) | MMD, HSIC V-statistics | wild bootstrap, W_t = e^{−1/l}W_{t−1}+…; centred W̃ | **centred ("artificially degenerated") multipliers keep Type I error above 5% (0.07–0.145), dropping with n; uncentred Vb1 accurate**; two-sample uses two independent multiplier chains | **H5 +** (cause 1, also 3), H4d ~ |
| 9 | chwialkowski2014random | A | T | φ/β-mixing | HSIC | shift test | short memory only | H6 · |
| 10 | shao2009portmanteau | F | T (ET) + Sim | short / long memory / anti-persistent linear processes | frequency-domain residual cross-correlation (Hong-type) | prewhitening by estimated spectral density (parametric or FAR(p,d) sieve) → N(0,1) | **prewhitening restores a standard null under long memory for the linear cross-correlation test** | **H10 + (linear)**, H2 ~ |
| 11 | tsay2000spurious | S | T (JoE) | independent stationary/nonstationary I(d) | OLS t-ratio | standard | "as long as [the] orders of integration sum up to a value greater than 1/2, the t-ratios become divergent" | **H2 + (linear sum rule is classical)** |
| 12 | ventosa2022spurious | F | T + Sim | multivariate I(d) | OLS t / F | standard | d_y+d_x<1/2: t not divergent but not N(0,1) in finite samples; =1/2: diverges at (log T)^{1/2}; >1/2 diverges | H2 +, H4 ~ (below-line distortion) |
| 13 | robinson1997regression | A | T (AoS) | LRD errors and regressors | LS / GLS | — | LS "need no longer be n^{1/2}-consistent" when singularities coincide; GLS-type estimates stay √n | H2 + (linear), H10 ~ (whitening idea) |
| 14 | ernst2025fbmindep | A | Tp | fBm / BM paths | path-level independence tests | exact / asymptotic | independence testing for fBm paths (non-stationary) | H2 · (adjacent) |
| 15 | zhai2026coverage | F | Tp + Sim | Gaussian linear residuals, 1/2<β<1 | realised coverage of split conformal, score \|W\| | **block sampling with estimated normalization** | coverage indicator cov ~ h^{2−4β}: Gaussian if β>3/4, log-normalised at 3/4, **Rosenblatt if β<3/4**; MBB and HAC oversized (13–44%) at β≤0.75; BS 13.9% at n=2^{12}, 5.7% at 2^{20} (β=0.6) | **H8 + (symmetric score → 3/4 line, known)**, **H4a + (known for scalars)**, H4c + |
| 16 | zhai2026splitconformal | A | Tp | FDM incl. Gaussian linear LRD | CQR / median conformal | — | interval-length theory under long memory | H8 ~ |
| 17 | shao2010sn | A | T (JRSSB) | weak dependence | SN for functionals | SN | baseline SN | H4b · |
| 18 | shao2011lrdcp | A + S | T (JTSA) | LRD | SN CUSUM | SN | limit "depends on the unknown Hurst index" (as reported in #22, F) | **H4b + (H-dependence confirmed)** |
| 19 | shao2015review | A | Review | various incl. long memory | SN | — | covers "inference for long memory time series"; kernel/U-statistic SN under LRD not visible in abstract | H4b · |
| 20 | betken2016sn | F | T (JTSA) + Sim | subordinated Gaussian LRD | SN Wilcoxon change-point | SN with H-tabulated critical values | limit is a functional of fBm B_H; 5% critical value 7.28 (H=0.6), 8.19 (0.7), 9.50 (0.8), 10.33 (0.9); size 4.8–5.5% for n≥50 with true H | **H4b ~ (H-dependence is not "mild": +42% from 0.6 to 0.9)** |
| 21 | betken2018subsampling | F | T (Stat Sin) + Sim | subordinated Gaussian LRD | general statistics; SN Wilcoxon | subsampling, l = O(n^{(1+D)/2−ε}) | subsampling an SN statistic removes the need to estimate H; reviews Hall–Jing–Lahiri, Nordman–Lahiri, Jach et al., Lahiri 1993 (MBB "might fail") | **H4a ~ (SN+subsampling is the H-free alternative)**, H4c + |
| 22 | bai2016snbs | F | T (SPA) + Sim | non-linear functions of Gaussian, SRD or LRD, heavy tails; Hermite rank ≥1 incl. Rosenblatt | SN block sums of the mean | SN block sampling | unified, no a-priori estimate of H; quotes Shao (2011): SN alone does not remove dependence on H | **H4b ~/H4a ~ (H-free route exists for scalars)** |
| 23 | bai2017resampling | A | T (AoS) | LRD Gaussian-based | general statistics | subsampling, b=o(n) | consistency irrespective of memory strength | H4a ~ |
| 24 | tewes2018blockbootstrap | A | T (JTSA) | LRD | empirical process | block bootstrap | bootstrap limit is semi-degenerate with Gaussian random part → "might fail" when the original limit is non-central | **H4c +** |
| 25 | zhang2022optimalblock | A | T (AoS) + Sim | Gaussian-subordinated LRD | block variance estimators | block resampling | properties depend on non-linearity beyond Hermite rank; block order n^{1/2} often reasonable | H4a ~, H1b ~ (block rule choice) |
| 26 | jach2012subsampling | A | T (JTSA) + Sim | heavy-tailed LRD, rank-1 g | SN mean | subsampling without tail or memory knowledge | H-free subsampling for an SN mean | H4a ~ |
| 27 | bertail1999unknownrate | A | T (JASA) | general | any | subsampling with **estimated rate** | rate estimation by subsampling is classical | **H4a ~ (prior art for rate correction)** |
| 28 | betken2017lmsv | F | Tp | LMSV LRD | SN change-point (rank-type) | SN + subsampling | restates that Shao (2011) adopted SN for CUSUM; "better size, less power" with SN | H4b ~ |
| 29 | bai2013multivariate | A | T (JTSA) | vectors of LRD Gaussian functionals | normalised sums | — | mixed Gaussian / Hermite limits; mixed ranks 1 and 2 treated; general mixed case conjectured | H2b ~, H3 ~ |
| 30 | bai2018instability | A | T (Stat Sci) | LRD | sample variance, empirical process, Whittle | — | **ranks other than 1 are statistically unstable**; assuming rank ≥2 under-estimates fluctuations | **H1c − (caution)**, H8 − (caution) |
| 31 | dehling1991bivariate / hsing2004weighted / levyleduc2011uprocess / beutner2014cma | I | T | LRD | U- and V-statistics (incl. degenerate) | — | classical limit theory exists for U/V-statistics under LRD (content not re-read) | H0 + (R items now identified) |
| 32 | shi2026gtcm | A | Tp | weak dependence | GTCM residual covariance | — | no LRD | H7 · |
| 33 | shah2020gcm | A | T | i.i.d. | GCM | — | no dependence | H7 · |
| 34 | poskitt2007ar | I | T | fractionally integrated | AR(p) approximation | — | exists; content not read | H7b (reading gap) |
| 35 | sato2023lmf | A + S | T-empirical (PRL) | TSE 2012–2020 order flow | sign ACF exponent γ | — | γ ≈ α−1 with metaorder tail α ≈ 3/2 (S) → γ ≈ 0.5, H ≈ 0.75 | H9 ~ (signs on the line) |
| 36 | toth2015persistent | A | empirical (JEDC) | LSE order flow | sign ACF | — | persistence "out to tens of thousands of orders"; exponent not extracted | H9 · |
| 37 | tewes2017localalt | A | T (EJS) | subordinated Gaussian LRD | KS / CvM change-point | — | theory covers a change in Hermite rank | H1c ~ |
| 38 | lahiri1993mbb | I + S | T | LRD | sample mean | MBB | MBB may fail under LRD (as reported in #21) | H4c + |

(Row 31 bundles four identity-only theory references; it counts as one extracted record.)

## Convergence summary

| Theme | For | Against / qualifying | Net | Strength | Confidence |
|---|---|---|---|---|---|
| T1. Quadratic and cross statistics change limit at the rank-2 / sum line | #1, #4, #5, #11, #12, #15 (6; F/S, refereed and preprint) | none | +6 | Strong (for linear, covariance, dCov-Gaussian and coverage indicators) | High that the line is real for those statistics; **no direct evidence for HSIC** |
| T2. The linear (Pearson/OLS) sum rule d_X+d_Y = 1/2 is classical | #11, #12, #13 | — | +3 | Strong | High |
| T3. Above the line the cross-statistic null is a squared scalar | #1 (G=id, dCov, explicit) | none found for non-linear G or HSIC | +1 | Moderate | Medium |
| T4. Self-normalised statistics under LRD have H-dependent limits | #18, #20, #22, #28 | — | +4 | Strong | High |
| T5. H-free calibration of scalar LRD statistics exists (SN + block sampling / subsampling) | #21, #22, #26, #23 | finite-sample oversize at small n (#15 Table 2) | +4 | Strong (scalars) | High for scalars; none for kernel statistics |
| T6. Estimated-rate block sampling works across the 3/4 transition | #15 (coverage), #27 (general), #5 (below only) | — | +2 | Moderate | Medium (one 2026 preprint across regimes) |
| T7. Block bootstrap fails under non-central limits | #24, #38, #15 | none | +3 | Moderate | Medium |
| T8. Centred dependent-wild-bootstrap multipliers over-reject at moderate n | #8 (and Leucht–Neumann Sec 5 as cited there) | — | +1–2 | Moderate | Medium |
| T9. Prewhitening restores standard nulls for linear cross-correlation tests under LRD | #10 | #12 (below-line finite-sample distortion without whitening) | +1 | Moderate | Medium (linear only) |
| T10. Higher Hermite ranks (even transforms) are reliable in practice | — | #30 (instability of ranks) | −1 | Contested | Low |
| T11. Kernel / HSIC / MMD / KCI test results under LRD | #1 (dCov only) | — | — | **Gap (not found in this bounded scan)** | — |

## Knowledge gaps (template form)

| Gap | Type | Implication | Priority |
|---|---|---|---|
| No HSIC (general kernel) null limit under LRD found | Theoretical | H2/H3 for HSIC stay worker inference | High |
| Mixed case (one D<1/2, sum >1) unproved even for dCov | Theoretical | H2b open | High |
| Non-linear G above the line unproved even for dCov | Theoretical | H3 for G3 (signs) open | High |
| No SN kernel statistic under LRD found | Methodological | H4b combination open; ingredients public | High (scoop risk) |
| No kernel CI / Granger test under LRD found | Theoretical + empirical | H7 open | Medium |
| No size study of kernel tests across a grid of H found | Empirical | the track's simulations are the only map found | High |
| SSRN and several journals not effectively reachable | Search | reading-gap risk remains | Medium |
