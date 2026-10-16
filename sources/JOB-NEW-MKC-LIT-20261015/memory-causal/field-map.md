# Field map — causal inference with memory (workstream memory-causal)

Job JOB-NEW-MKC-LIT-20261015 · execution `claude-mkc-causal-20261015` · 2026-10-15.
Plan SHA-256 6f6631ad…c99e3 (verified). Bounded scan; no novelty claims.

**Access labels.** F = I read the theorem/assumption section in the paper's full text
(PDF converted locally). S = summarising fetch of a publisher or arXiv page. A = abstract or
bibliographic record only. D = my own derivation (elementary; not checked by anyone else).
F(reused) = read at theorem level in an earlier ResearchHarness job and not re-read here.

**What "survives long memory" means here.** Long memory = autocovariance decaying like
k^{2d−1}, 0 < d < 1/2 (Hurst H = d + 1/2), or a self-exciting process at or near criticality,
or with power-law kernels. Such processes are not geometrically mixing; when α-mixing holds at
all, the coefficients decay polynomially and are often not summable. A result "survives" if its
proof does not use summable mixing or geometric ergodicity, or if its estimand is conditional on
the realised path.

---

## 1. The baseline that breaks: DML and orthogonal learning on one series

| Work | Access | Dependence assumption | Survives long memory? |
|---|---|---|---|
| Semenova, Goldman, Chernozhukov, Taddy (QE 2023), arXiv:1712.09988 | F(reused) | Geometric β-mixing; neighbours-left-out cross-fitting via Strassen coupling (Lemma 3.1) | No as stated. Coupling error is β(gap); with polynomial β it must beat T^{-1/2}, which fails for slow decay |
| Ballinari & Wehrli, arXiv:2411.10009v2 | **F** (Ass. 3, 6; §3.2) | Score process α-mixing with Σα(s)^{(β−2)/β} < ∞; single-series cross-fitting needs α(k_T)^ψ = o(T^{-1}) | **No, and the authors say so.** Their example (ii): for α(s) = O(s^{2d−1}), 0 < d < 1/2, "Assumption 6.5 is indeed never satisfied." They suggest transforming the data instead. Note: Ass. 3 is on the *score* process, not the raw data — this matters (see §6) |
| Yu, Chen, Du, Li, arXiv:2609.27531 (Sep 2026) | **F** (abstract, Ass. 2–3) | Geometric β-mixing; buffered block cross-fitting; IV for unmeasured confounding | No as stated |
| Ciganovic, D'Amario, Tancioni (Econometrics J. 2026), arXiv:2603.10999 | A | Stationarity + time-reversibility (reverse cross-fitting) | Not determined at abstract level |
| Balkus, Laith, Hejazi, arXiv:2601.10899v3 | **F** (Ass. 1–3, Thm 1, Cor. 4) | Generic: a CLT at some rate r_n plus nuisance MSE o(1/r_n). Time-series corollary proved only for **fixed m-dependence** | Open. Thm 1 allows a slower r_n, so it is not excluded in principle, but the time-series proof uses a bounded count of correlated pairs (Bienaymé bound), which long memory violates |
| Lewis & Syrgkanis (g-estimation DML), arXiv:2002.07285 | A (from S3) | i.i.d. units, multi-period within | Not applicable to one long series |

**Reading.** Every single-series DML result found assumes geometric or summable mixing. One paper
(Ballinari–Wehrli) states the long-memory failure explicitly. That statement is the clearest
"baseline that breaks" in this workstream.

## 2. Potential-outcome and design-based time series

| Work | Access | Assumption | Survives long memory? |
|---|---|---|---|
| Bojinov & Shephard (JASA 2019), arXiv:1706.07840 | F(reused) | Design-based: potential outcomes fixed, known assignment probabilities, non-anticipation. No outcome dependence assumption | **Yes.** Validity comes from the randomisation, not from mixing. Memory affects power, not size |
| Bojinov, Rambachan, Shephard (QE 2021); Rambachan & Shephard arXiv:1903.01637; Carlson & Shephard arXiv:2603.20394 | A (from S3) | Identification frameworks; general dependence | Identification survives; any estimation layer must be re-checked |
| Hu & Wager switchbacks, arXiv:2209.00197 | A (from S3) | Markov, geometric mixing; carryover bias rate depends on it | No as stated |
| Montiel Olea & Plagborg-Møller (Econometrica 2021), arXiv:2007.13888 | **F** (§2, Ass. 1–3, §4) | Finite-order VAR(p) with known p; innovations mean-independent of past and future | **Robust to near unit roots, not to long memory.** A fractionally integrated process is an infinite-order AR, outside the VAR(p) model. "Persistence-robust" here means local-to-unity, h_T/T → 0 |
| Baillie, Diebold, Kapetanios, Kim, Mora, arXiv:2203.04080 | A | HAC fails under strong autocorrelation; recommends dynamic regression | Background on why HAC is fragile; long memory not the focus |

## 3. Proxy, negative-control and instrument methods under memory

| Work | Access | What it does | Memory |
|---|---|---|---|
| Hao (PLOS One 2026), doi:10.1371/journal.pone.0347847 | S (first 100k characters of the article page) | Single series, no untreated units: a "common proximal variable" (default: time) plus Bernstein expansions and a conditional-mean-symmetry condition identify the effect. Ass. 4 allows long-range dependence with Hurst H; Thm 5 claims an fBm limit; Thm 4 claims bootstrap validity with a Hermite-rank condition (non-normal limit for rank ≥ 2) | **Explicitly targets long memory.** Single author, PLOS One, proofs in supplements I did not read. Using time as the proxy is a strong and unusual identifying device. Treat as an unverified claim, not as settled |
| Shi, Li, Miao, Hu, Tchetgen Tchetgen (proximal SC), arXiv:2108.13935; Ying et al. arXiv:2109.07030; Park & Tchetgen Tchetgen (JCI 2025) | A (from search + S3) | Proximal identification for time series / longitudinal data | No dependence rate theory for long memory found |
| Track J (ResearchHarness) | catalog | Kernel proxy causal estimation under β-mixing | Assumes β-mixing → **triple-intersection flag** (see §7) |

## 4. Granger causality and transfer entropy under long memory

| Work | Access | Result | Status |
|---|---|---|---|
| Tsay & Chung (J. Econometrics 2000) 96(1):155–182 | A (record) | Spurious regression between independent stationary fractionally integrated series | Known hazard: t-ratios can diverge without any causal link |
| Bauer & Maynard (J. Econometrics 2012) 169(2):293–300 | A | Surplus-lag Granger test with one χ² null across I(0), I(1), local-to-unity, **fractional**, and mean breaks in the forcing variable | Linear Granger testing under long memory is **largely handled** for this design |
| Chen (J. Forecasting 2006) doi:10.1002/for.981 | A | Spectral Granger estimator allowing long memory | Older; supports "linear case handled" |
| Barnett, Barrett, Seth (PRL 2009), arXiv:0910.4514 | A | Transfer entropy = Granger causality for Gaussian variables | So Gaussian long-memory TE inherits the linear results |
| Faes et al. (Entropy 2021), VARFI transfer entropy | A (search snippet; id not verified) | TE inside a VAR fractionally integrated model | Gaussian, parametric. Nonlinear / nonparametric TE under long memory not found |
| Runge et al. PCMCI / PCMCI+ (arXiv:1702.07007, 2003.03685) | A | Autocorrelation inflates false positives in CI-based discovery | No long-memory treatment found |

## 5. Self-exciting processes: Granger causality, interventions, near-criticality

**5a. Granger structure and inference for Hawkes.**

| Work | Access | Assumption | Near-critical / power-law? |
|---|---|---|---|
| Eichler, Dahlhaus, Dueck (JTSA 2017), arXiv:1605.06759 | **F** (Def. 2.3, 3.1, Thm 3.5, Thm 4.1) | Granger non-causality ⇔ link function φ_ij ≡ 0; global Markov property. Estimator consistency needs ‖∫φ‖ < 1 and ‖φ(u)‖ ≤ C/u | Allows heavy tails (up to 1/u) but gives **consistency only**: no rate, no test, no confidence set. Norm < 1 fixed |
| Wang, Kolar, Shojaie (JASA 2024), arXiv:2007.07448 | **F** (§2, Ass. 1–4, Thm 1–3) | Tests/CIs for Granger edges on one long path. Λ_max(ΩᵀΩ) ≤ γ² < 1 with γ fixed; intensity bounded between λ_min and λ_max; transition kernel shape known | **Not uniform as ρ → 1.** Bounded intensity is restrictive for self-exciting data. Parametric kernel shape |
| Achab et al., arXiv:1607.06333 | A | Integrated-cumulant estimate of the kernel-norm matrix ("causality" = Granger influence) | No inference |
| Christgau, Petersen, Hansen (AoS 2023), arXiv:2203.13559; Thams & Hansen arXiv:2110.12709 | A (from S8 + search) | DML-style local-independence test; i.i.d. replicates of the process | One long path with memory: not covered |
| Cordi, Challet, Muni Toke (J. Stat. Mech. 2018), arXiv:1709.08516 | A | Univariate and symmetric multivariate Hawkes are only weakly "causal" under time reversal | Contrary note: directional claims from Hawkes fits can be weak |
| Cai et al. THP, arXiv:2105.10884; TNPAR arXiv:2306.14114; Jin & Huang arXiv:2508.11727 (latent confounders); ICLR 2026 cross-covariance identification (id not captured) | A / S | Causal-structure learning on event sequences | No error control; no memory or criticality analysis found |

**5b. Interventions and potential outcomes for point processes.**

| Work | Access | What it gives | Memory |
|---|---|---|---|
| Noorbakhsh & Gomez-Rodriguez (NeurIPS 2022), arXiv:2111.07603 | A (from S8) | Counterfactual sequences via Gumbel-Max thinning SCM | Model-based, cross-world |
| Leclère, Ouazzani Chahdi, Rosenbaum, Szymanski, arXiv:2607.03239 (Jul 2026) | **F** (§1.1, §2.3–2.4, Ass. A–B) | Market impact = price functional on factual minus counterfactual path; exact conditional simulation of the latent Poisson measure | Counterfactual is **not identified from data** (they say so); it is fixed by Assumption B, a shared-noise coupling, plus a fitted model. Inserted metaorder "assumed not to excite the background Hawkes flow" in the aggressive case |
| Jiang, Chen, Ding, arXiv:2301.03246 | **F** (§2, Ass. 1–5) | IV with point-process treatment and outcome; Wald after Fourier transform = deconvolution | i.i.d. units. Deconvolution of a slowly decaying (power-law) response is the ill-posed case; not studied there |
| Ryalen, Stensrud, Røysland, arXiv:2604.12977 (Apr 2026) | S | Potential outcomes and g-formulas for marked point processes | Asymptotic regime not stated on page |
| Kresin, Clark, Davis, Hazelton, arXiv:2604.12124 (Apr 2026) | S | Spatiotemporal PP causal effects with spillover **and carryover**, Hawkes allowed; nonasymptotic EM error bounds | Closest 2026 paper to "intervention on a self-exciting process". Stationarity / branching assumptions not on the abstract page |
| Rytgaard et al. (AoS 2022); Røysland et al. (JRSSB 2024) | A (from S8) | Continuous-time TMLE; graphical identification of intensity interventions | i.i.d. subjects |

**5c. Near-criticality: is it even there?** (contested premise)

| Work | Access | Claim |
|---|---|---|
| Jaisson & Rosenbaum (AAP 2015), arXiv:1310.2033; (heavy-tailed) arXiv:1504.03100 | A | Nearly unstable Hawkes → integrated CIR; heavy-tailed nearly unstable → rough fractional diffusions |
| Hardiman, Bercot, Bouchaud (EPJB 2013), arXiv:1302.1405 | A | Markets near criticality (branching ratio ≈ 1) with power-law kernels |
| Filimonov & Sornette, arXiv:1308.6756 | S (search summary) | Branching-ratio estimates biased upward by outliers, power-law kernels, edge effects; regime-switching Poisson can look critical |
| Wheatley, Wehrli, Sornette (QF 2019) doi:10.1080/14697688.2018.1550266 | A | With a flexible background (EM + BIC) criticality is **rejected**; spurious long memory from trends, shocks, shifts |
| Muhle-Karbe, Ouazzani Chahdi, Rosenbaum, Szymanski, arXiv:2601.23172 (Jan 2026) | **F** (§1, §2.1–2.3) | Core flow (autonomous, metaorder splitting) and reaction flow (endogenous), both Hawkes, in a nearly unstable heavy-tailed regime; one parameter H₀ ≈ 3/4 pins signed-flow persistence, rough volume (H₀ − ½), rough volatility (2H₀ − 3/2) and impact exponent (2 − 2H₀). They note core and reaction trades cannot be told apart in anonymous data |

## 6. Market impact as a causal estimand

- **Structural / model-based:** Leclère et al. 2026 (F) and Muhle-Karbe et al. 2026 (F). Both
  define impact through a model; neither identifies it from data without the model.
- **Identification difficulty with public data:** Naviglio, Bormetti, Campigli, Rodikov, Lillo,
  arXiv:2501.17096 (S): models on public data give linear, non-reverting average paths because
  they misattribute the source of order-flow autocorrelation (memory). This is, in effect, a
  confounding-by-memory argument, though they do not use causal language.
- **Mechanical-impact evidence by scrambling:** Maitrier, Loeper, Kanazawa, Bouchaud,
  arXiv:2502.16246 (A): synthetic metaorders built by scrambling trader IDs show the same
  square-root law. This is a placebo-style design, not randomisation.
- **Memory has an observed cause:** Sato & Kanazawa (PRL 2023, arXiv:2301.13505; PRR 2023,
  arXiv:2308.01112) (A): with nine years of TSE account-level data, metaorder-length power laws
  quantitatively predict the order-sign memory exponent (γ ≈ α − 1). Caveat from Kanazawa's own
  slides: "virtual server" IDs are not exactly trader IDs.
- **Only design-based trading evidence:** Bojinov & Shephard 2019 (futures execution
  experiments) — exact randomisation tests, no memory assumption.

## 7. Where the three themes meet (memory × kernels × causal) — for the coordinator

1. **Track J vs long memory.** Track J's kernel-proxy claims assume β-mixing. Hao (2026) claims
   proximal identification *under long-range dependence* with an fBm limit (non-kernel, S,
   unverified). A kernel proxy estimator under long memory is a direct triple intersection.
2. **Local-independence testing (Christgau et al.)** is a DML test with ML/kernel-type nuisances
   for counting processes, valid for i.i.d. replicates. One long self-exciting path with
   power-law kernels is the open case — memory + nonparametric nuisance + causal (Granger) test.
3. **Deconvolution IV for point processes (Jiang–Chen–Ding)** works in the Fourier domain; a
   power-law (long-memory) response kernel makes the inverse problem ill-posed. Spectral /
   RKHS regularisation would be the kernel route.
4. **Word clash.** "Kernel" means the Hawkes excitation kernel (φ) in §5 and an RKHS kernel in
   the memory-kernel workstream. EDD's "nonparametric link functions" are the former. Keep them
   apart in the synthesis.

## 8. A derivation that organises the field (D — mine, elementary, not novelty)

For a binary treatment at time t with sequential ignorability given the full past, the AIPW
score splits into two parts:

ψ_t − θ = [τ(X_t) − θ] + [IPW residual terms].

- The residual part is a martingale difference with respect to the past (at horizon 0). It does
  not carry memory, whatever the memory of X_t.
- The first part, τ(X_t) − θ, is a function of the state. If the effect is heterogeneous and
  X_t has long memory (for instance τ depends on volatility, d ≈ 0.4), this part inherits it.
  With Hermite rank 1, its sample mean converges at T^{d−1/2}, not T^{-1/2}, and HAC intervals
  under-cover.
- Consequences: (a) the **population ATE** under long-memory heterogeneity is estimable only at a
  slower rate; (b) the **path-conditional ATE** (average of τ over the realised states) drops the
  first part and keeps √T inference; (c) design-based estimands (Bojinov–Shephard) survive for the
  same reason; (d) Ballinari–Wehrli's Ass. 3 is on the score, so it can hold under long-memory
  data when effects are homogeneous; it fails when effects load on the long-memory state.

This split is folklore in i.i.d. settings (sample vs population ATE). I did not find it stated
for long-memory time series in this scan. A bounded scan cannot show it is absent.
