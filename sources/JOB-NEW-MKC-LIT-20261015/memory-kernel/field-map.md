# Field map — kernels with memory (workstream `memory-kernel`)

Job JOB-NEW-MKC-LIT-20261015 · execution `claude-mkc-kernel-20261015` · 2026-10-15.
Plan SHA-256 6f6631ad4cb13f3abc3238c513498db95abe992ebff700fba1aeea8e938c99e3 (verified before start).

Bounded scan. No novelty claim. A gap below means "not found in this scan", never "does not exist".

**Evidence labels.** F = I read the relevant full-text section myself (PDF text extracted locally).
S = summarising fetch (WebFetch output, incl. arXiv abstract pages). A = abstract or search-snippet
level only. R = taken from an existing programme review (Track E strands E1/E2, S2, S8, Track A jobs),
not re-read here. Keys in `[brackets]` are in `bibliography.bib`.

**Two meanings of "kernel".** RKHS = positive-definite kernels (HSIC, MMD, embeddings, signature
kernels). MEM = memory kernels (Hawkes excitation, propagator/impact, fractional, Langevin memory).

---

## 0. The organising fact: three dependence tiers, and a boundary at H = 3/4

Almost every kernel guarantee in this field is stated for one of three tiers. Short-memory results
are the baseline; they break in a predictable order.

| Tier | Typical condition | What survives | Source (class) |
|---|---|---|---|
| 1. Geometric / fast polynomial mixing | β(k) or τ(k) ≤ Ck^{−r}, r large; or geometric ergodicity | Everything: degenerate-V-statistic limits, wild bootstrap, kernel CI tests, conformal coverage rates | [chwialkowski2014wild] F; E1 R |
| 2. Slow polynomial mixing / summable covariances | Σ_k ρ(k) < ∞ but r small (e.g. power-law Hawkes counts with moderate tail) | Consistency and √n limits often survive; **published bootstrap theorems do not** — they need τ(r) = o(r^{−4}) (bootstrap V-statistic) and τ(r) = O(r^{−6−ε}) (HSIC test) | [chwialkowski2014wild] F |
| 3. Long-range dependence (LRD) | ρ(k) ~ k^{−D}, D ∈ (0,1), i.e. Hurst H = 1 − D/2 ∈ (1/2, 1) | Consistency of embeddings survives at slower rate. Distributional limits change: for quadratic (kernel-type) statistics the limit stays √n-Gaussian only while H < 3/4; above that, non-central (Rosenblatt-type) limits with unknown laws | [betken2021ldr] F; [zhai2026ggmlong] F; [alquier2026mmd] F |

Why H = 3/4 matters. Kernel tests are quadratic functionals of the data. For a quadratic statistic of a
Gaussian-subordinated LRD process, the cross-product's autocovariance is ~k^{−2D}, summable iff D > 1/2,
i.e. H < 3/4. Two independent 2025–26 results land on the same boundary:

- Betken & Dehling (Bernoulli, arXiv v3 Jan 2026): distance-covariance (≡ HSIC with a distance kernel)
  cross-covariance is √n-Gaussian when D_ξ, D_η ∈ (1/2, 1); when D_ξ + D_η < 1 the rate and limit
  change and only Gaussian (not subordinated) processes are covered; D_ξ + D_η > 1 with one D < 1/2 is
  not covered at all (Thm 2.2, Remark 2.1). [betken2021ldr] F.
- Zhai, Zhong & Wu (arXiv Sep 2026): Gaussian-graphical CI test for high-dimensional linear processes,
  Berry–Esseen + block-bootstrap validity for coefficient decay β > 3/4 (≈ d < 1/4, H < 3/4); "ultra-long
  memory" 1/2 < β ≤ 3/4 is left open and their simulation at β = 0.55 fails "expectedly".
  [zhai2026ggmlong] F.

Where finance sits relative to the boundary (contested — see gap-map C1):
- Order signs: autocorrelation exponent ≈ 0.6, H ≈ 0.7 (LSE) — just inside H < 3/4. [lillo2004longmemory] S.
- Realised volatility: long-memory parameter d ≈ 0.4 is the standard ABDL finding, i.e. H ≈ 0.9 — outside.
  [andersen2003rv] (metadata verified via CrossRef; d-value is background knowledge, not re-read here).
- Impact propagators: power-law decay exponent ≈ 0.22 (Bucci et al.), so propagator-filtered flows are
  strongly persistent. Track A P2-CALIB R.
- Power-law Hawkes counts: count autocovariance ~ lag^{p−1} (polynomial; tier 2 or 3 depending on tail).
  [lee2026powerlawmarkov] S; the statement about count covariance is from that paper's abstract/snippet.

---

## 1. Lines of work

Each line: what it does, key results, the dependence it needs, and whether it survives long memory.

### L1. Kernel marginal, two-sample and serial tests under dependence (RKHS)

- **Short-memory field is closed** (reuse E1 R): HSIC/dCov/MMD under β/τ/α-mixing and NED
  [chwialkowski2014random], [chwialkowski2014wild], Wang–Li–Zhu 2018, Diz-Castro et al. JMVA 2026,
  Ghoshal & Shao EJS 2026 (self-normalised RKHS, α-mixing), Li–Xu–Zhou 2026 AutoHSIC. Impossibility under
  mere ergodicity: Yao 2026 (E1 [37]).
- **Exact mixing price of the wild bootstrap** [chwialkowski2014wild] F: Thm 1 (bootstrapped V-statistic)
  needs Σ r² τ(r) < ∞; Thm 2 needs τ(r) = o(r^{−4}); the HSIC/multi-lag test uses τ(r) = O(r^{−6−ε}).
  Gaussian LRD (τ ~ r^{−D}, D < 1) violates all three. Slow-polynomial processes (tier 2) also violate them.
- **LRD results that exist:** distance covariance with subsampling [betken2021ldr] F (above);
  two-sample U-statistic processes under LRD via Hermite expansion, with Dm < 1 non-central limits
  [dehling2014twosample] F (the MMD analogue, written for change-point tests, not for MMD);
  subsampling validity for general statistics under LRD [betken2018subsampling] A.
- **Self-normalisation under LRD** is established for CUSUM/Wilcoxon change-point statistics (Shao 2011;
  Betken 2016 [betken2016sn] A) but **not found for any kernel/HSIC/MMD statistic** in this scan.
- Survives long memory? **Consistency yes** (Betken–Dehling Thm 2.6 needs only stationary ergodic).
  **Calibration: only subsampling, only for (subordinated) Gaussian processes, and only cleanly below H = 3/4.**

### L2. Kernel conditional-independence and Granger tests under dependence (RKHS)

- Reuse E1/E2/GAP_ANALYSIS R: single-path RKHS conditional-mean Granger test (Tian, Waterbury & Düker
  2026, geometric ergodicity, homoskedastic Gaussian innovations); GTCM (Shi & Shah 2026, block-CI
  decay); distributional Granger by OT (Wang JTSA 2026, "weak dependence"); dGCM nonstationary
  (Wieck-Sosa et al. 2025); neural DR Granger on price–volume (Hui et al. 2025). All tier 1.
- **Long memory:** no kernel CI test found. Closest is linear/Gaussian: Zhai–Zhong–Wu
  [zhai2026ggmlong] F (precision-matrix zero pattern, β > 3/4). E1 already flagged "no kernel CI under LRD"
  as Track E gap G5 (ranked low).
- Survives long memory? **No result either way.** The H = 3/4 boundary above suggests where a residual-
  embedding statistic would stop being √n-Gaussian, but nobody has shown it (worker inference).

### L3. Kernel embeddings, MMD estimation and kernel regression with dependent data (RKHS)

- **MMD minimum-distance estimation** [alquier2026mmd] F, building on [cherief2022mmd] S: the key bound is
  E[MMD(P̂_T, P_0)] ≤ √((1 + 2Σ_T)/T), with Σ_T = Σ_{s≤T} ϱ_s and ϱ_s a kernel-specific dependence
  coefficient bounded by β(s) for radial kernels (Prop 4.1). Consistency needs only **Σ_T = o(T)**.
  Asymptotic normality needs β-mixing with Σ β(t)^{δ/(2+δ)} < ∞. So: consistency **survives LRD** with a
  slower rate (if ϱ_s ~ s^{−D}, rate T^{−D/2} — worker inference); inference does **not** (tier 1 needed).
- **Kernel regression rates under dependence:** KRR with random features under τ-mixing — exponential
  decay gives minimax rates, polynomial decay gives suboptimal ones [wang2024krrdependent] A; kernel PLS
  under LRD — rates slow down with memory [singer2017kpls] A. (Matsuoka & Torrent 2026,
  arXiv:2603.04681, S, is about smoothing-kernel averages under strong mixing — not RKHS; excluded.)
- **Koopman / transfer operators in RKHS** (operator-valued "kernels" for dynamics): [klus2020kto] A,
  [kostic2022koopman] S — bounds for trajectory data "in terms of mixing coefficients", and the object is
  the Koopman operator of a **Markov** chain. Non-Markovian memory is outside the framework; no 2025–26
  non-Markovian extension found (A-level search only).
- Survives long memory? **Point estimation/consistency: yes, slower. Rates and confidence statements: no
  (not shown).**

### L4. Signature and path kernels (RKHS on path space; natural home for rough/fractional processes)

- Foundations: signature moments characterise laws of processes and give sig-MMD
  [chevyrev2022sigmoments] A; the signature kernel solves a Goursat PDE [salvi2021goursat] A; survey
  [cass2024sigkernel] A. fBm has a rough-path lift for H > 1/4, so signature kernels are defined for rough
  volatility paths [concentration2025sig] A.
- Tests: sig-MMD two-sample tests on **collections of i.i.d. paths** [alden2025sigmmd] A (E1 excluded list);
  signature-kernel CI test + causal discovery for SDEs [manten2025sigci] F. Manten et al. (ICLR 2025):
  * the CI test needs **n i.i.d. sample paths** (App. A.13); within-path dependence is allowed;
  * power on **fBm-driven** systems approaches 1 for all H tried (App., Fig. 8) — so the test itself
    is not broken by fractional noise when replicates exist;
  * but **the causal-discovery algorithm is "not applicable anymore"** once the driving noise lacks
    independent increments (fBm) — the Markov/past–future split it relies on fails. ⚠ TRIPLE
    INTERSECTION (memory × kernel × causal): see gap K7.
- Hawkes via signatures: [noguer2026paths] F — expected truncated signatures of **exponential/affine**
  Hawkes have finite linear closures; parameters identifiable; an antisymmetric level-2 cross-area
  detects the **direction** of cross-excitation to first order (validated on 50,000 independent
  simulations). Power-law kernels and single-path estimation are not treated. Single author, unrefereed.
  ⚠ TRIPLE-adjacent (direction of excitation ≈ Granger direction).
- Other 2025–26 signature work found is pricing/computation (e.g. [bayer2025sigamerican] A,
  [buresdesantiago2026sigvol] S: classifies volatility models incl. close Hurst values with signature
  features + boosting, no guarantees).
- Survives long memory? **Defined and usable for rough paths. Every inferential guarantee found assumes
  i.i.d. replicate paths.** Finance has one path. No single-path signature-kernel test was found.

### L5. Nonparametric / RKHS estimation of memory kernels (MEM, sometimes via RKHS)

Hawkes excitation kernels:
- Second-order (Wiener–Hopf) inversion — unique causal solution, nonparametric, works for **power-law**
  kernels; error discussed but no consistency theorem in abstract [bacry2016wienerhopf] S;
  earlier symmetric-kernel version on HF data [bacry2012symmetric] A.
- **RKHS Hawkes estimators assume bounded support [0, A]:**
  * Bonnet & Sangnier (AISTATS 2025) [bonnet2025rkhs] F: g = (h + b)·1_[0,A], h ∈ RKHS, ReLU link;
    representer theorems for approximated likelihood/least squares; a bound on the approximation error
    only — **no statistical rate**; experiments use A = 5.
  * Kim & Iwata (ICLR 2026) [kim2026representer] F: representer theorem for exact penalized least
    squares, dual coefficients fixed at one; finite support window A; synthetic data only; **no rate**.
  A power-law tail is cut off at A by construction. The truncation bias near criticality is not studied.
- Near-critical heavy-tailed Hawkes → rough volatility: Jaisson & Rosenbaum [jaisson2016rough] S
  (tail exponent α ∈ (1/2,1) gives fractional CIR with H = α − 1/2); Horst, Xu & Zhang [horst2023heavy] A
  (stronger convergence, rough Heston); supercritical case [xu2025supercritical] A. These are scaling
  limits, **not estimators**. Two further leads seen only in search snippets and **not verified**
  (so not cited): a report that Hawkes MLE is biased near criticality, and a 2025 working paper on
  branching-ratio inference with a near-criticality test. Both are logged in search-log.md for the
  memory-causal worker.
- A Markov surrogate with power-law intensity and long memory near the stability edge
  [lee2026powerlawmarkov] S — positive Harris recurrence, so mixing tools become available.

Propagator / impact kernels:
- Neuman & Zhang 2023 [neuman2023regret] S (online, regret) and Neuman, Stockinger & Zhang 2023
  [neuman2023offline] A: Tikhonov-regularised nonparametric propagator estimation.
- **Hey, Neuman & Tuschmann (Oct 2025)** [hey2025crossimpact] F: nonparametric concave multi-asset
  propagator; constrained regularised least squares on an M-point lag grid; high-probability confidence
  region (Thm 2.10) from a self-normalised bound. Assumption 2.3: noise is **conditionally sub-Gaussian
  given the past (martingale-difference type)**. No mixing condition on the order flow is needed.
  Empirics: CFM metaorders + public S&P flow; self-impact concave with **shifted power-law decay**;
  asymmetric cross-impact. So long-memory order flow does not by itself void the confidence region, but
  the lag grid truncates a power-law tail and noise autocorrelation would void Assumption 2.3
  (worker reading).
- Track A P2-CALIB R: propagator exponent ≈ 0.22 (Bucci et al.); asymptote poorly identified.

Langevin memory kernels (physics sense):
- Lang & Lu [lang2024gle] A: RKHS-regularised recovery of the GLE memory kernel (first-kind Volterra,
  ill-posed); error bounded by correlation-function error. Same inverse-problem structure as the Hawkes
  Wiener–Hopf and propagator deconvolution. No finance use found.

Survives long memory? **Estimation of power-law kernels exists (Wiener–Hopf, propagators) but has no
rate theory; RKHS methods with theory truncate the tail.**

### L6. Kernel-based conformal prediction under dependence (RKHS weights / kernel scores)

- KOWCPI (ICLR 2025) [lee2025kowcpi] F: kernel-reweighted Nadaraya–Watson quantile scores; conditional
  coverage under **α-mixing of residuals with α(τ) = O(τ^{−(2+δ)})** (Assumption 4.2) — tier 1/2 only.
- MDCP kernel-PIT (β-mixing) and Zhai–Cheng–Wu coverage under FDM incl. long-memory Gaussian linear
  processes (residual/CQR scores, not kernel) — reuse S2 R (#10–12).
- Survives long memory? **Kernel-score methods: not shown. Non-kernel scores: yes (Zhai–Cheng–Wu).**
  Track F-KERNEL owns this.

### L7. Kernel causal estimators (proxy, IV, effects) and their dependence assumptions

- All i.i.d.: kernel proxy [mastouri2021proxy] A, KIV minimax [meunier2024kiv] A, KIV uniform inference
  [lob2025kivuniform] S (dependence not stated), counterfactual embeddings (E2 R). Track J's CL-J1 is
  the β-mixing extension.
- Long memory: nothing found. Plausible route via L3 (Σ_T = o(T) consistency) for point estimation; no
  route for inference without new limit theory (worker inference).

### L8. Fading-memory kernels (system identification)

- [huo2025fading] F: RKHS of fading-memory functionals with a weighted norm w; deterministic
  input–output identification; "fading memory" here is a continuity property, not a stochastic
  dependence class. [grigoryeva2025reservoir] A similar (reservoir/Volterra kernels). Low relevance;
  listed so the term is not confused with long memory.

---

## 2. Summary: does each line survive long memory?

| Line | Meaning of kernel | Tier-1 status | Survives LRD? | Best evidence |
|---|---|---|---|---|
| L1 marginal / two-sample tests | RKHS | closed | consistency yes; calibration only by subsampling, subordinated Gaussian, H < 3/4 clean | [betken2021ldr] F, [chwialkowski2014wild] F |
| L2 kernel CI / Granger | RKHS | crowded (E1) | not shown (linear Gaussian analogue only, H < 3/4) | [zhai2026ggmlong] F |
| L3 embeddings / MMD estimation / KRR / Koopman | RKHS | active | consistency yes (Σ_T = o(T)); rates slower; inference no; Koopman needs Markov | [alquier2026mmd] F |
| L4 signature kernels | RKHS on paths | active, i.i.d. paths | defined for rough paths; tests need replicates; discovery breaks with fBm noise | [manten2025sigci] F |
| L5 memory-kernel estimation (Hawkes, propagator, GLE) | MEM (some via RKHS) | active | power-law estimation without rates; RKHS theory truncates tails | [bonnet2025rkhs] F, [kim2026representer] F, [hey2025crossimpact] F |
| L6 kernel conformal | RKHS | active (F-KERNEL) | kernel scores not shown | [lee2025kowcpi] F |
| L7 kernel causal estimators | RKHS | i.i.d. only | nothing | E2 R; [lob2025kivuniform] S |
| L8 fading-memory kernels | RKHS of functionals | n/a | different notion | [huo2025fading] F |

## 3. Where RKHS kernels meet memory kernels (the author's "especially")

1. **Estimating a memory kernel inside an RKHS** — Hawkes (Bonnet–Sangnier, Kim–Iwata), GLE (Lang–Lu),
   propagators (Tikhonov ≈ finite-dimensional RKHS ridge; Neuman et al., Hey et al.). Common structure:
   a first-kind Volterra / Wiener–Hopf inverse problem regularised by an RKHS norm. Common gap: heavy
   tails and near-critical behaviour are where the inverse problem is worst conditioned, and no rate theory
   covers them.
2. **Kernel tests whose validity depends on the memory kernel** — a power-law Hawkes or propagator kernel
   sets the mixing tier of the observed series, which decides whether wild-bootstrap HSIC/MMD is valid.
   The tail exponent of the memory kernel therefore maps directly onto the r^{−6} / H = 3/4 thresholds
   above. Nobody has written this map down (worker inference; not found in scan).
3. **Signatures of self-exciting paths** — expected signatures identify exponential Hawkes parameters and
   excitation direction [noguer2026paths] F. Extending to power-law kernels would put a memory kernel,
   an RKHS kernel and a direction-of-influence (causal) question in one object.

## 4. Limits

- WebSearch returned few 2025–26 records for long-memory queries; arXiv listing search is AND-of-terms.
  The 2025–26 update is therefore thinnest for SSRN/econometrics journals.
- Full text read for 12 papers, sections only (assumptions, main theorems). Proofs not checked.
- Several useful statements are worker inference (marked). They are leads, not results.
