# Memory-removing kernels: what is already known — 1 November 2026

This note answers Q1–Q6 from the team lead's brief. It supports the open question at `paper/long-memory-kernel-tests/paper.tex:436`. BibTeX is in `literature/memory-kernel-refs-2026-11-01.bib`.

## Bottom line

1. **The mechanism is known.** Rank and filtering decide the long-memory boundary for scalar and linear statistics. Betken and Dehling already print a rank-weighted line for Pearson. Bardet and Surgailis print the differencing result for a rank-two scalar statistic.
2. **The kernel forms are known.** The increment kernel is printed verbatim in signature, preference and slow-feature papers. The Gram form `Q K Q` is the residual form of the KCI test and the contrast form of conditionally positive definite kernels.
3. **The application was not found.** I found no kernel built to weaken long memory. I found no single-path independence test on feature increments. I found no size or power study of either under long memory.
4. **"Holding size" is not open in principle.** Harris (2020) proves a finite-sample size bound for any statistic and any stationary series. It needs no mixing. The memory-removing kernel must therefore win on power or simplicity.
5. **A Hilbert-valued Breuer-Major theorem exists.** Düker and Zoubouloglou, Theorem 3.1, covers the Hermite-deflated kernel. It does not cover the increment kernel as stated.

## 0. Scope, method and labels

**Scope.** I did not repeat what `LITERATURE-AUDIT.md`, `assumption-audit-2026-10-23.md` and `litreview.tex` settle. I cite their findings M1, M2 and M3 by name.

**Evidence labels.**
- **Full:** I read the stated sections of the text with `pdftotext`.
- **Abstract:** I read the abstract or a journal record that holds it.
- **Excerpt:** I saw a passage quoted in a search result, and did not open the source.
- **Record:** I confirmed the title and venue on a landing page or in Crossref or Semantic Scholar. I read no content. I use these only through what an opened paper says about them.
- **Derived here:** my own working. No source states it. A referee should check it.

**Tools.** WebSearch, the arXiv export API, Crossref, Semantic Scholar, `curl` on open PDFs, and `pdftotext`.

**Failures.** OpenAlex returned HTTP 429 for the whole session. Project Euclid PDFs returned an error page, so I could not read Arcones (1994) or Dehling and Taqqu (1989). ScienceDirect, Wiley, Oxford and Taylor & Francis pages were closed or returned 403. SSRN and Google Scholar were not searched. A missed paper remains possible, mainly in journals without arXiv versions.

**Notation.** `D = 2 - 2H`, so `rho(k) ~ k^{-D}`. `H_X + H_Y = 3/2` equals `D_X + D_Y = 1`. `k_D` is the increment kernel from the brief. `Delta phi_t = phi(x_{t+1}) - phi(x_t)`.

## Verdicts

| Q | Topic | Verdict |
|---|---|---|
| Q1 | Raising the Hermite rank of a kernel to weaken long memory | **Partly known.** Known for scalar and linear statistics, and for kernels with other purposes. Not found as a kernel test. |
| Q2 | Increment kernel or feature differencing in testing | **Partly known.** The form is known. Raw differencing before a test is known. A single-path feature-increment test was not found. |
| Q3 | Differencing as the cure for spurious regression | **Known** at the linear level, including its costs. Kernel analogue not found. |
| Q4 | Limit theorems that could carry a proof | **Known.** A Hilbert-valued central theorem and a Hilbert-valued reduction both exist. Mixed pairs and the increment kernel are not covered. |
| Q5 | Other remedies | **Known** under mixing. **Partly known** for shifts without mixing. Long-memory kernel remedies other than Betken–Dehling were not found. |
| Q6 | Impossibility or trade-off | **Partly known.** Close impossibility theorems exist for neighbouring problems. The layer-level trade-off is proved. A kernel-level power-cost theorem was not found. |

---

## Q1. Has anyone raised the Hermite rank on purpose?

**Verdict: partly known.** The effect of rank is proved for scalars and for Pearson. Rank two arises by construction in robust scale estimation. Kernels with the deflated form exist for other purposes. No kernel, feature map or U-statistic kernel was built to weaken long memory.

**Closest sources.**

- **Betken and Dehling, arXiv 2107.03041v3 (Full).** This is the closest statement of the mechanism.
  - Theorem 4.1(i): the Pearson cross-covariance of `G_1(xi)` and `G_2(eta)` has a root-n normal limit if `D_xi r_1 + D_eta r_2 > 1`, with `r_i` the Hermite rank of `G_i`.
  - That is the rank-weighted sum line, and it includes mixed pairs.
  - It concerns a linear statistic on transformed data. They do not choose `G` to raise the rank, and they build no kernel.
- **Lévy-Leduc, Boistard, Moulines, Taqqu, Reisen 2011, JTSA (Full, arXiv 0912.4686).** They study the Rousseeuw–Croux `Q_n` scale estimator.
  - Theorem 6: root-n normal if `D > 1/2`, and rate `n^D` with a second-chaos limit if `D < 1/2`.
  - They prove the rates match the classical autocovariance estimators. So the `H = 3/4` boundary appears by construction, because the estimator is built to measure scale.
  - The aim is robustness to outliers, not memory removal.
- **Ho and Hsing 1996, Ann. Statist. (Abstract).** For long-memory moving averages, they expand the empirical process in `p` terms and study the remainder. Each term is a derivative of the distribution function times an uncorrelated sum `Y_{N,r}`. This is the closest analogue of the population form `phi(x) - mu - a_1 T(x)`. It is a limit-theory device, not a test.
- **Lévy-Leduc et al. 2011, Ann. Statist. (Full, arXiv 0912.4688).** They state that under Gaussian long memory the Hoeffding decomposition "may not provide the complete asymptotic behavior of U-statistics because all terms of this decomposition may contribute to the limit". Their eq. (52) expands the degenerate part with `c_{0,p} = c_{p,0} = 0`.
  - So Hoeffding degeneracy removes the one-variable layers only. The `(1,1)` layer stays.
  - Degeneracy alone therefore does not remove memory. This matters for MMD.
- **Dehling, Rooch, Wendler 2017 (Full, arXiv 1404.0551).** Definition 2.2 sets rank `m` as the smallest `k + l` with a nonzero coefficient. Theorem 2.1 gives the limit for any rank `m` with `Dm < 1`, under coefficient condition (15). See the correction in the "Corrections" section.
- **Zhang, Peters, Janzing, Schölkopf 2011, KCI test (Full, arXiv 1202.3775).** Eqs. (10)–(11) give `R_Z = epsilon (K_Z + epsilon I)^{-1}` and the residual Gram matrix `R_Z K R_Z`.
  - *Derived here:* with a finite-rank polynomial kernel on `z`, `R_Z` tends to the projector `Q` as `epsilon` tends to zero. So `Q K Q` is a limit case of this form.
  - KCI conditions on a shared variable and assumes i.i.d. data. It does not target memory.
- **Davydov 2021 (Full, arXiv 2006.13543).** For a conditionally positive definite kernel, `M^T K_X M` is positive definite, where `M` spans the null space of the polynomial matrix (eq. 9 and following).
  - *Derived here:* `{1, He_1, ..., He_r}` spans the polynomials of degree at most `r`. So `Q K Q` in `z` has the same algebra as a conditionally positive definite kernel of order `r + 1`.
  - The theory there is numerical approximation. It has no probability and no long memory.
- **Plumlee and Joseph 2018, Stat. Sinica (Full on Lemma 1 and Theorem 1, arXiv 1611.00203).** Their eq. (5) gives `c*(x,x') = c(x,x') - h(x)^T H^{-1} h(x')`, with `h` and `H` built from integrals of `c` against known functions `g`. Lemma 1 says `c*` is positive semi-definite. Theorem 1 says the resulting Gaussian process is orthogonal to the mean functions.
  - *Derived here:* work in `u = F(x)`, so the integral is against the law of `X`, and take `g = Phi^{-1}(u)`. Then the feature of `c*` has zero first Hermite coefficient. This is a kernel of Hermite rank at least two.
  - It differs from the brief's population form. It removes an RKHS direction, so it also changes the higher layers. The brief's form leaves them unchanged.
  - Their theorem assumes a bounded domain and bounded `g`. This `g` is unbounded, so the case falls outside the stated conditions. The authors say those conditions are "sufficient, not necessary". The aim is identifiability in computer experiments.
- **Chwialkowski, Strathmann, Gretton 2016, Stein kernel test (Full on method).** The Stein operator gives zero-mean features under a known target. Dependence is handled by a wild bootstrap under weak-dependence conditions.
  - *Derived here:* for `T f = f' - x f` under `N(0,1)`, `T He_n = -He_{n+1}`. So a Stein feature raises the Hermite rank by exactly one.
  - It tests fit to a fixed model. It is not an independence test, and no source links it to long memory.
- **Bergsma 2006 (Full on Sec. 2, arXiv math/0604627).** He decomposes a kernel independence coefficient into a weighted sum of squared correlations between marginal eigenfunctions. For one special marginal, the eigenfunctions use Hermite polynomials of the normal score. The setting is i.i.d. A referee could point to this as a layer view of kernel dependence.
- **Bai and Taqqu 2018 and 2019 (Full).** These are warnings, not prior constructions. See Risks.

**Not found after these searches:** a kernel, feature map or U-statistic kernel chosen to raise Hermite rank in order to weaken long memory. A Stein kernel used under long memory. A kernel test whose Gram matrix removes low Hermite layers.

---

## Q2. Is the increment kernel or feature differencing already used?

**Verdict: partly known.** The formula appears verbatim in three literatures. Raw differencing before a test is standard practice. A single-path independence test on consecutive feature increments was not found.

**Verbatim forms.**

- **Király and Oberhauser 2019, JMLR (Full, arXiv 1601.08169).** Section 1 prints `k(x(s_{i+1}), x'(t_{j+1})) + k(x(s_i), x'(t_j)) - k(x(s_i), x'(t_{j+1})) - k(x(s_{i+1}), x'(t_j))`. Section 5 calls this "kernel trick two".
  - The data unit is a whole sequence. The form is summed over ordered index tuples to give one sequence-level kernel.
  - Our construction indexes the Gram matrix by time on one path. That is a different object.
  - Remark 5.10(i) reports numerical instability if the differences are dropped.
- **Houlsby, Hernández-Lobato, Huszár, Ghahramani 2012 (Full).** They print `k_pref((x_i,x_j),(x_k,x_l)) = k(x_i,x_k) + k(x_j,x_l) - k(x_i,x_l) - k(x_j,x_k)`. Pairs are item comparisons, not time steps.
- **Herbrich, Graepel, Obermayer 2000 (Record).** This is the large-margin ranking chapter. Houlsby et al. say the kernel also follows from a large-margin view. I did not read the chapter.
- **Böhmer, Grünewälder, Nickisch, Obermayer 2011 (Abstract, formula from an Excerpt).** Kernel slow feature analysis minimises `tr(A^T K D D^T K^T A)` with `D` a temporal difference matrix. This is feature-level differencing with the `D K` Gram form. The aim is extraction of slow features, which is the opposite goal. There is no test.

**Frameworks that contain it.**

- **Zhang, Song, Gretton, Smola 2008, NeurIPS (Full).** This is the closest structural prior art. Structured HSIC uses cliques `(x_t, ..., x_{t+tau}, y_t, ..., y_{t+tau})` with a free choice of joint kernel, and `K_ij = k(x_{i,tau}, x_{j,tau})`.
  - *Derived here:* `k_D` on `(x_t, x_{t+1})` is a valid clique kernel. So the increment HSIC is a special case of the framework.
  - They do not make that choice. Their Corollary 2 needs universal clique kernels, and `k_D` is not universal (*derived here*).
  - Their limit needs absolute regularity, which M1 excludes for fGn. They recover the lagged linear method TD-SEP as a special case.
- **Besserve, Logothetis, Schölkopf 2013, NeurIPS (Full).** The KCSD operator is a Fourier transform of feature cross-covariances. It tests each frequency by block shuffling and a max-t correction.
  - Theorem 5 needs `sum |k| ||C_12(k)|| < inf`, `sum |k| Tr C_ii(k) < inf` and a cumulant condition. Its results exclude frequency zero modulo 1/2.
  - *Derived here:* the trace condition fails for fGn features, because a rank-one layer gives `Tr C_ii(k)` of order `k^{-D}`.
  - *Derived here:* under stationarity, `E[Delta phi_t ⊗ Delta psi_t] = 2 C(0) - C(1) - C(-1)`. This equals the integral of `S_12(w)` against `|1 - e^{iw}|^2`. So the increment HSIC targets a high-pass weighted KCSD.
- **Hong 1999, JASA (Abstract and sections).** The generalized spectral density applies the same idea with empirical characteristic functions. It treats one series under mixing.

**Raw differencing and prewhitening before a kernel test.**

- **Flaxman, Neill, Smola 2016, TIST (Full on the method).** They run GP regression on time or space, and then HSIC on the residuals. They name "removing trends, and taking first differences" as standard prewhitening. They state that beyond linear Gaussian models "the conditions under which it holds are not known", and they give conditions in their Section 4.1. A text search of the PDF finds no mention of long memory.
- **Chwialkowski and Gretton 2014 (audit).** They differenced forex data before shift HSIC and checked decay.
- **Wang, Li, Zhu 2021 (litreview entry).** They run HSIC on the innovations of fitted parametric models.

**Scalar proof that differencing removes the boundary.**

- **Bardet and Surgailis 2011, Bernoulli (Full, arXiv 0802.0489).** In their Theorem 4.2, a rank-two function of first-order increments gets a central limit under `H < 3/4`. With second-order increments the limit holds for every `H < 1`. They say this is because "the second order increment process" of fBm is "a short memory stationary" process. They add that the same holds for generalized quadratic variations.
  - fGn is the first increment of fBm. So their second-order case is exactly differenced fGn.
  - The statistic is a function of raw increments, on Gaussian processes. It is not an increment of a nonlinear feature.

**Two different things are both called differencing.** Keep these apart in the paper.

| Variant | Statistic | Covered by |
|---|---|---|
| Raw-level | `k(x_{t+1} - x_t, ...)` | Bardet–Surgailis for Gaussian scalars. Theorem 3.1 of Düker–Zoubouloglou for Gaussian data (Q4). |
| Feature-level | `k_D`, i.e. `phi(x_{t+1}) - phi(x_t)` | Not covered by either theorem as stated (Q4). |

For non-Gaussian marginals the two differ. A raw increment of `g(xi)` is not a function of `Delta xi`.

**Filters beyond first differences.**

- **Clausel, Roueff, Taqqu, Tudor 2014 (Full on Theorem 3.1 and intro, arXiv 1105.1011).** This is a warning. For Hermite rank at least two, the wavelet scalogram has a Rosenblatt limit, under a stated condition on the vanishing moments `M`. This happens "even if the original process is defined through a Hermite polynomial of order higher than 2".
  - So a linear filter on a feature does not make statistics quadratic in that feature Gaussian.
  - It concerns one series with growing scale. HSIC is a cross product, so the effect should mainly hit its normaliser (*derived here, not proved*).
- **Podobnik and Stanley 2008 (Abstract) and Blythe 2013 (Abstract, arXiv 1309.4073).** Detrended cross-correlation applies a window detrending filter to a linear cross-covariance. Blythe gives a test with "conservative" p-values for zero correlation. This is the linear high-pass cousin of construction 3.
- **Müller and Watson 2018 and 2017 (Full on the framing).** Their cosine projection keeps a small number of low-frequency coordinates. DCT deflation removes those same coordinates. So they study the exact complement (*derived here*).
- **Hosking 1981 (Record) and Düker 2018 (Abstract).** Hosking is the origin of fractional differencing. Düker treats Hilbert-valued linear processes with an operator exponent. These are linear in the Hilbert space, so RKHS features of a Gaussian series are not directly covered.

**Signature tests.** **Chevyrev and Oberhauser 2022 (Abstract)** and **Manten et al. 2025 (Full on the setup)** test with signature kernels. Manten et al. write that "the n process samples are i.i.d. realizations". They include fBm paths, each path being one data point. Neither faces single-path long memory.

**Even kernels.** **Haasdonk and Burkhardt 2007 (Record)** is the kernel-invariance reference. The symmetrised kernel is an invariance under the reflection group (*derived here*). I found no use for long memory.

**Not found after these searches:** a single-path independence or two-sample test using `k_D` or `D K D'`. HSIC with a lag-L feature increment. Fractional differencing in an RKHS. DCT deflation of a Gram matrix in a test. A long-memory study of an even or symmetrised kernel test.

---

## Q3. What does the spurious-regression literature say?

**Verdict: known at the linear level, including the costs.** The kernel analogue was not found.

- **Granger and Newbold 1974 (Full).** "we recommend taking first differences of all variables that appear to be highly autocorrelated. Once more, this may not completely remove the problem". They also say "we are not advocating first differencing as a universal sure-fire solution". Their rule is "to build one's models both with levels and also with changes". Their reason is efficiency: differenced terms carry nearly uncorrelated information.
- **Phillips 1986 (Abstract).** For independent random walks, t ratios diverge and the Durbin–Watson statistic converges to zero.
- **Phillips 1998 (Excerpt of abstract).** Spurious regressions reproduce the coordinate system of an orthonormal representation of the process. This supports reading the problem as shared low-frequency structure.
- **Tsay and Chung 2000 (Abstract, first-hand).** If the orders of integration "sum up to a value greater than 0.5, the t ratios become divergent". The cause is "the long memory, instead of nonstationarity or lack of ergodicity". They also warn: "the usual first-differencing procedure may not completely eliminate spurious effects when data possess strong long memory".
  - *Derived here:* for stationary fGn, first differencing gives `d - 1 < 0`, so the sum is negative. Their warning bites when orders are large, not on our grid.
  - I did not read the full text.
- **Plosser and Schwert 1977 (Excerpt of abstract).** Excess differencing induces a noninvertible moving average in the errors. The search excerpt reports that the problem is not serious given care with the disturbances, and that the main cost is reduced efficiency. RePEc holds no abstract, so treat this as second-hand.
  - Differencing fGn is over-differencing in this sense (*derived here*).
  - The spectrum at zero behaves like `|w|^{3-2H}`, so it vanishes.
- **Engle and Granger 1987 (Abstract).** "A vector autoregression in differenced variables is incompatible with these representations." This is the direct statement that differencing loses cointegrating dependence.
- **Müller and Watson 2018 (Full on framing) and 2017 (Full on framing).** They study long-run covariability through about `q = 12` cosine-weighted averages. Their confidence sets cover I(0), I(1), near unit roots and fractional processes. They note "the scarcity of low-frequency information in the data".
  - This gives the clearest reading of the cost. A low-frequency common factor lives in a handful of coordinates.
  - A filter removes few coordinates but nearly all of that signal.

**What this means for us.** The classical literature treats differencing as a useful rule with known costs. It says to report both levels and changes. It says a cointegrating or common low-frequency link is lost by design. So we should present the increment test as answering a different question, and report it next to a levels test.

**Not found after these searches:** a power-cost formula for differencing in any kernel test. A statement of the `H_X + H_Y > 3/2` line after feature-level differencing.

---

## Q4. Which limit theorems could carry a proof?

**Verdict: known, with two named gaps.** Yes, a Hilbert-valued Breuer-Major theorem exists and is usable. Its status is a preprint.

### Theorem map

| Theorem | Regime | Conditions as printed | Gives | Use for us |
|---|---|---|---|---|
| Düker–Zoubouloglou arXiv 2405.11452v3, **Theorem 3.1 with Assumption 3.1** (Full) | Central | Stationary Gaussian `X_k` in separable `H_1`. `G` in `L^2(H_1, gamma_Q : H_2)` of Hermite rank `q >= 1` (Definition 4.1). `sum_v ||R(v)||_op^q < inf`. | Gaussian limit of `n^{-1/2} sum G[X_k]` in `H_2`. Covariance series converges in trace norm. | RKHS-valued sum of a deflated feature. |
| Same paper, **Theorem 3.2**, **Theorem 6.1** | Central | Same | `H_2`-valued Brownian partial-sum limit. A quantitative bound. | Split-sample and block statistics. |
| Düker–Zoubouloglou arXiv 2509.20342v3, **Theorem 3.2, 3.4** (Abstract and intro) | Central | Trace-norm covariance convergence and fourth moments, chaos by chaos | Fourth-moment theorem on Hilbert spaces | The tool for a direct proof for `k_D`. |
| Betken–Dehling arXiv 2107.03041v3, **Theorem 2.7** (Full) | Non-central | Gaussian `X`, `f: R -> S = L^2(R, w)`. Rank `r` of the class `f_t`. `0 < D < 1/r`. `E||f_X||^2 < inf`. | Hilbert-valued reduction to one direction. Limit `J_r(t)/r! Z_{r,H}`. | Non-central half for a deflated feature. |
| Betken–Dehling **Theorem 2.2** (Full) | Product | (i) both `D > 1/2`, subordinated. (ii) `D_xi + D_eta < 1`, `G = id` only. Mixed case not covered (Remark 2.1). | Joint characteristic-function limit | Template for HSIC. |
| Betken–Dehling **Theorem 4.1** (Full) | Product, Pearson | (i) `D_xi r_1 + D_eta r_2 > 1`. (ii) `G = id`, `D_xi + D_eta < 1`. | Pearson limit, mixed pairs included | Scalar benchmark for the sharp line. |
| Betken–Dehling **Theorem B.3** (Full) | Central, vector | Follows from Arcones Theorem 4 | CLT for a vector function | How part (i) is proved. |
| Breuer–Major 1983 **Theorem 1** (Record, number from Lévy-Leduc) | Central, scalar | `sum |rho|^k < inf`, rank `k` | Root-n normal | Classical base case. |
| Arcones 1994 **Theorem 4** (Record, number from Lévy-Leduc and Betken–Dehling) | Central, vector | Sum of the sup over entries, to the power `q` | Root-n normal | Vector base case. |
| Dobrushin–Major 1979, Taqqu 1979 (Record, theorem number not checked) | Non-central | Rank `k`, `kD < 1` | Hermite process limit | Classical base case. |
| Dehling–Taqqu 1989 **Theorem 1.1** (Record, number from Lévy-Leduc) | Non-central | Long-memory Gaussian subordination | Empirical process reduction | Handles normal scores and rank statistics. |
| Dehling–Rooch–Wendler 2017 **Theorem 2.1, 2.2** (Full) | Non-central, split | `Dm < 1`, rank `m`, condition (15) | Limit of the split-sample U-process | Cross term of split MMD. |
| Nourdin–Nualart 2020 **Theorem 1.1** (Full on the statement) | Central, scalar | `sum |rho|^d < inf`, `E|phi|^p < inf` for some `p > 2` | Tightness in Skorokhod space | Functional scalar limit only. |
| Major 2024 arXiv 1901.04086 (Abstract, theorem number not checked) | Non-central, vector | Vector-valued fields | Vector non-central limit | He calls Arcones's discussion "incomplete". Cite for vector non-central steps. |
| Düker 2018 (Abstract) | Linear | `X_k = sum (j+1)^{-N} eps_{k-j}`, operator `N` | Hilbert-valued CLT and non-CLT | Not for subordinated features. |
| Betken–Wendler 2018 **Theorem 1** with Assumptions 1–3 (Full on the statement) | Subsampling | A limit exists. `X = G(xi)`. Block rule. | Subsampling validity | Calibration once a limit is proved. |
| Betken–Dehling **Theorem 2.4** (Full) | Subsampling | Theirs for distance covariance | Subsampling validity | Template. |
| Bai–Taqqu 2019 **Theorem 3.1, 3.3, 3.4, Proposition 4.1** (Full) | Perturbation | See Q6 and Risks | Rank instability and power rank | Controls estimated centring. |

### What this gives each construction

**Construction 1, population form.** *Derived here.* Take `H_1 = R` and `H_2` the RKHS. A bounded continuous kernel makes `G` square-integrable and `H_2` separable. Then `R(v) = rho(v)` and `G = phi - mu - a_1 T` has rank `q = 2`.
- For `D > 1/2`, Theorem 3.1 gives a Gaussian RKHS limit.
- For `D < 1/2`, Theorem 2.7 gives a one-direction second-chaos limit, after a Bochner identification of the RKHS with `L^2(w)`.
- Removing `r` layers replaces 2 by `r + 1`. The boundary moves to `H = 1 - 1/(2(r+1))`.

**HSIC product.** *Derived here.* Take `H_1 = R^2` with `R(v) = diag(rho_xi(v), rho_eta(v))`. The feature product has rank `q = r_X + r_Y`. The operator norm is the larger of the two correlations.

| Layers removed on each side | `q` | Equal-memory boundary from Theorem 3.1 |
|---|---|---|
| none | 2 | `H < 3/4` |
| first layer on one side | 3 | `H < 5/6` |
| first layer on both | 4 | `H < 7/8` |
| first two layers on both | 6 | `H < 11/12` |

**Gap 1: mixed pairs.** Assumption 3.1 needs `q D_min > 1`. The sharp line suggested by Theorem 4.1(i) is `r_X D_X + r_Y D_Y > 1`. So Theorem 3.1 leaves out the mixed region that Betken–Dehling Remark 2.1 names as open. In finite dimensions the authors say their condition matches (1.2), which is the sup-over-entries form. The gap is therefore not an artefact of the Hilbert setting.

**Gap 2: the increment kernel.** *Derived here.* Write `G(xi_t, xi_{t+1}) = phi(g(xi_{t+1})) - phi(g(xi_t))` on `H_1 = R^2`. This has Hermite rank one, with term `a_1 (xi_{t+1} - xi_t)`.
- `||R(v)||_op` decays like `v^{-D}`. So for HSIC, `q = 2` needs `D > 1/2`. Theorem 3.1 then gives nothing beyond `H < 3/4`.
- The true autocovariance of a layer decays like `v^{-qD-2}`, because of cancellation. The theorem's condition does not see that cancellation.
- A proof must therefore go chaos by chaos with the fourth-moment theorem, arXiv 2509.20342. I found no paper that does this.

**Raw-level differencing on Gaussian data.** *Derived here.* Take `H_1 = R^2` with `X_t = (Delta xi_t, Delta eta_t)`. Then `R(v)` decays like `v^{-D-2}`, which is summable. So Theorem 3.1 applies for every `H < 1`. This is the Hilbert form of Bardet–Surgailis.

**Estimated scores.** *Derived here, first order only.* By Dehling–Taqqu reduction, `Phi^{-1}(F_n(x))` is approximately `x - xi-bar`. So normal scores act like centring by the sample mean. Bai–Taqqu 2019 Theorem 3.4 covers polynomial `G` with a subtracted sample mean. It says the limit can change. I found no theorem for an RKHS feature with estimated scores.

**Status warning.** Düker–Zoubouloglou is an arXiv preprint (v3, 22 September 2026). I found no journal version. Its v1 relied on a Hilbert-space fourth-moment result that a later note showed uses a metric that does not metrize weak convergence. Version 3 replaces it. Cite v3 only.

---

## Q5. Other remedies for dependence

**Verdict: known under mixing, partly known for shifts without mixing.** Audit items are not repeated.

**Mixing-based kernel and distance rules (audit, M1, M2).** Chwialkowski–Gretton 2014 shift HSIC. Chwialkowski–Sejdinovic–Gretton 2014 wild bootstrap. Davis et al. 2018 auto-distance covariance. All assume conditions that exclude fGn with `H > 1/2`.

**Matsui, Mikosch, Samorodnitsky 2017 (Abstract).** They define distance covariance for two stochastic processes on an interval, and suggest a test of independence. I did not read how the sample is formed, so I do not claim it excludes one long path.

**Shift tests with no mixing condition. This is the most important finding in Q5.**

- **Harris 2020, arXiv 2012.06862 (Full).** Let `Y` be stationary and independent of `X`. `X` may be arbitrary. Compute any association `V_s = V(X[N:T-N], Y[s+N:s+T-N])` for `s = -N..N`. Let `m` count the shifts with `V_s >= V_0`.
  - The Corollary gives `P(m < M) <= M/(N+1)` for every finite `N`.
  - Theorem 2 shows the bound cannot be improved.
  - Theorem 3 shows a twice-as-powerful rule is valid as `N -> inf` with `D` fixed. The proof handles non-ergodic stationary series too.
  - No mixing appears anywhere. So the conservative test holds size for fGn at every `H`, for HSIC or any other `V` (*my reading, direct from the Corollary*).
  - At level 5% the conservative rule needs `N >= 19` and rejects only if `V_0` beats all 38 shifted values.
  - Harris writes that decorrelating "may increase statistical power but is not required to avoid false positives".
- **Liang, Dyson, Barber, Holz 2025, arXiv 2512.17372v2 (Full on Sec. 3–4).** This has two results that matter.
  - Theorem 1: for a local window statistic and any stationary pair, `P(p_t <= alpha) <= alpha / mar(t)`, with inflation about 4 at mid-sample. They say this allows "arbitrarily strong dependence", and they show the factor is tight.
  - Section 4.2: for a whole-series statistic with wrap-around shifts, they say a guarantee "is not possible without further assumptions". Theorem 3 then needs beta-mixing and a stability coefficient.
  - By M1 the beta term does not vanish for fGn. So the paper's sentence "No long-memory theory covers circular shifts" stands. These authors now support it in print. Cite them.

**Frequency-domain independence tests.** **Haugh 1976 (Record)** and **Hong 1996 (Record)** are the classical residual cross-correlation tests. Hong is cited as the basis of Shao's extension. **Shao 2009 (audit)** allows long memory and is oversized at small `n`. It is linear. **Besserve 2013** is the kernel version (Q2).

**Self-normalisation.** **Bai, Taqqu, Zhang 2016 (audit)** is for the mean. **Ghoshal and Shao 2026 (audit)** is the nearest kernel design and says its condition "requires the underlying temporal dependence to be of short-range type". Bai and Taqqu 2018 Section 6 recommend combining resampling with self-normalisation as the way around rank instability. I found no kernel self-normalised test under long memory.

**Block and subsampling.** **Betken and Wendler 2018 Theorem 1** and Betken–Dehling Theorem 2.4 are the validity results. The audit records Betken–Dehling's oversize. Bai and Taqqu 2017, Zhang et al. 2022 and Tewes 2018 are in `litreview.tex`.

**2023–2026.** Liang et al. 2025 is the only addition beyond `litreview.tex`. The Boniece et al., Diz-Castro et al., Podkopaev et al. and Shen et al. entries there still stand. My searches found no 2023–2026 kernel test on long-memory data.

---

## Q6. Impossibility or trade-off results

**Verdict: partly known.** Hard theorems exist for close neighbours. The layer-level trade-off is proved. No theorem was found for kernel tests or for the two-series common-factor problem.

**Impossibility in neighbouring problems.**

- **Pötscher 2002, Econometrica (Abstract).** Spectral density at frequency zero is an ill-posed target. The minimax risk is infinite at every sample size, and confidence sets are close to uninformative. "Similar results" hold for the long-memory parameter.
- **Faust 1996, Econometric Theory (Abstract).** Under general models, the size of unit root tests "will converge to one". Solving this requires a model that is small in a topological sense.
- **Blough 1992, JAE (Abstract).** "the maximum power of a generic unit root test against any stationary alternative is equal to the true level of the test."
- **Shah and Peters 2020, Ann. Statist. (Abstract).** A valid conditional independence test with continuous `Z` has no power against any alternative.
- **Lundborg, Shah, Peters 2022, JRSS-B (Abstract).** For functional `Z`, even in the Gaussian case, "the power of any test cannot exceed its size".
  - *Derived here, by analogy only:* a common slow factor acts like an unobserved functional confounder. Independence given such a factor looks like this problem. No source makes this link for our setting.

**Long memory against level shifts.**

- **Diebold and Inoue 2001 (Abstract).** Stochastic regime switching "is easily confused with long memory, even asymptotically", provided the switching is small in a sense they make precise.
- **Berkes, Horváth, Kokoszka, Shao 2006 (Abstract).** They give a test that separates weakly dependent series with a mean change from long memory. So separation is possible inside restricted classes.
- So the literature shows hardness and confusion, plus tests under restrictions. I found no theorem that a single path cannot separate them in general.

**The layer-level trade-off.**

- **Bai and Taqqu 2019, SPA, Proposition 4.1 (Full).** In the Gaussian case, Hermite rank equals power rank. The power rank is the first nonzero derivative at zero of `x -> E G(Z + x)`. They attribute the statement to **Ho and Hsing 1997** (Abstract).
  - So a rank-`k` feature responds to a small mean shift only at order `mu^k`.
  - The formula at `paper.tex:436` is the `k = 1` case. Cite this result.
- **Tewes 2017 (Full on the setup, arXiv 1506.07296).** Local alternatives for distributional change-point tests, including changes of Hermite rank. For a mean shift, the rank is unchanged and relative efficiency to CUSUM is one.

**Power cost, kernel level. Not found after these searches.** *Everything in this block is derived here and is heuristic.*

1. **Mean shift, MMD-type.** A rank-one feature detects at `mu ~ n^{H-1}`. A rank-two feature has signal `mu^2`. Its noise is `n^{-1/2}` below `H = 3/4` and `n^{2H-2}` above it. So it detects at `mu ~ n^{-1/4}` below `3/4` and `mu ~ n^{H-1}` above. The rate cost exists only below `H = 3/4`. Above it, the cost is in the constant.
2. **Gaussian dependence, HSIC-type.** Mehler's identity gives `Cov(He_k(xi), He_l(eta)) = delta_{kl} k! rho^k`. Removing the first layer on both sides turns the leading HSIC term from `rho^2` into `rho^4`. In the central regime, detection moves from `rho ~ n^{-1/2}` to `rho ~ n^{-1/4}`.
3. **Common factor, increment kernel.** Take `X = F + e`, `Y = F + u`, with `F` unit-variance fGn and noise variance `s^2`. The correlation is `1/(1+s^2)` in levels and `(1 - rho_F(1))/(1 - rho_F(1) + s^2)` in increments. At `H = 0.9` and `s^2 = 1` it falls from `0.50` to `0.21`.
4. **Constant factor.** If `F` does not move inside the sample, the increment HSIC target is exactly zero. The test then has power equal to size in that direction.

---

## Table: construction, closest prior work, what we must cite

| Construction | Closest prior work | How close | Must cite |
|---|---|---|---|
| 1a. Hermite-deflated Gram `Q K Q` | KCI `R_Z K R_Z`. Conditionally positive definite contrast form `M^T K M`. | **Algebra identical** in limit (derived). Purpose differs: conditioning, numerical approximation. | `zhang2011kci`, `davydov2021cpd`, `betken2021ldr` (Thm 4.1), `bai2019sensitivity`, `bai2018instability`. |
| 1b. Population `phi - mu - a_1 T` | Plumlee–Joseph `c*`. Ho–Hsing expansion. Dehling–Taqqu reduction. | **Close.** `c*` also has zero first Hermite coefficient (derived). Ho–Hsing subtracts layers for theory. Neither targets a test. | `plumlee2018orthogonal`, `ho1996expansion`, `dehling1989empirical`, `betken2021ldr` (Thm 2.7), `duker2026breuermajor`. |
| 2a. Increment kernel `D K D'` | Signature kernel, preference kernel, kernel SFA. Structured HSIC framework. | **Form identical.** Structured HSIC contains it as a special case (derived). Single-path test under long memory not found. | `kiraly2019sequential`, `houlsby2012preference`, `herbrich2000rank`, `bohmer2011rsksfa`, `zhang2008noniid`, `bardet2011roughness`, `granger1974spurious`. |
| 2b. Increment kernel, lag `L` | Structured HSIC with lag `tau`. Sun 2003 for the seasonal sum rule (litreview). | **Moderate.** Lagged cliques exist. Seasonal feature differencing not found. | `zhang2008noniid`, `sun2003gegenbauer` (in `litreview.bib`). |
| 3. Linear filters on features | KCSD. Generalized spectrum. DCCA. Müller–Watson cosine projection (the complement). | **Moderate.** KCSD filters in the RKHS by Fourier. Fractional differencing and DCT deflation in an RKHS not found. | `besserve2013kcsd`, `hong1999generalized`, `muller2018covariability`, `clausel2014wavelet`, `podobnik2008dcca`, `hosking1981fractional`. |
| 4. Even or symmetrised kernels | Kernel invariance. Rank-two scalar statistics. Hermite rank theory for even `G`. | **Weak on kernels, strong on the rank fact.** The rank claim is classical. Instability under shift is proved. | `haasdonk2007invariant`, `levyleduc2011scale`, `bai2019sensitivity` (Thm 2.1, 2.2, 3.4), `tewes2017localalt`. |

---

## Risks to our claim

1. **Do not claim that size control was open.** Harris (2020) gives a finite-sample bound under stationarity alone. Liang et al. (2025) give another. A referee will ask why a kernel is needed. The defensible answer is power, simplicity or a sharper rule. Their tests also pay for size: the conservative rule needs at least 38 shifts and a shorter segment.
2. **Rank above one is unstable.** Bai–Taqqu 2019 Theorem 2.1: a shift makes the rank one at all small nonzero shifts. Theorem 2.2: scaling does the same unless `G` is symmetric. Bai–Taqqu 2018 Section 6: "one may adopt the assumption that the rank is always 1".
   - Any stationary-marginal deflation depends on `T` and `F` being right.
   - Estimated location, scale or a wrong marginal puts rank one back. Report simulations with a misspecified `T`.
   - Bai–Taqqu 2019 Theorem 3.1 gives the crossover for a shrinking deterministic shift, `x_N ~ N^{(1/2-H)/(k-1)}`.
3. **Differencing does not restore mixing.** *Derived here.* Take tent weights on `2n` points, with a second tent `n` steps later. Summing by parts turns the increments into differences of block means. By self-similarity the correlation does not depend on `n`. I computed it for `n = 50, 200, 800`: about `-0.006`, `-0.019`, `-0.035` at `H = 0.6, 0.75, 0.9`.
   - So differenced fGn is not completely regular. Through Kolmogorov–Rozanov (Bradley 2005) it is not strongly mixing either.
   - Also, `tau(r) >= c r^{H-2}`, so the CSG condition fails too.
   - Consequence: no existing mixing-based kernel theorem applies after differencing. We cannot cite them for the increment kernel.
4. **The increment test has a blind spot by design.** It answers a high-pass question. It is blind to dependence at frequency zero, and it cannot see a factor that does not move inside the sample. Lag `L` is blind at all harmonics of `1/L`, which removes shared seasonality and any real dependence there. Say this first, not as a footnote.
5. **The kernels are not characteristic.** `k_D` is not universal on pairs. Its HSIC target can vanish under dependence. The MMD with `k_D` has mean element zero for any stationary law. Deflated kernels lose all first-layer signal. Do not describe either as a general independence test without qualification.
6. **Finite-sample behaviour near `H = 1/2`.** Differencing i.i.d. data gives an MA(1) with lag-one correlation `-1/2` (*derived here*). This is over-differencing in the Plosser–Schwert sense. Positive-only multiplier schemes may then misbehave. The increment kernel needs its own baseline at `H = 0.5`.
7. **Quadratic statistics in one filtered feature.** Clausel et al. show a Rosenblatt limit for squared wavelet coefficients of Hermite polynomials. HSIC normalisers and any MMD-type statistic may inherit slow convergence. A central-looking leading term does not protect the finite-sample scale.
8. **Referee overlap.** Zhang et al. 2008 has Gretton as an author. Besserve et al. 2013 and Chwialkowski–Gretton 2014 come from the same community. An MMD or HSIC reviewer may know this framework well. Cite it up front.
9. **Preprint dependence.** The central theorem rests on arXiv 2405.11452v3 and 2509.20342v3. Betken–Dehling is also an arXiv version in our bibliography. If a proof in the paper rests on these, say so.

**Safe wording.** "We could not find a prior kernel test built to remove low Hermite layers or to difference features for independence testing under long memory, after the searches listed in the review." Do not use stronger wording.

---

## Corrections and additions to existing repo notes

1. **`litreview.tex`, entries for `dehling2017twosample`.** The entry and gap map say "Non-degenerate kernels only". Theorem 2.1 as printed has no such restriction. It covers any rank `m` with `Dm < 1`, including `a_{1,1} != 0`. So the "split MMD null limit is open" gap may be narrower than stated for `D < 1/2`, at least for the cross term. A referee could raise this. Someone should check it against the exact MMD statistic.
2. **`paper.tex`, Discussion.** The displayed identity should cite Bai–Taqqu 2019 Proposition 4.1 and Ho–Hsing 1997.
3. **`paper.tex`, line 61 and Conclusion.** "No long-memory theory covers circular shifts" stands. Add Harris 2020 and Liang et al. 2025. State that a non-circular variant has a finite-sample bound.
4. **Closest prior work.** The audit should add Betken–Dehling Theorem 4.1(i) and Theorem 2.7. It discussed only Theorem 2.2, Corollary 2.1 and Theorem 2.4.
5. **Do not cite arXiv 1610.00684.** The authors withdrew it. Use `bai2019sensitivity`, which contains the same result as Proposition 4.1.
6. **Tsay–Chung.** The audit read this second-hand. I have now read the abstract first-hand, and it includes the first-differencing warning.

---

## Leads dropped or left unverified

**Dropped from memory, not verified this session.** I name them only so that nobody cites them on my authority. Bray and Martinez 2002 (kernel slow features). Micchelli 1986, Wahba 1990 and Smola et al. 1998 (polynomial null spaces, semiparametric SVM). Bergsma and Dassios 2014. Schölkopf and Smola 2002. Ebisuzaki 1997. Prichard and Theiler 1994. Wu 2003. Beran 1991. Samarov and Taqqu 1988. Adenstedt 1974. Perron and Qu 2010. Davidson and Sibbertsen 2005. Yau and Davis 2012. Cochrane 1991. Bahadur and Savage 1956. Domowitz and El-Gamal 2001. Patterson and Thompson 1971 (error contrasts). Theiler et al. 1992 has a Crossref record but I did not read it, so I dropped it too.

**Existence confirmed only, so no BibTeX entry.** I read nothing in these, or they carry no claim here. Istas and Lang 1997 and Coeurjolly 2001: Bardet and Surgailis name them for generalized quadratic variations. Tudor and Viens 2009. Dehling and Taqqu 1991: Lévy-Leduc et al. and DRW cite it for one-sample U-statistics under long memory. Podobnik et al. 2011. Müller 2014, Dufour 1997 and Faust 1999: abstracts read, but they carry no separate claim here.

**One lead worth checking.** The classical universal-kriging and error-contrast literature probably holds the earliest form of `Q K Q` with a polynomial mean. Plumlee and Joseph discuss universal kriging directly. I verified only Davydov for the matrix form.

---

## Exact search strings

**WebSearch.**

```
Hilbert space valued Breuer-Major theorem Gaussian subordinated functional central limit Hermite rank long-range dependence
kernel independence test time series differencing increments signature kernel HSIC long memory fractional
removing Hermite rank one component U-statistic degenerate kernel long-range dependence raise Hermite rank projection robust test
spurious regression fractionally integrated differencing cure power loss cointegration low-frequency common factor Tsay Chung
functional Breuer-Major theorem Nourdin Nualart Hilbert-valued central limit Gaussian functional tightness conditions
signature maximum mean discrepancy two-sample test stochastic processes rough paths Chevyrev Oberhauser signature moments
kernel cross-spectral density operator independence test time series frequency domain Fourier HSIC
GP regression on time prewhitening kernel conditional independence confounder time RCoT PCMCI "time index" HSIC autocorrelated residuals nonstationary
impossibility long memory versus mean shift single path cannot distinguish change point long-range dependence testing Mikosch Stărică Berkes Horváth Kokoszka Shao
Harris "shift test" independence generic time series conservative stationary autocorrelated no mixing assumption
kernel test long-range dependence 2025 2026 HSIC MMD long memory time series bootstrap self-normalization arXiv
independence between two time series impossibility single realization no valid test hardness "time series" nonparametric no-free-lunch persistence
kernel slow feature analysis differences in feature space Gram matrix phi(x_{t+1}) - phi(x_t) regularized sparse kernel SFA
detrended cross-correlation analysis DCCA coefficient long memory spurious correlation test independence power law cross-correlated processes statistical test
preference kernel GP pairwise k((u,v),(u',v')) = k(u,u') + k(v,v') - k(u,v') - k(v,u') collaborative Gaussian processes Houlsby
wavelet coefficients whitening 1/f fractional noise independence test two series wavelet cross-correlation long memory prewhitening decorrelation Hermite non-Gaussian scalogram Clausel Roueff Taqqu Tudor
Granger Newbold 1974 "Spurious regressions in econometrics" pdf first differences recommendation Plosser Schwert overdifferencing costs
Tsay Chung 2000 spurious regression fractionally integrated processes abstract "long memory" t-ratio diverge sum d greater 1/2 pdf
Potscher 2002 ill-posed estimation long memory parameter confidence sets infinite expected length lower risk bound abstract Econometrica
Flaxman Neill Smola Gaussian processes for independence tests with non-iid data causal inference residual HSIC regress on time space spatial temporal confounding arXiv pdf
Plosser Schwert 1978 "Money, income, and sunspots" differencing overdifferencing cost regression integrated abstract conclusions
Engle Granger 1987 co-integration abstract "vector autoregression in differenced variables is incompatible"
Phillips 1998 "New tools for understanding spurious regressions" abstract Karhunen-Loeve deterministic trends stochastic coordinate system low frequency
Müller 2014 "HAC corrections for strongly autocorrelated time series" abstract size control power tradeoff impossible long-run variance persistence uniform
Gaussian copula HSIC distance covariance Hermite expansion Mehler formula first order correlation rho^k decomposition kernel dependence measure eigenfunction Hermite
HSIC applied to differenced time series first differences before kernel independence test returns levels spurious kernel dependence nonstationary trend
"long-range dependence" OR "long memory" "Hermite rank" reduce transformation choose function "higher Hermite rank" test statistic to avoid long memory effect robust inference
"increment" kernel OR "difference kernel" k(x,y)-k(x,y')-k(x',y)+k(x',y') two-sample independence test pairs consecutive observations time series RKHS feature differences stationarity
Hong 1999 "Hypothesis testing in time series via the empirical characteristic function" generalized spectral density approach JASA abstract
Müller Watson low-frequency fractional long memory independent cross-correlation low frequency common factor effective sample size q cosine "long-run correlation" spurious regression fractional inference
spurious correlation between independent long memory series removed by differencing or detrending kernel; "detrended" OR "filtering" "fractional cointegration" "narrow band" power loss differenced regression low frequency relationship lost
circular shift permutation test independence time series validity theory cyclic shift exact finite sample stationary long memory rotation test surrogate 2023 2024 2025
kernel orthogonal to polynomials low-order removal RKHS independence test "Hermite" feature map kernel dependence copula normal scores HSIC rank transform kernel projected "orthogonal complement"
Plumlee Joseph "Orthogonal Gaussian process models" Statistica Sinica abstract orthogonal to linear effects identifiability covariance function
conditionally positive definite kernel order m polynomial null space Gram matrix projection semiparametric support vector kernel universal kriging error contrasts REML annihilate polynomial
```

**arXiv export API, field queries.** Each returned the count shown.

```
abs:"Hermite rank" AND abs:kernel  -> 2
abs:"Breuer-Major" AND abs:Hilbert  -> 2
abs:"long-range dependence" AND abs:kernel AND abs:Hermite  -> 2
abs:"signature kernel" AND (abs:independence OR abs:"two-sample" OR abs:MMD)  -> 10
abs:"long memory" AND abs:RKHS  -> 1
abs:"long-range dependent" AND abs:"reproducing kernel"  -> 0
abs:"long-range dependence" AND abs:"Hilbert space" AND abs:"limit theorem"  -> 2
abs:"long memory" AND abs:"functional time series" AND abs:Hermite  -> 0
abs:"long-range dependent" AND abs:"non-central limit" AND abs:"Hilbert"  -> 2
abs:"fractional" AND abs:"Hilbert space" AND abs:"non-central limit theorem"  -> 0
abs:differencing AND abs:"kernel" AND abs:"independence test" AND abs:"time series"  -> 0
abs:prewhitening AND (abs:HSIC OR abs:"kernel" OR abs:"distance covariance")  -> 4
abs:"preference kernel" OR abs:"pairwise kernel"  -> 20
abs:"Hermite rank" AND (abs:"kernel" OR abs:"U-statistic" OR abs:"feature")  -> 3
abs:"Hermite rank" AND abs:test  -> 2
abs:"Hermite rank" AND abs:increase  -> 1
abs:"long memory" AND abs:"kernel" AND abs:test AND abs:independence  -> 0
abs:"fractional Gaussian noise" AND (abs:kernel OR abs:MMD OR abs:HSIC) AND abs:test  -> 0
abs:"second-order" AND abs:"quadratic variation" AND abs:"fractional Brownian" AND abs:"central limit"  -> 3
abs:"Hermite variations" AND abs:"filter"  -> 0
abs:"over-differencing" OR abs:overdifferencing  -> 1
abs:"long-range dependent" AND abs:"U-statistic" AND abs:degenerate  -> 0
abs:"shift test" AND abs:independence AND abs:"time series"  -> 0
abs:"circular shift" AND abs:independence AND abs:"time series"  -> 0
abs:"long-run covariability"  -> 20
abs:"long memory" AND abs:"independence" AND abs:"two" AND abs:"time series" AND abs:test  -> 2
abs:"long-range dependence" AND abs:independence AND abs:test AND abs:"time series"  -> 2
abs:"energy distance" AND (abs:"long-range" OR abs:"long memory")  -> 1
abs:"Hilbert-Schmidt" AND (abs:"long-range" OR abs:"long memory")  -> 5
abs:"maximum mean discrepancy" AND (abs:"long-range" OR abs:"long memory")  -> 1
abs:"Stein" AND abs:"kernel" AND (abs:"long-range dependence" OR abs:"long memory")  -> 2
abs:"circular" AND abs:"shift" AND abs:"Hurst"  -> 1
abs:"cyclic shift" AND abs:"independence" AND abs:"test"  -> 1
abs:"generalized spectral" AND abs:independence AND abs:"time series"  -> 2
abs:"Hermite" AND abs:"U-statistics" AND abs:"long-range" AND abs:degenerate  -> 0
abs:"Rosenblatt" AND abs:"U-statistic" AND abs:"long memory"  -> 0
abs:"kernel" AND abs:"first differences" AND abs:"independence" AND abs:"time series"  -> 0
abs:"feature space" AND abs:"differencing"  -> 5
abs:"RKHS" AND abs:"fractional differencing"  -> 0
abs:"frequency domain" AND abs:"independence test" AND abs:"kernel" AND abs:"time series"  -> 0
abs:"high-pass" AND abs:kernel AND abs:"independence test"  -> 0
abs:"detrended" AND abs:"distance correlation"  -> 0
abs:"detrended cross-correlation" AND abs:kernel  -> 0
```

**arXiv export API, title lookups.**

```
ti:"Kernels for sequentially ordered data"
ti:"instability of ranks"
ti:"Sensitivity of the Hermite rank"
ti:"Hermite rank" AND ti:"power rank"
ti:"wavelet" AND ti:"Hermite" AND abs:"long memory"
ti:"Gaussian Processes for Independence Tests"
ti:"new correlation coefficient" AND ti:"orthogonal decomposition"
ti:"Variations and estimators for self-similarity parameters"
ti:"Robust estimation of the scale" AND abs:"long-range"
ti:"Slow Feature Analysis" AND ti:kernel
ti:"Hilbert space-valued linear processes" OR ti:"long-range dependent curve time series"
ti:"fourth moment" AND ti:"Hilbert" AND au:Duker
ti:"Yule" AND ti:"nonsense correlation"
ti:"increment ratio" OR ti:"roughness of random paths"
ti:"quadratic variations" AND abs:"fractional Brownian" AND abs:"central limit"
ti:"Distance covariance for discretized stochastic processes"
ti:"Distance covariance for stochastic processes"
ti:"Independence Testing for Temporal Data"
ti:"Spectral analysis of long range dependence functional time series"
```

**Crossref `query.bibliographic`.** These confirmed records, and they found a few neighbours.

```
Robust estimation of the scale and of the autocovariance function of Gaussian short- and long-range dependent processes
Bivariate symmetric statistics of long-range dependent observations
Deriving the asymptotic distribution of U- and V-statistics of dependent data using weighted empirical processes
How the instability of ranks under long memory affects large-sample inference
Sensitivity of the Hermite rank
Hermite rank, power rank and the generalized Weierstrass transform
On the asymptotic expansion of the empirical process of long-memory moving averages
Limit theorems for functionals of moving averages Ho Hsing
Density estimation under long-range dependence Csorgo Mielniczuk
Kernels for sequentially ordered data Kiraly Oberhauser
Large margin rank boundaries for ordinal regression Herbrich
Preference learning with Gaussian processes Chu Ghahramani
Gaussian processes for independence tests with non-iid data in causal inference Flaxman Neill Smola
Kernel-based conditional independence test and application in causal discovery Zhang
Causal discovery from heterogeneous/nonstationary data Huang Zhang
Statistical analysis of coupled time series with kernel cross-spectral density operators Besserve
Orthogonal Gaussian process models Plumlee Joseph
ANOVA kernels and RKHS of zero mean functions for model-based sensitivity analysis
Conditional independence testing in Hilbert spaces with applications to functional data analysis Lundborg
The hardness of conditional independence testing and the generalised covariance measure Shah Peters
Invariant kernel functions for pattern analysis and machine learning Haasdonk
Spurious regressions in econometrics Granger Newbold 1974
New tools for understanding spurious regressions Phillips 1998
Spurious regressions with stationary series Granger Hyung Jeon
The spurious regression of fractionally integrated processes Tsay Chung
A convergent t-statistic in spurious regressions Sun 2004
Money, income, and sunspots: measuring economic relationships and the effects of differencing Plosser Schwert
Co-integration and error correction: representation, estimation, and testing Engle Granger
Long-run covariability Muller Watson Econometrica
Testing models of low-frequency variability Muller Watson
Low-frequency econometrics Muller Watson
Semiparametric fractional cointegration analysis Robinson Marinucci 2003
Semiparametric estimation of multivariate fractional cointegration Chen Hurvich
Lower risk bounds and properties of confidence sets for ill-posed estimation problems Potscher
Near observational equivalence and theoretical size problems with unit root tests Faust
Conventional confidence intervals for points on spectrum have confidence level zero Faust
The relationship between power and level for generic unit root tests in finite samples Blough
Some impossibility theorems in econometrics with applications to structural and dynamic models Dufour
HAC corrections for strongly autocorrelated time series Muller
On discriminating between long-range dependence and changes in mean Berkes Horvath Kokoszka Shao
Long memory and regime switching Diebold Inoue
A test against spurious long memory Qu
The Hurst effect under trends Bhattacharya Gupta Waymire
Estimation of the memory parameter for nonstationary or noninvertible fractionally integrated processes Hurvich Ray
An efficient taper for potentially overdifferenced long-memory time series Hurvich Chen
Discrimination between monotonic trends and long-range dependence Kunsch
Fractional differencing Hosking
An introduction to long-memory time series models and fractional differencing Granger Joyeux
Central limit theorems for non-linear functionals of Gaussian fields Breuer Major
Non-central limit theorems for non-linear functionals of Gaussian fields Dobrushin Major
Convergence of integrated processes of arbitrary Hermite rank Taqqu
Limit theorems for nonlinear functionals of a stationary Gaussian sequence of vectors
The empirical process of some long-range dependent sequences with an application to U-statistics Dehling Taqqu
Checking the independence of two covariance-stationary time series: a univariate residual cross-correlation approach
Testing for independence between two covariance stationary time series Hong
Kernel measures of independence for non-iid data
Measuring nonlinear dependence in time-series, a distance correlation approach Zhou
Testing for nonlinearity in time series: the method of surrogate data Theiler
Long-range dependent curve time series Li Robinson Shang
On the validity of resampling methods under long memory Bai Taqqu
Occasional structural breaks and long memory Granger Hyung
Nonstationarities in financial time series, the long-range dependence, and the IGARCH effects
Detrended cross-correlation analysis: a new method for analyzing two nonstationary time series
Measuring correlations between non-stationary series with DCCA coefficient
Statistical tests for power-law cross-correlated processes
Quadratic variations and estimation of the local Holder index of a Gaussian process
Estimating the parameters of a fractional Brownian motion by discrete variations of its sample paths
Variations and estimators for self-similarity parameters via Malliavin calculus Tudor Viens
Estimation of a non-invertible moving average process: the case of overdifferencing
Regularized sparse kernel slow feature analysis
Optimizing search engines using clickthrough data
The Hermite rank of a special subordinated long-memory Gaussian process
Measuring the roughness of random paths by increment ratios
Limit theorems for Hilbert space-valued linear processes under long range dependence
Signature moments to characterize laws of stochastic processes
Consistent testing for pairwise dependence in time series Fokianos Pitsillou
Testing for structural change in a long-memory environment
Semiparametric frequency domain analysis of fractional cointegration
Asymptotic properties of U-processes under long-range dependence
The functional Breuer-Major theorem Nourdin Nualart
Applications of distance correlation to time series Davis Matsui Mikosch Wan
Approximation with conditionally positive definite kernels on deficient sets Davydov
```

**Crossref queries that failed.**

- HTTP 429, never recovered: `Understanding spurious regressions in econometrics Phillips 1986`, `Distance covariance for stochastic processes Matsui Mikosch Samorodnitsky`, `Copula-based kernel dependency measures Poczos`, `Comparison of tests of the independence of two covariance-stationary time series Geweke`, and the Hong 1999 title.
- No matching record: `A kernel test of goodness of fit` and `A kernelized Stein discrepancy for goodness-of-fit tests`.
- Recovered another way: Phillips 1986 came back through the 1998 query. Matsui, Hong 1999 and Chwialkowski 2016 were opened directly on arXiv, the author's site and PMLR.
- Dropped as unverified: Poczos and Geweke.
