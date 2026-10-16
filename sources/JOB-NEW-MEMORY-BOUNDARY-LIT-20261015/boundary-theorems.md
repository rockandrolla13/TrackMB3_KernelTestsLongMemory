# Boundary theorems — what is actually proved about H = 3/4

Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 · execution `claude-memboundary-20261015` · 2026-10-15.
Plan item (2). Evidence: F = full-text section read (PDF text extracted locally with `pdftotext`),
S = summarising fetch, A = abstract/snippet, R = reused from JOB-NEW-MKC-LIT-20261015.
"Worker inference" marks every step that is mine and not stated in a source.

Notation. D = LRD exponent of the autocovariance, ρ(k) = k^{−D}L(k). H = 1 − D/2. d = H − 1/2.
So D = 1/2 ⇔ H = 3/4 ⇔ d = 1/4.

---

## 1. Betken & Dehling, arXiv:2107.03041 v3 (27 Jan 2026; header "Submitted to Bernoulli") — F

PDF sha256 `5bc4b26cc738798e641af76eeddb27e50e6b43d21aeb0cb703680d6ba52b50df` (scratch copy).

**Setting.** X_i = G₁(ξ_i), Y_i = G₂(η_i); ξ, η are *independent* stationary Gaussian LRD processes
with D_ξ, D_η ∈ (0,1). Statistic: Σ_h a_h dcov_n(X,Y;h), i.e. distance covariance (weight
w(s,t) = c|s|⁻²|t|⁻² on the empirical characteristic-function process). Null = independence of the
two series.

**Theorem 2.2 (the boundary).**
- (i) If D_ξ, D_η ∈ (1/2, 1) — both H < 3/4 — then √n(φ̂_{XY;h} − φ̂_X φ̂_Y) → complex Gaussian
  process, so n·dcov_n → ∫∫|Z|²w (Corollary 2.1(i)). Holds for **subordinated** Gaussian G₁, G₂ with
  E|X| < ∞.
- (ii) If G₁ = G₂ = id (**Gaussian only**) and D_ξ + D_η < 1, the normalisation is
  n^{(D_ξ+D_η)/2} and the limit is a double Wiener–Itô (Rosenblatt-type) integral.
- Remark 2.1(i), verbatim in substance: the theorem "does not cover the case D_ξ + D_η > 1, but
  requires either D_ξ < 1/2 or D_η < 1/2." So the mixed case (one series above 3/4, one below,
  with D-sum > 1) is not covered by either part.
- Remark 2.1(ii): "In real-life data, one typically encounters LRD coefficients that are larger
  than 1/2, which corresponds to Hurst coefficients smaller than 0.75." No source is given for that
  claim in the remark.

**Calibration.** Limits are "unknown", so critical values come from subsampling (Theorem 2.4),
block length l_n = O(n^{(1+min(D_X,D_Y))/2 − ε}), subordinated Gaussian, Assumptions 1–2.
No bootstrap and no self-normalisation.

**Simulation (Section 4, Supplement Tables 1–4).** fGn via R `longmemo::simFGN0`; 5,000 reps;
n ∈ {100, 300, 500, 1000}; l_n = n^γ, γ ∈ {0.4, 0.5, 0.6}; **H ∈ {0.6, 0.7} only**. No H above 0.75
was simulated. Empirical size at 5% (Table 1, "linear" design, r = 0): dCov 0.073–0.131 at H = 0.6
and 0.097–0.175 at H = 0.7; it falls with n. **Contrary point to keep:** the test is clearly
oversized *below* the boundary at n ≤ 1000 (up to 3.5× nominal at H = 0.7).

**What the H = 3/4 line really is (worker inference, checked by formula — §5).** The statistic is
built from the cross term f(X_j,Y_j) = (e^{isX_j} − φ_X)(e^{itY_j} − φ_Y). With ξ ⊥ η, its
autocovariance is the product c_X(k)·c_Y(k) ≲ |ρ_ξ(k)ρ_η(k)| = k^{−(D_ξ+D_η)}. The √n Gaussian
limit needs that product to be summable: **D_ξ + D_η > 1, i.e. H_X + H_Y < 3/2.** The single-number
line H = 3/4 is the special case H_X = H_Y. Theorem 2.2(i) asks for more (each D > 1/2), which is a
proof restriction, not the true boundary.

## 2. Zhai, Zhong & Wu, arXiv:2609.30565 v1 (24 Sep 2026, stat.ME) — F

PDF sha256 `662f1397541c0d5c1440bd936aa341a699696f3569edcb80beabf811cffbfc07` (scratch copy).

**Setting.** p-dimensional stationary **Gaussian linear** process X_t = Σ A_l ε_{t−l}, ε i.i.d. N(0,I).
**Condition 1:** max_j |(A_l)_{j·}| ≤ C₀(1∨l)^{−β}, **β > 3/4**. Coefficient decay β gives
d = 1 − β, so β > 3/4 ⇔ d < 1/4 ⇔ H < 3/4. β ≤ 1 is long memory.
**Statistic.** Built from the sample covariance Σ̂_n (local inverse; Le et al. 2023) — a
**quadratic (Hermite rank 2)** functional. Results: Berry–Esseen-type Gaussian approximation
(Theorem 1) and overlapping block bootstrap validity with l = n^{2/3} (Theorem 2), both under
Conditions 1–5; two-sample extension (Corollary 3); consistency-empowered correction (Theorem 4/6).
**Boundary statement (Section 6):** "The theory relies on temporal dependence β > 3/4. For
ultra-long-memory time series with 1/2 < β ≤ 3/4, the asymptotic distribution will be
Rosenblatt-type, and is beyond the scope of this paper."
**Simulation (Section 4.1).** n ∈ {200, 400, 800}, p ∈ {2000, 4000, 8000}, β ∈ {2, 0.9, 0.55}
(H ≈ short, 0.6, 0.95), 200 copies. β = 0.55 "fails expectedly" (QQ plots, Fig. 3). Only one
point above the boundary; no grid near 0.75.
**Dependency.** Machinery comes from Zhai, Kolar & Wu (2026), "Simultaneous inference for
covariance and precision matrices of long-range dependent time series", IEEE Trans. Inf. Theory
72(6):4246–4296 (as cited in 2609.30565; **not independently verified**, budget reached).

## 3. Why both papers land on 3/4: Hermite rank (worker inference, standard theory)

For a functional F(ξ_t) of a Gaussian LRD process with Hermite rank m, the autocovariance of F
decays like ρ(k)^m = k^{−mD}. The sample mean is √n-Gaussian (Breuer–Major) when mD > 1, and
has a non-central (Hermite / Rosenblatt for m = 2) limit at rate n^{−mD/2} when mD < 1
(Dobrushin–Major 1979, Taqqu 1979; both cited in Betken–Dehling and Dehling–Rooch–Wendler,
not re-read here). So:

| Effective rank of the statistic's leading term | √n-Gaussian iff | Boundary |
|---|---|---|
| 1 (linear: sample mean, MMD between segments, DML state term with monotone effect) | D > 1 | **H = 1/2 — any long memory breaks √n** |
| 2 (quadratic: sample covariance, ZZW; rank-2 effect such as τ(X)=X²−1) | D > 1/2 | **H = 3/4** |
| (1,1) cross term of two independent series (HSIC/dCov null, Pearson cross-correlation) | D_X + D_Y > 1 | **H_X + H_Y = 3/2** |

The 3/4 line is therefore a property of **rank-2 and cross-type statistics**, not of long memory in
general.

## 4. Does it extend to HSIC, MMD and DML scores?

**HSIC (independence of two series), translation-invariant characteristic kernel — close to direct,
not stated.** By Bochner, HSIC with k(x,x′)=ψ(x−x′), l(y,y′)=χ(y−y′) equals
∫∫|φ_{XY} − φ_Xφ_Y|² dΛ_k(s)dΛ_l(t), with finite spectral measures. That is Betken–Dehling's object
with a finite measure in place of w(s,t). Their Hilbert-space argument should go through, and the
integrability steps get easier because Λ is finite. Betken–Dehling's own text says the cross-
covariance/MMD/HSIC equivalence holds for a specific kernel only (search snippet, A; Sejdinovic et al.
2013). Worker inference: direct for part (i) and the subsampling theorem; part (ii) and non-Gaussian
G in the high-memory case remain open, as they are for dCov. **Not direct:** non-translation-
invariant kernels (e.g. signature kernels), conditional (KCI) tests, and the *wild bootstrap*
calibration (Chwialkowski–Sejdinovic–Gretton need τ(r) = O(r^{−6−ε}) for HSIC, R/F in prior job;
no LRD version found).

**MMD (two-sample / homogeneity) — not the same boundary.** Comparing two segments of one LRD series,
or two independent LRD series with the same marginal, the leading term of μ̂₁ − μ̂₂ in the RKHS is
rank 1 (it is the Hermite-1 coefficient E[k(·,G(ξ))ξ]·(ξ̄₁ − ξ̄₂)). Worker inference: the √n
limit fails for **every H > 1/2**, and n^{D}·MMD² tends to a low-rank weighted χ² rather than the
usual infinite sum. Consistent with Dehling–Rooch–Wendler (arXiv:1404.0551, R/F): two-sample
U-statistic processes are non-central when mD < 1, m = Hermite rank of the kernel.
Exception: rank-2 kernels/G (e.g. symmetric G with even kernel) get the 3/4 boundary back.

**DML / AIPW scores — the 3/4 line is the wrong line for C1.** Decompose ψ_t = (τ(X_t) − θ) + r_t.
Under sequential ignorability r_t is a martingale difference, so its mean is √n-Gaussian by the
martingale CLT whatever the memory of X (needs ergodicity and finite variance only). The state
term has the Hermite rank of τ∘G. Monotone heterogeneity (rank 1): the population-ATE estimator
converges at n^{H−1} = n^{d−1/2}, and HAC-Wald intervals fail for **any** d > 0, with no special role
for 1/4. Even heterogeneity (rank 2): boundary at d = 1/4. A second channel — cross-fitting
fold independence and nuisance error products — is separate. Ballinari–Wehrli (R/F) state the
cross-fitting conditions fail for persistent processes. I found no theorem giving its rank or rate.
Ciganovic–D'Amario–Tancioni (arXiv:2603.10999, Econometrics Journal 2026, A/S) assume an FCLT
via summable α-mixing or Maxwell–Woodroofe. Worker inference: Maxwell–Woodroofe fails for linear
processes with d > 0, so their conditions exclude long memory.

**Linear regression analogue (pre-emption for the cross rule).** Robinson & Hidalgo (1997, AoS 25(1)
77–104, doi:10.1214/aos/1034276622; A via LSE eprint) state that least squares "need no longer be
n^{1/2}-consistent" "when sufficiently strong spectral singularities in the error and a regressor
coincide". The product-memory condition d_x + d_u < 1/2 is the textbook form (not read in full
text; budget). It is the same rule as H_X + H_Y < 3/2 above.

**Hawkes Granger tests — no Hermite analogue.** Point processes are not Gaussian-subordinated. The
likelihood score for zero cross-excitation is a martingale. Validity hinges on the observed
information obeying a law of large numbers as ρ → 1 and with power-law kernels. The relevant limit
theory is nearly-unstable Hawkes (Jaisson–Rosenbaum, AAP 25(2) 2015, arXiv:1310.2033, A): rescaled
processes tend to integrated CIR, and with heavy tails to rough (fractional) limits. No H = 3/4
boundary is predicted. Expect degradation in ρ and in kernel tail index, not in H.

## 5. Formula check (no simulation)

`sanity_variance_scaling.py` computes n·Var(sample mean) exactly from fGn autocovariances for n = 10³…10⁶.
Output in `sanity_variance_scaling.out`. Growth per decade matches theory to three decimals:
rank-1 slope 2H − 1 (0.2 at H = 0.6); rank-2 slope 4H − 3 (0 below 0.75, 0.6 at 0.9);
cross slope 2H_X + 2H_Y − 3. Order-sign-like × RV-like pairing (H 0.7 × 0.9): cross slope **0.207**
(theory 0.2). That pair is **outside** the √n regime even though one series is below 3/4.

## 6. Answers for plan item (2)

1. Betken–Dehling prove the normal/non-normal switch for distance covariance, for subordinated
   Gaussian series, with the normal side requiring both H < 3/4 and the non-normal side only for
   Gaussian series. The mixed case is unproved.
2. Zhai–Zhong–Wu prove Gaussian approximation and block bootstrap for sample-covariance-based GGM
   tests in Gaussian linear processes with H < 3/4. Above 3/4 they say the limit is Rosenblatt-type
   and leave it.
3. HSIC with translation-invariant kernels: extension looks direct (worker inference, unproved).
   MMD two-sample and DML population-ATE: **not** the 3/4 boundary — rank 1 gives H = 1/2.
   Hawkes Granger: different mechanism.
4. For two series, the boundary is a **sum rule, H_X + H_Y = 3/2**, not a per-series line.
