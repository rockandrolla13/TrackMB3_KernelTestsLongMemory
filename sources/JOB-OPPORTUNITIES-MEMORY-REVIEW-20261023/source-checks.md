# Source checks — JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023

Read depth: **F** = full text inspected by this reviewer. **P** = partial full text (range stated). **A** = abstract or metadata only. **X** = inaccessible.
"Disc." is discovery's read depth for the same source.

| Source | Depth | Disc. | Locator checked | Finding | Bears on |
|---|---|---|---|---|---|
| Betken & Dehling, arXiv:2107.03041v3 (27 Jan 2026), PDF sha256 `5bc4b26c…52b50df` | F | F (first 100k chars) | Remark 2.1(i)–(ii); Theorem 2.4; Assumptions 1–2; Section 4; Section 5 normalisation | The Remark 2.1(i) quote is accurate. Subsampling Theorem 2.4 needs l_n = O(n^{(1+min D)/2−ε}). The applied √n normalisation is chosen "Under a corresponding assumption" D ≈ 0.6, so it is not H-free. Simulations only at H = 0.6, 0.7. No journal reference. | G1, G2 |
| Düker, arXiv:1704.08609v4 (11 Feb 2020), PDF sha256 `5d6cb544…a6d8890` | F (abstract, §1, §2, §4) | not cited | Theorem 4.1; Remarks 4.2–4.4; Proposition 2.1 | Joint limits of sample cross-covariances, sorted by d_k + d_l: Brownian below 1/2, matrix Rosenblatt in (1/2, 1). Covers the mixed region for linear processes. Linear processes only. No journal reference on arXiv. | G1, G2, MAP |
| Dehling, Rooch & Wendler, arXiv:1404.0551v1, PDF sha256 `16033fd0…e304aa33c` | F (§1–2) | not cited | Introduction, Section 2 | Two-sample U-process from segments of one LRD process. Cites Lévy-Leduc–Taqqu for "two-dimensional Hermite rank". | G1 |
| Bai & Taqqu, arXiv:1211.0576 (JTSA 2013) | P (first 100k of 115k) | not cited | Setting, Theorem 3.6, Conjecture 3.7 | Functions of one common Gaussian series. Mixed SRD/LRD case proved only for LRD ranks 1 or 2. | G1 |
| Lévy-Leduc & Taqqu, "Hermite ranks and U-statistics", Metrika 77 (2014) 105–136 | A (HAL API) / X (full text) | A, authors unknown | Abstract | U-statistics Σ h(X_i, X_j) of one stationary Gaussian sequence. Compares Hermite-rank notions. | G1 |
| Müller & Watson, "Low-Frequency Analysis of Economic Time Series", Handbook draft May 2024, PDF sha256 `5bf91170…f560eedcb` | F (§8.2) | not cited | §8.2, Theorem 2 | "The inference results of this chapter crucially depend on the large sample Gaussianity of the suitably scaled cosine transform". The CLT covers linear processes with m.d.s. innovations. | G2 |
| Preinerstorfer & Pötscher, arXiv:1304.1383v2 (ET 2016), PDF sha256 `8f1231a4…f8997` | F (abstract, §3, §3.2.2) | A | Assumption 1; Theorem 3.3; Remark 3.5(i); Pötscher 2002 discussion | Negative results "rest on Assumption 1, i.e., C ⊇ C_AR(1)", weakened to correlation matrices converging to e₊e₊′. Pötscher 2002 "only discusses ill-posedness of the problem of estimating f(0)". | UNIFORM |
| Pötscher 2002, Econometrica 70(3) | X | A | — | The working-paper URL returns HTML. CORE redirects. | UNIFORM |
| Kulik & Soulier, arXiv:1109.5298 (AAP 2012) | P (first 100k of 150k) | P (same range) | §2 eq. (4); §4 τ_p; §5.1 τ*_p(s); Examples 5.2–5.3 | §5.1 τ* is the rank of the sample-covariance kernel K*, not of |Y|^p. §4 τ_p = rank(σ^p). In LMSV the rank "depends only on the function σ". σ = x² gives τ* = 2. Leverage gives rank 1 unless E[Z₀η₀] = 0. | G3 |
| Jaisson & Rosenbaum, arXiv:1504.03100 (AAP 2016) | A | not cited | Abstract | Kernel tail x^−(1+α), α ∈ (1/2, 1), near criticality gives a limit with H = α − 1/2. | HAWKES |
| Ernst, Viens & Yan, arXiv:2510.23563 (Oct 2025) | A | known (TrackMB3) | Abstract | Only citing paper of Betken–Dehling in Semantic Scholar. Non-stationary BM/fBm paths, no HSIC. | G1, MAP |

## Quotes from discovery that I re-checked

| Discovery quote | Check |
|---|---|
| BD Remark 2.1(i) "does not cover the case Dξ+Dη>1, but requires either Dξ<1/2 or Dη<1/2" | Matches the PDF |
| BD Theorem 2.2(ii) only for G1 = G2 = id | Matches ("If G_1 = G_2 = id and D_ξ + D_η < 1") |
| K–S §5.1 "for σ(x)=exp(x), τ*_p(s)=1 for all s" | Matches, but τ* is the covariance-kernel rank. Discovery's gloss "every power |Y|^p has Hermite rank 1" mislabels it. |
| K–S Ex. 5.2 rank of K is 1 "except if E[Z_0 η_0]=0" | Matches. The example uses σ(x) = x². |
| Pötscher 2002 implies limits on calibrations that estimate H | Not supported. Preinerstorfer–Pötscher state the 2002 paper is about estimating f(0) only. |
