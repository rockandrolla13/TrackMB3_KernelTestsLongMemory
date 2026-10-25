# Theory notes: Propositions 1 (seasonal mean) and 2 (half-mean spread)

Files in this folder:

- `propositions.tex`: statements (for a new Section 2.4) and proofs (for a new appendix). Preamble lines are listed at its top.
- `refs-additions.bib`: new bib entries. Every entry needs checking.
- `check_propositions.py`: numeric check. It runs single-threaded in about 3.5 minutes. Output is in `check_output.txt`.

The file compiles as a standalone test with the paper's preamble plus the lines listed at its top.

## What is proved

**Proposition 1 (SEAS-M), for 0 < H < 1, P ≥ 3, a > 0.**

1. (i) The sample correlation r_n tends to c(D) **almost surely**, where c(d) = (a²/2)cos(2πd/P)/(1+a²/2). The limit has mean zero. P(c(D)=0) is 2/P when 4 divides P, and 0 otherwise. For P = 50 it is 0.
2. (ii) Any test that rejects when |Ĉ_n| > z_n (v̂_n/n)^{1/2}, with z_n bounded and v̂_n/n → 0 in probability, has liminf size ≥ P(c(D) ≠ 0). That is 1 for P = 50. The proof places both implemented Pearson rules (`kmem/legacy/mbcore.py:lind_tests`) in this class without extra assumptions:
   - NW: v̂ ≤ (2S+1)γ̂(0), and the Andrews plug-in S_n is O(n^{1/3}) a.s. (the proof shows the limiting AR(1) coefficient has modulus below 1).
   - EWC: v̂ ≤ nγ̂(0)/ν by DCT orthonormality, with ν = ⌊0.4 n^{2/3}⌋ → ∞.
3. (iii) Gaussian-kernel HSIC with the paper's median bandwidths: V_n → h(D) a.s. and min_d h(d) > 0 for P ≥ 5. Any test with critical value κ_n = o_p(n) has size → 1. The proof shows deterministically that the subsampling critical value is at most 4ℓ_n. It also shows the wild-multiplier mean E_w T* ≤ 4(1+2ℓ_n), which bounds the rejection probability. Both hold as implemented, centred or raw multipliers.

**Proposition 2 (split mean shift), for 0 < H < 1.**

1. (i) Var(Δ_m) = (4 − 2^{2H}) m^{2H−2}. This is **exact** for every m, not asymptotic. So the "~" in the brief is actually "=".
2. (ii) The power of the oracle linear half-mean z-test is a function of R_m alone.
3. (iii) This part is new and stronger than the brief. Every test that is invariant to adding a constant and has level α has power at most 1 − Φ(z_{1−α} − δ I_n^{1/2}), where I_n is the GLS information. The proof uses the Neyman–Pearson lemma on the maximal invariant. Gaussian-kernel split MMD, its median bandwidth, and all three critical-value rules are location-invariant. The size-adjusted oracle is too. So the bound covers every MMD row of Table 4 rigorously. It does not rely on the Hermite heuristic.

## Where the statement had to be weakened

- Prop 1(ii) gives **liminf size ≥ P(c ≠ 0)**, not "size → P(c ≠ 0)". On {c(D) = 0} nothing is claimed. This does not matter for P = 50, where the bound is 1.
- Prop 1(iii) needs **P ≥ 5**, so that the sinusoid takes at least 3 distinct values. It needs no condition on cos. HSIC sees the dependence even when the correlation limit is 0.
- Prop 2 rate: we did not prove that I_n is of order m^{2−2H}, the GLS analogue of Adenstedt (1974) for the mean. The paper should quote the computed values only. Those are δ I_n^{1/2} = 2.24, 0.71, and 0.46 at H = 0.5, 0.7, 0.8.

## What is heuristic or sketched (labelled as such in the .tex)

- **Circular shifts (Remark).** The leading-order argument is a sketch. When P divides n, a shift gives an exact SEAS-M pair with gap D − s. The leading p-value term is then super-uniform because D is uniform. Two gaps remain:
  - We have not proved uniform convergence of T_n(s)/n over about 0.8n shifts.
  - h has ties, for example h(d) = h(P/2 − d). Within a tie the decision is set by lower-order long-memory noise.

  So the remark explains "near α, not near 1". It does not prove "= α".
- **Distance covariance (Remark).** This is stated without proof. It needs uniform integrability in the continuity step and a window-moment bound for subsampling.
- **Hermite rank-one link for Gaussian-kernel MMD.** The signal direction equals the rank-one coefficient a_1 = E[ξ φ(ξ)], by Stein's identity. This is exact. The claim that higher-rank noise explains the H = 0.5 gap between kernel and linear tests is a heuristic, and the text says so.

## Numeric check (check_output.txt)

- r_n vs c(d): at n = 20000 the mean of r_n matches c(d) within 0.006 for every gap and both H. The spread shrinks from n = 2000 to n = 20000. At H = 0.8 it shrinks more slowly, as long memory predicts.
- Pearson by gap, H = 0.5. NW rejects 91.8% and EWC 92.4% at n = 2000, against 91.4% in EXP2. The misses are only at gaps 11–14 and 36–39. At n = 20000 the rates are 98.0% and 99.0%, and only gaps 13 and 37 still miss.
- HSIC limit h(d), closed form for the Gaussian mixture: all 50 values are positive. The minimum is 3.9e-5 at d = 12, 13 (and 37, 38). At n = 2000, n·h(12) = 0.08, the same size as a no-season null statistic. So at those four gaps the asymptotics have not started at n = 2000. The empirical V_n at d = 0 (0.0055–0.0063) agrees with h(0) = 0.0059.
- Prop 2 variance, 10,000 Monte Carlo paths: the formula and Monte Carlo agree within 1.9 standard errors at H = 0.5, 0.7, 0.8. sd = 0.0447 and 0.2472 match the paper's 0.045 and 0.25.
- Ratios at n = 2000, δ = 0.1: **2.236 at H = 0.5 and 0.405 at H = 0.8** (the brief's ≈2.2 and ≈0.40 are confirmed). At H = 0.7 the ratio is 0.681.
- Linear half-mean test power: 60.9%, 10.5%, 6.9%. The invariant-test cap is 72.3%, 17.6%, 11.7%. Observed oracle MMD: 45.7% at H = 0.5 and 6.4% at H = 0.8.

## Citations to verify ([verify citation] in the .tex)

1. Maruyama (1949): Gaussian stationary sequence with covariance → 0 is mixing. Check the theorem number.
2. Cornfeld, Fomin & Sinai (1982), Ch. 14 (Gaussian dynamical systems): the vector-valued version.
3. Walters (1982): the product of mixing systems is mixing. The cited "Section 1.7" is a guess. Find the exact theorem.
4. Billingsley (1999), Thm 2.8: weak convergence of product measures.
5. Gretton et al. (2005, ALT): biased HSIC = n^{-2} tr(KHLH), and HSIC = MMD² to the product of marginals.
6. Sriperumbudur et al. (2010, JMLR): the Gaussian kernel on R^d is characteristic.
7. Lehmann & Romano (2005), Thm 6.2.1: invariant tests are functions of a maximal invariant. Check the number in the 3rd edition.
8. Andrews (1991, Econometrica): the AR(1) plug-in bandwidth. It is not in refs.bib now.

## Suggested placement in paper.tex

- **New §2.4 "Two exact results"**: after the "Independence statistics" subsection, at the end of its last paragraph ("…Table~\ref{tab:scope} separates…"), and before Table 1. Add two rows to Table 1: "Prop. 1: SEAS-M limit, any o(n)-critical-value rule" and "Prop. 2: invariant-test power cap for a split mean shift".
- **§4.3 (`sec:exp2`), first finding paragraph** ("First, a seasonal mean misleads…"): after "does not shrink as n grows", cite Proposition 1(i)–(iii). Add the 4/50 → 92% prediction. After "Circular shifts stay at 5.8%", cite Remark (shift).
- **§4.4 (`sec:power`), "Third, memory removes almost all the power…"**: replace the two sentences about 0.045/0.25 with a citation of Proposition 2(i). Add the cap: no location-invariant 5% test exceeds 12% at H = 0.8.
- **Proofs**: new appendix before "Reproducibility and audit trail".
