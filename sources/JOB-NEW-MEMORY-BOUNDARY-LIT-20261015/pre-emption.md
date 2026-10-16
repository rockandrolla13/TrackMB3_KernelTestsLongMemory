# Pre-emption check — has anyone already mapped this?

Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 · plan items (1) and (3) · 2026-10-15.
**Bounded scan.** Not finding something does not show it is absent. No novelty claim.

## Stopping rule

The plan says to stop early if a source already maps size/coverage across H for these methods.
**No such source was found, so the rule did not trigger.** The closest sources each cover one
method at two or three memory levels:

| Source | Method | Memory grid | What it maps | Class |
|---|---|---|---|---|
| Betken & Dehling 2107.03041 | dCov independence test, subsampling | H ∈ {0.6, 0.7} | size 7–17.5% at 5%, n ≤ 1000; never above 3/4 | F |
| Zhai, Zhong & Wu 2609.30565 | high-dim GGM test, block bootstrap | β ∈ {2, 0.9, 0.55} (H ≈ short, 0.6, 0.95) | QQ/Kolmogorov distances; fails at H = 0.95 | F |
| Betken, Dehling & Kroll 2112.14091 | dCov bootstrap independence test | α-mixing only | size under short memory | S (abstract) |
| Bauer & Maynard, J. Econometrics 2012 | surplus-lag linear Granger test | stationary, unit root, local-to-unity, long memory, breaks | one χ² null across persistence types | A (+R) |
| Lazarus, Lewis, Stock & Watson 2018; Lazarus, Lewis & Stock, Econometrica 2021 | HAR/HAC t-tests, fixed-b | persistent AR, not fractional | size–power frontier | A |
| Casini & Perron, arXiv:2103.02981 | HAC under strong dependence | "extreme" dependence, not fractional | size distortion | A |
| Achab et al., arXiv:1607.06333 | Hawkes Granger graph via integrated cumulants | rectangular, exponential, power-law kernels | estimation robustness, not test size | A |

## (1) Gap closure: SSRN, econometrics and probability journals

What was searched (see search-log.md for strings and counts):
- SSRN: one domain-restricted query. **The search tool returned arXiv records, not SSRN records**,
  so SSRN was not effectively searched. Google Scholar is not reachable with the available tools.
  Both remain reading-gap risks.
- J. Econometrics / JTSA / Econometric Theory / Bernoulli / AoS via publisher domains
  (sciencedirect, wiley, cambridge, projecteuclid) and open queries. Hits: Bauer–Maynard 2012,
  Shao 2011 (JTSA, change in mean under LRD), Robinson–Hidalgo 1997 (AoS), Tsay 2007 (JTSA),
  Robinson–Velasco 2015, spurious-long-memory tests. None maps kernel-test, DML or Hawkes-Granger
  size across H.

Per method:
- **Kernel tests under LRD.** Only Betken–Dehling (dCov) and Dehling–Rooch–Wendler (two-sample
  U-statistics, theory) found. No HSIC/MMD/KCI size study across H. Contrary evidence for
  "safe below 3/4": Betken–Dehling's own test is oversized at H = 0.6–0.7 for n ≤ 1000.
- **DML / AIPW intervals under LRD.** None found. The nearest is Ciganovic–D'Amario–Tancioni,
  "Double Machine Learning for Time Series" (arXiv:2603.10999; Econometrics Journal 2026 per arXiv
  page, A/S). It uses mixing or Maxwell–Woodroofe conditions, which (worker inference) exclude
  d > 0. Ballinari–Wehrli (R) remains the only source stating that cross-fitting conditions fail
  for persistent processes. Simulation studies of "treatment persistence" (arXiv:2609.04940, A)
  concern persistent *treatment* in panels, not long-memory states.
- **Hawkes Granger tests near criticality / power-law kernels.** None found with size results.
  A Politecnico di Milano thesis (A, snippet) reports "spurious cross-excitation in controlled
  synthetic experiments" for parametric multivariate exponential and power-law fits. That is
  direct qualitative support for C2(ii), from an unrefereed source. Nearly-unstable limit theory
  (Jaisson–Rosenbaum 2015; supercritical case arXiv:2407.19741) gives the probability side but not
  test size.

## (3) Is the C1 decomposition already stated?

C1: memory enters a DML/AIPW score only through effect heterogeneity in a long-memory state; the
residual term is a martingale difference; so population-ATE intervals fail and path-conditional
ones do not.

Found, each partial:
1. **Path-conditional estimands** are established. Bojinov & Shephard (arXiv:1706.07840, R) and
   Rambachan & Shephard (arXiv:1903.01637, A) define effects conditional on the realised path. The
   reason is identification, not memory.
2. **The product-memory rule** in regression is classical (Robinson & Hidalgo 1997, A). Least
   squares loses √n when error and regressor singularities coincide. That is the cross-term
   half of C1's logic.
3. **Sample vs population ATE** in i.i.d. data is textbook (not searched here; noted as a
   pre-emption risk in the prior job).
4. **Hermite-rank reasoning** (boundary-theorems.md §3) makes the C1 rate n^{d−1/2} an immediate
   corollary of Taqqu / Dobrushin–Major for a rank-1 effect function.

Not found: a statement of the decomposition for DML/AIPW scores on long-memory series, or a
coverage study of population vs path-conditional effects under fractional states.

**Assessment.** C1 is very likely "folklore-adjacent". Each ingredient is classical, and a
referee could call the combination a corollary. Its value would be as a clearly stated and
simulated warning for DML users, not as a theorem. The scan also **corrects C1**: the relevant
boundary is d = 0 for monotone effect heterogeneity (rank 1), not d = 1/4. Only even (rank-2)
heterogeneity gives the 3/4 line.

## Contrary evidence kept

- Betken–Dehling: oversized below 3/4 at realistic n (Table 1). The finite-sample problem is not
  confined to H > 3/4.
- Betken–Dehling Remark 2.1(ii): real data "typically" has H < 0.75. Not supported by the finance
  table (finance-memory-estimates.md): most activity series sit above 0.75.
- Gatheral–Jaisson–Rosenbaum: volatility may be rough, not long memory. If so, the "volatility is
  above 3/4" premise becomes a finite-sample mimicry question, not an asymptotic one.
- Bauer–Maynard: a linear Granger test with one χ² null across persistence types exists. For
  linear questions the breakdown map may be moot.

## Scoop watch

Betken/Dehling/Kroll (dCov, bootstrap, LRD), Zhai/Zhong/Wu/Kolar (LRD high-dim inference, Sep 2026),
and the Rosenbaum and Muhle-Karbe groups. ZZW's Section 6 names ultra-long memory (H > 3/4) as
their next step. An HSIC or GGM result above 3/4 from that group is plausible within a year.
