# Three-way scan (WHY / HOW / WHAT) — closest papers
Track NEW-MB-KERNEL · Agent 2 · execution claude-mb3-lit-20261016 · 2026-10-16. Class F unless stated.

1. Betken–Dehling 2107.03041v3. WHY: no dCov theory for long memory. HOW: Hilbert-space (non-)central limit theorem for the empirical characteristic-function process, plus subsampling. WHAT: √n limit when both H < 3/4; for Gaussian data above the line, a squared scalar limit; mixed case uncovered; size 7–17.5% at H = 0.6–0.7, n ≤ 1000.
2. ZZW 2609.30565. WHY: graphical-model CI testing for long-memory high-dimensional series. HOW: Gaussian approximation plus overlapping block bootstrap. WHAT: valid for β > 3/4 only; the Rosenblatt case is "beyond scope"; fails at β = 0.55.
3. ZKW IEEE TIT 2026. WHY: finite-sample covariance and precision inference under long memory. HOW: martingale and m-dependent approximation plus block sampling. WHAT: β > 3/4 only; states the Rosenblatt transition.
4. Zhai–Cheng–Wu 2609.33868. WHY: inference on realised conformal coverage. HOW: functional dependence measure, Bahadur representation, block sampling with estimated normalisation. WHAT: Gaussian, then log-normalised, then Rosenblatt at H = 3/4 for |W| scores; block sampling valid in all three; MBB and HAC oversized.
5. Ghoshal–Shao (EJS 2026). WHY: avoid bandwidth-dependent bootstrap. HOW: sample split, RKHS projection, self-normalisation, Brownian FCLT. WHAT: pivotal limits under α-mixing; no long memory.
6. Bai–Taqqu–Zhang (SPA 2016). WHY: SN alone still depends on H. HOW: SN block sums plus block sampling. WHAT: valid across short and long memory and heavy tails, including Rosenblatt limits, with no Ĥ.
7. Betken–Wendler (Stat. Sinica 2018). WHY: long-memory limits have unknown parameters. HOW: subsampling for general statistics, applied to an SN Wilcoxon statistic. WHAT: H-free change-point test; this is the subsampling theorem Betken–Dehling reuse.
8. Shao (Econometric Theory 2009). WHY: Haugh and Hong tests exclude long memory. HOW: frequency-domain prewhitening with an estimated spectrum. WHAT: N(0,1) null under long memory.
9. Chwialkowski–Sejdinovic–Gretton (NeurIPS 2014). WHY: permutation fails on dependent data. HOW: Markov-multiplier wild bootstrap under τ-mixing. WHAT: excludes every H > 1/2; centred multipliers over-reject at moderate n.

## Synthesis
- Common WHY: get valid critical values without knowing the strength of dependence.
- Divergent HOW: three routes. (a) Prove the limit, then subsample; this breaks where the normalisation changes. (b) Gaussian approximation plus blocks, now extended across 3/4 by estimating the normalisation. (c) Make the statistic pivotal, through SN, SN plus blocks, or prewhitening.
- Strongest WHAT: the lines are proved for covariance, indicator and Gaussian-dCov statistics. The above-line dCov law is a squared scalar. H-free calibration exists for scalar statistics.
- Unresolved gap: nobody joins the kernel route (HSIC/MMD limits, kernel SN) to the long-memory route (Hermite limits, H-free calibration). That leaves HSIC with non-linear G, the mixed case, kernel CI under long memory, and a size map across H.
