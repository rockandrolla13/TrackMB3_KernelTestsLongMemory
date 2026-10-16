# Candidate directions — kernels with memory (workstream `memory-kernel`)

Job JOB-NEW-MKC-LIT-20261015 · execution `claude-mkc-kernel-20261015` · 2026-10-15.

Candidates, not recommendations. The author decides. No novelty claim: each rests on a bounded scan.
"Fit" = fit to the author's stated interests (memory; both meanings of kernel; causal; finance).

---

## D1. Where do kernel tests break? A size map across memory levels

- **Question.** At what memory level do HSIC / MMD / kernel-CI tests lose size control under each
  calibration (permutation, wild bootstrap, block bootstrap, subsampling, self-normalisation)?
- **Assumptions.** None to prove; a simulation design. Grid: ARFIMA d ∈ {0, 0.1, 0.2, 0.25, 0.3, 0.4}
  (H 0.5–0.9), power-law Hawkes counts at several tail exponents and branching ratios, an order-sign-like
  series (exponent ≈ 0.6). Independent pairs for the null.
- **Closest prior work.** [chwialkowski2014wild] (wild bootstrap; τ(r) = O(r^{−6}) needed);
  [betken2021ldr] (dCov + subsampling under LRD); Shen et al. block permutation (E1 R); Ghoshal & Shao
  self-normalised RKHS (E1 R). None reports a size map across H.
- **Cheapest decisive check.** One simulation script, no data: n ∈ {500, 2000}, 1000 reps per cell,
  record empirical size at α = 0.05. Decisive if wild-bootstrap HSIC stays within ±0.02 of nominal for
  d < 0.25 and breaks above — that confirms the H = 3/4 boundary for kernel tests. If it holds everywhere,
  K2 is a scope note and C1 is moot.
- **Fit.** High on memory × kernels. Feeds Track E (CL-E4 motivation, N3 audit) and settles contested
  point C1. Low theory risk, low scoop value on its own; a JTSA/EJS-style note or a section of a bigger paper.
- **Overlap.** E (N3 size audit of i.i.d. KCI; G5). Not a new track.

## D2. A self-normalised kernel independence / CI statistic valid under long memory

- **Question.** Can a self-normalised RKHS statistic (Ghoshal–Shao style) have a pivotal limit that does
  not depend on H, at least for H < 3/4, and a stated non-central limit above it?
- **Assumptions.** Stationary subordinated-Gaussian or linear processes with memory D; bounded
  characteristic kernel; recursive-estimate self-normaliser.
- **Closest prior work.** Ghoshal & Shao (EJS 2026, α-mixing; E1 R); Shao's self-normalised CUSUM under
  LRD and [betken2016sn] (Wilcoxon, LRD); [betken2021ldr] (non-central Hilbert-space limit);
  [dehling2014twosample] (Hermite route for two-sample U-processes).
- **Cheapest decisive check.** Apply the Ghoshal–Shao statistic (code or re-implementation) inside the
  D1 grid. If size holds for d < 0.25 without tuning, the theory route is plausible; if it fails even there,
  drop D2.
- **Fit.** High (memory × RKHS). Theory-heavy; scoop risk from the Shao and Betken/Dehling groups.
- **Overlap.** E (G5). If pursued, it would extend CL-E4 beyond β-mixing rather than create a track.

## D3. RKHS estimation of power-law memory kernels near criticality (Hawkes and propagators)

- **Question.** Can a heavy-tail-weighted RKHS (e.g. a weighted Sobolev space whose norm penalises
  t^{a}·g(t)) estimate a power-law Hawkes kernel or propagator without truncating its tail, with a rate
  that degrades in a stated way as the branching ratio → 1?
- **Assumptions.** Linear Hawkes, kernel g(t) ~ t^{−1−α}, branching ratio ρ < 1; stationarity; or a
  propagator model with martingale-difference noise ([hey2025crossimpact] Assumption 2.3).
- **Closest prior work.** [bonnet2025rkhs], [kim2026representer] (RKHS, bounded support, no rates);
  [bacry2016wienerhopf] (Wiener–Hopf, power law, no rate theorem in abstract); [hey2025crossimpact]
  (propagators, confidence region, finite grid); [lang2024gle] (same inverse-problem structure);
  [jaisson2016rough] (why near-critical power law matters: rough volatility).
- **Cheapest decisive check.** Simulate near-critical power-law Hawkes (ρ ∈ {0.7, 0.9, 0.97},
  α ∈ {0.3, 0.6}); fit Kim–Iwata / Bonnet–Sangnier with support A ∈ {5, 50, 500}; measure bias in the
  estimated branching ratio and tail. If bias is small at feasible A, K4 is engineering only; if it is large
  and grows with ρ, the weighted-RKHS theory is worth it.
- **Fit.** Highest on "where RKHS kernels meet memory kernels" and on microstructure data. Links to
  rough volatility. Not causal by itself.
- **Overlap.** A (Hawkes check item 2: power-law kernels; Paper 2 calibration uses the propagator);
  G (CL-G1/G3 mixing lemma). Ownership of the Hawkes mixing lemma is already an open author decision
  in Track A's handoff.

## D4. ⚠ Single-path signature-kernel CI tests and causal discovery under fractional noise

- **Question.** Can a signature-kernel CI test be calibrated from blocks of **one** long path of a
  rough/fractional system, and what replaces the independent-increments Markov property that Manten
  et al.'s discovery algorithm needs?
- **Assumptions.** SDEs driven by fBm with H in a known range (> 1/4 for the rough lift); stationarity of
  path segments; a block-calibration scheme whose validity tracks the H = 3/4 boundary.
- **Closest prior work.** [manten2025sigci] (i.i.d. paths; fBm test works, discovery fails);
  [chevyrev2022sigmoments], [salvi2021goursat]; single-path kernel CI tests in E1 (Tian et al., Shi &
  Shah — tier-1 mixing only); [noguer2026paths] (excitation direction via signatures, exponential Hawkes).
- **Cheapest decisive check.** Re-run Manten et al.'s fBm example (App., eq. 22) with one long path cut
  into blocks instead of independent copies; measure null size and power across H. If size holds for
  H ≤ 0.7 with simple blocking, a method paper is plausible; if not, the obstacle is the calibration
  theory (D2's), not the signature.
- **Fit.** The only candidate that touches all three of memory, kernels and causal. Rough-volatility
  finance application is natural (volatility paths are rough; one path per asset).
- **Overlap.** Partly E (single-path CI), otherwise none. Highest intersection value for the coordinator.

## D5. Long-memory rates for kernel embeddings, feeding kernel conformal scores and kernel proxy methods

- **Question.** Under LRD with kernel dependence coefficients ϱ_t ~ t^{−D}, what rate does MMD /
  conditional-mean-embedding estimation achieve, and what does that imply for kernel conformity scores
  (F-KERNEL) and RKHS proxy bridges (J)?
- **Assumptions.** Stationary process; Σ_T = Σ ϱ_s with ϱ_s ≤ β(s) for radial kernels ([alquier2026mmd]
  Prop 4.1); LRD makes β-mixing fail, so ϱ_s must be bounded directly (e.g. via Hermite rank for
  Gaussian-subordinated data).
- **Closest prior work.** [alquier2026mmd], [cherief2022mmd]; [singer2017kpls] (kernel PLS slows under
  LRD); [lee2025kowcpi] (kernel conformal, α-mixing); S2 R (Zhai–Cheng–Wu long-memory coverage,
  non-kernel scores).
- **Cheapest decisive check.** Compute ϱ_t numerically for a Gaussian kernel on ARFIMA(0,d,0) and check
  whether ϱ_t decays like ρ(t)² (Hermite rank 2) or ρ(t). If ρ(t)², the MMD rate stays √T for d < 1/4 —
  meaning many kernel-embedding results transfer to moderate long memory almost for free.
- **Fit.** Medium: supporting lemma rather than a paper; useful to F-KERNEL and J.
- **Overlap.** F-KERNEL (conformal), J (proxy). Not a new track.

---

## Ranking against the author's interests (worker view, for the coordinator)

1. **D4** — only one at the triple intersection; cheap first check; one-path problem is the real
   obstacle for finance.
2. **D3** — best "RKHS meets memory kernel" fit and the strongest microstructure link; overlaps Track A.
3. **D1** — cheapest and decides contested point C1 for everything else; should probably run first
   whichever direction is chosen.
4. **D2** — theory payoff if D1 shows the H = 3/4 break.
5. **D5** — a lemma that serves F-KERNEL and J.
