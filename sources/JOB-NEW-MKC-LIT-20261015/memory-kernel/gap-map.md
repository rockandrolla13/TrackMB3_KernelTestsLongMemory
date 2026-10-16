# Gap map — kernels with memory (workstream `memory-kernel`)

Job JOB-NEW-MKC-LIT-20261015 · execution `claude-mkc-kernel-20261015` · 2026-10-15.
Method: research-gap-finder template (coverage → routes → bottlenecks → controversies → gaps →
entry points → risks), applied to field-map.md plus the reused Track E / S2 / S8 / Track A reviews.

**This is a bounded scan.** "Open" means *not found in this scan, and a named obstacle explains why it
might be open*. It never means "nobody has done it". Each gap says whether it looks like a true gap,
an engineering composition, or a reading gap (likely exists where we could not search).

Status: **open** / **contested** / **closed**. Evidence class of the deciding sources: F / S / A / R
(see field-map.md). Track overlap: programme/catalog.yaml ids, or "none".

⚠ = touches the triple intersection (memory × kernels × causal). For the coordinator.

---

## Input coverage

- Topic: kernel methods (RKHS and memory kernels) under dependence beyond geometric mixing.
- Sources: ~70 records screened, 25 extracted, 12 full-text section reads (see search-log.md), plus
  Track E E1/E2/GAP_ANALYSIS (56 + 42 records), S2, S8, Track A Hawkes and P2-CALIB jobs (reused, R).
- Boundary: 2025-01-01 → 2026-10-15 emphasised; SSRN and econometrics journals thinly covered.

## Gaps

| ID | Gap | Status | Kind | Deciding sources (class) | Track overlap |
|---|---|---|---|---|---|
| K1 | Valid calibration of HSIC/MMD/KCI tests in the **slow-polynomial** tier: dependence summable but slower than the τ(r) = O(r^{−6}) / Σr²τ(r) < ∞ the wild-bootstrap theorems need — e.g. power-law Hawkes counts, propagator-filtered order flow | open | true gap (theory) + engineering (audit) | [chwialkowski2014wild] F; [betken2021ldr] F; E1 R | E (CL-E4 assumes polynomial β-mixing r > 2 — **overlaps**, E owns the short-memory end); A / G for the Hawkes mixing ingredient |
| K2 | Kernel independence / two-sample tests under **LRD** with feasible calibration for **general (non-distance) kernels and non-Gaussian processes**; no self-normalised kernel statistic for LRD found | open | true gap | [betken2021ldr] F (dCov, subordinated Gaussian, subsampling); [dehling2014twosample] F; [betken2016sn] A (SN works for Wilcoxon) | E (GAP_ANALYSIS G5, ranked low) |
| K3 ⚠ | Nonlinear (kernel) **conditional-independence / Granger** test under LRD | open | true gap | [zhai2026ggmlong] F (Gaussian linear only, β > 3/4); E1 R (no kernel CI under LRD) | E (G5); intersects memory-causal |
| K4 | **RKHS / nonparametric Hawkes-kernel estimation with power-law tails and statistical rates**, esp. near criticality (branching ratio → 1) | open | true gap | [bonnet2025rkhs] F and [kim2026representer] F (support truncated to [0, A], no rates); [bacry2016wienerhopf] S (power law ok, no rate theorem in abstract); [jaisson2016rough] S (why it matters: power-law + near-critical ⇒ rough volatility) | A (Hawkes check item 2: power-law kernels), G (CL-G1/G3 need the mixing lemma) — neither owns kernel *estimation* |
| K5 | Inference for **propagator (impact) kernels** when the kernel is power-law and order flow is long-memory: truncation bias of a finite lag grid; validity of the martingale-noise assumption | contested | engineering + reading | [hey2025crossimpact] F (confidence region under conditionally sub-Gaussian MDS noise, finite grid M; no mixing on flow needed); Track A P2-CALIB R (exponent ≈ 0.22, asymptote poorly identified) | A (Paper 2 calibration) |
| K6 | **Single-path** signature-kernel tests (two-sample, independence, CI) for one long rough / fractional path | open | true gap | [manten2025sigci] F (needs i.i.d. paths, App. A.13); [alden2025sigmmd] A (i.i.d. paths); [chevyrev2022sigmoments] A | E partly (CL-E4 is single-path but not signature); otherwise none |
| K7 ⚠ | **Causal discovery for systems driven by fractional / rough noise** using signature-kernel CI tests | open | true gap | [manten2025sigci] F: test keeps power under fBm, but "Algorithm 1 is not applicable anymore" without independent increments | none |
| K8 | Rates and **inference** for kernel-embedding / MMD estimators under LRD; consistency is already implied by Σ_T = o(T) | open (rates), closed (consistency) | engineering → maybe theory | [alquier2026mmd] F; [cherief2022mmd] S; [singer2017kpls] A; [wang2024krrdependent] A | none directly; F-KERNEL and J use the same ingredient |
| K9 | **Kernel conformity scores** with coverage guarantees under long memory | open | true gap (possibly reading gap: Wu/Politis groups) | [lee2025kowcpi] F (α(τ) = O(τ^{−(2+δ)})); S2 R (#10–12: long memory handled only for residual/CQR scores) | **F-KERNEL owns** |
| K10 ⚠ | Kernel causal estimators (proxy, IV, counterfactual embeddings) on **one dependent, long-memory series** | open | true gap at LRD; short-memory version is Track J's | E2 R; [mastouri2021proxy] A; [lob2025kivuniform] S (dependence not stated) | **J owns** (CL-J1, β-mixing); LRD extension: none |
| K11 | Kernel **Koopman / conditional-embedding operator learning for non-Markovian (long-memory) dynamics** | open | likely reading gap (dynamical-systems venues not searched in depth) | [kostic2022koopman] S (Markov + mixing); [klus2020kto] A | none |
| K12 ⚠ | Signature-based identification of **excitation direction** for **power-law** Hawkes from one path | open | true gap; single unrefereed source | [noguer2026paths] F (exponential/affine kernels only; Monte Carlo over 50,000 independent paths) | none (adjacent to A/G Hawkes work and to memory-causal Hawkes Granger) |

### Closed (do not reopen)

| ID | What | Why closed | Source |
|---|---|---|---|
| X1 | Marginal HSIC/MMD/dCov tests under short-memory mixing | ≥5 published solutions | E1 R; [chwialkowski2014wild] F |
| X2 | Representer theorems for RKHS Hawkes kernels (bounded support) | two 2025–26 papers | [bonnet2025rkhs] F; [kim2026representer] F |
| X3 | Exponential-kernel Hawkes mixing → coverage gap | routine corollary per Track A | Track A HAWKES-CHECK R |
| X4 | Nonparametric multi-asset propagator estimation with confidence bounds (exponential-scale question) | done, with power-law finding | [hey2025crossimpact] F |
| X5 | Consistency of MMD estimation when Σ_T = o(T) | follows from stated bound | [alquier2026mmd] F |

### Contested

- **C1. Is financial data inside the H < 3/4 regime where kernel statistics stay √n-Gaussian?**
  Betken & Dehling (Remark 2.1) say real data typically has D > 1/2 (H < 0.75). Against: order-sign
  memory ≈ H 0.7 [lillo2004longmemory] S sits at the edge; realised volatility d ≈ 0.4 (H ≈ 0.9) in
  [andersen2003rv] is outside; propagator exponent ≈ 0.22 (Track A R) implies strong persistence of
  impact-filtered flow. Answer depends on which series is tested. This decides whether K2 is a
  "scope note" (Track E's view, G5) or the main obstacle for volatility-channel kernel tests.
- **C2. Is K5 a gap at all?** The Hey et al. bound needs no mixing of order flow, which is good news.
  But it needs martingale-difference noise and truncates the kernel; whether those bite on real data is
  untested. Could be routine.

## Bottlenecks (why these gaps persist)

- Theory: degenerate (quadratic) statistics under LRD have non-central limits with unknown laws above
  H = 3/4; subsampling is the only proven calibration, and it needs block lengths
  l_n = O(n^{(1+min D)/2 − ε}) [betken2021ldr] F — slow, wasteful at finance sample sizes.
- Theory: power-law memory kernels make Hawkes non-Markov; mixing tools need a Markov embedding
  ([lee2026powerlawmarkov] S offers one surrogate).
- Data: signature-kernel tests and kernel causal estimators assume replicate paths / units; markets give
  one path.
- Evaluation: no benchmark of kernel-test size across memory levels was found.

## Small entry points (cheapest first)

1. Size audit across a memory grid (no data). Decides K1/K2/C1 at once. → directions D1.
2. Truncation-bias study for RKHS Hawkes estimators on simulated near-critical power-law Hawkes. → D3.
3. Blocked single-path version of Manten et al.'s fBm experiment. → D4.

## Risks

- Reading-gap risk is highest for K9 and K11 (Wu/Politis and dynamical-systems groups), and for K4
  (point-process statistics venues, e.g. Bernoulli/SPA, not searched directly).
- Scoop risk: K1/K2 sit next to Betken/Dehling and Gretton/Sejdinovic groups; K4 next to
  Bonnet/Sangnier and Kim/Iwata (both posted in the last 12 months).
- Several bridge statements (the H = 3/4 map for HSIC; T^{−D/2} MMD rate) are worker inference.
