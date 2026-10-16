# Synthesis — memory × kernels × causal inference (coordinator)

Job JOB-NEW-MKC-LIT-20261015, 2026-10-15. Written by the coordinator from the two closed
workstreams: `memory-causal/` (handoff `55313475…`) and `memory-kernel/` (handoff `16c970ae…`).
Bounded scans (74 and 70 records screened; SSRN and econometrics journals not searched directly).
No novelty claim. Every item below is a candidate for the author's choice, not a finding.

## 1. The common picture

Both workstreams found the same structure, from opposite sides.

- **Estimation mostly survives memory; inference mostly does not.** Kernel estimators stay
  consistent at slower rates; DML point estimates stay consistent. What fails is calibration:
  wild-bootstrap HSIC/MMD needs τ(r) = O(r⁻⁶) (Chwialkowski–Sejdinovic–Gretton, F); DML
  cross-fitting conditions are "never satisfied" for long-memory processes (Ballinari–Wehrli, F);
  single-path Hawkes Granger tests need a fixed margin from criticality (Wang–Kolar–Shojaie, F).
- **There is a sharp boundary at H = 3/4 (d = 1/4).** Below it, kernel-type statistics keep a
  normal limit; above it the limit is non-central with unknown law (Betken–Dehling 2107.03041;
  Zhai–Zhong–Wu 2609.30565). Applying this to HSIC is a worker inference, not a stated theorem.
- **Finance straddles the boundary.** Order signs sit just below it (H ≈ 0.7, Lillo–Farmer);
  realised volatility sits above it (H ≈ 0.9). So the same method can be valid for order-flow
  questions and invalid for volatility-channel questions. That is the clearest single fact for
  the author's interests.
- **What survives regardless:** randomised designs (Bojinov–Shephard), linear Granger tests
  robust to long memory (Bauer–Maynard 2012), and propagator estimation with martingale-difference
  noise (Hey–Neuman–Tuschmann 2025), because their validity does not come from mixing.
- **The word "kernel" means two things here**: RKHS kernels (tests, embeddings, signatures) and
  memory kernels (Hawkes excitation, price-impact propagators). Their meeting point — estimating a
  power-law memory kernel with RKHS methods near criticality — is itself open (D3).

## 2. Candidate directions, ranked against the author's interests

1. **A memory breakdown map (first step for everything else).** One no-data simulation family
   combining kernel D1 (size of HSIC / MMD / kernel-CI tests across d = 0…0.4 and power-law
   Hawkes), causal C1 (coverage of population vs path-averaged DML/AIPW effects when the effect
   varies with a long-memory state), and causal C2 (Hawkes Granger test size as ρ → 1, plus a
   shared slow baseline). Shared simulators (ARFIMA, power-law Hawkes). It decides which theory
   directions are real and is publishable on its own as a short "where inference breaks under
   market memory" note if the breaks line up with H = 3/4. Fit: high. Risk: low; scoop value modest.
   C1's mechanism is a worker derivation that may be folklore.
2. **Single-path causal discovery for rough/fractional systems with signature kernels (D4).**
   The only direction at the full triple intersection. Manten et al. (ICLR 2025, F) show a
   signature CI test keeps power under fBm but their discovery algorithm is "not applicable"
   without independent increments, and all guarantees need i.i.d. replicate paths. Finance has
   one path per asset. Cheapest check: rerun their fBm example on blocks of one long path.
   Fit: very high. Risk: the calibration problem may reduce to D2's theory.
3. **Core vs reaction order flow as a causal decomposition of impact (C3).** Best fit on memory,
   near-critical Hawkes and causal impact. Muhle-Karbe et al. 2601.23172 (F) note the two flows
   cannot be told apart in anonymous data. Cheapest check is paper-only: is the reaction kernel
   identified from aggregate signed-flow spectra? Highest data risk (may need account labels).
4. **RKHS estimation of power-law Hawkes / propagator kernels near criticality (D3).** Current RKHS
   Hawkes estimators truncate support and give no rates. Cheapest check: simulate near-critical
   power-law Hawkes and measure truncation bias. Not causal by itself.
5. **Proxy causal identification under long memory (C4, with D5 as its lemma).** Track J assumes
   β-mixing; one unverified 2026 paper (Hao, PLOS One) claims long-memory proximal identification.
   Cheapest check: read Hao's proofs beside CL-J1's mixing step; compute the Gaussian-kernel
   dependence decay on ARFIMA (D5) — if it decays like ρ(t)², kernel embedding results transfer to
   d < 1/4 almost for free.

## 3. Overlap with existing tracks

| Direction | Owner or overlap |
|---|---|
| Breakdown map: kernel tests (D1, D2) | Track E (CL-E4, G5) |
| Breakdown map: DML coverage (C1) | Tracks K, F-CAUSAL |
| Breakdown map: Hawkes Granger (C2) | Track B2-COMP; closed Track A ρ item |
| Signature single-path discovery (D4) | Partly E; otherwise none |
| Core vs reaction flow (C3) | Tracks B2-HTE, B2-QUEUE |
| RKHS power-law kernels (D3) | Closed Track A; Track G (Hawkes lemma) |
| Proxy under long memory (C4, D5) | Track J; F-KERNEL for kernel conformal (K9) |

Only D4 and the cross-cutting breakdown map are not already inside an existing track's scope.

## 4. Limits

SSRN, Google Scholar and econometrics/probability journals were not searched directly (main
reading-gap risk for D3, D5, C1). Proofs were not checked. C1's decomposition and the H = 3/4
mapping for HSIC are worker inferences. Several groups (Betken/Dehling, Shao, Sornette,
Muhle-Karbe, Bacry/Muzy) publish monthly in this area; pre-emption risk is real for every item.
