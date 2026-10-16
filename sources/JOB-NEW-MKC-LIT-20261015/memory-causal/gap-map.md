# Gap map — memory-causal workstream

Execution `claude-mkc-causal-20261015`, 2026-10-15. Template: research-gap-finder
`gap_analysis_template.md` (sha256 a5c625dd…). Labels as in `field-map.md`: F, S, A, D, F(reused).

**This is a preliminary gap discussion from a bounded scan** (74 records screened — about 4 over budget — 24 extracted,
9 full-text reads). A gap marked "open" means "not found in this scan and the closest papers
state assumptions that exclude it". It is never a novelty claim. Status words:
**open** (closest work excludes the case), **contested** (sources disagree or the premise is
disputed), **closed** (a source handles it), **reading gap** (likely covered somewhere not searched).

Track labels come from `programme/catalog.yaml` (read-only). "none" = no catalog track found.

## Input coverage

- Topic: identification and inference for causal effects in time series and event data when
  dependence is beyond geometric mixing (long memory, near-critical or power-law self-excitation).
- Reused: research_assistant `conformal_causal_fin_ts` S3 and S8 (read), other files grepped for
  memory terms (they cover memory least: one row, R8, and the T6 Hawkes-mixing item); prior jobs
  JOB-A-HAWKES-CHECK and JOB-NEW-CAUSAL-RISK-EXCEED.
- Not covered: Google Scholar / SSRN full search; SSRN only through web search; the ICLR 2026
  Hawkes identification paper was seen as a proceedings PDF without an id.

## Gaps

| # | Gap | Status | Evidence | Closest work | Track overlap |
|---|---|---|---|---|---|
| G1 | Orthogonal / DML inference for causal effects on **one long-memory series** | **open** | F: Ballinari–Wehrli state their single-series cross-fitting condition "is indeed never satisfied" for α(s)=O(s^{2d−1}); Yu et al. and Semenova et al. assume geometric β-mixing | arXiv:2411.10009, 2609.27531, 1712.09988 | **K** (DML reliability), **F-CAUSAL** (horizon-h DR), JOB-NEW-CAUSAL-RISK-EXCEED baseline B1; B1 historical |
| G2 | Which part of the orthogonal score carries memory; population vs path-conditional estimand under long-memory effect heterogeneity | **open** (narrow) | D (field-map §8), supported by F: BW Ass. 3 is on the score process | Bojinov–Shephard (design-based, survives); BW; classical LRD sample-mean rate | **K**, **F-CAUSAL**, F2 (design-based successor) |
| G3 | Does "as-i.i.d." cross-fitting still remove the empirical-process bias under long memory? | **contested** | F: Balkus et al. prove it for fixed m-dependence only; Semenova / Yu / BW use buffers | arXiv:2601.10899 | **K**, **F-CAUSAL** (same question as S3 gap 3 in the reuse review) |
| G4 | Local-projection / IRF causal inference robust to **fractional** (not near-unit-root) persistence | **open** (narrow) | F: Montiel Olea–Plagborg-Møller require a finite VAR(p) with known p; ARFIMA is infinite-order | arXiv:2007.13888; Baillie–Kapetanios–Papailias (ER 2017, A) for fractional IRFs without LP | **F-CAUSAL** (horizon-h pivot) |
| G5 | Linear Granger tests under long memory and spurious-regression risk | **closed for linear designs** | A: Bauer–Maynard surplus-lag test covers fractional regressors; Tsay–Chung document the hazard; Gaussian TE = Granger | J. Econom. 169(2) 2012; J. Econom. 96(1) 2000; arXiv:0910.4514 | **E** (only for nonlinear/kernel versions) |
| G6 | Nonlinear / nonparametric Granger or transfer-entropy tests under long memory | **open / reading gap** | A: only Gaussian VARFI-TE found; PCMCI papers note autocorrelation inflates false positives | Faes et al. 2021 (id not verified); arXiv:2003.03685 | **E**; memory-kernel workstream |
| G7 | Valid Granger / local-independence inference for Hawkes **uniformly as branching ratio → 1** or with power-law kernels | **open** | F: Wang–Kolar–Shojaie fix γ < 1 and bound the intensity; EDD give consistency only (heavy tails allowed up to C/u); A: Christgau et al. need i.i.d. replicates | arXiv:2007.07448, 1605.06759, 2203.13559 | **B2-COMP** (coverage vs branching ratio), Track A/G Hawkes-mixing items (JOB-A-HAWKES-CHECK ρ-dependence heuristic), **E** partly |
| G8 | Endo–exo confounding: a time-varying exogenous baseline masquerading as self-excitation (and as long memory) | **contested** | A/S: Hardiman et al. (near-critical) vs Filimonov–Sornette and Wheatley–Wehrli–Sornette (criticality rejected with flexible baseline) | arXiv:1302.1405, 1308.6756; QF 2019 doi:10.1080/14697688.2018.1550266 | **B2-COMP**, B2-HTE (compensator-time blocking depends on the fitted intensity) |
| G9 | Market impact as an **identified** causal estimand when order flow has long memory | **open / contested** | F: Leclère et al. counterfactual fixed by a shared-noise assumption and a fitted model; F: Muhle-Karbe et al. say core and reaction trades are indistinguishable in anonymous data; S: Naviglio et al. on misattributed autocorrelation; A: Maitrier et al. scrambling; A: Sato–Kanazawa show memory has an observed cause (splitting) | arXiv:2607.03239, 2601.23172, 2501.17096, 2502.16246, 2301.13505 | **B2-HTE** ("when does order flow move prices"), **B2-QUEUE** (instrument) |
| G10 | Interventions on a single self-exciting path: identification plus single-path asymptotics | **open, high pre-emption risk** | S: 2026 papers on MPP potential outcomes (Ryalen et al.) and Hawkes with carryover (Kresin et al.) are moving fast; F: IV-for-point-processes uses i.i.d. units | arXiv:2604.12977, 2604.12124, 2301.03246, 2111.07603 | none directly; B2-HTE nearest |
| G11 | Proxy / negative-control identification under long-range dependence | **contested** | S: Hao (PLOS One 2026) claims identification with an fBm limit, using time as the proxy (proofs unread); Track J assumes β-mixing | doi:10.1371/journal.pone.0347847 | **J** (kernel proxy) — triple intersection |
| G12 | Rough volatility as a treatment, confounder or effect modifier | **open, premise contested** | A: Cont–Das argue roughness may be a microstructure artefact; Fukasawa et al. find it real; no causal paper found | arXiv:2203.13820, 1905.04852 | none (F-KERNEL uses volatility regimes but is not causal) |
| G13 | Causal discovery on event sequences (THP, TNPAR, latent-confounder Hawkes) with error control | **open, low priority** | A/S: methods have no error control and no memory analysis | arXiv:2105.10884, 2508.11727 | **E** |

## Existing coverage

- Well studied: identification of dynamic causal effects (Rambachan–Shephard, Carlson–Shephard);
  design-based time-series experiments; DML under geometric mixing; Hawkes Granger structure
  (support of the kernel) and its consistent estimation; linear Granger testing robust to
  persistence; scaling limits of nearly unstable Hawkes.
- Thin: anything that states a rate or a coverage guarantee under polynomially decaying
  dependence; single-path asymptotics for point-process causal estimands; designs that identify
  impact rather than model it.

## Main method routes

1. Mixing + blocked or buffered cross-fitting + HAC (Semenova, BW, Yu). Breaks under long memory.
2. Design-based randomisation (Bojinov–Shephard; switchbacks). Survives memory; needs known
   assignment.
3. Structural point-process models with model-based counterfactuals (Leclère, Muhle-Karbe,
   Noorbakhsh). Survives memory by assumption; identification rests on the model.
4. Fractional / self-normalised econometrics (Bauer–Maynard; self-normalisation reviews; Hao).
   Handles memory for means and linear tests; little causal machinery.

## Controversies and untested assumptions

- Are markets near-critical? Hardiman et al. yes; Wheatley et al. no (flexible baseline).
- Is rough volatility real? Fukasawa et al. yes; Cont–Das possibly an artefact.
- Is impact mechanical or informational? Maitrier et al. argue mechanical.
- Do bespoke folds matter? Balkus et al. say not (m-dependence); the econometric line uses buffers.
- Unchecked by any source read: whether the memory of financial series is in the *score* of a
  causal estimator, or only in the raw series. This decides whether G1 bites in practice.

## Risks

- Several 2026 preprints (Ryalen et al., Kresin et al., Leclère et al., Muhle-Karbe et al., Yu
  et al.) show active groups; pre-emption risk is real for G9 and G10.
- G2 may be folklore that an econometrician would call obvious. Its value is as a diagnostic
  and a framing, not as a theorem.
- G11 depends on a paper whose proofs I did not read and whose identifying device (time as proxy)
  may not survive scrutiny.
- Spurious memory (breaks, regime shifts) is itself a confounder. Any long-memory causal claim
  must rule it out first (Haldrup–Kruse; Qu 2011 tests; Wheatley et al.).
