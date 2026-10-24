# Research opportunities — 2026-10-23 22:19 UTC

**Search scope:** topic: Where does inference for kernel tests, double-machine-learning causal intervals and Hawkes Granger tests break as memory grows, and does the break follow the H = 3/4 boundary in financial data? Sub-questions from proposals NEW-MB-KERNEL-G1-HSIC (HSIC null limit with general bounded kernels on each side of H_X+H_Y=3/2), G2-HFREE (projected or self-normalised kernel statistic valid across the 3/4 line without estimating H) and G3-RANK (asymmetric |x|-type or squared transforms and Hermite rank-2 protection on financial data)..
**Catalogue SHA-256:** `83ea8ca007df48c188a31f8bac7c9b3df02ee60fc07986d331a04267c13df4af`. **Registry SHA-256:** `90fe403f680c8f30be0bbd3051b83eba9e6f2d9f3ac754a8c4991f039e84866d`.

Recommendations only. This report changes no track decision, hold, protocol or accepted claim.

## Ranked ideas

Criteria fixed before ranking: importance of the unanswered question; extent of what remains unresolved; confidence the gap survives contrary search; feasibility of the next check; and fit with recorded skills and resources. Gap strength and project attractiveness can differ. Ties are permitted.

| Rank | Idea and precise question | What remains unresolved | Why it matters | Confidence | Feasibility | Next investigation |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | **MB-G1-HSIC-MIXED-NULL** For two independent Gaussian-subordinated long-memory series with one memory exponent below 1/2 and D_X+D_Y>1 (the mixed region), and for non-Gaussian marginals above the H_X+H_Y=3/2 line, what is the null limit of HSIC with a fixed bounded characteristic kernel? | Arcones 1994 and the newer U-process literature (e.g. 'Hermite ranks and U-statistics', HAL hal-01197638, abstract only, authors unverified) were not read in full; a general theorem may already imply the case. | It decides whether the empirical sum-rule knot seen in the October 2026 simulation (H_X+H_Y≈1.5) has a theorem behind it for HSIC, not just for distance covariance. | moderate | Pen-and-paper for one kernel and Gaussian inputs; a specialist proof check is needed. No data required. | Compute the centred feature Hermite coefficients for a Gaussian kernel with identity marginals; check whether Dehling–Taqqu-type reduction plus Major's vector theorem already gives the mixed-region limit. |
| 2 | **MB-G2-LOWFREQ-KERNEL** Can low-frequency (fixed number of cosine averages) inference in the style of Müller–Watson be applied to kernel feature maps to give an HSIC or MMD calibration that stays valid on both sides of the memory boundary without estimating H? | Worker inference, not in a source: above the line, low-frequency averages of rank-2 or cross features converge to non-Gaussian Hermite/Rosenblatt limits, so the Gaussian low-frequency framework may break exactly where it is needed. | It is the only route found that has proven uniform coverage across fractional persistence, and it was missing from every earlier G2 review. | low | A derivation for one kernel, then a small simulation if the derivation is viable. Needs a separate author decision. | For a Gaussian kernel and fGn inputs, derive the limit of q low-frequency cosine averages of the centred cross feature at H_X+H_Y below, at and above 3/2. |
| 3 | **MB-UNIFORM-SIZE-LIMITS** Do known impossibility and fragility results (Pötscher 2002; Preinerstorfer–Pötscher 2016; Müller 2007) imply that no memory-agnostic kernel, DML or Granger calibration can control size uniformly near and above the boundary, or only that consistent variance estimation must be given up? | Pötscher 2002 and Preinerstorfer–Pötscher were read at abstract level only; whether their model classes include fractional processes near d=1/4 is unchecked. | Prevents an over-claim: a referee familiar with Pötscher or Müller will ask what uniformity is claimed. | low | Reading task: two full texts. | Read Pötscher 2002 Sections 2-3 and Preinerstorfer–Pötscher 2016 Section 3 for the covariance model classes; check if fractional d∈(0,1/2) is inside them. |

## Coordinator assessment and coverage

Both workers loaded the registered `research-gap-finder` SKILL.md (SHA-256 `e81cb1de…1d25919`) at matching hashes; the coordinator checked every `finished`-event output hash for both jobs on 2026-10-23. Retrieval was live WebSearch plus WebFetch primary-source opens. This is a bounded scan, not a systematic review; no-paper-found is a search result, not a novelty claim. Where the two workers disagree, the reviewer verdict is adopted only where the reviewer opened the primary text the discovery worker did not; a verdict resting on either worker's own algebra stays marked unresolved. No track state, hold, claim, protocol or experiment changed.

**Adjudicated.**
1. MB-G1-HSIC-MIXED-NULL: reviewer's **narrow** adopted with a qualification. Coordinator opened arXiv 1704.08609 (Düker): it proves functional limit theorems for sample means and cross-autocovariances of multivariate linear processes with mixed short/long-range components (limits include operator fBm and matrix Rosenblatt). It does not mention HSIC. The "linear case is done" reading therefore holds through the reduction of linear-kernel HSIC to squared cross-covariance, which the coordinator accepts but did not check against Theorem 4.1's exact conditions. The non-linear bounded-kernel mixed case remains the open part. Lévy-Leduc & Taqqu (Metrika 2014) is univariate by abstract; full text unread.
2. MB-G2-LOWFREQ-KERNEL: **narrow** adopted; the reviewer sourced discovery's "may fail above the line" (Müller–Watson need Gaussian averages; Düker gives non-Gaussian limits there). Priority lowered.
3. MB-UNIFORM-SIZE-LIMITS: **provisional, unranked** adopted; it is a reading task, and discovery's wording on Pötscher 2002 overreached per Preinerstorfer–Pötscher (read in full by the reviewer).
4. MB-G3-LMSV-RANK1: **reject** (both agree). Qualification from the reviewer kept: rank 2 is possible without leverage when the volatility function is symmetric, so G3's premise fails in the standard LMSV model but not in every model.
5. MB-HAWKES-NO-H-LINE: **reject** (both agree); near-critical Hawkes gives H below 1/2, so no 3/4 line arises.
6. MB-CROSS-CLASS-MAP: **provisional** (both agree); forward-citation sweeps found no cross-class map.

**Best next investigation.** G1 restricted to non-linear bounded kernels in the mixed-memory case, after confirming Düker Thm 4.1's conditions. Blocked: Pötscher 2002 and Lévy-Leduc–Taqqu full texts; Semantic Scholar rate limits.



## Idea cards

### 1. MB-G1-HSIC-MIXED-NULL: For two independent Gaussian-subordinated long-memory series with one memory exponent below 1/2 and D_X+D_Y>1 (the mixed region), and for non-Gaussian marginals above the H_X+H_Y=3/2 line, what is the null limit of HSIC with a fixed bounded characteristic kernel?

**Plain English and setting.** Kernel independence tests (HSIC) compare two series. With long memory their usual null theory can fail. The best existing theorem (Betken–Dehling, for distance covariance) covers two cases and openly leaves a third, 'mixed' case and the nonlinear high-memory case unproved.

**Closest work and support.**
- [https://arxiv.org/html/2107.03041v3](https://arxiv.org/html/2107.03041v3) — Remark 2.1(i): 'Theorem 2.2 does not cover the case Dξ+Dη>1, but requires either Dξ<1/2 or Dη<1/2'; Theorem 2.2(ii) only for G1=G2=id (full_text, retrieved 2026-10-23T21:58Z): The closest independence theorem explicitly leaves the mixed region and the non-Gaussian above-line case open. Re-inspected this job.

**Contrary evidence.**
- [https://arxiv.org/abs/2107.03041](https://arxiv.org/abs/2107.03041) — Search hits for v3 (27 Jan 2026) and author listing (abstract, retrieved 2026-10-23T21:53Z): No newer Betken–Dehling version closing the mixed case was found; a Betken–Dehling–Kroll bootstrap paper is mixing-only (prior jobs). Absence is a search result.

**Precise remaining gap.** Arcones 1994 and the newer U-process literature (e.g. 'Hermite ranks and U-statistics', HAL hal-01197638, abstract only, authors unverified) were not read in full; a general theorem may already imply the case.

**Proposed contribution.** A theorem, or a counterexample, for one fixed kernel and one specified pair of transforms, in the mixed region; with the kernel's leading Hermite coefficients written out.

**Importance, confidence and fit.** It decides whether the empirical sum-rule knot seen in the October 2026 simulation (H_X+H_Y≈1.5) has a theorem behind it for HSIC, not just for distance covariance. Confidence: moderate. Fit: Directly continues TrackMB3 theory and the G1 proposal; uses the existing simulation as a check..

**Small next investigation.** Compute the centred feature Hermite coefficients for a Gaussian kernel with identity marginals; check whether Dehling–Taqqu-type reduction plus Major's vector theorem already gives the mixed-region limit.

**Difficulty and resources.** Tightness of the RKHS-valued cross process after empirical centring. 

**Stopping criterion.** Stop if an existing U/V-process or vector non-central theorem is shown to imply the mixed-region HSIC limit directly.

**Track relationship.** NEW-MB-KERNEL-G1-HSIC, NEW-MB-KERNEL, NEW-MEMORY-BOUNDARY; parent sources: ResearchHarness@e79a7f30f82592bd19f1fad96493c034ece1de03, catalog.yaml sha256 83ea8ca007df48c188a31f8bac7c9b3df02ee60fc07986d331a04267c13df4af, TrackMB3_KernelTestsLongMemory@61cc8d9b5ec3473bc9179d10fb6ad79116abf437.

**Coverage limits and contrary searches.** 22 live searches across all items this job; no forward-citation index (OpenAlex/Scholar) used; Arcones and Major originals not re-read. Queries: Betken 2026 arXiv independence test long-range dependent kernel OR distance correlation OR HSIC new; empirical characteristic function independence test long-range dependent non-central limit Rosenblatt 2025 2026; arXiv 2026 independence test fractional Gaussian noise kernel OR distance covariance strongly dependent size simulation Hurst 0.9; Hermite rank kernel mean embedding OR reproducing kernel long-range dependence U-statistic degenerate limit.

### 2. MB-G2-LOWFREQ-KERNEL: Can low-frequency (fixed number of cosine averages) inference in the style of Müller–Watson be applied to kernel feature maps to give an HSIC or MMD calibration that stays valid on both sides of the memory boundary without estimating H?

**Plain English and setting.** G2 asks for a kernel test that works whatever the memory level. Econometrics already has a method that keeps coverage across many persistence types, but only for linear quantities. Whether it carries over to kernel tests is unknown.

**Closest work and support.**
- [https://www.nber.org/papers/w23186.pdf](https://www.nber.org/papers/w23186.pdf) — Abstract and Section 'Coverage properties of restricted versions of the model' (Table 2 discussion): the general model 'controls coverage uniformly across these models'; I(0) sets have 'coverage of just 1%' under other persistence models (full_text, retrieved 2026-10-23T21:50Z): Uniform coverage over I(0), I(1), local-to-unity and fractional persistence is achievable for a linear second-moment (covariability) parameter using q≈12 low-frequency cosine averages. Linear only; relies on Gaussian low-frequency limits.
- [https://www.princeton.edu/~umueller/robLRVP.pdf](https://www.princeton.edu/~umueller/robLRVP.pdf) — Theorem 1 and discussion p.1335-1338: consistent LRV estimators 'necessarily lack robustness'; Theorem 1 is 'Rather than being a statement about the impossibility of valid inference'; Theorem 2(iv) fixed-b kernel estimators not robust to C2T(d) (full_text, retrieved 2026-10-23T21:47Z): Any long-run-variance estimator consistent in a benchmark model is fragile to small low-frequency contamination; inconsistent (low-frequency) estimators can be robust. Not an impossibility theorem for inference.

**Contrary evidence.**
- [https://arxiv.org/abs/0810.2276](https://arxiv.org/abs/0810.2276) — Abstract (Shao 2009, Econometric Theory 25(1):195-210), via search snippet (abstract, retrieved 2026-10-23T21:52Z): Prewhitened frequency-domain independence statistics are N(0,1) while each series may have short memory, long memory or anti-persistence. Linear dependence only. TrackMB3 records an F read; this job saw the abstract only.
- [https://arxiv.org/abs/2601.13468](https://arxiv.org/abs/2601.13468) — Ghoshal & Shao (EJS 2026) Theorems 3.1-3.5, inherited F read from GAPLIT G2 (abstract, retrieved 2026-10-23T21:52Z): RKHS split-project-self-normalise exists under strong mixing only. Anticipates the method idea.

**Precise remaining gap.** Worker inference, not in a source: above the line, low-frequency averages of rank-2 or cross features converge to non-Gaussian Hermite/Rosenblatt limits, so the Gaussian low-frequency framework may break exactly where it is needed.

**Proposed contribution.** A check of whether low-frequency averages of centred kernel features have a usable (Gaussian or tabulable) joint limit across the line; if yes, a calibrated H-free kernel test; if no, a stated obstruction.

**Importance, confidence and fit.** It is the only route found that has proven uniform coverage across fractional persistence, and it was missing from every earlier G2 review. Confidence: low. Fit: Fits TrackMB3 and the author's econometrics strengths; reuses the simulation code base..

**Small next investigation.** For a Gaussian kernel and fGn inputs, derive the limit of q low-frequency cosine averages of the centred cross feature at H_X+H_Y below, at and above 3/2.

**Difficulty and resources.** Non-Gaussian limits above the line; choice of q. 

**Stopping criterion.** Stop if the derivation shows a non-Gaussian, H-dependent limit that cannot be tabulated without H.

**Track relationship.** NEW-MB-KERNEL-G2-HFREE, NEW-MB-KERNEL, NEW-MEMORY-BOUNDARY; parent sources: ResearchHarness@e79a7f30f82592bd19f1fad96493c034ece1de03, catalog.yaml sha256 83ea8ca007df48c188a31f8bac7c9b3df02ee60fc07986d331a04267c13df4af, TrackMB3_KernelTestsLongMemory@61cc8d9b5ec3473bc9179d10fb6ad79116abf437.

**Coverage limits and contrary searches.** Müller–Watson 2016 'Low-frequency econometrics' and 2008 'Testing models of low-frequency variability' not read; no search of Müller–Watson forward citations. Queries: low-frequency cosine transforms independence test nonlinear kernel long memory robust Müller Watson extension HSIC; kernel two-sample test long memory time series unknown Hurst exponent calibration self-normalized maximum mean discrepancy 2026; long-range dependence Hilbert space self-normalization OR fixed-b functional data change point Hurst-free.

### 3. MB-UNIFORM-SIZE-LIMITS: Do known impossibility and fragility results (Pötscher 2002; Preinerstorfer–Pötscher 2016; Müller 2007) imply that no memory-agnostic kernel, DML or Granger calibration can control size uniformly near and above the boundary, or only that consistent variance estimation must be given up?

**Plain English and setting.** Before promising an 'H-free' test it helps to know whether such a test can exist. Econometrics has several 'it cannot be done uniformly' results for related problems. None of the earlier jobs cited them.

**Closest work and support.**
- [https://www.princeton.edu/~umueller/robLRVP.pdf](https://www.princeton.edu/~umueller/robLRVP.pdf) — Theorem 1 and discussion p.1335-1338: consistent LRV estimators 'necessarily lack robustness'; Theorem 1 is 'Rather than being a statement about the impossibility of valid inference'; Theorem 2(iv) fixed-b kernel estimators not robust to C2T(d) (full_text, retrieved 2026-10-23T21:47Z): Any long-run-variance estimator consistent in a benchmark model is fragile to small low-frequency contamination; inconsistent (low-frequency) estimators can be robust. Not an impossibility theorem for inference.
- [https://www.econometricsociety.org/publications/econometrica/browse/2002/05/01/lower-risk-bounds-and-properties-confidence-sets-ill‐posed](https://www.econometricsociety.org/publications/econometrica/browse/2002/05/01/lower-risk-bounds-and-properties-confidence-sets-ill‐posed) — Abstract (Pötscher 2002, Econometrica 70(3):1035-1065) (abstract, retrieved 2026-10-23T21:45Z): Minimax risk for spectral density at zero is infinite; confidence sets nearly uninformative; similar results for persistence and long-memory parameters. Abstract only; full text unavailable (working-paper URL returned HTML).
- [https://arxiv.org/abs/1304.1383](https://arxiv.org/abs/1304.1383) — Abstract (Preinerstorfer & Pötscher, Econometric Theory 32 (2016) 261-358) (abstract, retrieved 2026-10-23T21:46Z): HAR F-type tests have size one or nuisance-infimal power zero under weak assumptions on stationary disturbances. Abstract does not mention long memory.

**Contrary evidence.**
- [https://www.nber.org/papers/w23186.pdf](https://www.nber.org/papers/w23186.pdf) — Abstract and Section 'Coverage properties of restricted versions of the model' (Table 2 discussion): the general model 'controls coverage uniformly across these models'; I(0) sets have 'coverage of just 1%' under other persistence models (full_text, retrieved 2026-10-23T21:50Z): Uniform coverage over I(0), I(1), local-to-unity and fractional persistence is achievable for a linear second-moment (covariability) parameter using q≈12 low-frequency cosine averages. Linear only; relies on Gaussian low-frequency limits.

**Precise remaining gap.** Pötscher 2002 and Preinerstorfer–Pötscher were read at abstract level only; whether their model classes include fractional processes near d=1/4 is unchecked.

**Proposed contribution.** A short framing note that states which uniformity is impossible and which is achievable, with the exact premises, to bound the G2 and DML repair claims.

**Importance, confidence and fit.** Prevents an over-claim: a referee familiar with Pötscher or Müller will ask what uniformity is claimed. Confidence: low. Fit: Framing for TrackMB1 (H5 repair) and G2; no compute..

**Small next investigation.** Read Pötscher 2002 Sections 2-3 and Preinerstorfer–Pötscher 2016 Section 3 for the covariance model classes; check if fractional d∈(0,1/2) is inside them.

**Difficulty and resources.** Full texts paywalled; working-paper copy failed. 

**Stopping criterion.** Stop once each result's model class is recorded against the d grid used in the simulation.

**Track relationship.** NEW-MEMORY-BOUNDARY, NEW-MB-KERNEL-G2-HFREE; parent sources: ResearchHarness@e79a7f30f82592bd19f1fad96493c034ece1de03, catalog.yaml sha256 83ea8ca007df48c188a31f8bac7c9b3df02ee60fc07986d331a04267c13df4af, TrackMB3_KernelTestsLongMemory@61cc8d9b5ec3473bc9179d10fb6ad79116abf437.

**Coverage limits and contrary searches.** Abstract-level for two of three sources. Queries: uniformly valid inference long memory parameter d uniform over [0,1/2) confidence interval mean or regression subsampling Andrews Guggenberger Mikusheva fractional; Müller 2007 A theory of robust long-run variance estimation Journal of Econometrics long memory contamination.


## Best next investigations

- **1. MB-G2-LOWFREQ-KERNEL:** For a Gaussian kernel and fGn inputs, derive the limit of q low-frequency cosine averages of the centred cross feature at H_X+H_Y below, at and above 3/2. Cheapest decisive check: one derivation tells whether the route is open or blocked.
- **2. MB-G1-HSIC-MIXED-NULL:** Compute the centred feature Hermite coefficients for a Gaussian kernel with identity marginals; check whether Dehling–Taqqu-type reduction plus Major's vector theorem already gives the mixed-region limit. Strongest supported gap: the closest source states the hole in its own words.
- **3. MB-UNIFORM-SIZE-LIMITS:** Read Pötscher 2002 Sections 2-3 and Preinerstorfer–Pötscher 2016 Section 3 for the covariance model classes; check if fractional d∈(0,1/2) is inside them. Cheap, and it constrains how the G2 and TrackMB1 H5 repairs can be worded.

## Appendix: provisional and rejected candidates

- **MB-CROSS-CLASS-MAP (provisional, uncertain):** Has any paper mapped size or coverage across a common memory grid for kernel tests, DML intervals and Hawkes Granger tests together? — Forward citations of Betken–Dehling, Zhai–Zhong–Wu and Zhai–Cheng–Wu still not enumerated. Müller–Watson Table 2 is a precedent format (coverage across persistence models) for a linear parameter, not a pre-emption. Next: Run the forward-citation sweep through a working scholarly index, as the 2026-10-22 job recommended.
- **MB-G3-LMSV-RANK1 (rejected, answered):** Does an |x|-type or squared transform of financial returns enjoy Hermite rank-2 protection under standard long-memory stochastic volatility models, with or without leverage? — Model-level answer only. Whether real returns follow exp-LMSV, and the empirical rank, remain open (rank tests: Beran–Möhrle–Ghosh, abstract only in GAPLIT). The step 'E[|r| | X] = σ(X)E|Z|, so rank = rank of σ' is worker derivation for the no-leverage case. Next: If G3 continues, state the volatility function explicitly and compute its first Hermite coefficient; cite Kulik–Soulier for the leverage case.
- **MB-HAWKES-NO-H-LINE (rejected, narrowed):** Does the size of Hawkes likelihood-ratio Granger tests break at an H=3/4 memory line? — Exploratory simulation (JOB-NEW-MEMORY-BOUNDARY-SIM-20261015 hawkes) shows size depends on ρ, kernel tail and baseline; the size-inflation mechanism under misspecification is not proved. Next: None from this job; TrackMB2 literature is the reference.

## Change since previous assessment

- MB-G1-HSIC-MIXED-NULL: new candidate
- MB-G2-LOWFREQ-KERNEL: new candidate
- MB-UNIFORM-SIZE-LIMITS: new candidate
- MB-G3-LMSV-RANK1: new candidate
- MB-CROSS-CLASS-MAP: new candidate
- MB-HAWKES-NO-H-LINE: new candidate
- HAWKES-EVENT-MIXING: absent from new assessment; prior finding retained
- KERNEL-DISTRIBUTIONAL-CI: absent from new assessment; prior finding retained
- DML-RESIDUAL-INVARIANCE: absent from new assessment; prior finding retained
- BOND-FORCED-FLOW: absent from new assessment; prior finding retained
- COST-AWARE-OFFLINE-RL: absent from new assessment; prior finding retained
- QUEUE-POSITION-IV: absent from new assessment; prior finding retained

## Separate reviewer challenge

The independent [review](../JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023/review.md) records its own source checks and disagreements. The coordinator assessment above states what was accepted and what remains open.


## Search and review record

Read the job-owned [search log](search-log.md), [evidence matrix](evidence-matrix.md) and [candidate records](candidates.json) for queries, read depth, exclusions and source revisions.
