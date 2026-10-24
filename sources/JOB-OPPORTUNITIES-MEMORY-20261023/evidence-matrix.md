# Evidence matrix — JOB-OPPORTUNITIES-MEMORY-20261023

Read depth: **F** = full-text passage inspected by this job. **A** = abstract or search snippet only.
**Inh-F** = a full-text read by an earlier job, not re-read here.
"New" means the source appears in none of the earlier job or track bibliographies listed in `search-log.md`.

| Source | Depth | New? | Exact passage or locator | Bears on | For / against |
|---|---|---|---|---|---|
| Kulik & Soulier, "Limit theorems for LMSV models with infinite variance", Adv. Appl. Probab. 44(4):1113–1141 (2012), arXiv:1109.5298 | F | New | §5.1: "for σ(x)=exp(x), τ*_p(s)=1 for all s". Introduction: partial sums "entirely similar" under leverage, but "sample autocovariances may be very different in the presence of leverage". Ex. 5.2: rank of K is 1 "except if E[Z_0 η_0]=0". | G3 | **Against** the G3 finance premise. In exp-volatility models, powers of returns have rank 1 in the latent Gaussian. Leverage changes ranks and rates of sample covariances. Infinite-variance setting. The finite-variance case was not re-read. |
| Müller, "A theory of robust long-run variance estimation", J. Econometrics 141(2):1331–1352 (2007) | F | New | Theorem 1: consistent scale-equivariant LRV estimators "necessarily lack robustness" in contamination neighbourhoods. The text says Theorem 1 is not "a statement about the impossibility of valid inference". Theorem 2(iv): fixed-b kernel estimators are not robust to C2T(d). | MB, G2, TrackMB1 H5 | Constrains H-free calibration. Plug-in HAC-type repairs are fragile. Robust inconsistent estimators exist. Not an impossibility result. |
| Müller & Watson, "Long-run covariability", Econometrica 86(3):775–804 (2018); NBER w23186 | F | New | Abstract: coverage is controlled over "arbitrary linear combinations of I(0), I(1), near unit roots and fractionally integrated processes". Table 2 text: I(0) sets have "coverage of just 1%" under other models. The general model "controls coverage uniformly". | G2, MB map | **For** feasibility of uniform-over-persistence inference, for a linear second-moment parameter only. It is a precedent format for a coverage map across persistence models. |
| Pötscher, "Lower risk bounds and properties of confidence sets for ill-posed estimation problems…", Econometrica 70(3):1035–1065 (2002) | A | New | Abstract: minimax risk infinite; confidence sets nearly uninformative; similar results for persistence and long-memory parameters. | MB, G2 | Suggests that any calibration needing an estimate of H, or of the zero-frequency spectrum, cannot be uniformly informative. Abstract only. |
| Preinerstorfer & Pötscher, "On size and power of HAR tests", Econometric Theory 32:261–358 (2016), arXiv:1304.1383 | A | New | Abstract: size equal to one or nuisance-infimal power zero under weak assumptions. | MB, TrackMB1 | A finite-sample impossibility for HAR tests in regression. The abstract does not state long memory. |
| Betken & Dehling, arXiv:2107.03041v3 (27 Jan 2026) | F (re-inspected) | Known | Remark 2.1(i): "Theorem 2.2 does not cover the case Dξ+Dη>1, but requires either Dξ<1/2 or Dη<1/2". Theorem 2.2(ii) needs G1=G2=id. Intro: HSIC, MMD and dCov "coincide" for a specific kernel. | G1 | **For** the G1 gap. The mixed and nonlinear above-line cases stay open in this source. |
| Shao, "A generalized portmanteau test for independence between two stationary time series", Econometric Theory 25(1):195–210 (2009), arXiv:0810.2276 | A (Inh-F in TrackMB3) | Known to TrackMB3; absent from G1/G2 reviews | Abstract: N(0,1) null limits with "short memory, long memory or anti-persistence" allowed. | G1, G2 | **Against** a broad "no H-robust independence test exists" claim. Linear dependence only. |
| Ghoshal & Shao, EJS 2026, arXiv:2601.13468 | Inh-F | Known | Theorems 3.1–3.5 (GAPLIT G2). | G2 | Against novelty of the split-project-self-normalise idea. Strong mixing only. |
| Bai & Taqqu 2019, arXiv:1710.01612 | Inh-F | Known | Theorems 2.1, 2.5, 3.1 (GAPLIT G3, coordinator-checked). | G3 | The generic mechanism is closed. Consistent with the Kulik–Soulier finding. |
| "Hermite ranks and U-statistics", HAL hal-01197638 | A (snippet) | New lead | Compares notions of Hermite rank for U-statistics under LRD; conditions for negligible remainder. | G1 | Possible prior art for the kernel-rank step of G1. Authors and content unverified. |
| Filimonov & Sornette, arXiv:1308.6756 | A | Known (TrackMB2) | Kernel misspecification biases the branching ratio upward. | Hawkes arm | Supports reframing the Hawkes arm away from an H line. |

## No-paper-found results (search results only)

- No HSIC or kernel independence theorem covering the mixed region (queries 4, 10, 15, 21).
- No H-free kernel two-sample or independence calibration under LRD (queries 7, 8, 19).
- No DML, orthogonal-score or cross-fitting coverage study under LRD beyond TrackMB1's corpus (queries 13, 24).
- No size study of Hawkes or fractional Granger tests organised by H (queries 14, 17, 18).
- No cross-statistic survey of where root-n inference fails by Hermite rank (query 25).
