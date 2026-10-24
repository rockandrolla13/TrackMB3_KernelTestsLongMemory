# Research Gap Analysis — memory boundary (NEW-MEMORY-BOUNDARY, G1–G3)

Sections follow `research-gap-finder/references/gap_analysis_template.md` (SHA-256 `a5c625dd…f34448`).
The original Chinese headings appear in brackets.
This is a preliminary gap discussion from a bounded scan. It makes no novelty claim.

## Input coverage (输入覆盖范围)

- Topic: where inference for kernel tests, DML intervals and Hawkes Granger tests breaks as memory grows. Does the break follow H = 3/4? G1–G3 are sub-questions.
- Papers: six earlier jobs and three track literature folders (well over 100 records), plus 25 new queries and 10 opens in this job. The job made four full-text reads, one of them a re-inspection.
- User ideas included: yes. These are the three proposals and the earlier C1 and C2 ideas.
- Boundary: no forward-citation index, no SSRN, two key impossibility papers at abstract level only.

## Topic summary (主题概括)

Earlier jobs settled the main frame. H = 3/4 is the boundary for **rank-2** statistics only.
Rank-1 statistics, such as the MMD halves test and the population ATE under monotone effects, lose root-n inference at any H > 1/2.
Two-series cross terms follow the sum rule H_X + H_Y = 3/2.
Hawkes tests degrade with branching ratio and kernel tail, not with H.
The October 2026 exploratory simulation agrees with each of these statements.

## Existing coverage (已有工作覆盖面)

- Well studied:
  - Hermite-rank limit theory (Taqqu, Dobrushin–Major, Dehling–Taqqu).
  - Distance-covariance theory on two sides of the line (Betken–Dehling).
  - High-dimensional covariance under LRD (Zhai–Zhong–Wu).
  - Rank instability (Bai–Taqqu).
  - H-free scalar calibration (Betken–Wendler, Bai–Taqqu–Zhang).
  - Prewhitened linear independence tests under LRD (Shao 2009).
- Newly added by this job:
  - Econometric fragility and impossibility results (Müller 2007; Pötscher 2002; Preinerstorfer–Pötscher 2016).
  - Uniform-over-persistence low-frequency inference (Müller–Watson 2018).
  - LMSV Hermite ranks with and without leverage (Kulik–Soulier 2012).
- Still weak:
  - HSIC in the mixed region.
  - Any H-free calibration for a kernel statistic.
  - Forward citations of the 2026 anchors.

## Common method routes (常见方法路线)

- Route 1, Hermite expansion plus a non-central limit theorem. This gives exact rates but needs H to calibrate.
- Route 2, subsampling or self-normalised block sampling. This is H-free for scalar statistics that converge, but slow near the line.
- Route 3, low-frequency (fixed-q) inference. It controls coverage uniformly across I(d) for linear parameters. It has not been tried for kernel features.
- Route 4, prewhitening by fractional filtering. It works for linear tests, but it estimates d.

## Current bottlenecks (当前瓶颈)

- Data or benchmark: no common memory grid across the three method classes outside this programme's own simulation.
- Method:
  - Tightness of RKHS-valued cross processes in the mixed region.
  - Non-Gaussian low-frequency limits above the line (worker inference).
- Evaluation: near-line finite-sample size is poor even where theory holds. Betken–Dehling report 10–17% size at H = 0.7.
- Deployment: finance proxies have uncertain and clock-dependent H. Under exp-LMSV the relevant rank is 1, not 2.

## Controversies (争议点)

- Inconsistent conclusions: whether finance series sit above 3/4 depends on proxy, clock and break treatment (DEEPLIT 2026-10-22).
- Assumptions often left untested:
  - "Symmetric transforms of returns give rank-2 protection." This is false in exp-volatility LMSV (Kulik–Soulier §5.1).
  - "An H-free test is possible." Müller 2007 shows consistent variance estimators are fragile. Pötscher 2002 (abstract) suggests uniform informativeness fails for long-memory parameters.

## Research gaps (research gap)

- **Gap 1, MB-G1-HSIC-MIXED-NULL (ranked 1, narrowed, moderate confidence).**
  - Why it may hold: Betken–Dehling Remark 2.1(i) states the mixed case is not covered. No newer version or HSIC paper surfaced.
  - Why existing work misses it: the proofs need each D > 1/2, or Gaussian inputs.
  - Confidence: moderate. An unread general U-process or vector theorem may already imply it.
- **Gap 2, MB-G2-LOWFREQ-KERNEL (ranked 2, open, low confidence; best next investigation).**
  - Why it may hold: Müller–Watson achieve uniform coverage across fractional persistence for linear parameters. No kernel extension was found.
  - Why existing work misses it: the econometrics and kernel-testing literatures do not cite each other.
  - Confidence: low. Above the line, feature averages may have non-Gaussian limits.
- **Gap 3, MB-UNIFORM-SIZE-LIMITS (ranked 3, narrowed, low confidence).** The question is which uniformity an "H-free" claim can honestly promise. This is a framing gap. Its two key sources were read at abstract level only.
- **Not gaps:**
  - MB-G3-LMSV-RANK1 is answered at model level.
  - MB-HAWKES-NO-H-LINE is reframed: Hawkes has no H line.
  - MB-CROSS-CLASS-MAP stays provisional until a citation sweep is done.

## Small entry points (可切入的小问题)

1. Compute the leading Hermite coefficients of the centred Gaussian-kernel cross feature for fGn inputs in the mixed region.
2. Derive the limit of q low-frequency cosine averages of that feature below, at and above H_X + H_Y = 3/2.
3. Record each impossibility result's model class against the simulation's d grid.

## Possible experiments (潜在实验验证思路)

These would need a separate author decision. This job ran none.
- Re-use the existing kernel-linear simulation cells. Add a low-frequency calibrated HSIC arm at H_X + H_Y ∈ {1.4, 1.5, 1.6}.
- An exp-LMSV arm with and without leverage, testing whether |r|-feature kernel tests break near H = 1/2 rather than 3/4.

## Risks (风险提示)

- Gap 1 may be a corollary of older vector non-central theory. A specialist check is needed before any contribution claim.
- Gap 2 may fail for the very reason that motivates it: non-Gaussian limits above the line.
- The search has no forward-citation coverage. Each gap's novelty remains unverified.
