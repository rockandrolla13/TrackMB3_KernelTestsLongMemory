# Gap analysis — kernel tests under long memory, against hypotheses H0–H10

Execution `claude-mb3-lit-20261016` · Agent 2 · 2026-10-16. Method: research-gap-finder SKILL.md
(SHA-256 `e81cb1de…5919`) and its `gap_analysis_template.md` (SHA-256 `a5c625dd…4448`), sections
translated to English. Inputs: `review.md`, `systematic-review.md`, `evidence-matrix.md`, plus the prior
maps in `sources/` (reused, not redone).

**Status words.** *Known* = a source states it (class given). *Partly known* = stated for a narrower
case. *Open* = not found in this bounded scan **and** a named obstacle explains why it might be open.
*Contested* = sources disagree, or a source qualifies the claim. "Open" never means "nobody has done it".
Each gap says whether it looks like a **true gap**, an **engineering/evaluation gap**, or a **reading
gap** (likely exists where I could not search).

---

## 1. Input coverage
- Topic: calibration of kernel two-sample, independence and CI tests under long-range dependence; the
  Hermite-rank map; calibration (SN, subsampling, block/wild bootstrap, prewhitening); finance memory.
- Papers: 103 screened, 38 extracted, 15 full-text section reads (14 usable). Plus prior field maps.
- Includes the author's own ideas: yes — the 11 hypotheses (H0–H10) and sub-items from Agent 1.
- Boundary: SSRN unreachable (HTTP 403); Google Scholar partial (4/10 queries returned no results);
  Semantic Scholar rate-limited; key classics (Hsing–Wu 2004, Dehling–Taqqu 1991, Arcones 1994,
  Shao 2011 full text, Poskitt 2007) identified but not read.

## 2. Topic summary
Long memory changes the null law of a test statistic when the autocovariance of the statistic's
leading term stops being summable. For linear statistics this is classical. For quadratic and cross
statistics the switch is at H = 3/4 or H_X + H_Y = 3/2, and has been proved for sample covariances
(ZKW, ZZW), distance covariance with Gaussian or subordinated-Gaussian data on each side separately
(Betken–Dehling), OLS t-ratios (Tsay–Chung) and conformal coverage indicators (Zhai–Cheng–Wu). No paper
found does this for HSIC, MMD or kernel CI tests. Calibration that survives the switch exists for
scalar statistics (SN + block sampling; estimated-normalisation block sampling) but not for kernel
statistics.

## 3. Coverage of existing work
- **Well covered:** Hermite limit theory (classical); SN under LRD for scalar change-point statistics;
  subsampling validity under LRD; linear spurious-regression sum rule; linear prewhitened independence
  tests under LRD; kernel tests under mixing (closed, prior job X1).
- **Thin:** dCov under LRD (one preprint, never simulated above H = 0.7); covariance / GGM inference
  across 3/4 (one group, 2026); conformal coverage across 3/4 (same group, Sep 2026).
- **Not found:** HSIC / MMD / energy-distance null laws under LRD; kernel CI under LRD; SN kernel tests
  under LRD; a size map of kernel tests across a grid of H.

## 4. Common method routes
1. **Hermite expansion + (non-)central limit theorem**, then subsampling (Betken–Dehling; DRW).
2. **Gaussian approximation + block bootstrap / block sampling** with finite-sample bounds, restricted to
   H < 3/4 (ZKW, ZZW), now extended across 3/4 with estimated normalisation (Zhai–Cheng–Wu).
3. **Self-normalisation**, with H-tabulated critical values (Betken 2016) or with subsampling to avoid
   Ĥ (Betken–Wendler 2018; Bai–Taqqu–Zhang 2016); kernel SN under mixing (Ghoshal–Shao 2026).
4. **Prewhitening**, linear (Haugh 1976; Hong 1996; Shao 2009 under LRD).

## 5. Bottlenecks
- *Theory:* Hilbert-space non-central limits for non-linear subordination above the line (Betken–Dehling
  cover only G = id there); tightness of the empirical characteristic-function process for a series with
  D < 1/2 (the mixed case).
- *Method:* every calibration that is valid on both sides of the line needs either the rate (Ĥ) or a
  self-normaliser plus subsampling; finite-sample oversize persists at n ≈ 10³–10⁴ (ZCW Table 2:
  13.9% at n = 4,096; Betken–Dehling 7–17.5% at n ≤ 1000).
- *Evaluation:* no published grid in H near the lines for any kernel test; Betken–Dehling stop at 0.7,
  ZZW/ZKW use β ∈ {0.9, 0.55}.
- *Data:* finance memory estimates move 0.1–0.3 with the estimator; order signs sit near 0.75
  (Sato–Kanazawa, S) — exactly where the asymptotics change.

## 6. Contested points
- **"Real data typically has H < 0.75"** (Betken–Dehling Rem 2.1(ii), unsourced) vs. finance estimates
  (activity series 0.75–0.96; order signs 0.65–0.9, TSE ≈ 0.75). Unchanged from the prior job.
- **"Safe below the line"**: Betken–Dehling (7–17.5% size at H ≤ 0.7, n ≤ 1000) and VSVL (bounded but
  non-normal t-ratios below 1/2) say finite-sample trouble starts below the line.
- **Rank-2 reliability**: H1c and H8's symmetric-score claims assume an exactly even G. Bai–Taqqu (2018)
  argue ranks > 1 are statistically unstable.
- **SN "mild" H-dependence** (H4b): contradicted by Betken (2016) Table 1.

## 7. Per-hypothesis gap table

| ID | Claim (short) | Status | Evidence class and sources | Kind of gap | Comment |
|---|---|---|---|---|---|
| H0 | Failure line set by effective Hermite rank of the leading null term | **Known** as a framework | Breuer–Major, Dobrushin–Major, Taqqu (I); instances F: ZKW, ZZW, BD, TC (S), ZCW | — | Organising lens, not a contribution. Referee "known" risk high (confirmed). |
| H1 | MMD homogeneity is rank 1; fails at every H > 1/2 | **Partly known** | DRW 2017 (F): two-sample U-processes, rank-m non-central limits; MMD itself not found | engineering → corollary | Rank-one χ² law and κ_H: not found (true gap, low theory value) |
| H1b | Size excess depends on (n, H) only via (n/l)^{2H−1} | **Open** | none found; Zhang–Lahiri–Nordman 2022 (A) say block behaviour depends on non-linearity beyond rank — a possible complication | true gap (empirical regularity) | Cheap to test; worth stating if it holds |
| H1c | Even G restores 3/4 for MMD | **Partly known / contested** | mechanism known for scalar and coverage statistics (ZCW F: \|W\| → h^{2−4β}); Bai–Taqqu 2018 (A) rank instability | engineering + caution | Will likely hold for exact \|ξ\|; fragile for real \|r_t\| |
| H2 | HSIC/dCov √n-Gaussian iff H_X + H_Y < 3/2 | **Partly known** | linear: TC 2000 (S), VSVL (F), RH 1997 (A); dCov: BD (F) for both H < 3/4 and, G = id, sum > 3/2; HSIC: not found | **true gap (HSIC, general G)** | Linear sum rule is classical — headline must be the kernel version |
| H2b | Mixed case (one H > 3/4, sum < 3/2) still √n | **Open** | BD Rem 2.1(i) (F) says not covered; no later work found (1 citer) | **true gap (theory)** | Highest-value theory gap found; obstacle named (ECF tightness for D < 1/2) |
| H3 | Above the line, the null collapses to a squared latent cross-correlation | **Known for dCov, G = id**; open otherwise | BD Cor 2.1(ii) + Thm 4.1 (F) | true gap for non-linear G and for HSIC | Lowers H3's novelty; raises its credibility. The non-linear-G case (signs) is open even for dCov |
| H4a | n-rate subsampling fails above the line; rate-corrected subsampling works | **Known for scalars; open for kernels** | Bertail–Politis–Romano 1999 (A); Zhang–Ho–Wendler–Wu 2013 (I); ZCW 2026 Thm 7 + Table 2 (F) | engineering gap for kernels | Pre-emption risk **high** (Wu group already does this across 3/4) |
| H4b | SN: pivotal below, H-dependent above; fails gracefully | **Partly known / contested** | SN limits depend on H (Betken 2016 F; Shao 2011 via BTZ, S); H-free SN+block sampling known (BTZ 2016 F; BW 2018 F); kernel SN under mixing (GS 2026 F); kernel SN under LRD not found | true gap (combination) | "Mild H-dependence" contradicted; "graceful failure" untested; risk **high** |
| H4c | Block bootstrap valid below the line, fails above | **Partly known** | fails for non-central limits: Lahiri 1993 (S), Tewes 2018 (A); valid below 3/4 for covariances: ZZW/ZKW (F); MBB oversized at β ≤ 0.75 (ZCW F) | engineering gap for degenerate kernel statistics | Low priority (as Agent 1 said) |
| H4d | Wild bootstrap invalid for rank-1 at any memory; may hold for cross statistics below the line | **Open** | CSG 2014 conditions exclude LRD (F, prior job); no LRD wild-bootstrap theory found | true gap (theory), cheap empirically | KL simulation will speak to it |
| H5 | Wild-bootstrap oversize at H = 0.5 is a multiplier artefact | **Known (hygiene)** | CSG 2014 Table 1 + text (F): centred multipliers over-reject, decreasing with n; Leucht–Neumann Sec 5 (cited there) | engineering | Fix before further wild-bootstrap reading; cite CSG |
| H6 | Shift test most robust below the line | **Open** | shift test theory under mixing only (Chwialkowski–Gretton 2014, A) | evaluation gap | Low theory value |
| H7 | Kernel CI: sum rule on residual memories + nuisance channel; memory only in Z is harmless | **Open** | GCM/GKCM/GTCM: i.i.d. or weak dependence (A); ZZW linear Gaussian β > 3/4 (F); Arcones 1994 (I) is the natural tool | **true gap**, with a **reading-gap** risk (SSRN, econometrics) | "Memory in Z only" half is clean and testable |
| H7b | Kernel Granger with p lags: small-constant sum rule | **Open (reading gap likely)** | Poskitt 2007, 2008 (I); FAR sieve (I) | reading gap | Read Poskitt before any simulation |
| H8 | Coverage fluctuation follows the Hermite rank of the score indicator | **Known for symmetric \|W\| scores (rank 2)**; open for signed scores and kernel weights | ZCW 2609.33868 Thm 6 (F) | pre-empted (symmetric); true gap (signed, kernel) | Hand to F-KERNEL with ZCW; risk **high** |
| H9 | Most finance pairs sit outside the √n regime | **Contested / partly supported** | prior `FME` (F/S/A); Sato–Kanazawa 2023 (A+S): signs H ≈ 0.75; BD Rem 2.1(ii) contrary (unsourced) | data gap | Needs own estimates with CIs |
| H10 | Fractional prewhitening restores HSIC size | **Known for linear tests; open for HSIC** | Shao 2009 (F): prewhitened portmanteau N(0,1) under LRD; Tsay 2000 (A) | engineering | Strong baseline with a citable linear theorem; "obvious" risk high |

## 8. Research gaps (template form)

**Gap 1 — the mixed case for cross statistics (H2b).**
- Why it may hold: the cross term's autocovariance is ρ_ξ ρ_η-type and summable when D_ξ + D_η > 1; the
  marginal centring terms are of order n^{H_X+H_Y−2} (D, Agent 1).
- Why not covered: Betken–Dehling's proof needs both D > 1/2 (Rem 2.1(i), F); the obstacle is tightness
  of the ECF process of the series with D < 1/2.
- Credibility: medium-high as a true theory gap; not found in forward citations (1 citer).

**Gap 2 — HSIC with general bounded kernels under LRD, both sides of the line (H2, H3 non-linear G).**
- Why it may hold: Bochner representation turns HSIC into BD's object with a finite spectral measure (D,
  prior job).
- Why not covered: BD do dCov only, and the above-line law only for G = id; HSIC papers are mixing-only.
- Credibility: medium. Risk of being judged a corollary: medium-high for the below-line half; lower for
  the non-linear-G above-line law.

**Gap 3 — an H-free calibration for kernel tests valid on both sides of the line (H4a/H4b).**
- Why it may hold: SN + block sampling is H-free for scalar statistics including Rosenblatt limits (BTZ
  2016); Ghoshal–Shao reduce kernel tests to a 1-D SN process under mixing.
- Why not covered: the projected process under LRD is not Brownian; nobody has combined the two (found).
- Credibility: medium. **Scoop risk high** (Shao and Wu groups hold both halves).

**Gap 4 — kernel CI under LRD with a residual-vs-conditioning memory split (H7).**
- Why it may hold: CI statistics are residual products; memory in Z enters only through nuisance error
  terms (D).
- Why not covered: kernel CI theory uses mixing; the only LRD CI result is linear Gaussian (ZZW).
- Credibility: medium; reading-gap risk (econometrics journals, SSRN) not closed.

**Gap 5 — a size map of kernel tests across an H grid near the lines (evaluation).**
- Why not covered: published simulations stop at H = 0.7 (BD) or jump from 0.6 to 0.95 (ZZW/ZKW).
- Credibility: high as an evaluation gap; the track's own MMD run already fills part of it (track S).

## 9. Small entry points (cheapest first)
1. H5 hygiene fix (cite CSG 2014; switch to uncentred or per-sample chains).
2. H1c with |ξ| and with a slightly skewed G (tests the Bai–Taqqu caution at the same time).
3. Read KL results against the **linear** sum rule first (Pearson arm = Tsay–Chung territory), then
   HSIC.
4. Re-analyse KL per-rep outputs for H3 (dCov/HSIC vs squared Pearson), knowing BD already prove the
   G = id case: the new information is G2/G3.
5. SN-MMD with H-free block sampling (BTZ-style) before any H-tabulated version.
6. HSIC on Shao-style prewhitened series (H10 baseline).

## 10. Experiment ideas (for the author; not run)
- Add a "skewed |ξ|" process (e.g. |ξ + 0.2|) to the H1c rerun: rank 1 re-enters with a small
  coefficient; the break should appear later in n, not disappear.
- In any calibration rerun, add Betken-style H-tabulated SN critical values and BTZ-style SN block
  sampling side by side.

## 11. Risk notes
- **Competitor risk:** Zhai/Wu group (ZZW, ZKW, two ZCW papers in 2026) — they name the 3/4 regime as the
  next limit and already handle it for coverage. Shao group — kernel SN (GS 2026) plus SN under LRD
  (Shao 2011; BTZ coauthor T. Zhang worked with Wu). Betken/Dehling — dCov under LRD, open mixed case.
- **"Known" risk:** H0 (high), H1 as a statement (medium-high), H2 linear side (known), H3 for G = id
  (known), H4a for scalars (known), H8 symmetric (known), H10 linear (known).
- **Reading-gap risk:** SSRN (unsearched), Scholar (partial), unread classics (Hsing–Wu, Dehling–Taqqu
  1991, Arcones, Poskitt). A U/V-statistic result covering HSIC under LRD could sit in one of these.
