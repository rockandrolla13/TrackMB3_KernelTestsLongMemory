# Hypotheses and claims — kernel tests under long memory (TrackMB3, NEW-MB-KERNEL)

Execution `claude-mb3-hyp-20261016` · Agent 1 · 2026-10-16.
Status: **proposals for the author.** Nothing here is an accepted claim, a locked protocol or a
Gate decision. The author decides which items go forward.

## How to read this file

**Evidence classes.**
- **S** = simulation in this programme (exploratory; cannot support an accepted claim).
- **F / S / A** for literature = full-text section read / summarising fetch / abstract or snippet.
  To avoid a clash with simulation S, literature classes are written **lit-F, lit-S, lit-A**.
- **D** = derivation or inference by a worker (including this one). Not proved, not refereed.
- **R** = recollection of standard theory or a known paper that has **not** been re-read in this
  track. Agent 2 must verify every R item before anyone relies on it.

**Notation.** ρ(k) ~ k^{−D}, H = 1 − D/2, d = H − 1/2. D = 1/2 ⇔ H = 3/4 ⇔ d = 1/4.
For a Gaussian-subordinated series X_t = G(ξ_t), the Hermite rank m is the lowest Hermite order with
a nonzero coefficient. For statistics, "effective rank" means the Hermite rank of the leading term
of the statistic's null expansion, not the rank of G alone.

**Paths.** All paths are relative to the track root
`/media/ak/10E1026C4FA6006E/GitRepos/TrackMB3_KernelTestsLongMemory/` unless they start with `/`.
Short names used below:
- `BT` = `sources/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/boundary-theorems.md`
- `FME` = `sources/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/finance-memory-estimates.md`
- `PRE` = `sources/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/pre-emption.md`
- `SD` = `sources/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/simulation-design.md` (SHA-256 `7bb1d6b9…559f`)
- `FM`, `GM`, `DIR` = `sources/JOB-NEW-MKC-LIT-20261015/memory-kernel/{field-map,gap-map,directions}.md`
- `MMD-R` = `experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_mmd-dml/results.md`
- `MMD-T` = `experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_mmd-dml/main/tables.md`
- `MMD-C` = `experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_mmd-dml/main/summary_khom.csv`
- `KL` = `/media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness/runs/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015/kernel-linear/`

**State of the HSIC / dCov / Pearson run (KL).** Still running when this file was written
(2026-10-16, about 13:15 UTC). There is no `results.md` or `handoff.md`. `main_part2.log` showed
860 of 4,900 chunks done. The author's amendment D-a28e4cb5fa dropped n = 4000 for K-ind, so
P1, KC1 and KC2 will be judged at n = 2000 (L-ind still goes to T = 8000). I did **not** read the
partial cell files or the 20-rep pilot: the pilot is not evidence, and partial cells of an unfinished
run should not steer hypotheses. Every item that depends on KL is marked **pending (KL)**.

**No novelty claims.** "Not found" below always means "not found in the bounded scans in `sources/`".

---

## 0. Summary of the hypotheses

| ID | One line | Type | Evidence now | Priority |
|---|---|---|---|---|
| H0 | The failure line is set by the effective Hermite rank of the statistic's leading null term | theorem (framework) | D, lit-F for two instances | High (organising) |
| H1 | MMD homogeneity is rank 1, so it fails at every H > 1/2 | theorem + empirical | **S (strong)**, D | High |
| H1b | Size excess of both MMD calibrations depends on (n, H) only through (n/l)^{2H−1} | empirical → theorem | D (post-hoc re-tabulation of S) | High |
| H1c | Rank-2 data (even G) restores the 3/4 line for MMD | theorem + empirical | D only | High (cheap, decisive for H0) |
| H2 | HSIC / dCov independence: √n-Gaussian iff H_X + H_Y < 3/2 | theorem + empirical | lit-F (dCov, part), D; S **pending (KL)** | High |
| H2b | The "mixed" case (one H > 3/4, sum < 3/2) is still √n-Gaussian | theorem | D only; S pending (KL, G4 pairs) | High |
| H3 | Above the line, the HSIC null collapses to a squared latent cross-correlation | theorem | D only | High |
| H4 | Calibration above the line: rate-corrected subsampling, H-dependent self-normalisation or prewhitening can restore size; i.i.d.-rate subsampling and the wild bootstrap cannot | method | S (negative side only), D | High |
| H5 | The wild bootstrap's oversize at H = 0.5 is a multiplier-scheme artefact, not memory | empirical | S | High (cheap; contaminates every other wild-bootstrap reading) |
| H6 | The shift test is the most robust of the three HSIC calibrations below the line | empirical | D; S pending (KL) | Medium |
| H7 | Kernel CI tests obey a sum rule on the *residual* memories plus a nuisance channel | theorem + empirical | D only; lit-F for linear Gaussian analogue | Medium–High |
| H7b | Kernel Granger tests with finite conditioning sets inherit the sum rule with small constants | theorem + empirical | D only | Medium |
| H8 | Kernel conformity scores: coverage fluctuation follows the Hermite rank of the score indicator | theorem + empirical | D only | Medium (F-KERNEL owns) |
| H9 | Most finance pairs of interest sit outside the √n regime; single-series homogeneity tests are broken for all of them | empirical | lit-F/S/A estimates, D | High (decides relevance) |
| H10 | Fractional prewhitening with a plug-in d̂ restores HSIC size across H | method | D only | Medium–High |

---

## H0. The Hermite-rank map (organising claim)

**Statement.** Let a test statistic T_n have a null expansion whose leading term is a sample mean
(1/n)Σ_t F(ξ_t, η_t) of a function of independent stationary Gaussian LRD processes ξ (memory D_ξ)
and η (D_η). Write m_ξ, m_η for the Hermite orders of F's leading component in each argument.
Then the leading term is √n-Gaussian iff the autocovariance of F is summable:
- single series: m·D > 1;
- cross term of two independent series: m_ξ·D_ξ + m_η·D_η > 1.

Otherwise T_n has a non-central (Hermite / Wiener–Itô) limit at rate n^{−(Σ m D)/2}. In
particular: rank 1 → line at H = 1/2; rank 2 → H = 3/4; rank (1,1) cross → H_X + H_Y = 3/2.

**Type.** Theorem (framework). Each instance needs its own proof.

**Current evidence.**
- lit-F: Betken–Dehling 2107.03041 Thm 2.2 (dCov, subordinated Gaussian, both D > 1/2);
  Zhai–Zhong–Wu 2609.30565 (sample covariance, Gaussian linear, β > 3/4). `BT` §1–2.
- D: `BT` §3–4 and the exact variance-scaling check `sources/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/sanity_variance_scaling.out`
  (slopes 2H − 1, 4H − 3, 2H_X + 2H_Y − 3 match to three decimals).
- S: DML arm confirms rank 1 → 1/2 and rank 2 → ≈ 3/4 for the population-ATE interval
  (`MMD-R`, DML section; knots 0.538 and 0.767). MMD confirms rank 1 → 1/2 (H1).
- R: Breuer–Major, Taqqu (1979), Dobrushin–Major (1979) — standard; not re-read here.

**Confirm if.** Every family in the memory-boundary run breaks at the line its leading-term rank
predicts (P1–P3 of `SD`), and a deliberately rank-shifted variant moves the line as predicted (H1c).

**Falsify if.** A statistic with known effective rank breaks at a clearly different line, with
growth in n. Example: HSIC on equal-H pairs breaks clearly below H = 0.7 with excess rising in n
(that would be KC2 of `SD`).

**Cheapest decisive check.** H1c (one small rerun of the existing K-hom code on |ξ|). It is the
only test of H0 that *moves* the line by design rather than reading it off.

**Closest prior work and pre-emption risk.** The map itself is classical probability. Applying it
to rank-1 vs rank-2 vs cross statistics is textbook for sample means and covariances (Robinson &
Hidalgo 1997 for regression, lit-A, `PRE`). Risk that a referee calls the map "known": **high**.
Its value is as the organising lens and as a size map for kernel tests, not as a theorem.

**Priority.** High, as the frame. Low as a stand-alone contribution.

---

## H1. MMD homogeneity is a rank-1 statistic and fails at every H > 1/2

**Statement.** For one stationary series X_t = G(ξ_t), split into halves, with a bounded
translation-invariant characteristic kernel k and a G whose Hermite-1 coefficient
c₁ = E[k(·, G(ξ))ξ] ∈ RKHS is nonzero:
1. μ̂₁ − μ̂₂ = c₁ (ξ̄₁ − ξ̄₂) + o_P(n^{H−1}) in the RKHS norm (D).
2. Hence n^{2−2H}·MMD²_n → κ_H ||c₁||² χ²₁ for a constant κ_H > 0, a **rank-one** weighted χ²,
   not the infinite weighted χ² sum of the i.i.d. case (D).
3. So n·MMD²_n diverges at rate n^{2H−1}, and any calibration that assumes the i.i.d. n-rate has
   size → 1 for every H > 1/2.

Here κ_H comes from Var(ξ̄₁ − ξ̄₂) for adjacent halves of fGn, which is positive for H < 1
(the halves do not cancel).

**Type.** Theorem (parts 1–2) + empirical (part 3 at finite n).

**Current evidence.**
- **S (strong, exploratory):** `MMD-R`, `MMD-T`. At n = 4000, size at nominal 5% is 11.7–19.2% at
  H = 0.6 and 44–71% at H = 0.9 for G1 (fGn), G2 (2Φ(ξ) − 1), G3 (sign ξ), both calibrations.
  Subsampling excess grows with n for all H > 0.5. Knots 0.50–0.55, CI contains 0.5 in 3 of 6 fits
  (the grid has no point in (0.5, 0.6)). Power check fine. 0 errors in 255,000 reps; analysis rerun
  byte for byte by the coordinator (`experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015.md`).
- D: `BT` §4 ("MMD … not the same boundary"); this file, part 1–2 above.
- lit-F: Dehling–Rooch–Wendler 1404.0551 (two-sample U-statistic processes under LRD;
  non-central when m·D < 1, m = Hermite rank of the kernel). Written for change-point tests.

**Confirm if.** (a) The rank-one χ² law in part 2 fits the simulated n^{2−2H}·MMD² (e.g. ratio of
variance to squared mean ≈ 2, the χ²₁ value) at H ≥ 0.7. (b) H1b holds.

**Falsify if.** n^{2−2H}·MMD² has a clearly non-χ²₁ shape that does not approach χ²₁ with n, or if
a G with c₁ ≠ 0 holds size at some H > 1/2 for large n.

**Cheapest decisive check.** Re-analysis only, no new run, if per-rep statistics are kept:
the `mmd2` column of `experiments/…_mmd-dml/main/khom/*.csv` gives n·MMD² per rep. Compute the
ratio Var/mean² of n^{2−2H}·MMD² per cell. χ²₁ gives 2; the i.i.d. weighted-sum law gives a value
well below 2 for a Gaussian kernel. One afternoon. (I have not done this; it belongs to Agent 3 or a
later job.)

**Closest prior work and pre-emption.** Dehling–Rooch–Wendler (lit-F) covers the U-process
machinery; the MMD wording is not there. CUSUM/Wilcoxon change-point tests under LRD (Shao 2011,
Betken 2016, lit-A) are the same phenomenon for scalar statistics. Pre-emption risk for "MMD
breaks under LRD": **medium–high** (likely seen as a corollary). Risk for the rank-one χ² law with
explicit κ_H for MMD: medium; not found.

**Priority.** High. Strongest evidence in the track. It is also the clearest warning for practice:
kernel change-point and homogeneity tests on market series.

---

## H1b. One index explains both MMD calibrations' size excess: (n/l)^{2H−1}

**Statement.** In the K-hom design both calibrations compare the full statistic with a reference
built at block scale l = ⌊√n⌋ and normalised as if i.i.d. (subsampling: window statistic
l·MMD²(window halves) vs n·MMD²; wild bootstrap: Markov multipliers with correlation e^{−k/l}).
Under H1, the true scale of n·MMD² relative to the reference scale grows like (n/l)^{2H−1}.
Hypothesis: the rejection rate of each calibration is, to first order, a single increasing function
of x = (n/l)^{2H−1}, the same across H and n.

**Type.** Empirical now; theorem if the ratio-of-scales argument is made rigorous.

**Current evidence — D (post-hoc, this worker, 2026-10-16).** I re-tabulated `MMD-C`
(20 null cells per process with H > 0.5, n ∈ {500, …, 4000}) against log x. Pure arithmetic on the
existing summary; no new simulation. A two-parameter linear fit of rejection rate on log x gives
R² = 0.978–0.988 (wild) and 0.991–0.992 (subsampling) for G1–G3, with max absolute residual
2.5–5.3 pp. A model with a separate mean for each H (five parameters, ignores n) gives R² =
0.926–0.954 (wild) and 0.909–0.927 (subsampling). Example, G1 wild: H = 0.7, n = 4000 (x ≈ 5.3)
rejects 25.8%; H = 0.75, n = 500 (x ≈ 4.8) rejects 25.9%. H = 0.8, n = 4000 (x ≈ 12.1) rejects
35.6%; H = 0.9, n = 500 (x ≈ 12.2) rejects 36.8%. **This was not pre-registered**, so it is a lead,
not a test.

**Why it matters.** It explains the odd pattern in `MMD-R` (wild-bootstrap excess "flat" at H = 0.6
over n = 500–4000): x only moves from 1.87 to 2.29 there. It also predicts that **no choice of
l = n^a with a < 1 fixes either calibration**, and it gives a direct correction: rescale the
reference by (n/l)^{2Ĥ−1} with an estimated Ĥ (see H4).

**Confirm if.** A pre-registered rerun with a different block rule (l = n^{1/3} and l = n^{2/3})
falls on the same curve in x within MC error.

**Falsify if.** Cells with equal x but different (H, n, l) differ by more than ~3 pp beyond MC
error, or a rescaled reference (l·MMD² multiplied by (n/l)^{2H−1} with true H) still over-rejects
at large n.

**Cheapest decisive check.** Rerun K-hom G1 at n ∈ {1000, 4000} with l ∈ {n^{1/3}, n^{2/3}} and
H ∈ {0.6, 0.8}: 8 cells × 1,000 reps, a few CPU-hours with the existing code. Add one column:
subsampling with the oracle rescale (true H). If the oracle rescale fixes size, H4(a) is
supported.

**Prior work and pre-emption.** Subsampling under LRD with estimated rate is studied for scalar
statistics (Hall, Horowitz & Jing 1998; Jach, McElroy & Politis 2012 — **R**, verify). Applying it to
MMD and showing the collapse: not found. Risk: low–medium.

**Priority.** High. Cheap, and it turns a negative result into a calibration recipe.

---

## H1c. Rank-2 data restores the 3/4 line for MMD

**Statement.** If G is even (e.g. X = |ξ|, X = ξ², or a symmetric volatility-like transform) and k is
translation-invariant, then c₁ = 0 and the leading term of μ̂₁ − μ̂₂ is rank 2. Prediction:
MMD homogeneity is √n-calibratable (wild bootstrap or i.i.d.-rate subsampling) for H < 3/4 and
fails above, with excess growing in n only for H > 3/4.

**Type.** Theorem (corollary of H0) + empirical.

**Current evidence.** D only (`BT` §4 notes the exception in one line). No simulation.

**Confirm if.** On X = |ξ| (fGn ξ), size stays within ±2 pp of the H = 0.5 baseline for H ≤ 0.7
and breaks for H ≥ 0.8, with knot CI containing 0.75.

**Falsify if.** |ξ| breaks at H = 0.6 like G1–G3 did. That would mean the effective rank of the
MMD leading term is not the rank of G, and H0 needs revision.

**Cheapest decisive check.** Add process G5: X = |ξ| to the existing K-hom script (one line);
H ∈ {0.5, 0.6, 0.7, 0.75, 0.8, 0.9}, n ∈ {1000, 4000}, 1,000 reps. Under 5 CPU-hours.

**Finance relevance.** Absolute returns and many volatility proxies are close to even functions of a
latent driver. If H1c holds, homogeneity tests on |r_t| behave very differently from tests on
signed series (order signs, signed flow).

**Prior work and pre-emption.** Classical for scalar statistics (rank 2 → Rosenblatt above 3/4).
For MMD: not found. Risk: medium (corollary).

**Priority.** High, because it is the cheapest **falsification** test of the whole Hermite-rank
framing.

---

## H2. HSIC / dCov independence: the cross rule H_X + H_Y = 3/2

**Statement.** Let X_t = G₁(ξ_t), Y_t = G₂(η_t), ξ ⊥ η stationary Gaussian LRD, G₁ and G₂ of
Hermite rank 1. Let HSIC use bounded translation-invariant characteristic kernels (Gaussian, median
bandwidth), lag 0. Then:
1. If H_X + H_Y < 3/2: n·HSIC_n → Σ_j λ_j Z_j² (an infinite weighted χ² with weights set by the
   long-run covariance operator of the cross term), and calibrations that estimate that law
   consistently are valid.
2. If H_X + H_Y > 3/2: n·HSIC_n diverges at rate n^{2(H_X+H_Y)−3}, and n^{D_ξ+D_η}·HSIC_n has a
   non-central limit (see H3).
3. If H_X = H_Y = H this reduces to the per-series line H = 3/4.

More generally, with Hermite ranks m_X, m_Y of the relevant feature maps, the line is
m_X D_X + m_Y D_Y = 1.

**Type.** Theorem (extension of Betken–Dehling to HSIC is D) + empirical.

**Current evidence.**
- lit-F: Betken–Dehling Thm 2.2(i) for dCov with both D > 1/2 (both H < 3/4), subordinated
  Gaussian, subsampling calibration (Thm 2.4). Their own simulation shows 7–17.5% size at 5% for
  H ∈ {0.6, 0.7}, n ≤ 1000 (`BT` §1). Contrary evidence to keep: trouble starts below the line in
  finite samples.
- D: the sum rule (`BT` §1 end, §3); HSIC extension via Bochner (`BT` §4).
- D (formula check): cross slope 2H_X + 2H_Y − 3 reproduced exactly; the (0.7, 0.9) pair has slope
  0.207 (`sanity_variance_scaling.out`).
- **S: pending (KL).** Predictions P1, KC1, KC2 in `SD`. HSIC (wild, shift, subsampling), dCov
  (wild, subsampling), Pearson (NW, EWC). G1–G3 equal-H, G4 unequal pairs
  (0.6, 0.85), (0.7, 0.8), (0.7, 0.9), (0.6, 0.9), R1.

**Confirm if (pre-registered in `SD`, judged at n = 2000 after D-a28e4cb5fa).**
P1: excess ≤ 2 pp for H_X + H_Y < 1.5 and growing with n above. Knot CI contains 1.5.

**Falsify if.** KC1 (no bite: wild-bootstrap HSIC in [3.5%, 6.5%] at every H ≤ 0.9) → memory is a
scope note for HSIC. KC2 (knot CI excludes 1.5, or > 5 pp excess already at H = 0.6 with no slope
change) → drop the cross-rule headline.

**Caution for reading KL.** With n_max = 2000, the excess just above the line may be small: by the
H1b logic, the scale mis-match at H = 0.8 is only (2000/44)^{0.2} ≈ 2.1 for a cross statistic with
sum 1.6. A weak slope above 1.5 is therefore not by itself evidence against P1. Pearson at T = 8000
is the sharper test of the rule; HSIC at 2000 is the weaker one.

**Cheapest decisive check.** Wait for KL. Then check the Pearson arm (same cross rule, larger T)
first, then HSIC.

**Prior work and pre-emption.** Betken–Dehling (dCov), Betken–Dehling–Kroll 2112.14091
(dCov bootstrap, α-mixing, lit-S); Zhai–Zhong–Wu (covariance-based, Gaussian linear). ZZW Section 6
names the ultra-long-memory case as next. Risk that an HSIC version appears within a year:
**medium–high** (`PRE`, scoop watch).

**Priority.** High. Central claim of the track's question.

---

## H2b. The mixed case: one series above 3/4, sum below 3/2

**Statement.** If D_ξ < 1/2 < D_η with D_ξ + D_η > 1 (e.g. H_X = 0.85, H_Y = 0.6, sum 1.45), the
cross term still has summable autocovariance, so n·HSIC_n (and n·dcov_n) keeps the √n-type
weighted χ² limit. Betken–Dehling's restriction to both D > 1/2 is a proof artefact.

**Type.** Theorem.

**Current evidence.** D (`BT` §1: Remark 2.1(i) says the case is not covered). lit-F for the gap.
S pending (KL, G4 pair (0.6, 0.85)).

**Confirm if.** G4 (0.6, 0.85) behaves like an equal-H pair with the same sum, not like H = 0.85.

**Falsify if.** (0.6, 0.85) breaks like an H = 0.85 equal pair. That would mean the marginal
(single-series) terms of the statistic matter, not only the cross term. Possible mechanism (D):
the centring by φ̂_X φ̂_Y involves single-series empirical characteristic functions whose
fluctuations are rank 1, of order n^{H_X−1} and n^{H_Y−1}. In the cross-covariance they enter as a
product, of order n^{H_X+H_Y−2}. That is negligible against n^{−1/2} exactly when H_X + H_Y < 3/2,
so this term obeys the same sum rule (D, this worker, not checked against a written proof). The
remaining risk is technical: tightness of the empirical characteristic-function process of the
series with D < 1/2 in the weighted L² space, which is the step Betken–Dehling avoid by assuming
both D > 1/2. That step is why the mixed case is a theorem question rather than a formality.

**Cheapest decisive check.** KL G4 (0.6, 0.85) cell when it lands. Then a targeted 1-cell rerun at
n = 4000 if the n = 2000 result is ambiguous.

**Prior work and pre-emption.** Not covered by Betken–Dehling (lit-F). Robinson–Hidalgo (1997,
lit-A) for linear regression has the product-memory form. Risk: medium.

**Priority.** High. It is the case finance actually presents (signs at 0.7 with volatility at 0.9 is
outside; signs at 0.6 with spreads at 0.8 is mixed and inside, if H2b holds).

---

## H3. Above the line, the HSIC null collapses to a squared latent cross-correlation

**Statement.** Above the cross line, the leading term of the empirical cross-covariance operator is
the rank-(1,1) projection a ⊗ b · (1/n)Σ_t ξ_t η_t, with a = E[k(·, X)ξ], b = E[l(·, Y)η]. So
n^{D_ξ+D_η}·HSIC_n → ||a||² ||b||² · W², where W is the limit of n^{(D_ξ+D_η)/2 − 1}Σ ξ_t η_t
(a product-type Wiener–Itô variable; for ξ ⊥ η it is not Rosenblatt but a "two independent
Gaussian measures" double integral; Betken–Dehling Thm 2.2(ii) gives this for G = id).
Consequences:
1. Above the line, HSIC's null is a one-dimensional object: a squared latent correlation.
2. HSIC and the Pearson test (on the data, for rank-1 G) become nearly the same test under the null.
   Their per-replication statistics should have rank correlation → 1 as n grows.
3. Calibrating HSIC above the line reduces to calibrating one LRD cross-correlation. That problem
   has a self-normalised or plug-in-H answer in the linear literature (H4).

**Type.** Theorem (G non-linear is open even for dCov per Betken–Dehling) + empirical.

**Current evidence.** D only (this worker, following `BT` §1 Thm 2.2(ii)). The Hilbert-space NCLT
in Betken–Dehling covers G = id only above the line (lit-F).

**Confirm if.** In KL per-rep outputs (if the HSIC and Pearson statistics are stored per rep),
Spearman correlation between n·HSIC and (√n ρ̂)² rises with H_X + H_Y and is high (> 0.8) at
sums ≥ 1.6. Rejection decisions agree in most reps at H = 0.9.

**Falsify if.** The correlation stays low above the line, or does not rise with the sum.

**Cheapest decisive check.** Re-analysis of KL raw cells after the run finishes (no new sim),
provided the per-rep CSVs carry both statistics. If not, a 2-cell rerun (H = 0.6, 0.9, n = 2000).

**Prior work and pre-emption.** Betken–Dehling (G = id). Not found for HSIC or for non-linear G.
Risk: medium (natural next step for Betken/Dehling).

**Priority.** High. If true, it gives a clean theorem statement *and* a practical calibration route
above the line, which is where finance pairs sit (H9).

---

## H4. Calibration under long memory — what can and cannot work

Four sub-hypotheses. Each is a method claim.

### H4a. i.i.d.-rate subsampling fails above the line; rate-corrected subsampling works

**Statement.** Betken–Dehling-style subsampling with the n-rate is valid only where the true rate
is n (below the line). Above it, the size excess is a function of (n/l)^{(effective excess
exponent)} (H1b). Rescaling the window statistics by (n/l)^{2Ĥ−1} (rank 1) or
(n/l)^{2(Ĥ_X+Ĥ_Y)−3} (cross) with consistent Ĥ restores size asymptotically, provided
(Ĥ − H)·log n → 0 in probability.

**Evidence.** S (negative side): MMD subsampling with the n-rate fails at every H > 1/2 and
worsens with n (`MMD-R`). D (positive side). lit-F: Betken–Dehling's block length
l_n = O(n^{(1+min D)/2 − ε}) assumes the n-rate below the line.

**Confirm / falsify.** Oracle-H rescale fixes size (confirm); plug-in local-Whittle Ĥ leaves
> 2 pp excess at n = 4000 because of log n·(Ĥ − H) (falsify the "practical" half, not the oracle half).

**Cheapest check.** The oracle-rescale column in the H1b rerun; then a plug-in column.

**Prior work.** Hall, Horowitz & Jing (1998) and Jach, McElroy & Politis (2012) on sampling
windows / subsampling with estimated rate under LRD — **R, verify**. Bertail, Politis & Romano
(1999) rate estimation for subsampling — R. Risk that the rescaled-subsampling recipe is "known":
high for scalars, medium for kernel statistics.

### H4b. Self-normalisation: pivotal below the line, H-dependent above

**Statement.** Use the Ghoshal–Shao construction (sample split; project the test half onto an
RKHS direction estimated from the training half; self-normalise the resulting 1-D partial-sum
process). Then:
1. Below the line (cross) and for rank ≥ 2 statistics with mD > 1, the projected process is a
   short-memory-like sum and the usual SN limit (Brownian functional, Shao 2010) applies. **Same
   critical values as under β-mixing.**
2. For MMD homogeneity (rank 1), the projected process converges to fBm with index H, and the SN
   ratio has a limit that depends on H only. Critical values can be tabulated in H (as in Shao 2011
   for the change-in-mean statistic) and used with a plug-in Ĥ. The dependence on H is mild, so
   size error from Ĥ error is small.
3. Above the cross line, the SN limit depends on (H_X + H_Y) through a Hermite-type process with
   self-similarity 2 − (D_ξ + D_η)/2. Using β-mixing critical values gives a **bounded** size
   distortion that does not grow with n (unlike wild bootstrap and n-rate subsampling).

Point 3 is the practically important one: SN fails *gracefully*.

**Evidence.** lit-S/A: Ghoshal & Shao, arXiv:2601.13468 ("Resampling-free inference for time
series via RKHS embedding", Jan 2026; abstract: pivotal limits under absolute regularity — β-mixing —
and mild moments). Found by this worker's WebSearch, abstract level only (lit-A).
lit-A: Shao (2011, JTSA) SN change-in-mean under LRD; Betken (2016) SN Wilcoxon under LRD (`FM`
L1; `BT`). **R:** the claim that Shao's LRD SN limit depends on H, with tabulated critical values, is
my recollection; verify. D: points 1–3.

**Confirm if.** Inside the D1 grid, SN-HSIC holds 5% ± 1.5 pp for H_X + H_Y < 1.5 and its
excess above the line is flat in n. SN-MMD with H-tabulated critical values holds size at all H.

**Falsify if.** SN-HSIC is oversized below the line (point 1 false), or its excess grows with n
above it (point 3 false: the SN denominator does not absorb the rate).

**Cheapest check.** Implement SN-MMD homogeneity (simplest case: one direction, 1-D process) on G1
at H ∈ {0.5, 0.7, 0.9}, n ∈ {1000, 4000}; compare β-mixing critical values with simulated
fBm-SN critical values at the true H. Few CPU-hours. Do this before any SN-HSIC work (`DIR` D2).

**Prior work and pre-emption.** Ghoshal–Shao (β-mixing). Shao's LRD SN (scalar). Not found for a
kernel statistic under LRD (`FM` L1). **Pre-emption risk: high** — Shao's group has both halves.

### H4c. Block bootstrap: valid for the cross term below the line, not above

**Statement.** Moving-block bootstrap of pairs (X_t, Y_t) with l → ∞, l/n → 0 consistently
estimates the null law of n·HSIC below the cross line (summable cross covariance), and fails above
it (it cannot reproduce a non-central limit's scaling). For rank-1 statistics (MMD) it fails for all
H > 1/2. Note the bootstrap must break the X–Y link to impose the null (resample blocks of X and of
Y independently), which keeps each series' within-block memory.

**Evidence.** lit-F: Zhai–Zhong–Wu Thm 2 (overlapping block bootstrap for covariance-based
statistics, β > 3/4). R: Lahiri (1993) and Kim & Nordman (2011) on MBB under LRD — I recall MBB
can be valid for Gaussian LRD sample means with proper normalisation but fails for non-Gaussian
(rank > 1) limits. **Verify.** D otherwise.

**Cheapest check.** Not in the current design. Add an independent-block-resampling calibration to a
small HSIC rerun (H ∈ {0.6, 0.8}, n = 2000).

**Priority.** Medium (likely subsumed by H4a/H4b).

### H4d. The wild bootstrap is invalid at any long memory for rank-1 statistics, and above the line for cross statistics; below the line it may hold despite failing the published mixing condition

**Statement.** Chwialkowski–Sejdinovic–Gretton's conditions (τ(r) = O(r^{−6−ε}) for HSIC) fail for
every H > 1/2. They are sufficient, not necessary. For a cross statistic with H_X + H_Y < 3/2 the
degenerate kernel's cross-covariance is summable, and a block-Markov multiplier with l → ∞ can still
capture its long-run covariance operator. So wild-bootstrap HSIC may hold size below the line. For
rank-1 MMD it cannot (S: `MMD-R`).

**Evidence.** lit-F: `FM` L1 (exact τ conditions). S: MMD wild bootstrap fails (rank 1).
**Pending (KL)** for HSIC.

**Confirm / falsify.** KL HSIC-wild within ±2 pp of its H = 0.5 level for H_X + H_Y < 1.5
(confirm); clear growth in n below 1.5 (falsify).

**Priority.** High; it is the most-used calibration in practice.

---

## H5. The wild bootstrap's memory-free oversize is a multiplier-scheme artefact

**Statement.** The K-hom wild bootstrap is oversized by about +2 to +4.5 pp at H = 0.5 (i.i.d. G1,
G2) at every n from 500 to 4000, and the oversize does not shrink with n. Hypothesis: this is a
finite-sample defect of the specific scheme — one Markov-Rademacher chain run across both halves
with correlation e^{−k/l}, l = ⌊√n⌋, then centred within each half, applied to the biased
(V-statistic) MMD² with the median bandwidth computed on the pooled sample — not a property of the
data. Candidate causes, in my order of suspicion:
1. Within-half centring of strongly autocorrelated multipliers removes a large share of their
   variance (a chain with l = √n has few effective independent values per half), so the bootstrap
   distribution is too narrow.
2. The V-statistic diagonal terms are treated differently in the observed and bootstrap
   statistics.
3. The chain crosses the half boundary, so the two halves' multipliers are correlated.

**Type.** Empirical (implementation and finite-sample behaviour).

**Current evidence.** S: `MMD-T` wild rows at H = 0.5: G1 +2.9, +2.5, +3.1, +2.2 pp; G2 +4.5,
+3.4, +1.3, +2.8 pp; G3 +4.2, +2.0, +1.9, +1.2 pp (n = 500…4000). Subsampling at H = 0.5 is on size.
D for the causes. The kernel-linear worker's HSIC wild bootstrap uses the same multiplier law
(`MMD-R` implementation notes) — **pending (KL)** whether HSIC shows the same oversize at H = 0.5.

**Why it matters.** Every wild-bootstrap number in this track is the sum of a memory effect and this
artefact. Until the artefact is removed, the memory effect for the wild bootstrap is over-stated by
2–4 pp, and the "2 pp" bands in P1/KC1 are partly consumed by it.

**Confirm if.** At H = 0.5, one of these variants brings size to 5% ± 1 pp: (i) l = 1 (i.i.d.
Rademacher); (ii) separate chains per half; (iii) no within-half centring (use the
Chwialkowski et al. two-sample form exactly as published); (iv) U-statistic instead of V-statistic.

**Falsify if.** All four variants keep the oversize. Then the cause is elsewhere (bandwidth
selection on the pooled sample, or a coding error), and Agent 3 should audit the code.

**Cheapest decisive check.** 4 variants × G1 H = 0.5 × n ∈ {500, 2000} × 2,000 reps. Under one
CPU-hour. **Should run before any further wild-bootstrap comparisons.**

**Prior work.** Chwialkowski, Sejdinovic & Gretton (2014, lit-F) report their own null sizes;
Leucht & Neumann (2013) dependent wild bootstrap for degenerate U/V-statistics — **R, verify**.
Not a novelty item; a hygiene item.

**Priority.** High (cheap and it contaminates comparisons).

---

## H6. The shift test is the most robust HSIC calibration below the line

**Statement.** The circular-shift test (Chwialkowski–Gretton 2014) keeps each series' full
dependence structure intact and only breaks the X–Y alignment. Prediction: below the cross line it
holds size at least as well as the wild bootstrap and subsampling. Above the line it fails because
shifted statistics are strongly correlated with each other and with the observed statistic: the
effective number of independent shifts falls like (n/A)^{…} as memory rises, so the empirical
quantile is too narrow.

**Type.** Empirical.

**Evidence.** D only. KL runs HSIC-shift with shift range [0.1n, 0.9n] — **pending (KL)**.

**Confirm / falsify.** Shift-test excess ≤ wild-bootstrap excess in every cell with sum < 1.5
(confirm); shift test worse than wild below the line (falsify).

**Priority.** Medium. Useful for practice; low theory value.

---

## H7. Kernel conditional-independence tests: a sum rule on residual memories, plus a nuisance channel

**Statement.** Consider a residual-type kernel CI test (GCM-style or KCI with regression
residuals) of X ⊥ Y | Z with X_t = f(Z_t) + ε^X_t, Y_t = g(Z_t) + ε^Y_t. The null statistic
decomposes into
- (i) a cross term in (ε^X, ε^Y) → sum rule on the **residual** memories: H_{εX} + H_{εY} < 3/2;
- (ii) nuisance-error products (f̂ − f)(Z_t)(ĝ − g)(Z_t), whose fluctuation inherits the memory of
  Z through the Hermite rank of the error functions; and
- (iii) mixed terms (f̂ − f)(Z_t)·ε^Y_t, a cross term in (Z, ε^Y) → sum rule on (H_Z, H_{εY}).

Prediction: memory that sits only in the conditioning variable Z (a common long-memory driver)
does **not** break the test if the residuals are short-memory and the nuisance rates are fast enough.
Memory that survives in the residuals does, by the sum rule.

**Finance reading (D).** Conditioning on a long-memory common factor (volatility, activity) can
move a pair from outside the line to inside it. Example: order signs vs trade durations conditional
on a volatility proxy.

**Type.** Theorem + empirical.

**Evidence.** D only. lit-F: Zhai–Zhong–Wu (Gaussian graphical CI, linear, β > 3/4) is the only
long-memory CI result found (`FM` L2, `GM` K3). E1 (R) found no kernel CI test under LRD.

**Confirm if.** In a simulation with Z = fGn(H = 0.9) and ε^X, ε^Y i.i.d., a GCM/KCI test holds
size; with ε^X, ε^Y fGn of H = 0.8 each, it breaks.

**Falsify if.** Memory in Z alone breaks the test even with oracle nuisances (then term (ii)/(iii) is
not the mechanism, or the residual is not what the statistic sees).

**Cheapest decisive check.** GCM with oracle and kernel-ridge nuisances, 4 cells
(H_Z ∈ {0.5, 0.9} × H_ε ∈ {0.5, 0.8}), n = 2000, 1,000 reps. A few CPU-hours.

**Prior work and pre-emption.** Shah & Peters GCM (i.i.d.); Shi & Shah GTCM, Tian–Waterbury–Düker
(time series, tier-1 mixing; E1, R); Zhai–Zhong–Wu and Zhai–Kolar–Wu (IEEE TIT 2026, cited in ZZW,
not verified) for linear LRD. Risk: medium–high (ZZW group).

**Priority.** Medium–High. Closest to finance practice (questions are conditional), and the "memory
in Z is harmless" half is a clear, testable message.

## H7b. Kernel Granger tests with finite conditioning sets inherit the sum rule, with small constants

**Statement.** Conditioning on p lags of a long-memory series removes most short-range dependence
but leaves residuals with the same power-law exponent at large lags and a small amplitude
(an AR(p) approximation to ARFIMA(0, d, 0) leaves residual autocorrelations ~ c_p k^{−(1−2d)} with
c_p small and decreasing in p). Prediction: a kernel Granger test conditioning on p lags looks
valid at moderate n and breaks only at large n, with the break n growing in p.

**Type.** Theorem (rate) + empirical (finite-n).

**Evidence.** D only. R: AR(p) approximation of long memory (standard; not re-read).

**Cheapest decisive check.** Exact computation (no simulation): residual autocovariance of the
best AR(p) predictor of ARFIMA(0, d, 0) for p ∈ {1, 5, 20}, d ∈ {0.2, 0.4}, then plug into the
cross-term variance formula to find the n at which the excess passes 2 pp.

**Priority.** Medium. Explains why short-memory Granger tests can look fine on finance data.

---

## H8. Kernel conformity scores under long memory

**Statement.** For split-conformal intervals with score s_t = S(X_t, Y_t) on a stationary
ergodic long-memory process:
1. Marginal coverage holds asymptotically (ergodic theorem). (Zhai–Cheng–Wu show this for residual
   / CQR scores under FDM, including long-memory Gaussian linear processes; S2, R.)
2. The fluctuation of realised coverage over a calibration window of size n (and over a test window
   of size W) is governed by the Hermite rank m of the indicator 1{s_t ≤ q}: scale n^{−1/2} if
   m·D > 1, scale n^{−mD/2} otherwise.
3. Symmetric scores of symmetric subordinated processes (|residual|, a symmetric kernel-density
   score) have m = 2 → the 3/4 line. Signed or one-sided scores (one-sided CQR, quantile scores on
   skewed series) have m = 1 → any memory widens the coverage fluctuation.
4. Kernel-weighted scores (KOWCPI-type localisation) add the memory of the weights; their effective
   rank is the rank of the weighted indicator.

**Type.** Theorem + empirical.

**Evidence.** D only. lit-F: KOWCPI needs α(τ) = O(τ^{−(2+δ)}) (`FM` L6). `GM` K9 (open,
possibly a reading gap: Wu / Politis groups).

**Confirm if.** On fGn with known H, the SD of realised coverage over W = 500 test points scales as
predicted, with a 3/4 break for |residual| scores and a 1/2 break for one-sided scores.

**Falsify if.** Both score types show the same scaling.

**Cheapest decisive check.** Pure computation for Gaussian fGn: the variance of the coverage
indicator average is an explicit function of the bivariate normal orthant probabilities and the
fGn autocovariance. No simulation needed for the Gaussian case.

**Prior work and pre-emption.** Zhai–Cheng–Wu (non-kernel scores, long memory), Chernozhukov–
Wüthrich–Zhu (ergodic conformal), Oliveira et al. (β-mixing split conformal) — R. Risk: medium.
**F-KERNEL owns** this line; this track should only hand it over.

**Priority.** Medium (ownership elsewhere).

---

## H9. Where financial series sit — what they need

**Statement (empirical, with uncertainty).** Using published estimates (`FME`):
- Order signs: H ≈ 0.70 (LSE), 0.65–0.90 (Paris, weakest row), ≈ 0.75 in metaorder theory.
  Binary sign data are rank 1 (G3 in the simulation).
- Volatility (log RV): d ≈ 0.3–0.45, H ≈ 0.8–0.95; or rough with H ≈ 0.1 and slow mean reversion
  (Gatheral–Jaisson–Rosenbaum). Close to an even function of a latent driver if one uses |r|.
- Durations, trade counts, volume, spreads: H ≈ 0.75–0.96 (durations above 0.75 in all ten
  Deo–Hsieh–Hurvich stocks, lit-F).

Claims:
1. **Every single-series MMD homogeneity or kernel change-point test on these series is outside
   the √n regime** (rank 1; H1). Calibrations at the i.i.d. rate will over-reject, by 8–15 pp at
   H ≈ 0.6–0.7 and 20–40 pp at H ≈ 0.8–0.9 at n = 500–4000 (`MMD-T`).
2. For independence tests, **pairs** matter (H2):
   - signs × signs (0.7 + 0.7 = 1.4): inside, but close; finite-sample oversize like
     Betken–Dehling's (10–17% at n ≤ 1000) is to be expected;
   - signs × volatility (0.7 + 0.9 = 1.6): outside;
   - durations × volatility, durations × volume (≥ 1.6): outside;
   - signs (0.6) × spread (0.8) = 1.4: inside only if H2b holds.
3. Because estimates move by 0.1–0.3 with the estimator (`FME` point 5), an "inside" verdict
   should use an **upper confidence bound** on H_X + H_Y, not a point estimate. With that rule,
   almost no activity-related pair qualifies as inside.
4. The rough-volatility arm (R1, λ = 0.01) mimics strong memory at n = 500–1000 (wild +26.9 pp and
   subsampling +51.1 pp at n = 500) and fades by n = 4000 (`MMD-T`). Daily RV samples
   (250–5,000 days) are in the danger zone whether volatility is long-memory or rough.

**Type.** Empirical.

**Evidence.** lit-F/S/A estimates (`FME`); S for MMD; D for pair arithmetic; HSIC **pending (KL)**.

**Confirm if.** Semiparametric H estimates with confidence intervals on the programme's own data
(marketdb order book, TRACE) place most pairs of interest above 1.5.

**Falsify if.** Own-data estimates put order signs and activity series robustly below 0.7 and
volatility below 0.75 after adjusting for intraday seasonality and level shifts (spurious long
memory). Then Betken–Dehling's Remark 2.1(ii) ("real data typically H < 0.75") would stand.

**Cheapest decisive check.** Needs data (not allowed in this job). A later job: local-Whittle
H with CIs for signs, |r|, log durations and spreads on 5–10 marketdb symbols, plus a
level-shift-robust estimator to separate true from spurious memory.

**Prior work.** `FME` table. Tóth et al. 2015 order-sign exponents (located, not read). Risk: low;
this is measurement.

**Priority.** High, because it decides whether the theory items matter for finance or only for
the boundary literature (`GM` C1).

---

## H10. Fractional prewhitening restores HSIC size

**Statement.** Estimate d for each series by local Whittle (bandwidth m), apply the fractional
difference (1 − L)^{d̂} to each series separately, and run HSIC (any calibration valid under
short memory) on the filtered pair. Because the filters are applied separately and are invertible,
whole-process independence of X and Y is equivalent to independence of the filtered series. The
residual memory d − d̂ = O_P(m^{−1/2}); its effect on the statistic is negligible if
m^{−1/2} log n → 0. Prediction: size within ±2 pp across H ∈ [0.5, 0.9] at n = 2000.

Caveats:
1. Lag-0 HSIC on filtered series tests a different **pairwise** relation than lag-0 HSIC on the raw
   series. The whole-process null is the same; power against specific alternatives changes.
2. For MMD **homogeneity** the trick does not transfer cleanly. Filtering changes the marginals, so
   homogeneity of filtered halves is a different hypothesis except in linear Gaussian models.
3. Non-linear subordination (G not linear) means (1 − L)^d G(ξ) is not short-memory in general.
   The residual is short-memory only in its rank-1 part. For rank-1 G that is enough for the
   leading term.

**Type.** Method.

**Evidence.** D only. R: prewhitening before cross-correlation tests (Haugh 1976) is classical for
ARMA; I do not know a long-memory kernel version.

**Confirm if.** Prewhitened HSIC (wild bootstrap) holds 5% ± 2 pp on G1 and G3 for H up to 0.9 and
on G4 pairs.

**Falsify if.** Prewhitened HSIC remains oversized at H ≥ 0.8 (d̂ error too large at n = 2000), or
is fine for G1 but not for G3 (non-linear subordination defeats the filter).

**Cheapest decisive check.** Add a prewhitening step to the KL code on a copy (not in the running
directory): G1 and G3 at H ∈ {0.6, 0.9}, n = 2000, 1,000 reps. A few CPU-hours.

**Prior work and pre-emption.** Prewhitening in linear cross-correlation testing (R). Fractional
filtering before testing is common in long-memory econometrics (R). Risk that a referee calls it
obvious: high. Value: the most practical fix, and a baseline any new method must beat.

**Priority.** Medium–High. It is the "dumb baseline" that H4b must beat.

---

## Priority order (worker view, for the author)

1. **H5** (wild-bootstrap artefact). Under one CPU-hour. Fix before reading any more wild numbers.
2. **H1c** (|ξ| restores 3/4 for MMD). Under 5 CPU-hours. Cheapest falsification of the whole
   Hermite-rank frame.
3. **H2 / H2b / H4d** — wait for KL; read Pearson (T = 8000) first, then HSIC at n = 2000.
4. **H1b + H4a** (collapse and rate-corrected subsampling). A few CPU-hours; turns the negative MMD
   result into a recipe.
5. **H3** — re-analysis of KL per-rep outputs if stored.
6. **H4b** (SN-MMD first, then SN-HSIC). The likely theory contribution; highest pre-emption risk.
7. **H10** (prewhitening baseline), **H7** (kernel CI), **H9** (own-data H estimates; needs a data job).
8. **H7b**, **H8** (hand H8 to F-KERNEL).

---

## Hand-off for Agent 2 (literature review)

Leads named by the coordinator: Betken–Dehling, Zhai–Zhong–Wu, Shao self-normalisation,
Ghoshal–Shao. Each question below says what answer would change a hypothesis. Please label every
source lit-F / lit-S / lit-A and verify identifiers. Items I marked **R** in this file are
unverified recollections; please confirm or correct each.

### For H1, H1b, H1c (MMD / two-sample / homogeneity under LRD)
1. Is there any paper that derives the null limit of MMD (or an RKHS two-sample / kernel
   change-point statistic) under Gaussian-subordinated LRD, and states a rank-one χ² limit or the
   n^{2H−1} divergence? Search Dehling–Rooch–Wendler (1404.0551) citers; "kernel change-point long
   memory"; "energy distance long-range dependence"; Betken's later papers.
2. Has subsampling with an estimated rate been applied to any U/V-statistic under LRD? Verify Hall,
   Horowitz & Jing (1998), Jach, McElroy & Politis (2012), Bertail–Politis–Romano (1999), and
   Betken–Wendler (subsampling for LRD, `betken2018subsampling`, lit-A so far).
3. Does any paper show that an even transform (|X|, X²) changes the size behaviour of a kernel
   two-sample or change-point test under LRD (H1c pre-emption)?

### For H2, H2b, H3 (HSIC / dCov cross rule)
4. **Betken–Dehling.** Read 2107.03041 v3 in full for: (a) whether the mixed case (Remark 2.1(i))
   is addressed anywhere in the supplement or later work; (b) the exact form of the non-central limit
   in Thm 2.2(ii) and whether it reduces to a squared latent cross-correlation (H3); (c) whether
   they discuss HSIC with general kernels. Also Betken–Dehling–Kroll 2112.14091 (bootstrap dCov) and
   `betken2024bootstrap` in the bibliography: do any of them cover LRD?
5. Search for HSIC / kernel independence tests under long memory published 2024–2026, including
   arXiv:2604.28104 ("Kernel-based independence and mean independence tests for weakly dependent
   data", seen at snippet level by this worker; check whether it treats LRD).
6. Bai & Taqqu on multivariate Hermite rank and the "instability" of Hermite rank under linear
   transforms (R): does it affect the rank-(1,1) claim when X and Y are built from dependent latent
   Gaussians (relevant if the null is conditional, H7)?
7. Robinson & Hidalgo (1997, AoS) and later work on the product-memory rule d_x + d_u < 1/2: get the
   exact statement (full text) as the linear prior for H2.

### For H4 (calibration)
8. **Shao self-normalisation.** Verify: Shao (2010, JRSSB) SN for time series; Shao (2011, JTSA)
   SN change-in-mean under LRD — does its limit depend on H, and are critical values tabulated in H?
   Shao (2015, JASA review) — any kernel or U-statistic SN under LRD? Betken (2016) SN Wilcoxon under
   LRD: what exactly is pivotal?
9. **Ghoshal–Shao.** arXiv:2601.13468 (Jan 2026) and the EJS 2026 paper cited in Track E (E1, R):
   are these the same work? Exact assumptions (β-mixing rate?), whether the projected process is
   handled by a Brownian FCLT, and whether they discuss long memory or fractional limits at all.
   This decides H4b's pre-emption risk.
10. Block bootstrap under LRD: verify Lahiri (1993), Kim & Nordman (2011) and any 2020–2026 work
    on block bootstrap for degenerate U/V-statistics under LRD.
11. Wild / dependent wild bootstrap: Leucht & Neumann (2013); Shao (2010, dependent wild bootstrap);
    any reported null oversize at independence for block-Markov multipliers on V-statistics (H5).
12. **Zhai–Zhong–Wu.** 2609.30565 Section 6 and Zhai–Kolar–Wu (IEEE TIT 2026, 72(6):4246–4296,
    cited in ZZW; verify): is the "ultra-long memory" Rosenblatt case announced as forthcoming?
    Any kernel (non-linear) extension? This is the biggest scoop risk for H2 and H7.

### For H7, H7b (kernel CI / Granger)
13. Any conditional-independence or Granger test (kernel or not) with a stated validity condition in
    terms of the memory of residuals vs the memory of the conditioning variable. Search GCM / KCI /
    "conditional independence long memory", "Granger causality fractional cointegration test size".
14. Results on residual autocorrelation of finite-order AR approximations to ARFIMA (H7b),
    e.g. Poskitt; Hurvich; "AR sieve long memory".

### For H8 (conformity scores; hand to F-KERNEL)
15. Verify Zhai–Cheng–Wu (conformal under FDM including long memory): what scores, what rate of
    coverage error, does the rate depend on the score's Hermite rank?

### For H9 (finance)
16. Order-sign exponents after 2010 (Tóth et al. 2015; crypto and futures), duration and spread
    memory estimates with CIs; spurious-long-memory tests applied to these series.
17. SSRN and Google Scholar were not effectively searched in earlier jobs (`PRE`). Please close that
    gap for questions 1, 5, 8, 9 and 12 first.

### For H10 (prewhitening)
18. Fractional prewhitening before independence / cross-correlation tests under long memory:
    any size study? Search "prewhitening long memory cross-correlation test", "fractionally
    differenced independence test".

---

## Provenance

- Inputs read (read-only): `README.md`, `FINDINGS.md`, `literature/literature-plan.yaml`, all files
  named in the path list above, `experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015.md`,
  `experiments/…_mmd-dml/khom.py` (calibration definitions), and the KL directory listing, job
  events and `mbcore.py` header (to confirm calibration definitions). KL partial cells and pilot were
  not read for results.
- Skills read for structure (read-only): `proof-methodology` SKILL.md, SHA-256
  `232534b51cb338399bd44cedda5bbd37e5fc4762b6810ddd6419139681146a6b` (matches expected);
  `quasi-theory-assistant` SKILL.md, SHA-256
  `f8460b76100a70678fe8926b9ff87fe28142f44e4d73df997437898e14aeba58` (no expected hash supplied).
- Web lookups by this worker (3 WebSearch queries, results lit-A): Ghoshal & Shao arXiv:2601.13468;
  arXiv:2604.28104 (title only); Dehling–Rooch–Wendler 1404.0551 and Betken's SN Wilcoxon paper
  (snippets). No full text was fetched.
- Computation: one arithmetic re-tabulation of `MMD-C` for H1b (printed, not saved). No simulation,
  no data, no commits.
