# First exploratory experiments — TrackMB3 (proposal, for the author)

Execution `claude-mb3-scope-20261016` · 2026-10-16. Nothing here has been run. Exploratory only.

## Proposal (200 words or fewer)

None can support an accepted claim. The HSIC/dCov run is unfinished; experiment 3 waits for it.

1. **Placebo pairs on real data.** Purpose: does memory break kernel tests on real flow? Setting: pair LQD trade signs from one stretch of days with HYG or JNK signs from a different stretch. Use data up to 2025-05-29 only. Such pairs are independent but keep real memory. Run HSIC and dCov on many pairs; estimate H, with a band, per series. Changes the plan if: rejection stays near 5%. Then memory is a footnote here, not the headline.

2. **Can MMD be rescued?** Purpose: test the Hermite-rank mechanism, and test a fix. Setting: simulated series, H from 0.5 to 0.9. Add a squared series, which should break only above H = 0.75. Add memory-aware calibrations. Changes the plan if: the squared series holds size to 0.75 (mechanism confirmed), or a calibration holds size to H = 0.8 (a method paper).

3. **Finer grid at the sum line.** Purpose: locate the HSIC break. Setting: H_X + H_Y from 1.3 to 1.7. Run only if pending results leave the break unclear. Changes the plan if: excess appears below 1.4. Then drop the 3/4 headline (KC2).

---

## Appendix — details (not counted)

### State of inputs at writing (2026-10-16, ~14:15 UTC)

- MMD (K-hom) and DML: finished. MMD breaks for every H > 0.5 on G1–G3. The wild bootstrap is
  also oversized by +2 to +4.5 pp at H = 0.5. R1 (rough vol, λ = 0.01) mimics strong memory at
  n = 500 and fades by n = 4000.
- HSIC/dCov/Pearson (kernel-linear): **pending.** Part 1 finished (12,250 chunks). Part 2 had
  logged about 860 of 4,900 chunks. Per decision D-a28e4cb5fa, K-ind n = 4000 was dropped, and
  P1, KC1 and KC2 are judged at n = 2000. I did not read the partial cell files. Partial cells
  are not a result.
- Real data: local credit ETFs LQD/HYG/JNK, 502 contiguous days 2023-05-30 → 2025-05-29
  (`db_book` family, Arca only, aggressor side missing on 18–31% of trades). **The hold-out
  (2025-05-30 onward, including MIDAS rows) is not touched by anything below.**

### E1 — Placebo pairs on real data (credit ETFs, pre-hold-out)

- **Series.** (a) Trade signs, in event time (trade index, not clock time). Use the aggressor side
  where present. Otherwise sign by the quote rule, and report the share signed each way.
  (b) Log trade durations, diurnally adjusted. These are the "activity" series that the literature
  puts at H ≈ 0.76–0.96.
- **Placebo construction.** Take a window of n consecutive events from ETF A starting in day
  block i. Take a window from ETF B starting in a day block j, with |i − j| ≥ 20 trading days.
  Random start offsets in event time. Aligning in event time with offsets avoids a shared
  time-of-day pattern, which would be a fake dependence. Pairings: LQD–HYG, LQD–JNK, HYG–JNK.
  Also a same-ETF, distant-block control (LQD–LQD).
- **n.** {500, 1000, 2000}, capped by the events available. First step is a feasibility count of
  Arca trades per day for each ETF. If n = 2000 needs many days per window, drop it and say so.
- **Tests.** HSIC (Gaussian kernel, median heuristic) with the wild bootstrap and with subsampling.
  dCov with subsampling. Pearson with HAC NW and EWC as a linear reference. Use the same code as
  kernel-linear (`mbcore.py`) so real and simulated numbers compare.
- **Memory estimates.** Local Whittle and GPH at m = n^{0.5}, n^{0.65}, n^{0.8}. Report the
  spread across bandwidths as the band. Place each pairing on the H_X + H_Y map.
- **Replications.** 500 placebo pairs per (pairing, series, n). MC s.e. ≈ 1 pp at 5%. Pairs
  reuse days, so they are not fully independent. Report a day-clustered interval.
- **Reads.** Rejection near 5% across pairings → memory is a scope note (close to KC1). Rejection
  well above 5%, rising with n and with estimated H_X + H_Y → the real-data hook for the paper.
- **Extension (not proposed now).** Tick Size Pilot control-arm stocks. That needs agreement with
  B2-HTE and with the NEW-MKC reframe question 1, which plans a control-only check of its own.
- **Compute.** HSIC wild, from the kernel-linear pilot under load: about 0.16 / 0.71 / 2.6 CPU-s
  per test at n = 500 / 1000 / 2000. 4 pairings × 2 series × 500 pairs × 3.5 s ≈ 4 CPU-h. With
  subsampling, dCov and Pearson, about 8 CPU-h. Data extraction is under 1 h. **Cap: 12 CPU-h.**
- **Governance.** Read-only marketdb queries, each bounded by `<= 2025-05-29`, logged in a
  queries audit as in the reframe job.

### E2 — Rank-2 MMD and memory-aware calibration (simulation)

- **Processes.** G1 (fGn) and G3 (sign) as before, as rank-1 references. New: G5, X = ξ² − 1
  (an even transform, so Hermite rank 2; squared-return-like). Prediction: MMD on G5 stays on size
  for H ≤ 0.7 and breaks above 0.75.
- **H grid.** {0.5, 0.55, 0.6, 0.7, 0.75, 0.8, 0.9}. 0.55 is added to resolve the rank-1 knot,
  which the first run put at 0.50–0.56 between grid points.
- **Calibrations.** (a) Wild bootstrap as before, plus multiplier lengths l = ⌊n^{0.3}⌋ and
  ⌊n^{0.4}⌋. This targets the +2–4.5 pp oversize at H = 0.5. (b) Subsampling with block length
  l = ⌊n^{(1+D̂)/2 − 0.05}⌋, where D̂ = 2 − 2Ĥ comes from local Whittle on the series. This is the
  Betken–Dehling rate. (c) Existing subsampling, l = ⌊√n⌋, as a reference.
- **n.** {500, 1000, 2000}. 2,000 reps per cell. MC s.e. ≈ 0.5 pp.
- **Reads.** G5 holds size to 0.75 → the Hermite-rank map is supported for MMD as well as DML.
  If a calibration keeps size within 2 pp up to H = 0.8 → a "fix" contribution. If none does →
  the track is a breakdown map only.
- **Compute.** The first MMD run cost 18.9 CPU-h for 84 cells up to n = 4000. Dropping n = 4000
  and adding the cells and calibrations above: about 3 processes × 7 H × 3 n = 63 cells, roughly
  3× the bootstrap work per cell. Estimate 12–18 CPU-h. **Cap: 24 CPU-h.** Run a 20-rep timing
  pilot first (not evidence).

### E3 — Fine grid at the cross-rule line (conditional)

- **Trigger.** Run only if the kernel-linear analysis gives a knot CI wider than ±0.1 around 1.5,
  or if P1 is mixed at n = 2000. If kernel-linear is decisive either way, skip E3.
- **Pairs (H_X, H_Y).** Equal: (0.65, 0.65), (0.7, 0.7), (0.725, 0.725), (0.75, 0.75),
  (0.775, 0.775), (0.8, 0.8), (0.85, 0.85). Sums 1.3–1.7. Unequal at the same sums:
  (0.55, 0.85), (0.6, 0.9), (0.65, 0.9). These separate the sum rule from Betken–Dehling's
  per-series condition (both H < 3/4).
- **Tests.** HSIC (wild bootstrap, subsampling), dCov (subsampling). n ∈ {500, 1000, 2000}.
  2,000 reps per cell.
- **Reads.** Excess > 2 pp at sums ≤ 1.4 → KC2: the headline becomes "fails at any memory".
  Unequal pairs at sum < 1.5 on size, with one H > 0.75 → the sum rule, not the per-series rule,
  is the true boundary. That would be a new result beyond Betken–Dehling.
- **Compute.** 10 pairs × 2,000 reps × ≈3.5 CPU-s ≈ 19 CPU-h for HSIC wild. About 30 CPU-h with
  all tests. **Cap: 48 CPU-h** (the design's per-family cap).

### Total

E1 + E2 about 20–26 CPU-h. E3 adds about 30 CPU-h only if triggered. The kernel-linear pilot
showed n = 4000 too costly under shared load, so no experiment here uses n > 2000.
