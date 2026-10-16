# Handoff: MMD and DML under long memory

Exploratory simulation. It is not evidence for an accepted claim. Full detail is in results.md,
and every cell is in main/tables.md.

**MMD breaks as soon as there is any long memory.** We compared the first half of a series with
its second half. At H = 0.6 the 5% test already rejects 12–19% of the time. At H = 0.9 it rejects
44–71% of the time. This holds for Gaussian, bounded and sign series, and for both calibrations.
For subsampling the problem gets worse as n grows. For the wild bootstrap it gets worse with n from
H = 0.7 upward. There is no 3/4 line for this test. One caveat: the wild bootstrap is oversized by
about 2–4 points even with no memory. That is a calibration flaw, not a memory effect.

**The DML population-ATE interval breaks for any long memory when the effect is monotone.**
- With τ = 1 + X, 95% coverage at T = 8000 falls as H rises: 84% at H = 0.6, 54% at 0.7 and
  12% at 0.9. More data makes it worse.
- With an even effect (τ = X²), the break comes only near H = 3/4. Coverage is fine up to
  H = 0.7. It is 88% at H = 0.75, 76% at 0.8 and 26% at 0.9.
- With a constant effect, coverage is fine at every H.
- Learning the nuisances, and blocked versus random cross-fitting folds, make almost no
  difference. The whole failure comes from the slowly moving state.

**The path-conditional ATE stays valid.** The same interval covers the in-sample average effect
at least 93.45% of the time in every cell. For heterogeneous effects it is 96–100%, so the
interval is wider than that target needs.

**How bad at realistic sizes.** Finance estimates put H at about 0.6–0.8. With a monotone effect
at T = 2000, coverage is 40–88%. The EWC fixed-b variance helps a little (50–90%) but does not
fix it.

**Kill criteria.**
- KC3 is not triggered: coverage is 12%, far below 93%. **C1 survives.**
- KC4 is not triggered: path coverage at H = 0.9 with known propensity is at least 94.4%. No stop
  is needed.

**Predictions.**
- P2 is confirmed in substance. The strict knot test is mixed: the fitted knot sits at 0.50–0.56,
  between the grid points 0.5 and 0.6.
- P3 holds below 0.75 and above it. Exactly at 0.75 there is moderate under-coverage (about 7
  points with NW).
- P4 is confirmed.
