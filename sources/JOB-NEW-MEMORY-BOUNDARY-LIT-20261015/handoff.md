# Handoff — memory boundary

Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 · execution `claude-memboundary-20261015` · 2026-10-15.
Bounded scan. No novelty claim. Nothing was run except an exact formula check.

## Has anyone already mapped this?

No, as far as this scan reached. No paper maps test size or interval coverage across memory
levels for kernel tests, DML intervals or Hawkes Granger tests. The two closest papers test
only two or three memory levels each, and never near the line. So the stopping rule did not fire.

SSRN was not really searched: the tool returned arXiv results instead. Google Scholar was not
reachable. Treat this as a real gap in the search.

## Is the H = 3/4 line real for these methods?

Partly. It is real, but it belongs to some statistics and not others.

- **It holds for "square" statistics.** Both papers prove it for statistics built from products or
  squares: distance covariance and sample covariance. A formula check here confirms the rates.
- **For two series it is a sum rule.** The real line is H_X + H_Y = 3/2. Order signs (0.7) paired
  with volatility (0.9) are outside it, even though one series is below 3/4.
- **For "linear" statistics it is the wrong line.** This covers the MMD two-sample test and the
  DML population average effect with effect heterogeneity. They break as soon as there is any long
  memory (H > 1/2). So the C1 idea needs correcting. Its boundary is d = 0, not d = 1/4.
- **Extending to HSIC looks direct** for the usual kernels. That is my inference; nobody has
  stated it.
- **Hawkes Granger tests follow a different mechanism.** They should degrade with the branching
  ratio, not with an H line.
- **Contrary evidence.** The one proven test (Betken–Dehling) already rejects 10–17% of the time at
  a 5% level at H = 0.7 with n ≤ 1000. Trouble starts below the line in finite samples.

## Where financial series sit

Mostly above the line.
- Durations, trade counts, volume and spreads are reported at H ≈ 0.75–0.96.
- Realised volatility: d ≈ 0.3–0.45, unless volatility is rough rather than long memory, which is
  disputed.
- Order signs straddle the line: about 0.7 in London, 0.65–0.9 in Paris, and about 0.75 in
  metaorder theory.
- Estimates move by 0.1–0.3 with estimator settings.

## The proposed simulation (not run)

- Simulate memory levels H from 0.5 to 0.9, plus unequal pairs, a rough-volatility arm and
  near-critical Hawkes.
- Measure the size of HSIC, distance-covariance, MMD and Pearson tests. Use three calibrations.
- Measure AIPW coverage for the population effect vs the path-conditional effect, with three
  effect shapes.
- Measure the size of a Hawkes likelihood-ratio Granger test.
- Settings: 2,000 replications, n up to 4,000–8,000, capped at 48 CPU-hours per family.
- Six kill criteria are fixed in advance. Examples: no size distortion at all; a break at the
  wrong place; C1 coverage staying nominal.

## Decisions for the author

1. Approve, change or reject the simulation design (simulation-design.md). It is a separate job.
2. Accept the reframing: the headline becomes "a Hermite-rank map of where inference breaks", not
   "everything breaks at 3/4"?
3. Who owns each arm: one shared job here, or split to Tracks E, K/F-CAUSAL and B2-COMP?
4. Should C1 continue, given it is likely seen as a corollary of classical results? Its value would
   be as a stated warning, not a theorem.
5. Fund a follow-up search of SSRN and Google Scholar, which this job could not reach?
