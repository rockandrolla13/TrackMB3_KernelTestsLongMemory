# Handoff — kernel-linear arm (K-ind, L-ind)

Worker `claude-mb-kernel-20261015` · exploratory simulation, not evidence for any claim.
**n = 4000 was dropped for the kernel tests by author decision D-a28e4cb5fa.** So the kernel verdicts are at n = 2000.

**Does the H = 3/4 (cross-rule) line show up?** Yes, in where the line sits. Size gets sharply worse once
H_X + H_Y passes about 1.5, which is H ≈ 0.75 per series. All seven tests put the kink between about 1.4 and 1.6, and
1.5 is inside every interval. But the full prediction P1 fails for every kernel test. P1 also says size keeps growing
with n beyond the line. Up to n = 2000 that only shows clearly at H = 0.9. Pearson with the EWC correction follows the
rule exactly: nominal below the line, growing with sample size above it.

**How bad is it at realistic n (2,000 points)?** Nominal is 5%.
- **Mild memory** (H ≤ 0.7): most tests reject 5–9% of the time. Subsampling reaches 11–14% just below the line (G4, sum
  1.45). Newey–West reaches 8–9% there.
- **At or just past the line** (H 0.75–0.8): wild bootstrap 6–10%, subsampling 10–25%, Pearson-NW 10–15%.
- **H = 0.9**: wild bootstrap 13–18%, subsampling 34–52%, Pearson-NW 21–34% (T = 2000–8000).

**Which calibration holds best.** For kernel tests, the shift test. It stays within 3.9–7.6% everywhere, even at
H = 0.9. For the linear test, EWC fixed-b, but only below the line. Subsampling is the worst, even though it is the only
calibration with a proof under long memory. Wild bootstrap is up to 1.6 pp oversized even with no memory at n = 2000.

**Rough-volatility arm.** A strongly persistent short-memory series (λ = 0.01) fakes long memory. It drives subsampling
to 43–46% and NW to 18%. Its distortion shrinks as n grows, which separates it from true long memory.

**Kill criteria.** KC1 (no bite): **not triggered**. Most wild-bootstrap cells at n = 2000 are outside [3.5%, 6.5%].
KC2 (wrong line): **not triggered** for any kernel calibration.

Power check (r = 0.25, H = 0.7): every test rejects 99–100% at every n.

Run facts:
- 2,000 reps per null cell; 0 errors.
- 29.6 of 48 CPU-hours used.
- My Newey–West uses Andrews' weights 1 − j/S and no clip on ρ̂. The mmd-dml arm uses 1 − j/(S+1) with ρ̂ clipped
  to ±0.97.

Details and every cell are in `results.md`.
