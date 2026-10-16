# Handoff — literature review for kernel tests under long memory

Execution `claude-mb3-lit-20261016` · Agent 2 · 2026-10-16. Bounded review. No novelty claim.
Details: `review.md` (answers to all 18 questions), `systematic-review.md`, `gap-analysis.md`,
`evidence-matrix.md`, `search-log.md`, `bibliography.bib`.

## What is already known

- **The linear sum rule is old.** For Pearson-type and regression tests, the line d_X + d_Y = 1/2
  (H_X + H_Y = 3/2) was shown by Tsay and Chung in 2000. So the cross rule is only new for kernel tests.
- **H3 is already proved for one case.** Betken and Dehling show that, for Gaussian series above the
  line, the distance-covariance null is a squared cross-correlation variable times a constant. The open
  part is non-linear series (like order signs) and HSIC.
- **Self-normalisation depends on H.** This is confirmed. The dependence is not mild: one published
  5% critical value rises 42% from H = 0.6 to 0.9. A known fix avoids estimating H: self-normalise,
  then block-sample (Bai, Taqqu and Zhang 2016; Betken and Wendler 2018).
- **Rate-corrected subsampling is known for simple statistics.** The Wu group used it in September 2026
  to cross the 3/4 line for conformal coverage.
- **H8 is largely pre-empted.** The same group proved the 3/4 line for coverage with |residual| scores.
- **Prewhitening works for linear independence tests under long memory** (Shao 2009).
- **The wild-bootstrap oversize (H5) is a known effect of centred multipliers.** The method's own authors
  report it (Chwialkowski, Sejdinovic and Gretton 2014).

## What is still open (not found)

- An HSIC or MMD null law under long memory, for any kernel.
- **The mixed case (H2b):** one series above 3/4, sum below 3/2. Betken and Dehling say their proof does
  not cover it. Nobody citing them has closed it. This is the clearest theory gap.
- A self-normalised kernel test under long memory (H4b). Both halves exist in print; the join does not.
- Any kernel conditional-independence test under long memory (H7).
- A size map of kernel tests across H near the lines.

## Main competitor risk

The Zhai / Wu group. In 2026 they posted four papers on long-memory inference. They name the 3/4
regime as their next limit, and in one paper they already cross it. They show no kernel work yet.
Second risk: Shao's group. They hold the kernel self-normalisation method and the long-memory
self-normalisation results.

## What the author must decide

1. Should the headline drop the linear sum rule and claim only the kernel version and the mixed case?
2. Should H4b move to the H-free route (self-normalise, then block-sample) instead of tabulating in H?
3. Should H8 go to F-KERNEL now, with the two Zhai–Cheng–Wu papers attached?
4. Should H5 be fixed by following the published two-sample bootstrap before reading more wild numbers?
5. Should someone with library access search SSRN and read the unread classics? SSRN blocked us
   (HTTP 403). The unread classics are Hsing–Wu 2004, Dehling–Taqqu 1991, Arcones 1994 and Shao 2011.
   A U-statistic result covering HSIC could sit there.

## Corrections to the hypotheses file

- "Hall, Horowitz & Jing (1998)" should be Hall, Jing & Lahiri (1998), Statistica Sinica.
- Ghoshal–Shao assume strong (α-) mixing, not β-mixing, in the published version. It is the same work
  as the EJS 2026 paper.
- Zhai–Kolar–Wu (IEEE TIT 2026) is verified.
- One caution for H1c: Bai and Taqqu argue that Hermite ranks above 1 are unstable in practice.
  A real |r_t| may behave like rank 1 at large n.
