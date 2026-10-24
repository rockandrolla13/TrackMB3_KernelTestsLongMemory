# Handoff — JOB-OPPORTUNITIES-MEMORY-20261023

Discovery worker `opp-memory-discovery-20261023`, 2026-10-23. Read-only opportunity assessment.
It contains no decision, no experiment and no novelty claim. A separate reviewer should check it next.

## Plain-English status

**NEW-MEMORY-BOUNDARY — narrowed.** The literature does not support one universal 3/4 line.
The line holds for "square-type" statistics, which have Hermite rank 2.
Linear-type statistics break as soon as any long memory appears (H > 1/2).
For two series the line is a sum, H_X + H_Y = 3/2. Hawkes tests break with the branching ratio, not with H.
This job found no paper that maps all three method classes on one memory grid. That is a search result, not proof that no such paper exists.
New this time: econometrics already has "fragility" results. Müller (2007) shows that any long-run variance estimator that is consistent can be thrown off by tiny low-frequency contamination.
Müller and Watson (2018) show the other side: coverage can be kept across many persistence types, for linear quantities.
These results bound what any "works for every H" claim can promise.

**G1 (HSIC null limit) — open but narrowed.** HSIC is a kernel measure of dependence between two series.
The closest theorem (Betken–Dehling, distance covariance) says in its own words that it does not cover the "mixed" case. In that case one series has strong memory, the other weak, and their sum is still below the line.
It also covers the high-memory side only for Gaussian data. This job found no newer version or other paper closing either case.
This is the best-supported gap. An older general theorem may still imply it, so a specialist must check before anyone calls it new.

**G2 (H-free kernel calibration) — open.** "H-free" means the test needs no estimate of the memory level.
Earlier reviews found only scalar H-free methods and kernel methods for weak dependence.
This job adds a route none of them cited: low-frequency inference (Müller–Watson). It averages the data over a few slow cosine waves and keeps coverage across many memory types. So far it works only for linear quantities.
Whether it works for kernel features is unknown. It may fail above the line, because there the averages may not be normally distributed (my inference, not a source).

**G3 (asymmetry and rank-2 protection) — answered at model level.** "Rank 2" is the property that would delay the break to 3/4.
The generic mechanism was already closed (Bai–Taqqu 2019).
New here: in the standard stochastic-volatility model of returns, with exponential volatility, every power of returns has rank 1 (Kulik–Soulier 2012, §5.1). That holds even with no asymmetry.
So "rank-2 protection on financial data" fails in that model, and the break for |r| or r² is expected at H = 1/2.
Leverage, the link between return shocks and volatility shocks, changes the ranks again. Real-data rank stays an empirical question.

## Shortlist (candidates.json, validated)

| ID | Status | Assessment |
|---|---|---|
| MB-G1-HSIC-MIXED-NULL | ranked 1, narrowed, moderate | Strongest supported gap |
| MB-G2-LOWFREQ-KERNEL | ranked 2, open, low; next_priority 1 | Best next investigation (one derivation decides it) |
| MB-UNIFORM-SIZE-LIMITS | ranked 3, narrowed, low | Framing; two sources abstract-only |
| MB-G3-LMSV-RANK1 | rejected, answered | Correction to G3 premise |
| MB-CROSS-CLASS-MAP | provisional, uncertain | Waits for a citation sweep |
| MB-HAWKES-NO-H-LINE | rejected, narrowed | Reframe Hawkes by branching ratio |

## Closest primary sources

- Betken & Dehling, arXiv:2107.03041v3, Remark 2.1(i) (full text).
- Müller & Watson, Econometrica 2018 / NBER w23186, Table 2 discussion (full text).
- Müller, J. Econometrics 2007, Theorem 1 (full text).
- Kulik & Soulier, Adv. Appl. Probab. 2012, arXiv:1109.5298, §5 (full text, first 100k characters).

## Contrary results

- Shao 2009: an H-robust linear independence test already exists, so a broad "no robust independence test" claim is wrong.
- Ghoshal–Shao 2026: kernel split, project and self-normalise exists under strong mixing.
- Müller 2007 is explicitly *not* an impossibility of valid inference.
- Bai–Taqqu 2019 already closed G3's generic mechanism.

## Uncertainties

- Pötscher 2002 and Preinerstorfer–Pötscher 2016 were read at abstract level. Their model classes may exclude fractional d.
- "Hermite ranks and U-statistics" (HAL hal-01197638) was blocked. It may pre-empt the kernel-rank step of G1.
- No forward-citation index, Google Scholar or SSRN was used.
- Kulik–Soulier treats infinite-variance LMSV. The finite-variance rank statement is my derivation.

## Suggested challenge searches for the reviewer

1. Forward citations of Betken–Dehling 2107.03041 and Zhai–Cheng–Wu 2609.33868 in a working index.
2. "Arcones 1994" and "Major 2019 vector non-central" applied to products H_r(ξ)H_s(η). Does that settle G1's mixed case?
3. Müller–Watson forward citations with "nonlinear", "kernel" or "distribution" (G2).
4. Pötscher 2002 §§2–3 full text: is fractional d ∈ (0, 1/2) inside the model class?
5. Finite-variance LMSV with leverage: Surgailis–Viano 2002, Deo–Hurvich, for |r| Hermite rank (G3).

## Output receipts

SHA-256 values are recorded in the `finished` event of `events.jsonl`.
