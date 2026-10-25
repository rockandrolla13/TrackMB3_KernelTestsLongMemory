# Framing drafts for Major concern 2 (decision rule)

Drafted 2026-10-25. Nothing in `paper.tex` was edited. These files are drop-in replacements.

## Files and where they go

| File | Replaces | Location in paper.tex |
|---|---|---|
| `abstract.tex` | whole abstract | lines 24–26 |
| `intro.tex` | `\section{Introduction}` plus paragraphs 1 and 2 | lines 30–33 |
| `guidance.tex` | new Box 1 | after line 43 (end of Section 1) |
| `conclusion.tex` | Conclusion body only | lines 363–364. Keep the data and ethics paragraphs. |

Add this line to the preamble: `\newcommand{\todo}[1]{\textbf{[TODO: #1]}}`.
All four files compile together against the paper's preamble (scratch test, one pass, no errors).

## What changed and why

1. **The abstract leads with the decision rule.** It names the rule, gives its evidence (3.4–7.0% null rates at n=2000 and a power shortfall of at most 4.9 points), and states two limits: no long-memory theory, and 9.7% under slow drift. Four contributions became two findings plus a short stock sentence. Without the PROP-A sentence it has 173 words, leaving 20 for that sentence under the 200 cap.
2. **The intro now opens with the question and the answer.** The old paragraph 2 (test definitions) moves to paragraph 1, which ends with the three-step rule and defines "circular shift" at first use. The contributions paragraph comes second and has three items, not four. The "different nulls, not a ranking" caveat stays.
3. **Box 1 is new.** It has five numbered items: adjust, test with shifts, expect drift risk, what the rule does not promise, and split MMD. Each claim carries a number from the paper and a section reference.
4. **The conclusion leads with the recommendation.** Its order is rule, then three limits, then the two findings behind the rule, then split MMD and the stock arm.

## Limits stated honestly (where each appears)

- No long-memory theory for shifts: abstract, Box item 4, conclusion.
- 9.7% under slow drift (6.2% at n=500): abstract, Box item 3, conclusion.
- Adjustment raises size under memory (HSIC sub 17.5→23.8, shift 5.3→6.4 at H=0.8): Box item 4, conclusion.
- Shared fixed phase defeats shifts (28–34%). This is why adjustment comes first: Box item 1.
- MMD power loss at a fixed δ (half-mean SD 0.045→0.25): Box item 5, conclusion.
- Stock fractions are descriptive: abstract, intro, conclusion.

## How the 3.4%–7.0% range was derived (please verify)

It is the min and max of the HSIC shift column at n=2000 across Table 3 (rules), Table 4 (EXP2, rows with independent phases) and Table 6 (EXP4, A0/A1/A2 rows). It excludes SEAS-V-slow (9.7, 10.3) and the shared-fixed-phase diagnostic (27.9–34.3). The minimum 3.4 is SEAS-V A1 at H=0.8. The maximum 7.0 is NL-TAR at H=0.8. The paper tables report no n=500 shift sizes, so the range is stated for n=2000 only. These are uncorrected descriptive comparisons (review minor point 6). Box 1 says so in a comment, and the main text could say it too.

## Placeholders

| Tag | Files | Waits on |
|---|---|---|
| `PROP-A` | abstract, intro, conclusion | theory agent: random-phase seasonal-mean proposition (one sentence each, ≤20 words in the abstract) |
| `PROP-B` | intro, guidance, conclusion | theory agent: MMD half-mean spread proposition |
| `STOCK-D2` | intro, guidance, conclusion | stock-diurnal agent: stock rates after the smooth diurnal factor |
| `MMD-SCALED` | guidance, conclusion | mmd-scale-n agent: power at a shift scaled to null-SD units |
| `N4000` | guidance, conclusion | mmd-scale-n agent: shift null rates at n=4000, sums 1.4 and 1.6 |
| `BW result` | guidance, conclusion | bandwidth agent: does shift HSIC keep size and power at 0.5× and 2× the median? |

If BW or N4000 shows that shifts fail, items 2 and 4 of the box need to change. A narrower rule would also need a new abstract sentence.

## New citations (placeholder keys, not in refs.bib)

- `engle1998acd`: Engle, R. F. and Russell, J. R. (1998). Autoregressive conditional duration: a new model for irregularly spaced transaction data. *Econometrica* 66(5), 1127–1162.
- `wu2012intraday`: Wu, Z. (2012). On the intraday periodicity duration adjustment of high-frequency data. *Journal of Empirical Finance* 19(2), 282–291.

Both come from the peer review's missing-citation list. Check the bibliographic details against the publisher before adding them.

## Knock-on edits for the editor (not made)

- The old intro paragraph 3 onward repeats material. Paragraph 5 (9.2/7.35/5.45%) and the last two sentences of paragraph 6 ("Two checks bound this finding...") overlap with the new contributions paragraph and Box 1. Consider trimming them.
- The Discussion paragraph at line 352 restates the shift and adjustment remedies. It could point to Box 1 instead.
- The keyword list could add "circular shift".
