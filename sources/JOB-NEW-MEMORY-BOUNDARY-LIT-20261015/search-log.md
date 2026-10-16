# Search log

Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 · execution `claude-memboundary-20261015` · 2026-10-15
(searches run 22:10–22:40 UTC). Tools: WebSearch (standard / extended), WebFetch, local `pdftotext`.
Plan sha256 9c1d4ffa4a5ef386cfb5b854fe9bf127769f9d4522a066443c8d02452918065f (verified at start).

## Budget

| Item | Cap | Used | Note |
|---|---|---|---|
| Records screened | ≤ 60 | **79** | **Overrun by 19.** Each search returns ~5–10 records. I counted per query but crossed 60 during queries 14–15 and ran 6 more to finish item (4). Those final queries are marked ⚠. Records already in the prior job count again here. |
| Records extracted | ≤ 20 | 19 | entries in bibliography.bib with an `extracted` note |
| Full-text reads | ≤ 12 | 3 | Betken–Dehling 2107.03041v3; Zhai–Zhong–Wu 2609.30565v1; Deo–Hsieh–Hurvich 2005 |

Full-text PDFs were fetched once via WebFetch, then text was extracted locally (scratch copies,
not in the job directory). sha256: Betken–Dehling `5bc4b26c…52b50df`; ZZW `662f1397…fbfc07`;
DHH `bab85a8a…59319f7` (full hashes in boundary-theorems.md and below).

## Queries (in order) — records returned / relevant

| # | Query (abbreviated) | Mode / domains | Returned | Kept |
|---|---|---|---|---|
| 1 | HSIC kernel independence test long-range dependence size simulation fractional | extended | 8 | 2107.03041, Chwialkowski–Gretton ICML14 |
| 2 | double machine learning long memory time series coverage simulation ARFIMA | extended | 4 | none relevant |
| 3 | Hawkes Granger causality test size near criticality branching ratio | extended | 5 | 1605.06759 (prior), 2605.24595 (context) |
| 4 | long memory causal inference treatment effect … persistence | allowed: ssrn.com → **returned arXiv only** | 9 | 2609.27531 (prior), 2609.04940, 2602.00836 |
| 5 | long memory Granger causality test size robust | sciencedirect, wiley, cambridge, projecteuclid | 9 | Bauer–Maynard 2012; Shao 2011 JTSA |
| 6 | distance covariance/HSIC LRD Bernoulli/EJS/JTSA | standard | 2 new | 2112.14091 |
| 7 | "Distance correlation for long-range dependent time series" Betken | standard | 0 new | — |
| 8 | kernel two-sample MMD LRD Hermite rank non-central | standard | 3 | 1404.0551 (prior) |
| 9 | debiased/double ML "long memory"/"fractionally integrated" | standard | 5 | 2603.10999 |
| 10 | nearly unstable Hawkes inference branching ratio close to one | standard | 2 | 1310.2033, 2407.19741 |
| 11 | time series causal estimand path average vs population | extended | 4 | 1706.07840 (prior), 2509.12985 |
| 12 | Granger Hawkes power-law kernel size simulation order flow | standard | 5 | 1607.06333; PoliMi thesis (spurious cross-excitation) |
| 13 | Robinson Hidalgo 1997 time series regression LRD | standard | 1 | Robinson–Hidalgo 1997 |
| 14 | Müller HAC strongly autocorrelated; LLSW HAR | standard | 3 | LLSW 2018; LLS 2021; Casini–Perron |
| 15 | Rambachan Shephard nonparametric causal meaning | standard | 2 | 1903.01637 |
| 16 ⚠ | Deo Hsieh Hurvich long memory durations counts RV | standard | 3 | DHH 2005 (full text) |
| 17 ⚠ | Lillo Farmer 2004 order signs Hurst | standard | 4 | cond-mat/0311053 |
| 18 ⚠ | bid-ask spread long memory Hurst | standard | 3 | physics/0701017 (Mike–Farmer secondhand) |
| 19 ⚠ | Lobato Velasco 2000 trading volume | standard | 2 | Lobato–Velasco (secondhand d) |
| 20 ⚠ | Gatheral Jaisson Rosenbaum rough vs long memory | standard | 1 | 1410.3394 |
| 21 ⚠ | Tóth Palit Lillo Farmer order flow persistent | standard | 4 | 1108.1632 (number not read) |

## Fetches

| URL | Purpose | Class |
|---|---|---|
| arxiv.org/abs/2112.14091 | identity and abstract | S |
| arxiv.org/abs/2107.03041 | version history (v1 7 Jul 2021; v2 4 Jan 2024; v3 27 Jan 2026) | S |
| arxiv.org/pdf/2107.03041v3 | full text (local pdftotext) | F |
| arxiv.org/pdf/2609.30565 | full text (local pdftotext) | F |
| arxiv.org/abs/2603.10999 | abstract, venue | S |
| researchonline.lse.ac.uk/id/eprint/35727 | Robinson–Hidalgo abstract; no full text | S |
| projecteuclid.org (R–H full) | failed: empty page | — |
| archive.nyu.edu/…/26309/2/05-02.pdf | DHH full text (local pdftotext) | F |
| arxiv.org/abs/cond-mat/0311053 | Lillo–Farmer abstract | S |
| arxiv.org/abs/physics/0701017 | Gu–Chen–Zhou abstract | S |

## Identifier checks

- 2107.03041: v3 date 27 Jan 2026 confirmed on the abs page. A search summary said "last revised
  4 Jan 2024" — that is v2; the prior job's bib (v3, 27 Jan 2026) is correct. No journal ref on arXiv;
  the PDF header says "Submitted to Bernoulli".
- 2609.30565: v1 24 Sep 2026, stat.ME (PDF stamp).
- Robinson–Hidalgo: AoS 25(1):77–104, doi 10.1214/aos/1034276622 (Project Euclid per search).
  Some aggregators list pages 2054–2083, which a search summary attributed to a different Robinson
  1997 paper. I take the Euclid record.
- 2603.10999: Econometrics Journal 2026 per the arXiv page (not checked on the publisher site).
- Zhai, Kolar & Wu 2026 IEEE TIT 72(6):4246–4296: as cited in 2609.30565 only. **Not verified.**

## Not done / limits

- **SSRN not effectively searched**: the domain-restricted query returned only arXiv records.
- **Google Scholar not reachable** with the tools available.
- Econometric Theory and EJS were not hit directly. Projecteuclid appeared only via the
  Robinson–Hidalgo and Jaisson–Rosenbaum records.
- Not read in full (budget): Robinson–Hidalgo, Dehling–Rooch–Wendler (reused F from prior job),
  Ciganovic et al., Tóth et al., Plerou et al., Lobato–Velasco, Andersen et al.
