# Search log — memory-causal

Execution `claude-mkc-causal-20261015`, 2026-10-15, 21:41–21:55 UTC.
Tools: WebSearch (standard; three "extended"), WebFetch (abstract and publisher pages),
arXiv export API (identifier check), local `pdftotext` on arXiv PDFs kept in the session
scratchpad (outside the repository). Window emphasised: 2025-01-01 to 2026-10-15.

## Budget

| Item | Limit | Used |
|---|---|---|
| Records screened | ≤ 70 | **74 — over the limit by about 4** (recount at write-up; see note below) |
| Records extracted (assumption-level notes in field-map) | ≤ 25 | 24 |
| Full-text reads (F) | ≤ 12 | 9 new (+2 reused F from JOB-NEW-CAUSAL-RISK-EXCEED, not re-read) |

## Reuse step (read before searching)

- research_assistant `docs/topics/conformal_causal_fin_ts/S3_causal_time_series.md` and
  `S8_point_processes.md`: read in full.
- `LIT_REVIEW.md`, `GAP_ANALYSIS.md`, `DEVILS_ADVOCATE.md`, `VERIFICATION.md`, `S2`, `S5`: grepped
  for long memory / long-range / fractional / rough / Hurst / near-critical / branching /
  power-law / non-mixing / ARFIMA. Hits: GAP_ANALYSIS T6 and R8, LIT_REVIEW D2 (FDM covers long
  memory for conformal), S2 rows 10–12. Confirms the reuse corpus treats long memory least.
- `programme/catalog.yaml`: tracks B1, B2-HTE, B2-COMP, B2-QUEUE, F-CAUSAL, F2, J, K, E, F-KERNEL.
- `runs/JOB-A-HAWKES-CHECK/handoff.md`; `runs/JOB-NEW-CAUSAL-RISK-EXCEED/handoff.md`, `baseline.md`.
- The `.bib` files of the reuse corpus were not re-parsed; entries reused here were re-verified
  through the arXiv API instead.

## Queries (in order) and yield

| # | Query (abbreviated) | Mode | New relevant records |
|---|---|---|---|
| 1 | double machine learning long memory time series causal non-mixing 2025 | std | 3 (2603.10999, 2411.10009, 2311.06012) — none on long memory |
| 2 | spurious Granger causality long memory fractional | std | 5 spurious-long-memory test papers (Haldrup–Kruse; Leccadito et al.; Kuswanto; Frederiksen–Nielsen; Qu R pkg) |
| 3 | causal discovery Hawkes THP SHTP 2025 2026 | std | 3 (2105.10884, 2306.14114, 2608.21815); "SHTP" not found |
| 4 | nearly unstable Hawkes branching ratio near one inference | std | 4 (1310.2033, 2407.19741, 1507.02887, 2605.24595) |
| 5 | "long-range dependence" causal effect time series 2025/2026 | ext | 4 (Hao PLOS One 2026; 2609.27531; 2504.11076; 2609.24153) |
| 6 | Granger causality fractionally integrated VAR realized volatility | std | 4 (Chen 2006; Bauer–Maynard; Hassler et al.; FCVAR) |
| 7 | causal effect intervention Hawkes market impact order flow 2025 | std | 4 (2607.03239, 2601.23172, 2511.02518, Bacry–Muzy 1301.1135) |
| 8 | transfer entropy long memory fractional Gaussian noise | std | 3 (0910.4514, VARFI-TE Entropy 2021, fractional TE) |
| 9 | Eichler Dahlhaus Dueck Hawkes Granger | std | 1 (1605.06759) |
| 10 | local projection persistent data Montiel Olea Plagborg-Møller | std | 1 (2007.13888) |
| 11 | β-mixing nonlinear Hawkes power-law polynomial rate 2026 | std | 2 (Cheysson talk abstract; 1904.06051) — no polynomial nonlinear result found |
| 12 | rough volatility causal inference fBm confounding | std | 4 (1905.04852, 2203.13820, 2312.01426, 2210.01214) — none causal |
| 13 | DML functional/physical dependence measure cross-fitting | std | 3 (2601.10899, 2602.11333, 2511.10995) — none use Wu's FDM |
| 14 | interventional Hawkes do-intervention counterfactual 2025 | std | 2 (2111.07603, 2507.10809) |
| 15 | Granger Hawkes HF data near-critical spurious | std | 4 (1308.6756, 1607.06333, 1709.08516, Wheatley et al. 2019) |
| 16 | synthetic control / DiD long-memory errors | std | 2 (Tsay difference-based; Ferman–Pinto) — none on fractional errors |
| 17 | causal discovery long memory PCMCI | std | 2 (1702.07007, 2003.03685) |
| 18 | Christgau Hansen local independence single path 2025 | std | 2 (2110.12709, 2512.20280) |
| 19 | negative control proximal time series 2025 | std | 4 (2109.07030, 2608.01352, 2402.00335, 2108.13935) |
| 20 | self-normalized inference long memory treatment effect | std | 3 (Shao review; McElroy–Politis; Bai–Taqqu–Zhang) |
| 21 | arXiv 2026 Hawkes Granger test single trajectory | ext | 4 (2508.11727, 2308.06106, 2309.02027, ICLR 2026 paper — id not captured) |
| 22 | causal market impact randomized metaorders 2025 | std | 3 (2502.16246, 2509.05065, 2501.17096) |
| 23 | Hao long-range dependency endogenous treatment PLOS One | std | 1 (confirmation) |
| 24 | Baillie Diebold Kapetanios Kim robust inference | std | 1 (2203.04080) |
| 25 | local projections long memory fractional | std | 1 (Baillie–Kapetanios–Papailias ER 2017) |
| 26 | order flow long memory LMF metaorder splitting | std | 2 (2301.13505, 2502.17906) |
| 27 | estimation near-critical Hawkes MLE | std | 0 new |
| 28 | Tsay Chung spurious regression fractional | std | 1 (J. Econom. 96(1) 2000) |
| 29 | Hawkes potential outcomes treatment effect single path 2025/2026 | std | 3 (2604.12124, 2604.12977, 2301.03246) |
| 30 | arXiv 2025/2026 long-memory Granger new test | ext | 0 new 2025–26 (Bauer–Maynard re-found) |
| 31 | Hawkes debiased / orthogonal / DML 2025/2026 | std | 0 |
| 32 | Sato Kanazawa trader ID TSE | std | 1 (2308.01112) |
| 33 | Wehrli Sornette endo-exo | std | 1 (QF 2019, plus SFI 20-39 follow-up) |

Records screened: the yield column sums to 78; removing 4 repeats (Hao, Wheatley et al.,
Bauer–Maynard, 2111.07603) gives **74 distinct records screened from searches — about 4 over the
≤ 70 limit.** I did not track the running count during searching; the overrun was found at
write-up. The reused-corpus anchors (1712.09988, 1706.07840, 1903.01637, 2603.20394, 2203.13559,
2007.07448, 2311.11730, 2003.04314, 1905.11053, 2512.08828, 1504.03100, 1609.05177, 1302.1405)
were only re-verified by identifier, not re-screened, and are not counted. Some yield
entries were judged irrelevant on screening (2311.06012, 2608.21815, 2511.02518, 2605.24595,
2609.24153, 2504.11076, 2511.10995, 2602.11333, 2512.20280, 2608.01352, 2402.00335, 2308.06106,
2309.02027, 2509.05065, 2502.17906, 2210.01214, 2312.01426).

## Fetches and full-text reads

| Record | Access | What was read |
|---|---|---|
| arXiv:2007.13888 (MOPM) | F | §1 intro, §2 Ass. 1, §4.2 Ass. 2–3 |
| arXiv:2411.10009v2 (BW) | F | §3 Ass. 3 + Thm 2 remark; Procedure 2, Ass. 6 and dependence examples (i)–(iv) |
| arXiv:2609.27531v1 (Yu et al.) | F | abstract; Ass. 1–3 headers; mixing statements |
| arXiv:2601.10899v3 (Balkus et al.) | F | Ass. 1–3, Lemma 3, §3.3 Cor. 4 scope (m-dependence) |
| arXiv:2007.07448 (WKS) | F | §2 model, §4 Ass. 1–4 |
| arXiv:1605.06759 (EDD) | F | Ass. 2.1, Def. 2.3, 3.1, Thm 3.5 statement, Thm 4.1 conditions |
| arXiv:2601.23172v2 (Muhle-Karbe et al.) | F | abstract, §1, §2.1–2.3 |
| arXiv:2607.03239v1 (Leclère et al.) | F | §1.1, §2.3–2.4, Ass. A–B |
| arXiv:2301.03246v1 (Jiang–Chen–Ding) | F | §2 Ass. 1–5 |
| doi:10.1371/journal.pone.0347847 (Hao) | S | publisher page, first 100k of 172k characters; supplements not read |
| arXiv:2501.17096, 2604.12124, 2604.12977, 2508.11727, 2507.10809 | S | arXiv abstract pages |

PMC mirror of Hao blocked by a browser check; publisher page used instead.

## Identifier verification

arXiv export API returned title, authors, first-version date, journal-ref and DOI for 41 ids
(all matched the cited titles and authors). One mistaken id (1507.01729, unrelated) was queried
and discarded. Non-arXiv records are taken from publisher or RePEc pages seen in search results
and labelled A. Not verified: the VARFI transfer-entropy paper's id and author list; the ICLR
2026 Hawkes paper's title/authors; the Filimonov–Sornette journal venue.

## Contrary evidence kept

- Wheatley–Wehrli–Sornette; Filimonov–Sornette (against near-criticality).
- Cont–Das (roughness may be an artefact).
- Balkus et al. (bespoke folds unnecessary — against the buffered-fold line).
- Cordi–Challet–Muni Toke (Hawkes directional causality is weak).
- Bauer–Maynard (linear Granger testing under long memory already handled — kills G5 as a gap).

## Limits

- No Google Scholar or SSRN direct search; arXiv full-text search not used.
- Unsuccessful searches are logged, not interpreted as absence.
- No author-page or citation chasing beyond the queries above.
