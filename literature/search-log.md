# Search log — TrackMB3 deep literature review

Execution `claude-mb3-lit-20261016` · Agent 2 · 2026-10-16, searches run about 13:22–13:45 UTC.
Hand-off verified first: `hypotheses/hypotheses-and-claims.md` SHA-256
`b937424a3c1b1595c06bbb529df603eb036019590eb2ce89e40733b4c93afb04` (matches).

Tools: arXiv export API (`export.arxiv.org/api/query`), Crossref REST API, OpenAlex API (abstracts),
Semantic Scholar Graph API (forward citations; often rate-limited), WebSearch (standard and extended,
some with `allowed_domains`), WebFetch, Google Scholar result pages fetched with curl and parsed, local
`curl` + `pdftotext` for full text. PDFs were kept in the session scratchpad, not in the track. Their
SHA-256 hashes are in `events.jsonl` (`primary_sources_checked`).

Evidence labels: **F** full-text section read by me; **S** summarising fetch or a verbatim quote read
inside another full text; **A** abstract or snippet; **I** identity verified only (title, venue, DOI),
content not read.

## Budget (running count, final)

| Item | Cap | Used | Note |
|---|---|---|---|
| Records screened | ≤ 120 | **103** | Unique records whose title and abstract (or snippet) I assessed. Off-field noise (hyperspectral imaging, LLM papers, etc. returned by keyword engines) was not counted. |
| Records extracted | ≤ 40 | **38** | Rows in `evidence-matrix.md` with fields filled. |
| Full-text reads | ≤ 20 | **15** | 14 useful + 1 failed numeric extraction (Tóth et al. 2015). List below. |

Full-text reads (sections, not proofs): Betken–Dehling 2107.03041v3; Zhai–Zhong–Wu 2609.30565v1;
Zhai–Kolar–Wu 2604.16219v1; Ghoshal–Shao 2601.13468v2; Dehling–Rooch–Wendler 1404.0551v1;
Diz-Castro et al. 2604.28104v2 (assumption scan); Betken 1403.0265v2; Betken–Wendler 1509.05720v3;
Bai–Taqqu–Zhang 1512.00820v2; Betken–Kulik 1706.06351v1 (Sec 1 and 3); Shao 0810.2276v1;
Zhai–Cheng–Wu 2609.33868v1; Ventosa-Santaulària et al. (AAU accepted manuscript); Chwialkowski–
Sejdinovic–Gretton 1408.5404; Tóth et al. 1108.1632 (extraction failed for the exponent).

## Queries in order

`raw` = records returned with plausibly relevant titles; `new` = records not already seen in this job.

### arXiv export API (field-restricted, AND of terms)

| # | Query | raw | new kept |
|---|---|---|---|
| X1 | id_list 2609.30565, 2601.13468, 2604.28104, 1404.0551, 2112.14091 (+2107.03041) | 6 | 6 (identity + abstracts) |
| X2 | ti:"long-range dependent" AND ti:precision | 1 | 2604.16219 (ZKW preprint) |
| X3 | au:Zhai AND au:Wu_W AND abs:memory | 0 | — |
| X4 | abs:"long-range dependence" AND abs:HSIC | 0 relevant (2 off-field) | — |
| X5 | abs:"long-range dependence" AND abs:"maximum mean discrepancy" | 0 | — |
| X6 | abs:"long memory" AND abs:kernel AND abs:"independence test" | 0 | — |
| X7 | abs:"long-range dependent" AND abs:"distance covariance" | 1 | (2107.03041 again) |
| X8 | abs:"long-range dependence" AND abs:"self-normalized" | 6 | 1509.05720, 1403.0265, 1706.06351, 1512.00820, 2004.06574 (2609.00717 excluded) |
| X9 | abs:"long memory" AND abs:"self-normalized" | 3 | 1006.1572 (1312.3120 excluded) |
| X10 | abs:"long-range dependent" AND abs:"self-normalization" | 6 | 0 new |
| X11 | abs:"long-range dependence" AND abs:"two-sample" AND abs:"U-statistic" | 0 | — |
| X12 | abs:"long-range dependent" AND abs:"energy distance" | 0 | — |
| X13 | abs:"long-range dependence" AND abs:"change-point" AND abs:kernel | 1 | (1404.0551 again) |
| X14 | au:Betken (all) | 28 | 2211.01108, 2601.14223, 2408.10570, 2511.04870, 1905.11033, 1612.08261 screened; rest off-topic (Carina Betken) not counted |
| X15 | abs:"long-range dependence" AND abs:"conditional independence" | 0 relevant | — |
| X16 | au:Zhai_P AND abs:conformal | 0 (author string mismatch) | — |
| X17 | ti:"instability of ranks" | 1 | 1610.00690 |
| X18 | au:Poskitt AND abs:fractionally | 6 | none kept (bias-correction papers; off-question) |
| X19 | abs:"Hermite rank" AND abs:kernel AND abs:test | 0 | — |
| X20 | abs:"long memory" AND abs:Granger AND abs:nonparametric | 0 | — |
| X21 | abs:"long-range dependence" AND abs:bootstrap AND abs:"U-statistic" | 0 | — |
| X22 | abs:"long memory" AND abs:"independence test" | 0 relevant | — |
| X23 | abs:"long-range dependent" AND abs:independence AND abs:test | 25 | 1611.06828 screened (excluded: α-dependent, not LRD); rest off-topic |
| X24 | abs:conformal AND abs:"long memory" | 5 | 2609.33866, 2609.33868 |
| X25 | abs:conformal AND abs:"functional dependence measure" | 2 | 0 new |
| X26 | id_list abstracts: 0810.2276, 2510.23563, 2605.22025, 1908.06486, 2208.01713, 1512.00819, 1601.01122, 1706.07237, 2401.11272, 2301.13505, 1108.1632, 2609.02525, 2609.20772, 2510.02471 | 14 | all screened |

### WebSearch (general and domain-restricted)

| # | Query (abbreviated) | Mode / domains | Result |
|---|---|---|---|
| W1 | Ghoshal Shao Electronic Journal of Statistics 2026 SN RKHS | standard | arXiv only; EJS not confirmed by search → confirmed via Crossref (C3) |
| W2 | Shao 2011 "A simple test of changes in mean…" critical values Hurst | standard | Wiley/RePEc/WUSTL records; no full text |
| W3 | HSIC independence test "long-range dependence" OR "long memory" | extended | 2107.03041; 1402.4501; 2605.22025; 1804.09866; 2212.09108 (JMLR) |
| W4 | kernel change-point MMD "long-range dependent" | extended | 1210.4903, 2206.01853, Li–Xie (NeurIPS 2015), 1903.01661, 1202.3878; DRT 2013 |
| W5 | independence test long memory kernel distance covariance | **allowed_domains ssrn.com** | **Returned arXiv records only** (2107.03041, 2310.08575, 1908.06486, 2604.28104, 2112.14091, 0810.2276, 2212.09108, 2510.23563). SSRN not effectively searched by this route (same failure as the 2026-10-15 job). |
| W6 | "long-range dependence" kernel test … Hilbert space | projecteuclid.org, bernoullisociety.org, imstat.org | AoS local permutation CI (2022), EJS 2021 kernel metrics; no LRD kernel test |
| W7 | long memory independence test nonparametric dCov / HSIC / ECF | sciencedirect, wiley, tandfonline, cambridge, OUP | Beering (JNPS 2025) title; JMVA 2025 projection HSIC; WIREs review; 1903.00708 |
| W8 | kernel two-sample / independence non-mixing or LRD | neurips, PMLR, JMLR, openreview | 1402.4501, 2010.00271, 1408.5404, 2604.28104, 1908.06486 |
| W9 | Granger causality fractionally integrated size nonparametric | standard | Chen (2006) spectrum Granger; Bauer–Maynard; Andersson–Gredenhoff; Lach |
| W10 | Zhai Cheng Wu conformal FDM long memory | standard | not found by WebSearch (found via X24) |
| W11 | "conditional independence" "long memory" residuals | standard | 2609.20772 (GTCM), 2504.21647 (dGCM), Beran (1992) GOF, 2110.08089 |
| W12 | change-point LRD Hermite rank absolute value / even function | standard | 1506.07296 (Tewes), DRT 2013, Rooch thesis |
| W13 | order sign Hurst crypto / futures 2020–2026 | standard | only pre-2020 equity/FX and one 2020 estimator paper; no crypto/futures numbers |
| W14 | prewhitening long memory cross-correlation independence | standard | Hong 1996; Qu (2011) spurious-LM test; no LRD kernel prewhitening |
| W15 | Tsay Chung 2000 "d1 + d2" | standard | Sinica repository; AAU accepted manuscript of Ventosa-Santaulària et al. (read F) |
| W16 | Sato Kanazawa LMF TSE exponent | standard | PRL 131 197401; PRR; JSP — summary says α ≈ 3/2 |
| W17 | wild bootstrap V-statistic MMD oversized centering | standard | 1408.5404 only → read F |
| W18 | LRD two-sample / GOF ECF / energy distance non-central | wiley, sciencedirect, projecteuclid, cambridge, springer | 0907.1787; CSDA 2008 ECF test; ECF fast test 2016; DRT 2013 |
| W19 | kernel CI / GCM long-range dependence fGn | standard | GKCM (PMLR 323, 2026); GCM (AoS 2020); nothing under LRD |

### Google Scholar (curl, result pages parsed)

| # | Query | Parsed hits | Note |
|---|---|---|---|
| G1 | "long-range dependence" HSIC independence test | 3 (all off-field) | results look recency-sorted and noisy |
| G2 | "long memory" "kernel" "two-sample" test "maximum mean discrepancy" | 7 | 2609.22545 (online kernel CPD) screened; rest off-topic |
| G3 | "long-range dependent" "distance covariance" | **0 parsed** | page returned (149 KB) but no result blocks: treat as a failed query, not as absence |
| G4 | "self-normalization" "long-range dependence" kernel | 2 | off-question |
| G5 | "conditional independence" test "long memory" time series | 10 | 2609.30565 (known); Greaves-Tunnell thesis (screened, excluded) |
| G11 | "distance covariance" "long-range dependent" Betken | **0 parsed** | failed query |
| G12 | source:ssrn long memory independence test | 3 | Hwu–Tsay SSRN 7546859 (adjacent); 2 off-topic |
| G13 | "Hilbert-Schmidt independence criterion" "fractional" "Hurst" | **0 parsed** | failed query |
| G14 | "kernel" "Granger" "long memory" test size | 7 | none on question |
| G15 | Betken Dehling distance covariance cited | **0 parsed** | failed query |

**Google Scholar status.** Reachable (HTTP 200, no captcha), but 4 of 10 queries returned pages with no
parsable result blocks, and returned results were dominated by 2026 items. I treat Scholar coverage as
partial and unreliable. Nothing below is inferred from an empty Scholar page.

**SSRN status.** Direct fetch of `papers.ssrn.com/sol3/results.cfm?...` returned **HTTP 403**.
WebSearch with `allowed_domains: ssrn.com` returned arXiv records only. Scholar `source:ssrn`
returned 3 SSRN records, one adjacent. **SSRN remains effectively unsearched for this question.**

### Citation chasing (Semantic Scholar Graph API)

| Seed | Result |
|---|---|
| 2107.03041 (Betken–Dehling) | 1 citer: 2510.23563 (Ernst–Viens–Yan) |
| 1404.0551 (DRW) | 8 citers, all change-point / U-statistic work; none on kernel tests under LRD |
| 1509.05720 (Betken–Wendler) | 20 citers; new and relevant: 1512.00819 (Bai–Taqqu AoS 2017), 1601.01122 (Tewes JTSA 2018), 2208.01713 (Zhang–Lahiri–Nordman AoS 2022), 2609.33868 (Zhai–Cheng–Wu) |
| 2609.30565 (ZZW) | 0 citers (posted 3 weeks ago) |
| 2601.13468, 1512.00820 | **HTTP 429 rate limit**, not retrieved |

### Identifier checks (Crossref unless stated)

All DOIs in `bibliography.bib` were resolved through Crossref, except: Hall–Jing–Lahiri (1998, Statistica
Sinica, no DOI; taken from the Betken–Wendler reference list), McElroy–Politis (2007, same),
the four arXiv-only preprints (arXiv API), and Hwu–Tsay (SSRN; Scholar byline only).

Notable checks:
- **Zhai, Kolar & Wu (2026)**: IEEE TIT 72(6):4246–4296, doi 10.1109/TIT.2026.3685246 — **verified**
  (was "not verified" in the prior job). Preprint arXiv:2604.16219 (17 Apr 2026).
- **Ghoshal & Shao**: arXiv:2601.13468 and Electronic Journal of Statistics 2026, doi
  10.1214/26-EJS2488, identical title — **same work** (answers hand-off Q9). v1 abstract says "absolute
  regularity"; v2 abstract and Proposition 1 say strong (α-) mixing with Σ α(m)^{δ/(2+δ)} < ∞.
- **Betken–Wendler**: Statistica Sinica 2018, doi 10.5705/ss.202015.0435 — verified (prior bib had "not verified").
- **Betken (2016)**: JTSA 37(6):785–809 — verified.
- **Dehling–Rooch–Wendler**: Statistics 51(1):84–104 (2017) — verified; prior key `dehling2014twosample`.
- **Shao (2009)** generalized portmanteau: Econometric Theory 25(1):195–210 — verified.
- **Hand-off error corrected:** "Hall, Horowitz & Jing (1998)" for subsampling under LRD should be
  **Hall, Jing & Lahiri (1998), Statistica Sinica 8(4):1189–1204**.
- Semantic Scholar DOI lookups returned HTTP 429 for 5 of 8 calls; OpenAlex used instead.

## Failures and limits (exact)

- `papers.ssrn.com` search: HTTP 403.
- Project Euclid full text of Robinson–Hidalgo (1997): WebFetch returned an empty page.
- Wiley (Shao 2011 JTSA), ScienceDirect (Tsay–Chung 2000), Taylor & Francis (Beering 2025): paywalled
  or HTTP 403; abstracts from OpenAlex / RePEc where available (RePEc had none for Shao 2011 and TC 2000).
- Tóth et al. 1108.1632: `pdftotext` reported "xref num 524 not found"; text partly extracted; no
  numeric sign-ACF exponent found in the extracted text.
- Shao (2015) JASA review: no open copy found (OpenAlex locations list publisher only).
- Semantic Scholar: HTTP 429 on several calls (citations of 2601.13468 and 1512.00820 not obtained).
- Lahiri (1993), Kim–Nordman (2011), Leucht–Neumann (2013), Zhang–Ho–Wendler–Wu (2013), Poskitt (2007):
  no abstract in OpenAlex; content known only second-hand or not at all (labelled I or S).
