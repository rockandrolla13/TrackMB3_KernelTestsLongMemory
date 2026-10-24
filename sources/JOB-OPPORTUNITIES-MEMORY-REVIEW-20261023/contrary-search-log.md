# Contrary search log — JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023

Execution `opp-memory-review-20261023`. All searches ran on 2026-10-23 between 22:03 and 22:20 UTC.
Tools: WebSearch (standard and extended), WebFetch, and `curl` against the Semantic Scholar Graph API, OpenAlex and the HAL search API.
Downloads went to the session scratchpad, not to this directory. A search that finds nothing is a search result, not a novelty claim.

## Forward-citation sweeps

| # | Index | Query | Candidate | Result |
|---|---|---|---|---|
| C1 | Semantic Scholar | `paper/arXiv:2107.03041/citations` | G1, MAP | 1 citing paper: Ernst–Viens–Yan 2510.23563, non-stationary paths (abstract opened). `citationCount` = 1. |
| C2 | OpenAlex | `works?search=Betken Dehling independence long-range dependent` | G1 | W4390632085, `cited_by_count` 0 |
| C3 | Semantic Scholar | `paper/DOI:10.3982/ECTA15047/citations` (Müller–Watson 2018) | G2, MAP | 86 records. Keyword filter on abstracts (nonlinear, kernel, independen, Hermite, Rosenblatt, long memory) found no kernel or nonlinear extension. |
| C4 | Semantic Scholar | `paper/arXiv:1109.5298/citations` (Kulik–Soulier) | G3 | 15 records, mostly extremes and SLLN. None revisits power ranks under leverage. |
| C5 | Semantic Scholar | `paper/arXiv:2609.33868/citations` (Zhai–Cheng–Wu) | MAP | 0 records. Metadata call rate-limited (429), so indexing is unverified. |

## Retrieval of the blocked item

| # | Route | Result |
|---|---|---|
| R1 | HAL API `halId_s:hal-01197638` | Identified it: Lévy-Leduc & Taqqu, Metrika 2014, DOI 10.1007/s00184-013-0474-4. The abstract covers a single Gaussian sequence. |
| R2 | hal.science page and `/document` | Anubis bot challenge |
| R3 | Springer article page | Login redirect |
| R4 | WebSearch "Hermite ranks and U-statistics" arXiv | No arXiv copy. Led to 1404.0551 and 0912.4688. |
| R5 | WebSearch Lévy-Leduc Taqqu Metrika pdf | No open PDF |

## Contrary and equivalent-formulation searches

| # | Query (exact) | Mode | Candidate | Result and decision |
|---|---|---|---|---|
| S1 | Bai Taqqu "Multivariate limit theorems in the context of long-range dependence" short-range long-range mixed components | std | G1 | Found Düker 1704.08609 (kept, decisive for the linear case) and Bai–Taqqu 1211.0576 (one common series, so excluded for G1) |
| S2 | HSIC OR "kernel independence test" long-range dependent time series asymptotic null distribution 2025 2026 arXiv Hurst | ext | G1 | No HSIC-under-LRD null theory. 2605.22025 and 2604.28104 are weak dependence, excluded. |
| S3 | low-frequency cosine weighted averages Hermite rank two non-Gaussian Rosenblatt limit long memory fixed-q inference nonlinear transformation | ext | G2 | No combination found. Located the Müller–Watson handbook chapter (kept). |
| S4 | Surgailis Viano 2002 long memory properties volatility absolute powers stochastic volatility leverage Hermite rank finite variance | std | G3 | Surgailis–Viano ESAIM P&S 2002 located, not opened |
| S5 | testing Hermite rank empirical long memory transformed series absolute returns rank estimation Beran | std | G3 | Beran–Möhrle–Ghosh rank test (abstract). Iowa State thesis on log returns only. |
| S6 | Pötscher 2002 "Lower risk bounds and properties of confidence sets for ill-posed estimation problems" pdf long memory parameter | std | UNIFORM | Only the vie0202 working paper, whose fetch returns HTML |
| S7 | test of the mean uniformly valid over fractional integration parameter d long memory size control impossibility OR "uniformly valid" confidence interval Hurst | ext | UNIFORM | No uniform-in-d result found |
| S8 | Hawkes process power-law kernel long-range dependence counts Hurst exponent heavy-tailed kernel Jaisson Rosenbaum nearly unstable fractional | std | HAWKES | Jaisson–Rosenbaum 1504.03100, H = α − 1/2 ∈ (0, 1/2) (kept) |
| S9 | Granger causality test Hawkes process size distortion long memory OR power-law kernel misspecification exponential kernel simulation false positive | std | HAWKES | Filimonov–Sornette (known). Diouane–Schoenberg–Mohler spurious cross-excitation near criticality (snippet). No H-indexed size study. |
| S10 | Monte Carlo study size of statistical tests under long-range dependence across Hurst grid kernel test causal inference point process benchmark "long memory" inference breakdown | std | MAP | Single-class size studies only (1506.07296, 1509.05720). No cross-class map. |

## Searches per candidate

- G1: C1, C2, R1–R5, S1, S2
- G2: C3, S3
- UNIFORM: S6, S7
- G3: C4, S4, S5
- HAWKES: S8, S9
- MAP: C1, C3, C5, S10

## Coverage limits

- No Google Scholar and no SSRN. Semantic Scholar rate-limited several calls.
- The arXiv export API returned empty responses from the shell. arXiv PDFs were fetched directly or through WebFetch.
