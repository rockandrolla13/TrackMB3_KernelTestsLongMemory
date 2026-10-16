# Published memory estimates for market series

Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 · plan item (4) · 2026-10-15.
Conversions: H = d + 1/2; for autocorrelation ∝ k^{−γ}, H = 1 − γ/2. Boundary: H = 3/4 (d = 1/4).
Evidence class per row. "Secondhand" means the number was read in a citing source, not the original.

## Table

| Series | Market, sample | Estimator | Reported | H equivalent | Side of 3/4 | Source (class) |
|---|---|---|---|---|---|---|
| Market-order signs | LSE (20 stocks and Lo R/S per a citing snippet) | ACF power law | γ ≈ 0.6 | ≈ 0.70 | below | Lillo & Farmer, SNDE 8(3) 2004, cond-mat/0311053 (S, abstract quoted) |
| Market-order signs | Paris Bourse | ACF power law | γ from 0.2 to 0.7 | 0.65–0.90 | **both** | Bouchaud et al. 2004, secondhand via a search-result snippet of a citing paper; citing source not pinned (A, weakest row) |
| Core (metaorder) order flow | model + data | theory/calibration | H₀ ≈ 3/4 | ≈ 0.75 | **on the line** | Muhle-Karbe et al., arXiv:2601.23172 (R/F, prior job) |
| Bid–ask spread | LSE | DFA / R/S | 0.75 < H < 0.85 | 0.75–0.85 | above | Mike & Farmer, secondhand via Gu–Chen–Zhou EPJB 57 (2007) (A) |
| Bid–ask spread | Shenzhen SE, one stock | DFA | long memory, no number in abstract | — | — | Gu, Chen & Zhou, EPJB 57:81–87, physics/0701017 (S) |
| Daily trading volume | 30 DJIA stocks, 1962–1994 | tapered two-step semiparametric | d = 0.30 ± 0.08 | ≈ 0.80 | above | Lobato & Velasco, JBES 18(4) 2000, secondhand via arXiv 0904.1042 snippet (A) |
| Log trade durations (diurnally adjusted) | NYSE TAQ 2002, 10 stocks | GPH, m = n^0.4 / n^0.5 | d = 0.26–0.46 | 0.76–0.96 | above (all 10) | Deo, Hsieh & Hurvich 2005, NYU SOR-2005-2, Table 2 (**F**) |
| Log durations, parametric | same | LMSD Whittle | d = 0.35–0.42 | 0.85–0.92 | above | same, Table 8 text (F) |
| Trade counts, 5 min / 30 min / 1 day | same | GPH, m₁ and m = n^0.8 | d = 0.22–0.54 | 0.72–1.04 | mostly above | same, Table 3 (F) |
| Squared returns, 5 min / 30 min / 1 day | same | GPH | d = 0.02–0.44; with m = n^0.8 mostly 0.05–0.29 | 0.52–0.94 | **both; bandwidth-dependent** | same, Table 4 (F) |
| Daily log realised volatility (5-min) | same, **n = 249 days** | GPH, m = n^0.6 / n^0.8 | d = 0.28–0.87 (s.e. 0.15 / 0.08) | 0.78–1.37 | above (non-stationary d > 0.5 in several) | same, Table 5 (F) |
| Realised volatility, FX | DM/$, ¥/$ | GPH-type | d ≈ 0.4 | ≈ 0.90 | above | Andersen, Bollerslev, Diebold & Labys, Econometrica 2003 (R/S, prior job) |
| Log-volatility increments | indices, many scales | scaling of moments | H ≈ 0.1 (rough, not long memory) | n/a | **contrary framing** | Gatheral, Jaisson & Rosenbaum, arXiv:1410.3394 (A) |

## What the numbers say

1. **Order signs straddle the line.** LSE estimates put them just below (H ≈ 0.7). The Paris range
   runs from 0.65 to 0.9. The metaorder ("core") flow theory puts them on it (H₀ ≈ 3/4). Calling
   order signs "safely below" is not supported.
2. **Everything activity-like sits above it.** Durations, counts, volume and spreads are all
   reported at H ≈ 0.75–0.96. The Deo–Hsieh–Hurvich durations are above 0.75 in all ten stocks.
3. **Volatility is above it, or the question is ill-posed.** Classical estimates give d ≈ 0.3–0.45.
   The rough-volatility view says log-volatility is not long memory at all (H ≈ 0.1 for an fBm-
   like log-vol). It adds that standard long-memory tests applied to such data still report long
   memory (abstract). For inference this matters as finite-sample behaviour: a rough series with
   slow mean reversion can produce the same size distortion over realistic horizons. The
   simulation includes a rough-volatility arm for this reason.
4. **For two-series tests, use the sum rule.** Order signs (0.7) paired with RV, volume or durations
   (≥ 0.8) gives H_X + H_Y ≥ 1.5. That is at or beyond the cross-statistic boundary (boundary-
   theorems.md §3). So most microstructure *pairs* of interest sit outside the proved √n regime.
5. **Estimates are fragile.** Deo et al. show d moving by 0.1–0.3 with GPH bandwidth for squared
   returns. Their RV sample is one year. Any boundary claim about real data needs an uncertainty
   band on H, not a point value.

## Gaps in this table (not searched or not reached within budget)

- Order-sign exponents from the Tóth, Palit, Lillo & Farmer JEDC 2015 paper (arXiv:1108.1632) —
  located, number not read.
- Plerou et al. PRE 71 (2005) spread exponents — only secondhand.
- Crypto and futures LOB data; post-2015 equity order-sign estimates.
- The programme's own Track A propagator exponent (≈ 0.22, R) implies strong persistence of
  impact-filtered flow; not converted here because the mapping to H depends on the model.
