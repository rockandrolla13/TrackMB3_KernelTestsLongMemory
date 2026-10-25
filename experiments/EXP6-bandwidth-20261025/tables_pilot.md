# EXP6 tables (pilot)

κ grid: [0.5, 1.0, 2.0, 4.0]. Rates in %. CPU-hours: 0.10. Errors: 0. Reps: {'G0__null__H0.5__n2000': 8, 'G0__null__H0.8__n2000': 8, 'SEAS-V__null__H0.5__n2000': 8, 'G0__lin__H0.5__n2000': 8, 'G0__sq__H0.5__n2000': 8, 'SEAS-V__lin__H0.5__n2000': 8}

## Raw null rejection rate (size), Wilson 95%

### G0__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 0.0 [0.0, 32.4] valid (McN p=1) | 0.0 [0.0, 32.4] valid | 0.0 [0.0, 32.4] valid (McN p=1) | 0.0 [0.0, 32.4] valid (McN p=1) |
| hsic_shift | 0.0 [0.0, 32.4] valid (McN p=1) | 0.0 [0.0, 32.4] valid | 0.0 [0.0, 32.4] valid (McN p=1) | 0.0 [0.0, 32.4] valid (McN p=1) |
| hsic_sub | 0.0 [0.0, 32.4] valid (McN p=1) | 0.0 [0.0, 32.4] valid | 0.0 [0.0, 32.4] valid (McN p=1) | 0.0 [0.0, 32.4] valid (McN p=1) |

κ-free references: dcov_wild 0.0, dcov_sub 0.0, pearson_nw 0.0, pearson_ewc 0.0.

### G0__H0.8__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 37.5 [13.7, 69.4] over (McN p=1) | 25.0 [7.1, 59.1] over | 25.0 [7.1, 59.1] over (McN p=1) | 25.0 [7.1, 59.1] over (McN p=1) |
| hsic_shift | 12.5 [2.2, 47.1] valid (McN p=1) | 25.0 [7.1, 59.1] over | 25.0 [7.1, 59.1] over (McN p=1) | 25.0 [7.1, 59.1] over (McN p=1) |
| hsic_sub | 37.5 [13.7, 69.4] over (McN p=1) | 37.5 [13.7, 69.4] over | 62.5 [30.6, 86.3] over (McN p=0.62) | 50.0 [21.5, 78.5] over (McN p=1) |

κ-free references: dcov_wild 25.0, dcov_sub 37.5, pearson_nw 25.0, pearson_ewc 25.0.

### SEAS-V__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 75.0 [40.9, 92.9] over (McN p=1) | 75.0 [40.9, 92.9] over | 62.5 [30.6, 86.3] over (McN p=1) | 12.5 [2.2, 47.1] valid (McN p=0.062) |
| hsic_shift | 12.5 [2.2, 47.1] valid (McN p=1) | 0.0 [0.0, 32.4] valid | 12.5 [2.2, 47.1] valid (McN p=1) | 12.5 [2.2, 47.1] valid (McN p=1) |
| hsic_sub | 75.0 [40.9, 92.9] over (McN p=1) | 75.0 [40.9, 92.9] over | 37.5 [13.7, 69.4] over (McN p=0.25) | 12.5 [2.2, 47.1] valid (McN p=0.062) |

κ-free references: dcov_wild 62.5, dcov_sub 37.5, pearson_nw 0.0, pearson_ewc 0.0.

## Size-adjusted power, bootstrap 95% (null and alternative resampled)

Oracle-only κ (bare statistic, no rule): [].

### G0__lin0.065__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 25.0 [0.0, 62.5] | 50.0 [25.0, 87.5] | 87.5 [62.5, 100.0] | 87.5 [62.5, 100.0] |
| hsic_shift | 25.0 [0.0, 62.5] | 50.0 [25.0, 87.5] | 75.0 [50.0, 100.0] | 87.5 [62.5, 100.0] |
| hsic_sub | 25.0 [0.0, 62.5] | 50.0 [12.5, 87.5] | 75.0 [37.5, 100.0] | 87.5 [62.5, 100.0] |
| oracle_hsic | 25.0 [0.0, 62.5] | 50.0 [25.0, 87.5] | 87.5 [62.5, 100.0] | 87.5 [62.5, 100.0] |

κ-free references: dcov_wild 75.0, dcov_sub 75.0, pearson_nw 75.0, pearson_ewc 87.5, oracle_dcov 75.0, oracle_pearson 87.5.

Paired gain vs κ = 1 (pp, bootstrap 95%):

- hsic_wild|0.5: -25.0 [-62.5, +0.0] small/unclear
- hsic_wild|2.0: +37.5 [+0.0, +62.5] small/unclear
- hsic_wild|4.0: +37.5 [+0.0, +62.5] small/unclear
- hsic_shift|0.5: -25.0 [-62.5, +0.0] small/unclear
- hsic_shift|2.0: +25.0 [+0.0, +62.5] small/unclear
- hsic_shift|4.0: +37.5 [+0.0, +62.5] small/unclear
- hsic_sub|0.5: -25.0 [-62.5, +0.0] small/unclear
- hsic_sub|2.0: +25.0 [+0.0, +62.5] small/unclear
- hsic_sub|4.0: +37.5 [+0.0, +75.0] small/unclear
- oracle_hsic|0.5: -25.0 [-62.5, +0.0] small/unclear
- oracle_hsic|2.0: +37.5 [+0.0, +62.5] small/unclear
- oracle_hsic|4.0: +37.5 [+0.0, +62.5] small/unclear

### G0__sq0.075__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 50.0 [12.5, 87.5] | 75.0 [49.7, 100.0] | 87.5 [62.5, 100.0] | 62.5 [25.0, 100.0] |
| hsic_shift | 25.0 [0.0, 75.0] | 37.5 [12.5, 100.0] | 87.5 [62.5, 100.0] | 62.5 [25.0, 100.0] |
| hsic_sub | 25.0 [0.0, 75.0] | 75.0 [37.5, 100.0] | 87.5 [62.5, 100.0] | 62.5 [25.0, 100.0] |
| oracle_hsic | 25.0 [0.0, 62.5] | 37.5 [12.5, 100.0] | 87.5 [62.5, 100.0] | 62.5 [25.0, 100.0] |

κ-free references: dcov_wild 75.0, dcov_sub 87.5, pearson_nw 25.0, pearson_ewc 25.0, oracle_dcov 75.0, oracle_pearson 25.0.

Paired gain vs κ = 1 (pp, bootstrap 95%):

- hsic_wild|0.5: -25.0 [-62.5, +0.0] small/unclear
- hsic_wild|2.0: +12.5 [+0.0, +50.0] small/unclear
- hsic_wild|4.0: -12.5 [-37.5, +37.5] small/unclear
- hsic_shift|0.5: -12.5 [-62.5, +0.0] small/unclear
- hsic_shift|2.0: +50.0 [+0.0, +87.5] small/unclear
- hsic_shift|4.0: +25.0 [-37.5, +75.0] small/unclear
- hsic_sub|0.5: -50.0 [-75.0, +0.0] small/unclear
- hsic_sub|2.0: +12.5 [+0.0, +37.5] small/unclear
- hsic_sub|4.0: -12.5 [-37.5, +25.0] small/unclear
- oracle_hsic|0.5: -12.5 [-75.0, +0.0] small/unclear
- oracle_hsic|2.0: +50.0 [+0.0, +87.5] small/unclear
- oracle_hsic|4.0: +25.0 [-37.5, +75.0] small/unclear

### SEAS-V__lin0.065__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 12.5 [0.0, 50.0] | 12.5 [0.0, 50.0] | 75.0 [50.0, 100.0] | 75.0 [50.0, 100.0] |
| hsic_shift | 12.5 [0.0, 37.5] | 12.5 [0.0, 50.0] | 37.5 [12.5, 100.0] | 50.0 [25.0, 100.0] |
| hsic_sub | 12.5 [0.0, 50.0] | 25.0 [0.0, 75.0] | 62.5 [37.5, 100.0] | 62.5 [37.5, 100.0] |
| oracle_hsic | 12.5 [0.0, 50.0] | 12.5 [0.0, 50.0] | 25.0 [0.0, 100.0] | 50.0 [25.0, 100.0] |

κ-free references: dcov_wild 50.0, dcov_sub 37.5, pearson_nw 75.0, pearson_ewc 75.0, oracle_dcov 25.0, oracle_pearson 62.5.

Paired gain vs κ = 1 (pp, bootstrap 95%):

- hsic_wild|0.5: +0.0 [+0.0, +0.0] small/unclear
- hsic_wild|2.0: +62.5 [+25.0, +100.0] full recovery (within 5 pp of dCov wild)
- hsic_wild|4.0: +62.5 [+12.5, +100.0] full recovery (within 5 pp of dCov wild)
- hsic_shift|0.5: +0.0 [-25.0, +0.0] small/unclear
- hsic_shift|2.0: +25.0 [-12.5, +87.5] small/unclear
- hsic_shift|4.0: +37.5 [-12.5, +100.0] small/unclear
- hsic_sub|0.5: -12.5 [-50.0, +0.0] small/unclear
- hsic_sub|2.0: +37.5 [+0.0, +87.5] small/unclear
- hsic_sub|4.0: +37.5 [-12.5, +87.5] small/unclear
- oracle_hsic|0.5: +0.0 [+0.0, +0.0] small/unclear
- oracle_hsic|2.0: +12.5 [-25.0, +87.5] small/unclear
- oracle_hsic|4.0: +37.5 [-12.5, +100.0] small/unclear

## Mechanism diagnostics

Spearman correlation of the HSIC statistic with corr(X²,Y²)² and with r², and R² of a linear fit HSIC ~ n r² + n corr(X²,Y²)², per cell.

| cell | median σ_X | ρS(·,c²) κ=0.5 | ρS(·,r²) κ=0.5 | R² κ=0.5 | ρS(·,c²) κ=1 | ρS(·,r²) κ=1 | R² κ=1 | ρS(·,c²) κ=2 | ρS(·,r²) κ=2 | R² κ=2 | ρS(·,c²) κ=4 | ρS(·,r²) κ=4 | R² κ=4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G0__null__H0.5__n2000 | 0.957 | -0.07 | 0.33 | 0.15 | -0.07 | 0.33 | 0.17 | -0.24 | 0.31 | 0.48 | -0.21 | 0.76 | 0.98 |
| G0__null__H0.8__n2000 | 0.928 | -0.05 | -0.17 | 0.02 | -0.48 | 0.21 | 0.54 | -0.71 | 0.74 | 0.97 | -0.38 | 0.95 | 1.00 |
| SEAS-V__null__H0.5__n2000 | 0.942 | 0.90 | 0.10 | 0.88 | 0.93 | 0.19 | 0.91 | 0.93 | 0.19 | 0.91 | 0.14 | 0.69 | 0.94 |
| G0__lin__H0.5__n2000 | 0.957 | 0.26 | 0.31 | 0.16 | -0.07 | 0.81 | 0.67 | -0.05 | 0.90 | 0.91 | -0.19 | 1.00 | 0.99 |
| G0__sq__H0.5__n2000 | 0.970 | -0.24 | 0.02 | 0.13 | -0.24 | 0.24 | 0.04 | -0.21 | 0.52 | 0.25 | 0.14 | 0.98 | 0.93 |
| SEAS-V__lin__H0.5__n2000 | 0.947 | 0.74 | -0.43 | 0.79 | 0.74 | -0.33 | 0.87 | 0.07 | 0.76 | 0.90 | -0.55 | 1.00 | 0.98 |
