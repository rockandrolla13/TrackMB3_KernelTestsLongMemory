# EXP6 tables (main)

κ grid: [0.5, 1.0, 2.0, 4.0]. Rates in %. CPU-hours: 4.34. Errors: 0. Reps: {'G0__null__H0.5__n2000': 500, 'G0__null__H0.8__n2000': 500, 'SEAS-V__null__H0.5__n2000': 500, 'G0__lin__H0.5__n2000': 500, 'G0__sq__H0.5__n2000': 500, 'SEAS-V__lin__H0.5__n2000': 500}

## Raw null rejection rate (size), Wilson 95%

### G0__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 7.2 [5.2, 9.8] mild (McN p=0.72) | 7.8 [5.8, 10.5] mild | 5.0 [3.4, 7.3] valid (McN p=0.013) | 5.4 [3.7, 7.7] valid (McN p=0.081) |
| hsic_shift | 5.4 [3.7, 7.7] valid (McN p=1) | 5.4 [3.7, 7.7] valid | 4.4 [2.9, 6.6] valid (McN p=0.38) | 5.0 [3.4, 7.3] valid (McN p=0.86) |
| hsic_sub | 6.0 [4.2, 8.4] valid (McN p=0.7) | 6.6 [4.7, 9.1] valid | 5.8 [4.1, 8.2] valid (McN p=0.56) | 5.4 [3.7, 7.7] valid (McN p=0.39) |

κ-free references: dcov_wild 5.2, dcov_sub 5.2, pearson_nw 5.4, pearson_ewc 4.4.

### G0__H0.8__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 6.4 [4.6, 8.9] valid (McN p=0.13) | 8.0 [5.9, 10.7] mild | 8.6 [6.4, 11.4] mild (McN p=0.61) | 8.6 [6.4, 11.4] mild (McN p=0.68) |
| hsic_shift | 4.0 [2.6, 6.1] valid (McN p=0.39) | 4.8 [3.2, 7.0] valid | 5.0 [3.4, 7.3] valid (McN p=1) | 5.2 [3.6, 7.5] valid (McN p=0.8) |
| hsic_sub | 12.2 [9.6, 15.4] over (McN p=1.6e-05) | 19.2 [16.0, 22.9] over | 24.2 [20.7, 28.1] over (McN p=0.00025) | 26.4 [22.7, 30.4] over (McN p=3.2e-06) |

κ-free references: dcov_wild 8.6, dcov_sub 23.2, pearson_nw 13.6, pearson_ewc 8.2.

### SEAS-V__H0.5__n2000

| rule | κ=0.5 | κ=1 | κ=2 | κ=4 |
|---|---|---|---|---|
| hsic_wild | 72.0 [67.9, 75.8] over (McN p=0.036) | 74.4 [70.4, 78.0] over | 56.8 [52.4, 61.1] over (McN p=1.5e-25) | 6.6 [4.7, 9.1] valid (McN p=1.9e-97) |
| hsic_shift | 5.0 [3.4, 7.3] valid (McN p=0.65) | 5.6 [3.9, 8.0] valid | 5.8 [4.1, 8.2] valid (McN p=1) | 4.8 [3.2, 7.0] valid (McN p=0.67) |
| hsic_sub | 71.8 [67.7, 75.6] over (McN p=0.16) | 73.6 [69.6, 77.3] over | 52.0 [47.6, 56.3] over (McN p=2.4e-30) | 6.8 [4.9, 9.4] valid (McN p=2.2e-93) |

κ-free references: dcov_wild 53.4, dcov_sub 47.4, pearson_nw 4.2, pearson_ewc 4.8.

## Size-adjusted power, bootstrap 95% (null and alternative resampled)

Oracle-only κ (bare statistic, no rule): [0.25, 8.0].

### G0__lin0.065__H0.5__n2000

| rule | κ=0.25 | κ=0.5 | κ=1 | κ=2 | κ=4 | κ=8 |
|---|---|---|---|---|---|---|
| hsic_wild | – | 23.6 [15.2, 32.2] | 53.8 [48.8, 60.6] | 80.4 [69.8, 85.8] | 81.8 [76.8, 86.6] | – |
| hsic_shift | – | 22.8 [12.4, 32.4] | 55.0 [45.2, 61.8] | 79.0 [70.6, 83.6] | 81.8 [74.2, 86.6] | – |
| hsic_sub | – | 24.6 [13.0, 32.4] | 54.4 [46.0, 60.4] | 76.0 [68.4, 82.2] | 80.4 [74.6, 85.4] | – |
| oracle_hsic | 12.0 [8.4, 16.8] | 22.8 [12.2, 33.4] | 55.0 [44.4, 61.8] | 78.8 [71.4, 83.4] | 81.2 [74.0, 86.4] | 82.0 [75.0, 86.6] |

κ-free references: dcov_wild 77.4, dcov_sub 74.2, pearson_nw 81.4, pearson_ewc 82.0, oracle_dcov 77.6, oracle_pearson 81.2.

Paired gain vs κ = 1 (pp, bootstrap 95%):

- hsic_wild|0.5: -30.2 [-39.0, -23.2] loss
- hsic_wild|2.0: +26.6 [+16.2, +31.2] full recovery (within 5 pp of dCov wild)
- hsic_wild|4.0: +28.0 [+21.2, +33.4] full recovery (within 5 pp of dCov wild)
- hsic_shift|0.5: -32.2 [-41.6, -21.0] loss
- hsic_shift|2.0: +24.0 [+17.0, +32.4] full recovery (within 5 pp of dCov wild)
- hsic_shift|4.0: +26.8 [+19.2, +36.2] full recovery (within 5 pp of dCov wild)
- hsic_sub|0.5: -29.8 [-40.6, -21.2] loss
- hsic_sub|2.0: +21.6 [+15.2, +27.8] full recovery (within 5 pp of dCov wild)
- hsic_sub|4.0: +26.0 [+19.8, +33.4] full recovery (within 5 pp of dCov wild)
- oracle_hsic|0.25: -43.0 [-50.0, -32.6] loss
- oracle_hsic|0.5: -32.2 [-41.0, -20.8] loss
- oracle_hsic|2.0: +23.8 [+17.6, +33.2] full recovery (within 5 pp of dCov wild)
- oracle_hsic|4.0: +26.2 [+19.6, +36.2] full recovery (within 5 pp of dCov wild)
- oracle_hsic|8.0: +27.0 [+19.0, +36.2] full recovery (within 5 pp of dCov wild)

### G0__sq0.075__H0.5__n2000

| rule | κ=0.25 | κ=0.5 | κ=1 | κ=2 | κ=4 | κ=8 |
|---|---|---|---|---|---|---|
| hsic_wild | – | 22.0 [13.0, 29.0] | 37.2 [31.8, 45.2] | 41.0 [25.2, 51.6] | 7.2 [4.0, 13.8] | – |
| hsic_shift | – | 20.0 [11.2, 29.8] | 39.8 [28.2, 48.6] | 36.4 [22.2, 47.2] | 6.8 [2.8, 13.2] | – |
| hsic_sub | – | 21.8 [10.4, 27.8] | 38.8 [29.2, 49.0] | 34.8 [21.2, 49.2] | 8.2 [4.8, 14.2] | – |
| oracle_hsic | 11.8 [7.8, 15.8] | 20.4 [10.8, 30.2] | 39.0 [28.4, 48.6] | 37.4 [25.0, 49.4] | 6.8 [3.0, 14.2] | 3.6 [1.6, 7.4] |

κ-free references: dcov_wild 35.4, dcov_sub 28.2, pearson_nw 2.8, pearson_ewc 3.8, oracle_dcov 32.0, oracle_pearson 3.2.

Paired gain vs κ = 1 (pp, bootstrap 95%):

- hsic_wild|0.5: -15.2 [-24.0, -9.8] loss
- hsic_wild|2.0: +3.8 [-10.8, +12.4] small/unclear
- hsic_wild|4.0: -30.0 [-37.4, -23.0] loss
- hsic_shift|0.5: -19.8 [-29.2, -8.2] loss
- hsic_shift|2.0: -3.4 [-15.8, +8.8] small/unclear
- hsic_shift|4.0: -33.0 [-41.8, -21.8] loss
- hsic_sub|0.5: -17.0 [-28.4, -8.8] loss
- hsic_sub|2.0: -4.0 [-16.0, +7.6] small/unclear
- hsic_sub|4.0: -30.6 [-39.6, -21.6] loss
- oracle_hsic|0.25: -27.2 [-37.0, -17.4] loss
- oracle_hsic|0.5: -18.6 [-29.4, -7.6] loss
- oracle_hsic|2.0: -1.6 [-13.0, +9.6] small/unclear
- oracle_hsic|4.0: -32.2 [-41.4, -21.2] loss
- oracle_hsic|8.0: -35.4 [-44.6, -24.4] loss

### SEAS-V__lin0.065__H0.5__n2000

| rule | κ=0.25 | κ=0.5 | κ=1 | κ=2 | κ=4 | κ=8 |
|---|---|---|---|---|---|---|
| hsic_wild | – | 10.2 [7.2, 13.4] | 15.6 [12.0, 19.2] | 42.6 [36.6, 50.0] | 72.8 [67.4, 81.0] | – |
| hsic_shift | – | 10.2 [6.4, 14.6] | 15.0 [9.8, 20.0] | 41.2 [34.6, 49.8] | 69.0 [63.6, 80.4] | – |
| hsic_sub | – | 10.8 [7.0, 14.8] | 11.0 [8.0, 23.0] | 44.2 [35.6, 48.8] | 71.0 [65.6, 80.8] | – |
| oracle_hsic | 6.0 [3.6, 11.0] | 10.2 [7.0, 13.6] | 13.8 [10.0, 18.2] | 44.6 [35.2, 51.0] | 70.0 [64.4, 80.2] | 69.2 [64.6, 80.4] |

κ-free references: dcov_wild 48.2, dcov_sub 49.2, pearson_nw 75.2, pearson_ewc 73.2, oracle_dcov 47.2, oracle_pearson 70.2.

Paired gain vs κ = 1 (pp, bootstrap 95%):

- hsic_wild|0.5: -5.4 [-8.4, -2.2] small/unclear
- hsic_wild|2.0: +27.0 [+22.0, +33.8] gain
- hsic_wild|4.0: +57.2 [+51.4, +65.8] full recovery (within 5 pp of dCov wild)
- hsic_shift|0.5: -4.8 [-8.6, -1.0] small/unclear
- hsic_shift|2.0: +26.2 [+20.2, +34.4] gain
- hsic_shift|4.0: +54.0 [+47.4, +66.8] full recovery (within 5 pp of dCov wild)
- hsic_sub|0.5: -0.2 [-10.8, +2.6] small/unclear
- hsic_sub|2.0: +33.2 [+21.6, +36.8] full recovery (within 5 pp of dCov wild)
- hsic_sub|4.0: +60.0 [+47.8, +69.8] full recovery (within 5 pp of dCov wild)
- oracle_hsic|0.25: -7.8 [-11.4, -3.0] small/unclear
- oracle_hsic|0.5: -3.6 [-6.8, -0.8] small/unclear
- oracle_hsic|2.0: +30.8 [+22.0, +36.6] full recovery (within 5 pp of dCov wild)
- oracle_hsic|4.0: +56.2 [+49.8, +66.8] full recovery (within 5 pp of dCov wild)
- oracle_hsic|8.0: +55.4 [+49.6, +67.6] full recovery (within 5 pp of dCov wild)

## Mechanism diagnostics

Spearman correlation of the HSIC statistic with corr(X²,Y²)² and with r², and R² of a linear fit HSIC ~ n r² + n corr(X²,Y²)², per cell.

| cell | median σ_X | ρS(·,c²) κ=0.25 | ρS(·,r²) κ=0.25 | R² κ=0.25 | ρS(·,c²) κ=0.5 | ρS(·,r²) κ=0.5 | R² κ=0.5 | ρS(·,c²) κ=1 | ρS(·,r²) κ=1 | R² κ=1 | ρS(·,c²) κ=2 | ρS(·,r²) κ=2 | R² κ=2 | ρS(·,c²) κ=4 | ρS(·,r²) κ=4 | R² κ=4 | ρS(·,c²) κ=8 | ρS(·,r²) κ=8 | R² κ=8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G0__null__H0.5__n2000 | 0.953 | 0.01 | 0.06 | 0.02 | 0.05 | 0.17 | 0.07 | 0.06 | 0.39 | 0.27 | 0.03 | 0.70 | 0.71 | -0.01 | 0.92 | 0.96 | -0.01 | 0.99 | 1.00 |
| G0__null__H0.8__n2000 | 0.930 | 0.08 | 0.16 | 0.09 | 0.14 | 0.35 | 0.30 | 0.17 | 0.63 | 0.64 | 0.13 | 0.86 | 0.91 | 0.08 | 0.97 | 0.99 | 0.06 | 0.99 | 1.00 |
| SEAS-V__null__H0.5__n2000 | 0.936 | 0.71 | -0.04 | 0.42 | 0.82 | -0.05 | 0.56 | 0.89 | -0.04 | 0.69 | 0.88 | 0.12 | 0.75 | 0.36 | 0.73 | 0.87 | -0.00 | 0.96 | 0.99 |
| G0__lin__H0.5__n2000 | 0.953 | 0.04 | 0.25 | 0.08 | 0.06 | 0.49 | 0.27 | 0.07 | 0.75 | 0.59 | 0.06 | 0.93 | 0.88 | 0.06 | 0.99 | 0.99 | 0.07 | 1.00 | 1.00 |
| G0__sq__H0.5__n2000 | 0.952 | 0.05 | 0.03 | 0.01 | 0.03 | 0.03 | 0.01 | 0.04 | 0.12 | 0.04 | 0.10 | 0.34 | 0.25 | 0.12 | 0.71 | 0.77 | 0.09 | 0.93 | 0.98 |
| SEAS-V__lin__H0.5__n2000 | 0.934 | 0.74 | 0.10 | 0.50 | 0.80 | 0.13 | 0.60 | 0.80 | 0.27 | 0.67 | 0.52 | 0.71 | 0.78 | 0.10 | 0.97 | 0.96 | 0.03 | 1.00 | 0.99 |
