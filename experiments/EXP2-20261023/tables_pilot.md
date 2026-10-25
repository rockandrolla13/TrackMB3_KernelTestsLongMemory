# EXP2 tables (pilot)

Rejection % at nominal 5% [Wilson 95%]. Rows = generator, columns = H (H_X = H_Y for independence).

## mmd_unc_wild, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.0 [2.8, 12.5] | 23.0 [15.8, 32.2] | 16.0 [10.1, 24.4] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 6.0 [2.8, 12.5] | 12.0 [7.0, 19.8] | 21.0 [14.2, 30.0] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [2.2, 11.2] | 16.0 [10.1, 24.4] | 24.0 [16.7, 33.2] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 7.0 [3.4, 13.7] | 17.0 [10.9, 25.5] | 23.0 [15.8, 32.2] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 10.0 [5.5, 17.4] | 19.0 [12.5, 27.8] | 19.0 [12.5, 27.8] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 1.0 [0.2, 5.4] | 2.0 [0.6, 7.0] | 12.0 [7.0, 19.8] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 2.0 [0.6, 7.0] | 9.0 [4.8, 16.2] | 20.0 [13.3, 28.9] | 100/100/100 |

## mmd_unc_wild, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 3.0 [1.0, 8.5] | 15.0 [9.3, 23.3] | 25.0 [17.5, 34.3] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 6.0 [2.8, 12.5] | 15.0 [9.3, 23.3] | 31.0 [22.8, 40.6] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [2.2, 11.2] | 16.0 [10.1, 24.4] | 41.0 [31.9, 50.8] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 11.0 [6.3, 18.6] | 18.0 [11.7, 26.7] | 22.0 [15.0, 31.1] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 11.0 [6.3, 18.6] | 18.0 [11.7, 26.7] | 22.0 [15.0, 31.1] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 3.7] | 17.0 [10.9, 25.5] | 23.0 [15.8, 32.2] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 3.0 [1.0, 8.5] | 15.0 [9.3, 23.3] | 24.0 [16.7, 33.2] | 100/100/100 |

## mmd_sub_sqrt, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.0 [2.2, 11.2] | 29.0 [21.0, 38.5] | 39.0 [30.0, 48.8] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 8.0 [4.1, 15.0] | 20.0 [13.3, 28.9] | 37.0 [28.2, 46.8] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 6.0 [2.8, 12.5] | 28.0 [20.1, 37.5] | 49.0 [39.4, 58.7] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 17.0 [10.9, 25.5] | 35.0 [26.4, 44.7] | 50.0 [40.4, 59.6] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 10.0 [5.5, 17.4] | 36.0 [27.3, 45.8] | 44.0 [34.7, 53.8] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 1.0 [0.2, 5.4] | 7.0 [3.4, 13.7] | 17.0 [10.9, 25.5] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 1.0 [0.2, 5.4] | 15.0 [9.3, 23.3] | 34.0 [25.5, 43.7] | 100/100/100 |

## mmd_sub_sqrt, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 8.0 [4.1, 15.0] | 31.0 [22.8, 40.6] | 51.0 [41.3, 60.6] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 5.0 [2.2, 11.2] | 31.0 [22.8, 40.6] | 59.0 [49.2, 68.1] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [1.6, 9.8] | 43.0 [33.7, 52.8] | 58.0 [48.2, 67.2] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 19.0 [12.5, 27.8] | 43.0 [33.7, 52.8] | 45.0 [35.6, 54.8] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 11.0 [6.3, 18.6] | 33.0 [24.6, 42.7] | 38.0 [29.1, 47.8] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 3.7] | 3.0 [1.0, 8.5] | 19.0 [12.5, 27.8] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 0.0 [0.0, 3.7] | 26.0 [18.4, 35.4] | 39.0 [30.0, 48.8] | 100/100/100 |

## mmd_sn_block, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.0 [2.2, 11.2] | 12.0 [7.0, 19.8] | 5.0 [2.2, 11.2] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 8.0 [4.1, 15.0] | 8.0 [4.1, 15.0] | 6.0 [2.8, 12.5] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 6.0 [2.8, 12.5] | 5.0 [2.2, 11.2] | 8.0 [4.1, 15.0] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.0 [1.6, 9.8] | 12.0 [7.0, 19.8] | 17.0 [10.9, 25.5] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 0.0 [0.0, 3.7] | 3.0 [1.0, 8.5] | 3.0 [1.0, 8.5] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 3.7] | 0.0 [0.0, 3.7] | 3.0 [1.0, 8.5] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 2.0 [0.6, 7.0] | 4.0 [1.6, 9.8] | 7.0 [3.4, 13.7] | 100/100/100 |

## mmd_sn_block, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.0 [2.2, 11.2] | 7.0 [3.4, 13.7] | 6.0 [2.8, 12.5] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 7.0 [3.4, 13.7] | 9.0 [4.8, 16.2] | 16.0 [10.1, 24.4] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [1.6, 9.8] | 7.0 [3.4, 13.7] | 14.0 [8.5, 22.1] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 7.0 [3.4, 13.7] | 11.0 [6.3, 18.6] | 12.0 [7.0, 19.8] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 6.0 [2.8, 12.5] | 5.0 [2.2, 11.2] | 0.0 [0.0, 3.7] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 3.7] | 0.0 [0.0, 3.7] | 0.0 [0.0, 3.7] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 0.0 [0.0, 3.7] | 2.0 [0.6, 7.0] | 3.0 [1.0, 8.5] | 100/100/100 |

## hsic_wild, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 10.0 [5.5, 17.4] | 4.0 [1.6, 9.8] | 11.0 [6.3, 18.6] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 5.0 [2.2, 11.2] | 2.0 [0.6, 7.0] | 8.0 [4.1, 15.0] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 7.0 [3.4, 13.7] | 7.0 [3.4, 13.7] | 7.0 [3.4, 13.7] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 6.0 [2.8, 12.5] | 7.0 [3.4, 13.7] | 12.0 [7.0, 19.8] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 14.0 [8.5, 22.1] | 15.0 [9.3, 23.3] | 11.0 [6.3, 18.6] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 84.0 [75.6, 89.9] | 71.0 [61.5, 79.0] | 71.0 [61.5, 79.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 45.0 [35.6, 54.8] | 43.0 [33.7, 52.8] | 31.0 [22.8, 40.6] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## hsic_wild, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 9.0 [4.8, 16.2] | 6.0 [2.8, 12.5] | 9.0 [4.8, 16.2] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 7.0 [3.4, 13.7] | 6.0 [2.8, 12.5] | 5.0 [2.2, 11.2] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 6.0 [2.8, 12.5] | 9.0 [4.8, 16.2] | 13.0 [7.8, 21.0] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 8.0 [4.1, 15.0] | 5.0 [2.2, 11.2] | 10.0 [5.5, 17.4] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 5.0 [2.2, 11.2] | 8.0 [4.1, 15.0] | 14.0 [8.5, 22.1] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 89.0 [81.4, 93.7] | 92.0 [85.0, 95.9] | 87.0 [79.0, 92.2] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 77.0 [67.8, 84.2] | 72.0 [62.5, 79.9] | 66.0 [56.3, 74.5] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## hsic_shift, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 7.0 [3.4, 13.7] | 4.0 [1.6, 9.8] | 5.0 [2.2, 11.2] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 3.0 [1.0, 8.5] | 0.0 [0.0, 3.7] | 5.0 [2.2, 11.2] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 3.0 [1.0, 8.5] | 1.0 [0.2, 5.4] | 2.0 [0.6, 7.0] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 3.0 [1.0, 8.5] | 4.0 [1.6, 9.8] | 5.0 [2.2, 11.2] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 7.0 [3.4, 13.7] | 7.0 [3.4, 13.7] | 6.0 [2.8, 12.5] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 9.0 [4.8, 16.2] | 5.0 [2.2, 11.2] | 1.0 [0.2, 5.4] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 3.0 [1.0, 8.5] | 10.0 [5.5, 17.4] | 8.0 [4.1, 15.0] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 18.0 [11.7, 26.7] | 16.0 [10.1, 24.4] | 21.0 [14.2, 30.0] | 100/100/100 |

## hsic_shift, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.0 [2.2, 11.2] | 6.0 [2.8, 12.5] | 1.0 [0.2, 5.4] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 8.0 [4.1, 15.0] | 1.0 [0.2, 5.4] | 1.0 [0.2, 5.4] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [1.6, 9.8] | 7.0 [3.4, 13.7] | 8.0 [4.1, 15.0] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.0 [1.6, 9.8] | 3.0 [1.0, 8.5] | 8.0 [4.1, 15.0] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 4.0 [1.6, 9.8] | 3.0 [1.0, 8.5] | 9.0 [4.8, 16.2] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 1.0 [0.2, 5.4] | 4.0 [1.6, 9.8] | 2.0 [0.6, 7.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 4.0 [1.6, 9.8] | 5.0 [2.2, 11.2] | 6.0 [2.8, 12.5] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 38.0 [29.1, 47.8] | 35.0 [26.4, 44.7] | 19.0 [12.5, 27.8] | 100/100/100 |

## hsic_sub, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 8.0 [4.1, 15.0] | 7.0 [3.4, 13.7] | 19.0 [12.5, 27.8] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 3.0 [1.0, 8.5] | 1.0 [0.2, 5.4] | 14.0 [8.5, 22.1] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [1.6, 9.8] | 9.0 [4.8, 16.2] | 21.0 [14.2, 30.0] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.0 [1.6, 9.8] | 7.0 [3.4, 13.7] | 10.0 [5.5, 17.4] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 14.0 [8.5, 22.1] | 23.0 [15.8, 32.2] | 29.0 [21.0, 38.5] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 84.0 [75.6, 89.9] | 72.0 [62.5, 79.9] | 78.0 [68.9, 85.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 37.0 [28.2, 46.8] | 39.0 [30.0, 48.8] | 46.0 [36.6, 55.7] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## hsic_sub, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.0 [2.8, 12.5] | 8.0 [4.1, 15.0] | 17.0 [10.9, 25.5] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 8.0 [4.1, 15.0] | 7.0 [3.4, 13.7] | 12.0 [7.0, 19.8] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [2.2, 11.2] | 11.0 [6.3, 18.6] | 36.0 [27.3, 45.8] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 1.0 [0.2, 5.4] | 5.0 [2.2, 11.2] | 11.0 [6.3, 18.6] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 6.0 [2.8, 12.5] | 15.0 [9.3, 23.3] | 27.0 [19.3, 36.4] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 90.0 [82.6, 94.5] | 93.0 [86.3, 96.6] | 87.0 [79.0, 92.2] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 75.0 [65.7, 82.5] | 72.0 [62.5, 79.9] | 70.0 [60.4, 78.1] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## dcov_wild, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 9.0 [4.8, 16.2] | 4.0 [1.6, 9.8] | 14.0 [8.5, 22.1] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 11.0 [6.3, 18.6] | 3.0 [1.0, 8.5] | 10.0 [5.5, 17.4] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [2.2, 11.2] | 11.0 [6.3, 18.6] | 9.0 [4.8, 16.2] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 11.0 [6.3, 18.6] | 5.0 [2.2, 11.2] | 6.0 [2.8, 12.5] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 10.0 [5.5, 17.4] | 14.0 [8.5, 22.1] | 9.0 [4.8, 16.2] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 83.0 [74.5, 89.1] | 71.0 [61.5, 79.0] | 70.0 [60.4, 78.1] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 18.0 [11.7, 26.7] | 20.0 [13.3, 28.9] | 18.0 [11.7, 26.7] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## dcov_wild, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 8.0 [4.1, 15.0] | 6.0 [2.8, 12.5] | 12.0 [7.0, 19.8] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 5.0 [2.2, 11.2] | 8.0 [4.1, 15.0] | 4.0 [1.6, 9.8] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 3.0 [1.0, 8.5] | 13.0 [7.8, 21.0] | 14.0 [8.5, 22.1] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.0 [1.6, 9.8] | 4.0 [1.6, 9.8] | 7.0 [3.4, 13.7] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 7.0 [3.4, 13.7] | 6.0 [2.8, 12.5] | 13.0 [7.8, 21.0] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 92.0 [85.0, 95.9] | 93.0 [86.3, 96.6] | 88.0 [80.2, 93.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 52.0 [42.3, 61.5] | 45.0 [35.6, 54.8] | 40.0 [30.9, 49.8] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## dcov_sub, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 10.0 [5.5, 17.4] | 8.0 [4.1, 15.0] | 21.0 [14.2, 30.0] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 2.0 [0.6, 7.0] | 2.0 [0.6, 7.0] | 13.0 [7.8, 21.0] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 1.0 [0.2, 5.4] | 9.0 [4.8, 16.2] | 18.0 [11.7, 26.7] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 3.0 [1.0, 8.5] | 2.0 [0.6, 7.0] | 3.0 [1.0, 8.5] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 13.0 [7.8, 21.0] | 27.0 [19.3, 36.4] | 37.0 [28.2, 46.8] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 84.0 [75.6, 89.9] | 74.0 [64.6, 81.6] | 78.0 [68.9, 85.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 8.0 [4.1, 15.0] | 22.0 [15.0, 31.1] | 30.0 [21.9, 39.6] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## dcov_sub, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.0 [2.2, 11.2] | 7.0 [3.4, 13.7] | 22.0 [15.0, 31.1] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 4.0 [1.6, 9.8] | 7.0 [3.4, 13.7] | 16.0 [10.1, 24.4] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 3.0 [1.0, 8.5] | 12.0 [7.0, 19.8] | 32.0 [23.7, 41.7] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 2.0 [0.6, 7.0] | 0.0 [0.0, 3.7] | 2.0 [0.6, 7.0] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 5.0 [2.2, 11.2] | 15.0 [9.3, 23.3] | 32.0 [23.7, 41.7] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.0 [83.8, 95.2] | 92.0 [85.0, 95.9] | 88.0 [80.2, 93.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 45.0 [35.6, 54.8] | 46.0 [36.6, 55.7] | 52.0 [42.3, 61.5] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## pearson_nw, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.0 [2.8, 12.5] | 5.0 [2.2, 11.2] | 16.0 [10.1, 24.4] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 8.0 [4.1, 15.0] | 1.0 [0.2, 5.4] | 13.0 [7.8, 21.0] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [1.6, 9.8] | 8.0 [4.1, 15.0] | 11.0 [6.3, 18.6] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 6.0 [2.8, 12.5] | 5.0 [2.2, 11.2] | 7.0 [3.4, 13.7] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 5.0 [2.2, 11.2] | 13.0 [7.8, 21.0] | 9.0 [4.8, 16.2] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 83.0 [74.5, 89.1] | 71.0 [61.5, 79.0] | 71.0 [61.5, 79.0] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 5.0 [2.2, 11.2] | 8.0 [4.1, 15.0] | 16.0 [10.1, 24.4] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## pearson_nw, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 3.0 [1.0, 8.5] | 8.0 [4.1, 15.0] | 17.0 [10.9, 25.5] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 1.0 [0.2, 5.4] | 9.0 [4.8, 16.2] | 7.0 [3.4, 13.7] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [2.2, 11.2] | 7.0 [3.4, 13.7] | 18.0 [11.7, 26.7] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.0 [1.6, 9.8] | 4.0 [1.6, 9.8] | 2.0 [0.6, 7.0] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 5.0 [2.2, 11.2] | 8.0 [4.1, 15.0] | 12.0 [7.0, 19.8] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 92.0 [85.0, 95.9] | 91.0 [83.8, 95.2] | 89.0 [81.4, 93.7] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 6.0 [2.8, 12.5] | 5.0 [2.2, 11.2] | 21.0 [14.2, 30.0] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## pearson_ewc, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.0 [2.8, 12.5] | 3.0 [1.0, 8.5] | 9.0 [4.8, 16.2] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 7.0 [3.4, 13.7] | 2.0 [0.6, 7.0] | 10.0 [5.5, 17.4] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 6.0 [2.8, 12.5] | 7.0 [3.4, 13.7] | 3.0 [1.0, 8.5] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.0 [1.6, 9.8] | 4.0 [1.6, 9.8] | 3.0 [1.0, 8.5] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 3.0 [1.0, 8.5] | 7.0 [3.4, 13.7] | 7.0 [3.4, 13.7] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 82.0 [73.3, 88.3] | 70.0 [60.4, 78.1] | 66.0 [56.3, 74.5] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 3.0 [1.0, 8.5] | 6.0 [2.8, 12.5] | 9.0 [4.8, 16.2] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

## pearson_ewc, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 3.0 [1.0, 8.5] | 2.0 [0.6, 7.0] | 11.0 [6.3, 18.6] | 100/100/100 |
| HT-S | yes (instantaneous, rank 1) | 2.0 [0.6, 7.0] | 4.0 [1.6, 9.8] | 2.0 [0.6, 7.0] | 100/100/100 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [2.2, 11.2] | 5.0 [2.2, 11.2] | 10.0 [5.5, 17.4] | 100/100/100 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.0 [2.2, 11.2] | 2.0 [0.6, 7.0] | 1.0 [0.2, 5.4] | 100/100/100 |
| NL-TAR | no (nonlinear path functional) | 5.0 [2.2, 11.2] | 7.0 [3.4, 13.7] | 9.0 [4.8, 16.2] | 100/100/100 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.0 [83.8, 95.2] | 92.0 [85.0, 95.9] | 89.0 [81.4, 93.7] | 100/100/100 |
| SEAS-V | no (deterministic periodic scale, random phase) | 8.0 [4.1, 15.0] | 4.0 [1.6, 9.8] | 14.0 [8.5, 22.1] | 100/100/100 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100.0 [96.3, 100.0] | 100/100/100 |

# Verdicts (pre-declared rule, DESIGN.md)

alone = F(0.5) − G0(0.5) in pp; I(h) = (F(h) − F(0.5)) − (G0(h) − G0(0.5)) in pp, z in brackets.

## n = 500

| feature | generator | rule | alone pp | I(0.7) | I(0.8) | verdict |
|---|---|---|---|---|---|---|
| heavy tails | HT-S | mmd_unc_wild | +0.0 | -11.0 (-1.8) | +5.0 (+0.8) | no effect |
| heavy tails | HT-S | mmd_sub_sqrt | +3.0 | -12.0 (-1.7) | -5.0 (-0.7) | no effect |
| heavy tails | HT-S | mmd_sn_block | +3.0 | -7.0 (-1.3) | -2.0 (-0.4) | no effect |
| heavy tails | HT-S | hsic_wild | -5.0 | +3.0 (+0.7) | +2.0 (+0.4) | no effect |
| heavy tails | HT-S | hsic_shift | -4.0 | +0.0 (+0.0) | +4.0 (+0.9) | no effect |
| heavy tails | HT-S | hsic_sub | -5.0 | -1.0 (-0.2) | +0.0 (+0.0) | no effect |
| heavy tails | HT-S | dcov_wild | +2.0 | -3.0 (-0.6) | -6.0 (-1.0) | no effect |
| heavy tails | HT-S | dcov_sub | -8.0 | +2.0 (+0.4) | +0.0 (+0.0) | no effect |
| heavy tails | HT-L | mmd_unc_wild | -1.0 | -6.0 (-0.9) | +9.0 (+1.4) | no effect |
| heavy tails | HT-L | mmd_sub_sqrt | +1.0 | -2.0 (-0.3) | +9.0 (+1.2) | no effect |
| heavy tails | HT-L | mmd_sn_block | +1.0 | -8.0 (-1.6) | +2.0 (+0.4) | no effect |
| heavy tails | HT-L | hsic_wild | -3.0 | +6.0 (+1.2) | -1.0 (-0.2) | no effect |
| heavy tails | HT-L | hsic_shift | -4.0 | +1.0 (+0.3) | +1.0 (+0.2) | no effect |
| heavy tails | HT-L | hsic_sub | -4.0 | +6.0 (+1.2) | +6.0 (+0.9) | no effect |
| heavy tails | HT-L | dcov_wild | -4.0 | +11.0 (+2.1) | -1.0 (-0.2) | only with memory |
| heavy tails | HT-L | dcov_sub | -9.0 (below G0) | +10.0 (+2.0) | +6.0 (+0.9) | only with memory |
| nonlinear dynamics | NL-SV | mmd_unc_wild | +1.0 | -7.0 (-1.1) | +6.0 (+0.9) | no effect |
| nonlinear dynamics | NL-SV | mmd_sub_sqrt | +12.0 | -6.0 (-0.8) | -1.0 (-0.1) | no effect |
| nonlinear dynamics | NL-SV | mmd_sn_block | -1.0 | +1.0 (+0.2) | +13.0 (+2.5) | only with memory |
| nonlinear dynamics | NL-SV | hsic_wild | -4.0 | +7.0 (+1.4) | +5.0 (+0.8) | no effect |
| nonlinear dynamics | NL-SV | hsic_shift | -4.0 | +4.0 (+1.0) | +4.0 (+0.9) | no effect |
| nonlinear dynamics | NL-SV | hsic_sub | -4.0 | +4.0 (+0.8) | -5.0 (-0.8) | no effect |
| nonlinear dynamics | NL-SV | dcov_wild | +2.0 | -1.0 (-0.2) | -10.0 (-1.7) | no effect |
| nonlinear dynamics | NL-SV | dcov_sub | -7.0 | +1.0 (+0.2) | -11.0 (-2.0) | dampens memory |
| nonlinear dynamics | NL-TAR | mmd_unc_wild | +4.0 | -8.0 (-1.2) | -1.0 (-0.1) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sub_sqrt | +5.0 | +2.0 (+0.3) | -0.0 (-0.0) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sn_block | -5.0 | -4.0 (-0.9) | +3.0 (+0.8) | no effect |
| nonlinear dynamics | NL-TAR | hsic_wild | +4.0 | +7.0 (+1.1) | -4.0 (-0.6) | no effect |
| nonlinear dynamics | NL-TAR | hsic_shift | +0.0 | +3.0 (+0.6) | +1.0 (+0.2) | no effect |
| nonlinear dynamics | NL-TAR | hsic_sub | +6.0 | +10.0 (+1.5) | +4.0 (+0.5) | no effect |
| nonlinear dynamics | NL-TAR | dcov_wild | +1.0 | +9.0 (+1.6) | -6.0 (-1.0) | no effect |
| nonlinear dynamics | NL-TAR | dcov_sub | +3.0 | +16.0 (+2.3) | +13.0 (+1.7) | only with memory |
| seasonality | SEAS-M | mmd_unc_wild | -5.0 | -16.0 (-3.1) | +1.0 (+0.2) | dampens memory |
| seasonality | SEAS-M | mmd_sub_sqrt | -4.0 | -18.0 (-3.1) | -18.0 (-2.7) | dampens memory |
| seasonality | SEAS-M | mmd_sn_block | -5.0 | -7.0 (-1.8) | +3.0 (+0.8) | no effect |
| seasonality | SEAS-M | hsic_wild | +74.0 | -7.0 (-1.0) | -14.0 (-1.9) | breaks alone |
| seasonality | SEAS-M | hsic_shift | +2.0 | -1.0 (-0.2) | -6.0 (-1.3) | no effect |
| seasonality | SEAS-M | hsic_sub | +76.0 | -11.0 (-1.6) | -17.0 (-2.3) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_wild | +74.0 | -7.0 (-1.0) | -18.0 (-2.4) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_sub | +74.0 | -8.0 (-1.1) | -17.0 (-2.3) | breaks alone; dampens memory |
| seasonality | SEAS-V | mmd_unc_wild | -4.0 | -10.0 (-1.7) | +8.0 (+1.3) | no effect |
| seasonality | SEAS-V | mmd_sub_sqrt | -4.0 | -10.0 (-1.6) | -1.0 (-0.1) | no effect |
| seasonality | SEAS-V | mmd_sn_block | -3.0 | -5.0 (-1.1) | +5.0 (+1.2) | no effect |
| seasonality | SEAS-V | hsic_wild | +35.0 | +4.0 (+0.5) | -15.0 (-1.9) | breaks alone |
| seasonality | SEAS-V | hsic_shift | -4.0 | +10.0 (+2.1) | +7.0 (+1.5) | only with memory |
| seasonality | SEAS-V | hsic_sub | +29.0 | +3.0 (+0.4) | -2.0 (-0.2) | breaks alone |
| seasonality | SEAS-V | dcov_wild | +9.0 | +7.0 (+1.1) | -5.0 (-0.7) | no effect |
| seasonality | SEAS-V | dcov_sub | -2.0 | +16.0 (+2.5) | +11.0 (+1.5) | only with memory |

## n = 2000

| feature | generator | rule | alone pp | I(0.7) | I(0.8) | verdict |
|---|---|---|---|---|---|---|
| heavy tails | HT-S | mmd_unc_wild | +3.0 | -3.0 (-0.5) | +3.0 (+0.4) | no effect |
| heavy tails | HT-S | mmd_sub_sqrt | -3.0 | +3.0 (+0.4) | +11.0 (+1.4) | no effect |
| heavy tails | HT-S | mmd_sn_block | +2.0 | -0.0 (-0.0) | +8.0 (+1.4) | no effect |
| heavy tails | HT-S | hsic_wild | -2.0 | +2.0 (+0.4) | -2.0 (-0.4) | no effect |
| heavy tails | HT-S | hsic_shift | +3.0 | -8.0 (-1.9) | -3.0 (-0.8) | no effect |
| heavy tails | HT-S | hsic_sub | +2.0 | -3.0 (-0.6) | -7.0 (-1.1) | no effect |
| heavy tails | HT-S | dcov_wild | -3.0 | +5.0 (+1.0) | -5.0 (-1.0) | no effect |
| heavy tails | HT-S | dcov_sub | -1.0 | +1.0 (+0.2) | -5.0 (-0.8) | no effect |
| heavy tails | HT-L | mmd_unc_wild | +2.0 | -1.0 (-0.2) | +14.0 (+2.0) | only with memory |
| heavy tails | HT-L | mmd_sub_sqrt | -4.0 | +16.0 (+2.1) | +11.0 (+1.4) | only with memory |
| heavy tails | HT-L | mmd_sn_block | -1.0 | +1.0 (+0.2) | +9.0 (+1.8) | no effect |
| heavy tails | HT-L | hsic_wild | -3.0 | +6.0 (+1.1) | +7.0 (+1.2) | no effect |
| heavy tails | HT-L | hsic_shift | -1.0 | +2.0 (+0.4) | +8.0 (+1.9) | no effect |
| heavy tails | HT-L | hsic_sub | -1.0 | +4.0 (+0.8) | +20.0 (+2.9) | only with memory |
| heavy tails | HT-L | dcov_wild | -5.0 | +12.0 (+2.3) | +7.0 (+1.2) | only with memory |
| heavy tails | HT-L | dcov_sub | -2.0 | +7.0 (+1.4) | +12.0 (+1.8) | no effect |
| nonlinear dynamics | NL-SV | mmd_unc_wild | +8.0 | -5.0 (-0.8) | -11.0 (-1.6) | no effect |
| nonlinear dynamics | NL-SV | mmd_sub_sqrt | +11.0 | +1.0 (+0.1) | -17.0 (-2.0) | dampens memory |
| nonlinear dynamics | NL-SV | mmd_sn_block | +2.0 | +2.0 (+0.4) | +4.0 (+0.8) | no effect |
| nonlinear dynamics | NL-SV | hsic_wild | -1.0 | +0.0 (+0.0) | +2.0 (+0.3) | no effect |
| nonlinear dynamics | NL-SV | hsic_shift | -1.0 | -2.0 (-0.5) | +8.0 (+1.9) | no effect |
| nonlinear dynamics | NL-SV | hsic_sub | -5.0 | +2.0 (+0.5) | -1.0 (-0.2) | no effect |
| nonlinear dynamics | NL-SV | dcov_wild | -4.0 | +2.0 (+0.4) | -1.0 (-0.2) | no effect |
| nonlinear dynamics | NL-SV | dcov_sub | -3.0 | -4.0 (-1.1) | -17.0 (-3.3) | dampens memory |
| nonlinear dynamics | NL-TAR | mmd_unc_wild | +8.0 | -5.0 (-0.8) | -11.0 (-1.6) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sub_sqrt | +3.0 | -1.0 (-0.1) | -16.0 (-2.0) | dampens memory |
| nonlinear dynamics | NL-TAR | mmd_sn_block | +1.0 | -3.0 (-0.6) | -7.0 (-1.8) | no effect |
| nonlinear dynamics | NL-TAR | hsic_wild | -4.0 | +6.0 (+1.2) | +9.0 (+1.6) | no effect |
| nonlinear dynamics | NL-TAR | hsic_shift | -1.0 | -2.0 (-0.5) | +9.0 (+2.1) | only with memory |
| nonlinear dynamics | NL-TAR | hsic_sub | +0.0 | +7.0 (+1.2) | +10.0 (+1.5) | no effect |
| nonlinear dynamics | NL-TAR | dcov_wild | -1.0 | +1.0 (+0.2) | +2.0 (+0.3) | no effect |
| nonlinear dynamics | NL-TAR | dcov_sub | +0.0 | +8.0 (+1.5) | +10.0 (+1.4) | no effect |
| seasonality | SEAS-M | mmd_unc_wild | -3.0 | +5.0 (+0.9) | +1.0 (+0.2) | no effect |
| seasonality | SEAS-M | mmd_sub_sqrt | -8.0 (below G0) | -20.0 (-3.5) | -24.0 (-3.5) | dampens memory |
| seasonality | SEAS-M | mmd_sn_block | -5.0 | -2.0 (-0.6) | -1.0 (-0.3) | no effect |
| seasonality | SEAS-M | hsic_wild | +80.0 | +6.0 (+1.1) | -2.0 (-0.3) | breaks alone |
| seasonality | SEAS-M | hsic_shift | -4.0 | +2.0 (+0.5) | +5.0 (+1.7) | no effect |
| seasonality | SEAS-M | hsic_sub | +84.0 | +1.0 (+0.2) | -14.0 (-2.2) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_wild | +84.0 | +3.0 (+0.6) | -8.0 (-1.3) | breaks alone |
| seasonality | SEAS-M | dcov_sub | +86.0 | -1.0 (-0.2) | -20.0 (-3.1) | breaks alone; dampens memory |
| seasonality | SEAS-V | mmd_unc_wild | +0.0 | +0.0 (+0.0) | -1.0 (-0.1) | no effect |
| seasonality | SEAS-V | mmd_sub_sqrt | -8.0 (below G0) | +3.0 (+0.4) | -4.0 (-0.5) | no effect |
| seasonality | SEAS-V | mmd_sn_block | -5.0 | -0.0 (-0.0) | +2.0 (+0.6) | no effect |
| seasonality | SEAS-V | hsic_wild | +68.0 | -2.0 (-0.3) | -11.0 (-1.5) | breaks alone |
| seasonality | SEAS-V | hsic_shift | -1.0 | +0.0 (+0.0) | +6.0 (+1.5) | no effect |
| seasonality | SEAS-V | hsic_sub | +69.0 | -5.0 (-0.7) | -16.0 (-2.1) | breaks alone; dampens memory |
| seasonality | SEAS-V | dcov_wild | +44.0 | -5.0 (-0.6) | -16.0 (-2.0) | breaks alone |
| seasonality | SEAS-V | dcov_sub | +40.0 | -1.0 (-0.1) | -10.0 (-1.2) | breaks alone |

Cells: 90. Errors: 0. CPU-hours (row sum): 2.80.
