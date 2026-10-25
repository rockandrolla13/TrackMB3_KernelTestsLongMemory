# EXP2 tables (main)

Rejection % at nominal 5% [Wilson 95%]. Rows = generator, columns = H (H_X = H_Y for independence).

## mmd_unc_wild, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 4.7 [3.6, 6.2] | 14.8 [12.7, 17.1] | 20.3 [17.9, 22.9] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.5 [3.4, 6.0] | 13.6 [11.6, 15.9] | 21.7 [19.3, 24.4] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 5.2 [4.0, 6.8] | 17.5 [15.3, 20.0] | 21.8 [19.4, 24.5] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 9.6 [7.9, 11.6] | 17.9 [15.6, 20.4] | 21.5 [19.1, 24.2] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 5.7 [4.4, 7.3] | 16.6 [14.4, 19.0] | 21.5 [19.1, 24.2] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 0.4] | 3.2 [2.3, 4.5] | 8.8 [7.2, 10.7] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 1.8 [1.1, 2.8] | 10.1 [8.4, 12.1] | 19.5 [17.2, 22.1] | 1000/1000/1000 |

## mmd_unc_wild, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 4.9 [3.7, 6.4] | 20.1 [17.7, 22.7] | 26.5 [23.9, 29.3] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.7 [3.6, 6.2] | 18.2 [15.9, 20.7] | 29.6 [26.9, 32.5] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [3.0, 5.4] | 21.9 [19.4, 24.6] | 30.6 [27.8, 33.5] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 6.6 [5.2, 8.3] | 18.6 [16.3, 21.1] | 26.3 [23.7, 29.1] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 6.1 [4.8, 7.8] | 20.6 [18.2, 23.2] | 28.6 [25.9, 31.5] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.4 [0.2, 1.0] | 10.8 [9.0, 12.9] | 22.1 [19.6, 24.8] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 2.2 [1.5, 3.3] | 16.4 [14.2, 18.8] | 27.6 [24.9, 30.5] | 1000/1000/1000 |

## mmd_sub_sqrt, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.0 [4.7, 7.6] | 26.2 [23.6, 29.0] | 42.2 [39.2, 45.3] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.4 [3.3, 5.9] | 26.3 [23.7, 29.1] | 39.4 [36.4, 42.5] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 5.4 [4.2, 7.0] | 33.2 [30.4, 36.2] | 47.3 [44.2, 50.4] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 21.6 [19.2, 24.3] | 36.5 [33.6, 39.5] | 41.1 [38.1, 44.2] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 8.4 [6.8, 10.3] | 28.9 [26.2, 31.8] | 43.8 [40.8, 46.9] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.1 [0.0, 0.6] | 5.9 [4.6, 7.5] | 16.5 [14.3, 18.9] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 1.8 [1.1, 2.8] | 16.4 [14.2, 18.8] | 35.0 [32.1, 38.0] | 1000/1000/1000 |

## mmd_sub_sqrt, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.1 [3.9, 6.6] | 33.7 [30.8, 36.7] | 49.1 [46.0, 52.2] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.9 [3.7, 6.4] | 32.7 [29.9, 35.7] | 52.0 [48.9, 55.1] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 4.0 [3.0, 5.4] | 39.2 [36.2, 42.3] | 57.3 [54.2, 60.3] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 13.3 [11.3, 15.5] | 38.3 [35.3, 41.4] | 52.6 [49.5, 55.7] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 6.9 [5.5, 8.6] | 35.6 [32.7, 38.6] | 49.9 [46.8, 53.0] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 0.4] | 3.9 [2.9, 5.3] | 19.5 [17.2, 22.1] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 0.2 [0.1, 0.7] | 20.9 [18.5, 23.5] | 43.2 [40.2, 46.3] | 1000/1000/1000 |

## mmd_sn_block, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 4.2 [3.1, 5.6] | 8.6 [7.0, 10.5] | 7.6 [6.1, 9.4] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 3.7 [2.7, 5.1] | 8.2 [6.7, 10.1] | 9.4 [7.7, 11.4] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 4.7 [3.6, 6.2] | 9.7 [8.0, 11.7] | 9.8 [8.1, 11.8] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.8 [4.5, 7.4] | 11.2 [9.4, 13.3] | 13.8 [11.8, 16.1] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 1.4 [0.8, 2.3] | 1.9 [1.2, 2.9] | 2.1 [1.4, 3.2] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 0.4] | 0.5 [0.2, 1.2] | 0.7 [0.3, 1.4] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 1.1 [0.6, 2.0] | 3.6 [2.6, 4.9] | 5.7 [4.4, 7.3] | 1000/1000/1000 |

## mmd_sn_block, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 3.7 [2.7, 5.1] | 8.6 [7.0, 10.5] | 7.9 [6.4, 9.7] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.8 [3.6, 6.3] | 9.0 [7.4, 10.9] | 10.7 [8.9, 12.8] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 3.6 [2.6, 4.9] | 10.0 [8.3, 12.0] | 9.6 [7.9, 11.6] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 3.2 [2.3, 4.5] | 9.4 [7.7, 11.4] | 13.0 [11.1, 15.2] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 3.1 [2.2, 4.4] | 5.0 [3.8, 6.5] | 3.1 [2.2, 4.4] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 0.0 [0.0, 0.4] | 0.1 [0.0, 0.6] | 1.8 [1.1, 2.8] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 0.4 [0.2, 1.0] | 4.6 [3.5, 6.1] | 7.6 [6.1, 9.4] | 1000/1000/1000 |

## hsic_wild, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 8.2 [6.7, 10.1] | 7.7 [6.2, 9.5] | 12.3 [10.4, 14.5] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 7.7 [6.2, 9.5] | 8.9 [7.3, 10.8] | 11.4 [9.6, 13.5] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 8.3 [6.7, 10.2] | 9.7 [8.0, 11.7] | 10.6 [8.8, 12.7] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 8.2 [6.7, 10.1] | 11.2 [9.4, 13.3] | 9.1 [7.5, 11.0] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 7.7 [6.2, 9.5] | 9.9 [8.2, 11.9] | 14.4 [12.4, 16.7] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 81.5 [79.0, 83.8] | 77.1 [74.4, 79.6] | 74.4 [71.6, 77.0] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 38.7 [35.7, 41.8] | 38.3 [35.3, 41.4] | 34.6 [31.7, 37.6] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 99.7 [99.1, 99.9] | 1000/1000/1000 |

## hsic_wild, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.8 [5.4, 8.5] | 7.4 [5.9, 9.2] | 8.1 [6.6, 10.0] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 6.1 [4.8, 7.8] | 7.4 [5.9, 9.2] | 7.0 [5.6, 8.8] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 8.9 [7.3, 10.8] | 9.5 [7.8, 11.5] | 10.6 [8.8, 12.7] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 7.9 [6.4, 9.7] | 7.8 [6.3, 9.6] | 6.7 [5.3, 8.4] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 7.1 [5.7, 8.9] | 8.2 [6.7, 10.1] | 11.9 [10.0, 14.1] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.2 [89.3, 92.8] | 90.3 [88.3, 92.0] | 85.2 [82.9, 87.3] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 73.3 [70.5, 75.9] | 72.4 [69.5, 75.1] | 68.0 [65.0, 70.8] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## hsic_shift, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.1 [3.9, 6.6] | 5.3 [4.1, 6.9] | 7.7 [6.2, 9.5] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.5 [3.4, 6.0] | 5.1 [3.9, 6.6] | 5.9 [4.6, 7.5] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 5.0 [3.8, 6.5] | 4.7 [3.6, 6.2] | 6.5 [5.1, 8.2] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.7 [3.6, 6.2] | 5.9 [4.6, 7.5] | 6.5 [5.1, 8.2] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 5.0 [3.8, 6.5] | 5.6 [4.3, 7.2] | 7.0 [5.6, 8.8] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 5.2 [4.0, 6.8] | 5.0 [3.8, 6.5] | 4.4 [3.3, 5.9] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 5.5 [4.2, 7.1] | 5.9 [4.6, 7.5] | 5.4 [4.2, 7.0] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 23.0 [20.5, 25.7] | 20.3 [17.9, 22.9] | 20.3 [17.9, 22.9] | 1000/1000/1000 |

## hsic_shift, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.3 [4.1, 6.9] | 4.5 [3.4, 6.0] | 5.3 [4.1, 6.9] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.5 [3.4, 6.0] | 6.3 [5.0, 8.0] | 3.8 [2.8, 5.2] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 6.4 [5.0, 8.1] | 6.3 [5.0, 8.0] | 5.3 [4.1, 6.9] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.8 [4.5, 7.4] | 5.9 [4.6, 7.5] | 3.9 [2.9, 5.3] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 4.7 [3.6, 6.2] | 5.8 [4.5, 7.4] | 7.0 [5.6, 8.8] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 5.8 [4.5, 7.4] | 4.6 [3.5, 6.1] | 5.9 [4.6, 7.5] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 5.0 [3.8, 6.5] | 5.3 [4.1, 6.9] | 5.1 [3.9, 6.6] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 34.3 [31.4, 37.3] | 29.4 [26.7, 32.3] | 27.9 [25.2, 30.8] | 1000/1000/1000 |

## hsic_sub, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.1 [4.8, 7.8] | 10.5 [8.7, 12.6] | 21.0 [18.6, 23.6] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 5.2 [4.0, 6.8] | 8.1 [6.6, 10.0] | 16.5 [14.3, 18.9] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 5.9 [4.6, 7.5] | 9.8 [8.1, 11.8] | 24.1 [21.6, 26.8] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.3 [4.1, 6.9] | 6.8 [5.4, 8.5] | 8.5 [6.9, 10.4] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 9.2 [7.6, 11.2] | 21.0 [18.6, 23.6] | 37.8 [34.8, 40.8] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 82.0 [79.5, 84.3] | 79.6 [77.0, 82.0] | 79.0 [76.4, 81.4] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 29.2 [26.5, 32.1] | 37.0 [34.1, 40.0] | 45.8 [42.7, 48.9] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## hsic_sub, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.8 [4.5, 7.4] | 8.2 [6.7, 10.1] | 17.5 [15.3, 20.0] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.3 [3.2, 5.7] | 8.0 [6.5, 9.8] | 15.7 [13.6, 18.1] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 6.8 [5.4, 8.5] | 11.4 [9.6, 13.5] | 25.9 [23.3, 28.7] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 4.9 [3.7, 6.4] | 6.2 [4.9, 7.9] | 6.2 [4.9, 7.9] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 7.2 [5.8, 9.0] | 14.6 [12.5, 16.9] | 29.8 [27.0, 32.7] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 90.9 [89.0, 92.5] | 89.7 [87.7, 91.4] | 85.4 [83.1, 87.5] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 71.5 [68.6, 74.2] | 71.4 [68.5, 74.1] | 73.0 [70.2, 75.7] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## dcov_wild, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 8.1 [6.6, 10.0] | 8.0 [6.5, 9.8] | 12.1 [10.2, 14.3] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 7.7 [6.2, 9.5] | 9.5 [7.8, 11.5] | 12.3 [10.4, 14.5] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 7.1 [5.7, 8.9] | 9.5 [7.8, 11.5] | 12.1 [10.2, 14.3] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 9.0 [7.4, 10.9] | 9.7 [8.0, 11.7] | 8.8 [7.2, 10.7] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 9.0 [7.4, 10.9] | 10.8 [9.0, 12.9] | 15.6 [13.5, 18.0] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 83.8 [81.4, 86.0] | 79.2 [76.6, 81.6] | 75.4 [72.6, 78.0] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 19.3 [17.0, 21.9] | 19.7 [17.4, 22.3] | 17.3 [15.1, 19.8] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 99.9 [99.4, 100.0] | 1000/1000/1000 |

## dcov_wild, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.7 [5.3, 8.4] | 8.9 [7.3, 10.8] | 9.1 [7.5, 11.0] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 6.4 [5.0, 8.1] | 7.6 [6.1, 9.4] | 9.7 [8.0, 11.7] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 7.0 [5.6, 8.8] | 7.6 [6.1, 9.4] | 10.9 [9.1, 13.0] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 8.4 [6.8, 10.3] | 7.6 [6.1, 9.4] | 7.5 [6.0, 9.3] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 5.9 [4.6, 7.5] | 7.7 [6.2, 9.5] | 11.8 [9.9, 13.9] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.7 [89.8, 93.3] | 90.6 [88.6, 92.3] | 86.1 [83.8, 88.1] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 51.7 [48.6, 54.8] | 44.7 [41.6, 47.8] | 34.2 [31.3, 37.2] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## dcov_sub, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.4 [4.2, 7.0] | 11.1 [9.3, 13.2] | 25.5 [22.9, 28.3] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 3.1 [2.2, 4.4] | 6.7 [5.3, 8.4] | 16.4 [14.2, 18.8] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 2.7 [1.9, 3.9] | 8.0 [6.5, 9.8] | 22.8 [20.3, 25.5] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 2.6 [1.8, 3.8] | 3.2 [2.3, 4.5] | 4.0 [3.0, 5.4] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 9.9 [8.2, 11.9] | 23.3 [20.8, 26.0] | 41.6 [38.6, 44.7] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 84.0 [81.6, 86.1] | 82.0 [79.5, 84.3] | 81.9 [79.4, 84.2] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 9.1 [7.5, 11.0] | 17.4 [15.2, 19.9] | 27.6 [24.9, 30.5] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## dcov_sub, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 6.8 [5.4, 8.5] | 10.9 [9.1, 13.0] | 22.8 [20.3, 25.5] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 3.2 [2.3, 4.5] | 6.4 [5.0, 8.1] | 16.5 [14.3, 18.9] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 3.2 [2.3, 4.5] | 8.4 [6.8, 10.3] | 25.9 [23.3, 28.7] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 2.2 [1.5, 3.3] | 1.9 [1.2, 2.9] | 1.9 [1.2, 2.9] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 6.4 [5.0, 8.1] | 16.0 [13.9, 18.4] | 33.5 [30.6, 36.5] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.3 [89.4, 92.9] | 89.7 [87.7, 91.4] | 86.6 [84.3, 88.6] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 45.6 [42.5, 48.7] | 46.4 [43.3, 49.5] | 53.3 [50.2, 56.4] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## pearson_nw, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.2 [4.0, 6.8] | 7.3 [5.8, 9.1] | 14.7 [12.6, 17.0] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 5.0 [3.8, 6.5] | 6.7 [5.3, 8.4] | 13.5 [11.5, 15.8] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 3.9 [2.9, 5.3] | 7.2 [5.8, 9.0] | 12.0 [10.1, 14.2] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.2 [4.0, 6.8] | 4.9 [3.7, 6.4] | 4.7 [3.6, 6.2] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 7.2 [5.8, 9.0] | 9.1 [7.5, 11.0] | 12.8 [10.9, 15.0] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 83.5 [81.1, 85.7] | 78.6 [76.0, 81.0] | 74.2 [71.4, 76.8] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 5.8 [4.5, 7.4] | 7.1 [5.7, 8.9] | 11.6 [9.8, 13.7] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 99.9 [99.4, 100.0] | 1000/1000/1000 |

## pearson_nw, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.0 [3.8, 6.5] | 9.1 [7.5, 11.0] | 14.1 [12.1, 16.4] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 5.7 [4.4, 7.3] | 7.3 [5.8, 9.1] | 11.5 [9.7, 13.6] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 4.2 [3.1, 5.6] | 7.1 [5.7, 8.9] | 16.2 [14.0, 18.6] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.2 [4.0, 6.8] | 4.1 [3.0, 5.5] | 4.6 [3.5, 6.1] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 6.6 [5.2, 8.3] | 8.2 [6.7, 10.1] | 12.1 [10.2, 14.3] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.4 [89.5, 93.0] | 89.9 [87.9, 91.6] | 84.5 [82.1, 86.6] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 4.7 [3.6, 6.2] | 8.4 [6.8, 10.3] | 13.4 [11.4, 15.7] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

## pearson_ewc, n = 500

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.5 [4.2, 7.1] | 5.1 [3.9, 6.6] | 8.4 [6.8, 10.3] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 5.1 [3.9, 6.6] | 4.8 [3.6, 6.3] | 7.6 [6.1, 9.4] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 3.8 [2.8, 5.2] | 5.3 [4.1, 6.9] | 6.6 [5.2, 8.3] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.3 [4.1, 6.9] | 4.2 [3.1, 5.6] | 5.0 [3.8, 6.5] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 5.3 [4.1, 6.9] | 5.9 [4.6, 7.5] | 7.6 [6.1, 9.4] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 84.1 [81.7, 86.2] | 75.7 [72.9, 78.3] | 70.7 [67.8, 73.4] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 5.1 [3.9, 6.6] | 5.4 [4.2, 7.0] | 6.9 [5.5, 8.6] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 99.7 [99.1, 99.9] | 1000/1000/1000 |

## pearson_ewc, n = 2000

| generator | subordinated | H = 0.5 | H = 0.7 | H = 0.8 | reps |
|---|---|---|---|---|---|
| G0 | yes (rank 1) | 5.3 [4.1, 6.9] | 6.9 [5.5, 8.6] | 8.5 [6.9, 10.4] | 1000/1000/1000 |
| HT-S | yes (instantaneous, rank 1) | 4.9 [3.7, 6.4] | 4.8 [3.6, 6.3] | 7.1 [5.7, 8.9] | 1000/1000/1000 |
| HT-L | no (linear filter, t3 innovations) | 3.8 [2.8, 5.2] | 5.1 [3.9, 6.6] | 8.6 [7.0, 10.5] | 1000/1000/1000 |
| NL-SV | multivariate only (function of two Gaussian processes) | 5.6 [4.3, 7.2] | 4.7 [3.6, 6.2] | 4.6 [3.5, 6.1] | 1000/1000/1000 |
| NL-TAR | no (nonlinear path functional) | 5.2 [4.0, 6.8] | 6.5 [5.1, 8.2] | 9.8 [8.1, 11.8] | 1000/1000/1000 |
| SEAS-M | no (deterministic periodic mean, random phase) | 91.7 [89.8, 93.3] | 89.8 [87.8, 91.5] | 85.1 [82.8, 87.2] | 1000/1000/1000 |
| SEAS-V | no (deterministic periodic scale, random phase) | 4.5 [3.4, 6.0] | 6.5 [5.1, 8.2] | 8.7 [7.1, 10.6] | 1000/1000/1000 |
| SEAS-M-shared | no (fixed common phase; diagnostic) | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 100.0 [99.6, 100.0] | 1000/1000/1000 |

# Verdicts (pre-declared rule, DESIGN.md)

alone = F(0.5) − G0(0.5) in pp; I(h) = (F(h) − F(0.5)) − (G0(h) − G0(0.5)) in pp, z in brackets.

## n = 500

| feature | generator | rule | alone pp | I(0.7) | I(0.8) | verdict |
|---|---|---|---|---|---|---|
| heavy tails | HT-S | mmd_unc_wild | -0.2 | -1.0 (-0.6) | +1.6 (+0.8) | no effect |
| heavy tails | HT-S | mmd_sub_sqrt | -1.6 | +1.7 (+0.8) | -1.2 (-0.5) | no effect |
| heavy tails | HT-S | mmd_sn_block | -0.5 | +0.1 (+0.1) | +2.3 (+1.5) | no effect |
| heavy tails | HT-S | hsic_wild | -0.5 | +1.7 (+1.0) | -0.4 (-0.2) | no effect |
| heavy tails | HT-S | hsic_shift | -0.6 | +0.4 (+0.3) | -1.2 (-0.8) | no effect |
| heavy tails | HT-S | hsic_sub | -0.9 | -1.5 (-0.9) | -3.6 (-1.8) | no effect |
| heavy tails | HT-S | dcov_wild | -0.4 | +1.9 (+1.1) | +0.6 (+0.3) | no effect |
| heavy tails | HT-S | dcov_sub | -2.3 | -2.1 (-1.4) | -6.8 (-3.4) | dampens memory |
| heavy tails | HT-L | mmd_unc_wild | +0.5 | +2.2 (+1.1) | +1.0 (+0.5) | no effect |
| heavy tails | HT-L | mmd_sub_sqrt | -0.6 | +7.6 (+3.3) | +5.7 (+2.3) | only with memory |
| heavy tails | HT-L | mmd_sn_block | +0.5 | +0.6 (+0.4) | +1.7 (+1.1) | no effect |
| heavy tails | HT-L | hsic_wild | +0.1 | +1.9 (+1.1) | -1.8 (-1.0) | no effect |
| heavy tails | HT-L | hsic_shift | -0.1 | -0.5 (-0.4) | -1.1 (-0.7) | no effect |
| heavy tails | HT-L | hsic_sub | -0.2 | -0.5 (-0.3) | +3.3 (+1.5) | no effect |
| heavy tails | HT-L | dcov_wild | -1.0 | +2.5 (+1.4) | +1.0 (+0.5) | no effect |
| heavy tails | HT-L | dcov_sub | -2.7 (below G0) | -0.4 (-0.2) | +0.0 (+0.0) | no effect |
| nonlinear dynamics | NL-SV | mmd_unc_wild | +4.9 | -1.8 (-0.9) | -3.7 (-1.7) | breaks alone |
| nonlinear dynamics | NL-SV | mmd_sub_sqrt | +15.6 | -5.3 (-2.1) | -16.7 (-6.3) | breaks alone; dampens memory |
| nonlinear dynamics | NL-SV | mmd_sn_block | +1.6 | +1.0 (+0.6) | +4.6 (+2.7) | only with memory |
| nonlinear dynamics | NL-SV | hsic_wild | +0.0 | +3.5 (+1.9) | -3.2 (-1.7) | no effect |
| nonlinear dynamics | NL-SV | hsic_shift | -0.4 | +1.0 (+0.7) | -0.8 (-0.5) | no effect |
| nonlinear dynamics | NL-SV | hsic_sub | -0.8 | -2.9 (-1.8) | -11.7 (-6.2) | dampens memory |
| nonlinear dynamics | NL-SV | dcov_wild | +0.9 | +0.8 (+0.5) | -4.2 (-2.3) | dampens memory |
| nonlinear dynamics | NL-SV | dcov_sub | -2.8 (below G0) | -5.1 (-3.5) | -18.7 (-10.7) | dampens memory |
| nonlinear dynamics | NL-TAR | mmd_unc_wild | +1.0 | +0.8 (+0.4) | +0.2 (+0.1) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sub_sqrt | +2.4 | +0.3 (+0.1) | -0.8 (-0.3) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sn_block | -2.8 (below G0) | -3.9 (-3.2) | -2.7 (-2.2) | dampens memory |
| nonlinear dynamics | NL-TAR | hsic_wild | -0.5 | +2.7 (+1.5) | +2.6 (+1.3) | no effect |
| nonlinear dynamics | NL-TAR | hsic_shift | -0.1 | +0.4 (+0.3) | -0.6 (-0.4) | no effect |
| nonlinear dynamics | NL-TAR | hsic_sub | +3.1 | +7.4 (+3.7) | +13.7 (+5.9) | only with memory |
| nonlinear dynamics | NL-TAR | dcov_wild | +0.9 | +1.9 (+1.1) | +2.6 (+1.3) | no effect |
| nonlinear dynamics | NL-TAR | dcov_sub | +4.5 | +7.7 (+3.8) | +11.6 (+4.8) | breaks alone and amplifies memory |
| seasonality | SEAS-M | mmd_unc_wild | -4.7 (below G0) | -6.9 (-4.9) | -6.8 (-4.0) | dampens memory |
| seasonality | SEAS-M | mmd_sub_sqrt | -5.9 (below G0) | -14.4 (-8.2) | -19.8 (-9.4) | dampens memory |
| seasonality | SEAS-M | mmd_sn_block | -4.2 (below G0) | -3.9 (-3.5) | -2.7 (-2.5) | dampens memory |
| seasonality | SEAS-M | hsic_wild | +73.3 | -3.9 (-1.8) | -11.2 (-4.9) | breaks alone; dampens memory |
| seasonality | SEAS-M | hsic_shift | +0.1 | -0.4 (-0.3) | -3.4 (-2.3) | dampens memory |
| seasonality | SEAS-M | hsic_sub | +75.9 | -6.8 (-3.2) | -17.9 (-7.7) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_wild | +75.7 | -4.5 (-2.1) | -12.4 (-5.5) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_sub | +78.6 | -7.7 (-3.7) | -22.2 (-9.7) | breaks alone; dampens memory |
| seasonality | SEAS-V | mmd_unc_wild | -2.9 (below G0) | -1.8 (-1.1) | +2.1 (+1.1) | no effect |
| seasonality | SEAS-V | mmd_sub_sqrt | -4.2 (below G0) | -5.6 (-2.8) | -3.0 (-1.3) | dampens memory |
| seasonality | SEAS-V | mmd_sn_block | -3.1 (below G0) | -1.9 (-1.5) | +1.2 (+0.9) | no effect |
| seasonality | SEAS-V | hsic_wild | +30.5 | +0.1 (+0.0) | -8.2 (-3.2) | breaks alone; dampens memory |
| seasonality | SEAS-V | hsic_shift | +0.4 | +0.2 (+0.1) | -2.7 (-1.8) | no effect |
| seasonality | SEAS-V | hsic_sub | +23.1 | +3.4 (+1.4) | +1.7 (+0.7) | breaks alone |
| seasonality | SEAS-V | dcov_wild | +11.2 | +0.5 (+0.2) | -6.0 (-2.7) | breaks alone; dampens memory |
| seasonality | SEAS-V | dcov_sub | +3.7 | +2.6 (+1.3) | -1.6 (-0.7) | breaks alone |

## n = 2000

| feature | generator | rule | alone pp | I(0.7) | I(0.8) | verdict |
|---|---|---|---|---|---|---|
| heavy tails | HT-S | mmd_unc_wild | -0.2 | -1.7 (-0.8) | +3.3 (+1.5) | no effect |
| heavy tails | HT-S | mmd_sub_sqrt | -0.2 | -0.8 (-0.3) | +3.1 (+1.3) | no effect |
| heavy tails | HT-S | mmd_sn_block | +1.1 | -0.7 (-0.5) | +1.7 (+1.1) | no effect |
| heavy tails | HT-S | hsic_wild | -0.7 | +0.7 (+0.4) | -0.4 (-0.2) | no effect |
| heavy tails | HT-S | hsic_shift | -0.8 | +2.6 (+1.9) | -0.7 (-0.5) | no effect |
| heavy tails | HT-S | hsic_sub | -1.5 | +1.3 (+0.8) | -0.3 (-0.2) | no effect |
| heavy tails | HT-S | dcov_wild | -0.3 | -1.0 (-0.6) | +0.9 (+0.5) | no effect |
| heavy tails | HT-S | dcov_sub | -3.6 (below G0) | -0.9 (-0.6) | -2.7 (-1.3) | no effect |
| heavy tails | HT-L | mmd_unc_wild | -0.9 | +2.7 (+1.3) | +5.0 (+2.2) | only with memory |
| heavy tails | HT-L | mmd_sub_sqrt | -1.1 | +6.6 (+2.8) | +9.3 (+3.9) | only with memory |
| heavy tails | HT-L | mmd_sn_block | -0.1 | +1.5 (+1.0) | +1.8 (+1.2) | no effect |
| heavy tails | HT-L | hsic_wild | +2.1 | +0.0 (+0.0) | +0.4 (+0.2) | no effect |
| heavy tails | HT-L | hsic_shift | +1.1 | +0.7 (+0.5) | -1.1 (-0.8) | no effect |
| heavy tails | HT-L | hsic_sub | +1.0 | +2.2 (+1.3) | +7.4 (+3.5) | only with memory |
| heavy tails | HT-L | dcov_wild | +0.3 | -1.6 (-1.0) | +1.5 (+0.9) | no effect |
| heavy tails | HT-L | dcov_sub | -3.6 (below G0) | +1.1 (+0.7) | +6.7 (+3.1) | only with memory |
| nonlinear dynamics | NL-SV | mmd_unc_wild | +1.7 | -3.2 (-1.6) | -1.9 (-0.8) | no effect |
| nonlinear dynamics | NL-SV | mmd_sub_sqrt | +8.2 | -3.6 (-1.4) | -4.7 (-1.8) | breaks alone |
| nonlinear dynamics | NL-SV | mmd_sn_block | -0.5 | +1.3 (+0.9) | +5.6 (+3.5) | only with memory |
| nonlinear dynamics | NL-SV | hsic_wild | +1.1 | -0.7 (-0.4) | -2.5 (-1.5) | no effect |
| nonlinear dynamics | NL-SV | hsic_shift | +0.5 | +0.9 (+0.6) | -1.9 (-1.4) | no effect |
| nonlinear dynamics | NL-SV | hsic_sub | -0.9 | -1.1 (-0.7) | -10.4 (-6.0) | dampens memory |
| nonlinear dynamics | NL-SV | dcov_wild | +1.7 | -3.0 (-1.8) | -3.3 (-1.9) | no effect |
| nonlinear dynamics | NL-SV | dcov_sub | -4.6 (below G0) | -4.4 (-3.1) | -16.3 (-9.8) | dampens memory |
| nonlinear dynamics | NL-TAR | mmd_unc_wild | +1.2 | -0.7 (-0.3) | +0.9 (+0.4) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sub_sqrt | +1.8 | +0.1 (+0.0) | -1.0 (-0.4) | no effect |
| nonlinear dynamics | NL-TAR | mmd_sn_block | -0.6 | -3.0 (-2.2) | -4.2 (-3.2) | dampens memory |
| nonlinear dynamics | NL-TAR | hsic_wild | +0.3 | +0.5 (+0.3) | +3.5 (+2.0) | only with memory |
| nonlinear dynamics | NL-TAR | hsic_shift | -0.6 | +1.9 (+1.4) | +2.3 (+1.6) | no effect |
| nonlinear dynamics | NL-TAR | hsic_sub | +1.4 | +5.0 (+2.8) | +10.9 (+5.0) | only with memory |
| nonlinear dynamics | NL-TAR | dcov_wild | -0.8 | -0.4 (-0.2) | +3.5 (+2.0) | only with memory |
| nonlinear dynamics | NL-TAR | dcov_sub | -0.4 | +5.5 (+2.9) | +11.1 (+4.9) | only with memory |
| seasonality | SEAS-M | mmd_unc_wild | -4.5 (below G0) | -4.8 (-2.7) | +0.1 (+0.1) | dampens memory |
| seasonality | SEAS-M | mmd_sub_sqrt | -5.1 (below G0) | -24.7 (-14.0) | -24.5 (-11.5) | dampens memory |
| seasonality | SEAS-M | mmd_sn_block | -3.7 (below G0) | -4.8 (-4.5) | -2.4 (-2.1) | dampens memory |
| seasonality | SEAS-M | hsic_wild | +84.4 | -1.5 (-0.9) | -7.3 (-3.9) | breaks alone; dampens memory |
| seasonality | SEAS-M | hsic_shift | +0.5 | -0.4 (-0.3) | +0.1 (+0.1) | no effect |
| seasonality | SEAS-M | hsic_sub | +85.1 | -3.6 (-2.1) | -17.2 (-8.5) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_wild | +85.0 | -3.3 (-1.9) | -8.0 (-4.3) | breaks alone; dampens memory |
| seasonality | SEAS-M | dcov_sub | +84.5 | -5.7 (-3.1) | -20.7 (-9.9) | breaks alone; dampens memory |
| seasonality | SEAS-V | mmd_unc_wild | -2.7 (below G0) | -1.0 (-0.5) | +3.8 (+1.8) | no effect |
| seasonality | SEAS-V | mmd_sub_sqrt | -4.9 (below G0) | -7.9 (-3.8) | -1.0 (-0.4) | dampens memory |
| seasonality | SEAS-V | mmd_sn_block | -3.3 (below G0) | -0.7 (-0.6) | +3.0 (+2.2) | only with memory |
| seasonality | SEAS-V | hsic_wild | +66.5 | -1.5 (-0.7) | -6.6 (-2.8) | breaks alone; dampens memory |
| seasonality | SEAS-V | hsic_shift | -0.3 | +1.1 (+0.8) | +0.1 (+0.1) | no effect |
| seasonality | SEAS-V | hsic_sub | +65.7 | -2.5 (-1.1) | -10.2 (-4.2) | breaks alone; dampens memory |
| seasonality | SEAS-V | dcov_wild | +45.0 | -9.2 (-3.6) | -19.9 (-8.0) | breaks alone; dampens memory |
| seasonality | SEAS-V | dcov_sub | +38.8 | -3.3 (-1.3) | -8.3 (-3.1) | breaks alone; dampens memory |

Cells: 90. Errors: 0. CPU-hours (row sum): 26.40.
