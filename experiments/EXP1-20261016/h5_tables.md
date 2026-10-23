# H5 — wild-bootstrap variants at H = 0.5 (MMD homogeneity, l = floor(sqrt n), B = 500)

Size at nominal 5% (%, Wilson 95% CI). Exploratory simulation.

| cell | reps ok / errors | old (mmd-dml) | csg (published) | iid | unc | ustat |
|---|---|---|---|---|---|---|
| h5__G1_H0.5__n500 | 2000 / 0 | 9.0 [7.8, 10.3] | 8.9 [7.7, 10.2] | 5.1 [4.2, 6.2] | 4.4 [3.6, 5.4] | 6.1 [5.1, 7.2] |
| h5__G1_H0.5__n2000 | 2000 / 0 | 7.0 [6.0, 8.2] | 7.2 [6.1, 8.4] | 5.3 [4.4, 6.4] | 5.0 [4.1, 6.0] | 5.9 [4.9, 7.0] |
| h5__G3_H0.5__n500 | 2000 / 0 | 9.0 [7.8, 10.3] | 9.1 [7.9, 10.4] | 5.3 [4.4, 6.4] | 4.6 [3.8, 5.6] | 5.3 [4.4, 6.4] |
| h5__G3_H0.5__n2000 | 2000 / 0 | 6.6 [5.5, 7.7] | 6.1 [5.1, 7.2] | 5.3 [4.4, 6.4] | 4.4 [3.6, 5.4] | 4.8 [3.9, 5.8] |

Rule (fixed before the run): a variant is on size if it is within 5 ± 1 pp in every H5 cell.

- old: not on size
- csg: not on size
- iid: on size
- unc: on size
- ustat: not on size
