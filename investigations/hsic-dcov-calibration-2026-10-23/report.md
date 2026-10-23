# HSIC and distance-covariance multiplier control

23 October 2026. This is a paired finite-sample check under independent `H=0.5` Gaussian fractional noise. It does not establish validity under long memory and does not rerun the trade-window placebos.

## Result

The uncentred Markov multiplier reduces the centred implementation's oversize at `n=500`, but is conservative in all six tested cells. It **fails the advancement rule** in `DESIGN.md`: at `n=1000` for dCov and `n=2000` for HSIC, its absolute distance from 5% is slightly larger than the centred rule's. These small differences are not asserted to be statistically decisive; the predeclared gate is deliberately conservative. No candidate is called a corrected calibration, and no E1 wild-bootstrap p-value is replaced.

| `n` | HSIC centred | HSIC raw | dCov centred | dCov raw |
| ---: | ---: | ---: | ---: | ---: |
| 500 | 9.0 [7.4, 10.9]% | 3.9 [2.9, 5.3]% | 8.0 [6.5, 9.8]% | 3.7 [2.7, 5.1]% |
| 1,000 | 6.4 [5.0, 8.1]% | 3.6 [2.6, 4.9]% | 6.5 [5.1, 8.2]% | 3.4 [2.4, 4.7]% |
| 2,000 | 5.9 [4.6, 7.5]% | 4.0 [3.0, 5.4]% | 7.0 [5.6, 8.8]% | 3.9 [2.9, 5.3]% |

Each rate uses 1,000 independently generated pairs, 500 multiplier draws per pair, and a nominal 5% test. Brackets are Wilson 95% counting intervals. Both calibrations use the same pair and the same raw multiplier draws; only column centering changes. The original centred draws were checked to match `mbcore.markov_rademacher` exactly at all three sample sizes under the same random seed. The earlier 100-pair pilot was not used to choose a winner.

## Decision for the manuscript

The observed E1 wild-bootstrap rejection rates remain descriptive placebo rates. The manuscript can report this H=0.5 calibration check as a limit on interpretation, but cannot claim that subtracting the baseline excess or switching to raw multipliers corrects E1. The predeclared branch in `PLAN-revise-manuscript.md` is to omit the matched E1 rerun when no candidate passes the known-null gate.

## Reproduction

```sh
python investigations/hsic-dcov-calibration-2026-10-23/run_check.py --n 500,1000,2000 --reps 1000 --workers 4
```

The script writes `control_pairs.csv` and `control_summary.json`. The summary records seed `2026102301`, 500 multiplier draws, Python/NumPy-independent input hashes, and SHA-256 `e039836f3560a251df51563b06ce1fa0c05805e9206ed43684ab88045539c54c` for the output CSV. The source `mbcore.py` hash is in the JSON. This is a new run; original EXP1 cell files were not changed.
