#!/bin/bash
# Amendment 1: D0 only at n = 2000 (regression check). D1 and D2 at all sizes. Resumable.
cd "$(dirname "$0")/../.."
python -I experiments/STOCK-DIURNAL-20261025/run_diurnal_placebo.py --pairs 500 --workers 2 --n 2000 --adj D0 D1 D2
python -I experiments/STOCK-DIURNAL-20261025/run_diurnal_placebo.py --pairs 500 --workers 2 --n 500 1000 --adj D1 D2
echo ALLDONE
