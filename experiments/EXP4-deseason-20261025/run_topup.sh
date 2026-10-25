#!/bin/bash
# EXP4 top-up (DESIGN.md amendment 5): restore n = 2000 cells to 1,000 replications. Resumable.
cd "$(dirname "$0")"
python -I run_exp4.py --mode main --reps 1000 --out main --workers 16 \
  --only "G0__H0.5__n2000,G0__H0.8__n2000,SEAS-V__H0.5__n2000,SEAS-V__H0.8__n2000" --adj2000 A1,A2,A3
python -I run_exp4.py --mode main --reps 1000 --out main --workers 16 \
  --only "SEAS-V-slow__H0.5__n2000,SEAS-V-slow__H0.8__n2000" --adj A0 --adj2000 A0
python -I analyze_exp4.py --in main --tag main
echo ALLDONE
