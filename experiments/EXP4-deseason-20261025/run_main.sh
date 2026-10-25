#!/bin/bash
# EXP4 main run, sequential stages, 5 workers. Resumable: re-run this script.
cd "$(dirname "$0")"
python -I run_exp4.py --mode main --reps 500 --out main --workers 5 \
  --only "G0__H0.5__n2000,G0__H0.8__n2000,SEAS-V__H0.5__n2000,SEAS-V__H0.8__n2000" --adj2000 A1,A2,A3
python -I run_exp4.py --mode main --reps 1000 --out main --workers 5 \
  --only "G0__H0.5__n500,G0__H0.7__n500,G0__H0.8__n500,SEAS-M__H0.5__n500,SEAS-M__H0.7__n500,SEAS-M__H0.8__n500,SEAS-V__H0.5__n500,SEAS-V__H0.7__n500,SEAS-V__H0.8__n500,SEAS-M-shared__H0.5__n500,SEAS-M-shared__H0.7__n500,SEAS-M-shared__H0.8__n500"
python -I run_exp4.py --mode main --reps 1000 --out main --workers 5 \
  --only "SEAS-V-slow__H0.5__n500,SEAS-V-slow__H0.8__n500" --adj A0
python -I run_exp4.py --mode main --reps 300 --out main --workers 5 \
  --only "SEAS-V-slow__H0.5__n2000,SEAS-V-slow__H0.8__n2000" --adj A0 --adj2000 A0
echo ALLDONE
