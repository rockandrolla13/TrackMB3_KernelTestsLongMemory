#!/bin/bash
# EXP5 main run (DESIGN.md amendment A1). Resumable. Worker count is re-checked before each
# stage: 3 while run_exp4.py is running, 6 after it finishes.
cd "$(dirname "$0")"
W() { if pgrep -f run_exp4.py > /dev/null; then echo 3; else echo 6; fi; }
run() { local reps=$1 only=$2; local w; w=$(W); echo "== $(date -u +%FT%TZ) stage reps=$reps workers=$w only=$only"
  python -I -B run_exp5.py --mode main --reps "$reps" --out main --workers "$w" --only "$only"; }
# Part 1: split MMD, all 18 cells, 1,000 reps (n = 2000 cells in two halves to re-check workers)
run 500 "mmd__"
run 1000 "mmd__"
# Part 2: sums 1.4 and 1.6 at n = 4000, in steps of 250 reps
for r in 250 500 750 1000; do run $r "ind__G4__H0.7__n4000,ind__G4__H0.8__n4000"; done
# H = 0.5 baseline at n = 4000
for r in 250 500; do run $r "ind__G4__H0.5__n4000"; done
echo MAINDONE
