#!/bin/bash
# EXP5 amendment A2: Part 2 in parallel with Part 1. 3 workers while a Part 1 pool runs, else 6.
cd "$(dirname "$0")"
W() { if pgrep -f '^python -I -B run_exp5.py.*--only mmd__' > /dev/null; then echo 3; else echo 6; fi; }
run() { local reps=$1 only=$2; local w; w=$(W); echo "== $(date -u +%FT%TZ) stage reps=$reps workers=$w only=$only"
  python -I -B run_exp5.py --mode main --reps "$reps" --out main --workers "$w" --only "$only"; }
for r in 125 250 375 500 625 750 875 1000; do run $r "ind__G4__H0.7__n4000,ind__G4__H0.8__n4000"; done
for r in 250 500; do run $r "ind__G4__H0.5__n4000"; done
echo PART2DONE
