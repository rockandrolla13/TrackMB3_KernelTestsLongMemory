#!/bin/bash
# EXP5 amendment A3: optional n = 8000 cells after Part 2, 6 workers, stages of 48 reps.
# Each stage starts only if the wall clock is before its cut-off (budget end about 21:50 UTC).
cd "$(dirname "$0")"
until grep -q PART2DONE part2.log; do sleep 15; done
before() { [ "$(date -u +%H%M)" \< "$1" ]; }
run() { local reps=$1 only=$2; echo "== $(date -u +%FT%TZ) stage reps=$reps workers=6 only=$only"
  python -I -B run_exp5.py --mode main --reps "$reps" --out main --workers 6 --only "$only"; }
for r in 48 96 144; do before 2125 && run $r "ind__G4__H0.8__n8000"; done
for r in 48 96 144; do before 2135 && run $r "ind__G4__H0.7__n8000"; done
echo PART3DONE
