#!/bin/bash
# EXP5 pilot (DESIGN.md "Pilot"). 3 workers while EXP4 runs.
cd "$(dirname "$0")"
python -I -B run_exp5.py --mode pilot --reps 2 --out pilot --workers 3 --only ind__G4__H0.8__n8000
python -I -B run_exp5.py --mode pilot --reps 4 --out pilot --workers 3 --only ind__G4__H0.5__n4000,ind__G4__H0.8__n4000
python -I -B run_exp5.py --mode pilot --reps 50 --out pilot --workers 3 --only mmd__G0__null__H0.5__n2000,mmd__G0__null__H0.8__n2000,mmd__G0__shiftc2.2361__H0.5__n2000,mmd__G0__shiftc2.2361__H0.8__n2000
echo PILOTDONE
