#!/bin/bash
# EXP5 amendment A2: Part 1 second half (reps 500..999) after the first stage finishes. 3 workers.
cd "$(dirname "$0")"
while pgrep -f '^python -I -B run_exp5.py.*--only mmd__' > /dev/null; do sleep 10; done
echo "== $(date -u +%FT%TZ) stage reps=1000 workers=3 only=mmd__"
python -I -B run_exp5.py --mode main --reps 1000 --out main --workers 3 --only mmd__
echo PART1DONE
