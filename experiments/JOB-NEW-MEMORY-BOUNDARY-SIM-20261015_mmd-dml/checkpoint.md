# Checkpoint — mmd-dml worker (claude-mb-mmddml-20261015)

State at 2026-10-15 ~22:58 UTC: validation passed (validation.out); pilot done (pilot/, not evidence;
pilot_v1_blas_multithread/ is a discarded mis-timed pilot); parameters_fixed and parameters_amended
logged. Main run started in background: DML (2000 reps/cell) then K-hom (2000 reps, 1000 at n=4000).

To resume after interruption: re-run the same commands; run.py skips reps already in the per-cell CSVs.

    source ~/.venvs/researchHarness/bin/activate
    python -I run.py dml main --reps 2000 --chunk 50
    python -I run.py khom main --reps 2000 --reps-largest 1000 --chunk 25
    python -I analyze.py main

Remaining after the run: analyze, results.md, handoff.md, finished event with hashes.

Update 2026-10-16: main runs complete (DML 00:58 UTC, K-hom 04:07 UTC per coordinator; 0 errors);
analysis, results.md, handoff.md and the finished event written. Nothing remains to resume.
