# Job JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023 — independent contrary review

Status: closed  
Track: read-only review of `JOB-OPPORTUNITIES-MEMORY-20261023`; no track state change  
Owner: one Claude Code reviewer, execution ID `opp-memory-review-20261023` (different worker from discovery)  
Date opened: 2026-10-23

## Question

Do the shortlisted candidates in `runs/JOB-OPPORTUNITIES-MEMORY-20261023/candidates.json` survive an independent check of their closest sources and a fresh contrary search?

## Inputs

- Closed discovery job `runs/JOB-OPPORTUNITIES-MEMORY-20261023/` (handoff.md, candidates.json, evidence-matrix.md, gap-analysis.md, search-log.md; hashes in its `finished` event) and the sources it cites.
- Workflow `programme/RESEARCH_OPPORTUNITIES.md` step 6; research-gap-finder skill as in the discovery job record (same hashes).

## Output

Only `runs/JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023/`: `events.jsonl`, `review.md`, `contrary-search-log.md`, `source-checks.md`, `handoff.md`. For each shortlisted ID: closest work with section/theorem, direct answer or partial result, contrary evidence, remaining question, read depth, and a verdict keep / narrow / provisional / reject. At least one contrary or equivalent-formulation search per candidate. Mark inaccessible material. Do not assume the discovery worker is right.

## Completion condition

Each candidate assessed against sources the reviewer opened; dated searches logged; disagreements with the discovery worker stated explicitly; `finished` event with output hashes. No decisions, source edits or experiments.

## Scope

Reviewer writes only `runs/JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023/`. Everything else is read-only, including the discovery directory.

## State, result and next action

Assigned 2026-10-23 to `opp-memory-review-20261023`.

Coordinator close 2026-10-23: `finished` event present; all output SHA-256 values match; `skill_loaded` hash matches. Disagreements adjudicated in `runs/JOB-OPPORTUNITIES-MEMORY-20261023/coordinator-adjudication.md`; report rendered at `runs/JOB-OPPORTUNITIES-MEMORY-20261023/report.md` and linked from programme-status.md.
