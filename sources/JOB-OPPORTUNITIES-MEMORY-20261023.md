# Job JOB-OPPORTUNITIES-MEMORY-20261023 — literature review and gap search (topic: NEW-MEMORY-BOUNDARY plus proposals NEW-MB-KERNEL-G1-HSIC, G2-HFREE, G3-RANK)

Status: closed — discovery complete; independent review assigned  
Track: read-only opportunity assessment for topic: NEW-MEMORY-BOUNDARY plus proposals NEW-MB-KERNEL-G1-HSIC, G2-HFREE, G3-RANK; no track state change  
Owner: one Claude Code discovery worker, execution ID `opp-memory-discovery-20261023`  
Date opened: 2026-10-23

## Question

See scope.json topic: where inference for kernel tests, DML intervals and Hawkes Granger tests breaks as memory grows, and the three G1-G3 sub-questions.

## Inputs

- Author decision 2026-10-23: "the rest investigate, do lit-review and gap-search" (programme/decisions.jsonl, `investigate` entries dated 2026-10-23 for the catalogue/started tracks in scope; proposals G1-G3 remain proposed and are scanned read-only).
- Harness main at `e79a7f30f82592bd19f1fad96493c034ece1de03`; `programme/catalog.yaml` SHA-256 `83ea8ca007df48c188a31f8bac7c9b3df02ee60fc07986d331a04267c13df4af`; `registry.yaml` SHA-256 `90fe403f680c8f30be0bbd3051b83eba9e6f2d9f3ac754a8c4991f039e84866d`. Frozen scope: `runs/JOB-OPPORTUNITIES-MEMORY-20261023/scope.json`.
- Track IDs in scope: NEW-MEMORY-BOUNDARY NEW-MB-KERNEL-G1-HSIC NEW-MB-KERNEL-G2-HFREE NEW-MB-KERNEL-G3-RANK NEW-MB-KERNEL. Their catalogue entries, proposals (`programme/proposals/`), track repositories, and every earlier literature, WIDEGAP, DEEPLIT, GAPLIT and review job under `runs/` for these IDs. Earlier jobs are the starting point: refresh and extend them, do not repeat them.
- Previous programme scan: `runs/JOB-PROGRAMME-OPPORTUNITIES-20261012/` (report.md, candidates.json). Reuse stable candidate IDs where an idea recurs.

## Required skill and version

- `research-gap-finder`: `/home/ak/skills-sources/jup-local/academic-skills/research-gap-finder/SKILL.md` SHA-256 `e81cb1de0f36e7378f315418a0e2386571c03a84d25ceec15cec6e18f1d25919`; `references/gap_analysis_template.md` `a5c625dd5437974ea4bad80b9f153767842f887f4ac9247ebc17571904f34448`; conditional `references/problem_framing_template.md` `80a2e390c3ac17898b62aae9f5463b6e19acdb86f12fdeca4036016a52da40f2`. Hashes verified by the coordinator 2026-10-23.
- Broad literature route: `/home/ak/academic-research-skills/deep-research/SKILL.md` SHA-256 `f903840f07c447abe51fc7f52c6fe2354a9b15af25bfd254728a199cd26985be` with its `references/source_quality_hierarchy.md` and `references/semantic_scholar_api_protocol.md`. Bounded retrieval guidance only; no full multi-agent pipeline.
- Retrieval tools: WebSearch and WebFetch (live search plus primary-source open). Record any other tool actually used.

## Dependencies

Earlier jobs for these tracks are closed and ready. A separate reviewer job starts after this job's handoff.

## Output

Only `runs/JOB-OPPORTUNITIES-MEMORY-20261023/`: `events.jsonl`, `search-log.md`, `evidence-matrix.md`, `candidates.json`, `gap-analysis.md`, `handoff.md`, optional `checkpoint.md`. Broad pass across every ID in scope, then at most six shortlisted questions examined with comparable depth. `candidates.json` follows `research_harness.programme.opportunities.Candidate` and must pass `PYTHONPATH=src python -m research_harness.programme opportunities validate runs/JOB-OPPORTUNITIES-MEMORY-20261023/candidates.json`.

## Completion condition

Live searches and primary-source opens logged with dates; every ranked candidate backed by an inspected primary passage and a contrary search; full-text versus abstract-only read depth marked; actual skill hashes logged in `skill_loaded`; candidates validate; no decisions, source edits or experiments. A limited scan cannot establish novelty: no-paper-found is a search result.

## Required handoff

To the coordinator and a separate reviewer: per-track one-paragraph plain-English status of the literature, shortlisted candidates, closest primary sources, contrary results, uncertainties and suggested challenge searches.

## Scope

Worker writes only `runs/JOB-OPPORTUNITIES-MEMORY-20261023/`. Everything else is read-only, including track repositories, catalogue, decisions, registry and programme status.

## Execution checkpoints

Worker-owned `runs/JOB-OPPORTUNITIES-MEMORY-20261023/events.jsonl`: `started`, `skill_loaded`, search/source events and `finished` (with output paths and SHA-256), each with `execution_id` and `timestamp_utc`.

## State, result and next action

Assigned 2026-10-23 to execution `opp-memory-discovery-20261023`.

Coordinator close 2026-10-23: `finished` at 2026-10-23T22:02:06.874449+00:00 under `opp-memory-discovery-20261023`; `skill_loaded` hashes match the job record; all 5 output SHA-256 values in the finished event match the files; `candidates.json`: Validated 6 candidates. Result: [handoff](JOB-OPPORTUNITIES-MEMORY-20261023/handoff.md). Next: separate reviewer job `JOB-OPPORTUNITIES-MEMORY-REVIEW-20261023`.
