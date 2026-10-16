# Job JOB-NEW-MKC-LIT-20261015 — initial literature plan

Status: closed — field map, gap maps and ranked directions delivered; author decision pending  
Track: NEW-MEMORY-KERNEL-CAUSAL  
Owner: Claude Code coordinator; two workers — `claude-mkc-causal-20261015` (writes `memory-causal/` only) and `claude-mkc-kernel-20261015` (writes `memory-kernel/` only); coordinator writes `synthesis.md`

## Question

Where do memory (long-range dependence, self-excitation, rough or fractional dynamics), kernel methods and causal inference for time series meet, and what is open there for financial data?

## Inputs

- /media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness@f35d38c61a2bc04b5a81d7678581a161a01f2600
- /media/ak/10E1026C4FA6006E/GitRepos/research_assistant@be4179e5940ac5b61b378a73a6d3810688aea128

## Required skill and version

Bound routes (`binding.json`): `broad_literature` (deep-research `f903840f…`), `recent_papers` (paper-lit-finder `9632a242…`), `gap_synthesis` (research-gap-finder `e81cb1de…`), `scholarly_metadata`. Plan approved and track started by the author's request (`D-e0e0a1852b` start). Workers log `skill_loaded` with actual hashes.

## Dependencies

Author review of this plan and a separate programme start decision.

## Output

/media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness/runs/JOB-NEW-MKC-LIT-20261015; plan: /media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness/runs/JOB-NEW-MKC-LIT-20261015/literature-plan.yaml (SHA-256 6f6631ad4cb13f3abc3238c513498db95abe992ebff700fba1aeea8e938c99e3).

## Completion condition

Approved scope and source versions recorded; route and tools checked; worker outputs retain searches, selection, contrary evidence and limitations.

## Scope

Coordinator may prepare this plan. No retrieval or worker dispatch is authorised.

## State, result and next action

Pending author review of the literature plan. The experiment protocol remains unapproved.

## Author literature-plan decision

Approved plan SHA-256 6f6631ad4cb13f3abc3238c513498db95abe992ebff700fba1aeea8e938c99e3 as ordinary Decision D-e0e0a1852b. Next action: record a separate programme start decision before assigning retrieval. No experiment protocol was approved.

## Programme start decision

Author start event at 2026-10-15T21:40:50.962199Z. Next action: check source pins, required tools, dependencies and assignment readiness before assigning one worker.


## Coordinator close, 2026-10-15

Both workstreams finished with `finished` events; output hashes checked (memory-causal 6/6, memory-kernel 6/6); plan hash `6f6631ad…` unchanged. Budget: memory-causal screened 74 of ≤70 (overrun logged by the worker), memory-kernel at its caps. Coordinator synthesis `synthesis.md` (`8bb461a17adf2225…`).

Result: estimation mostly survives long memory, inference mostly does not; kernel-type statistics have a boundary at H = 3/4, and finance straddles it (order signs ≈ 0.7, realised volatility ≈ 0.9). Ranked directions: (1) a no-data memory breakdown map combining kernel-test size, DML coverage and Hawkes Granger size; (2) single-path signature-kernel causal discovery under fractional noise (the only full triple-intersection item); (3) core vs reaction order flow as causal impact; (4) RKHS power-law Hawkes/propagator kernels near criticality; (5) proxy causal identification under long memory (Track J). Most items overlap existing tracks (E, K, F-CAUSAL, B2-*, J, G); only (1) and (2) are outside existing scopes. Evidence class: bounded literature scan; no novelty claim.

Next action: author chooses a direction; the breakdown-map simulation is the recommended first job.
