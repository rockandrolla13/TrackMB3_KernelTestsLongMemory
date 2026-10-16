# Job JOB-NEW-MEMORY-BOUNDARY-LIT-20261015 — initial literature plan

Status: closed — no prior map found; boundary depends on statistic type; simulation design proposed; author decision pending  
Track: NEW-MEMORY-BOUNDARY  
Owner: Claude Code coordinator; one assigned worker, execution `claude-memboundary-20261015`

## Question

Where does inference for kernel tests, double-machine-learning causal intervals and Hawkes Granger tests break as memory grows, and does the break follow the H = 3/4 boundary in financial data?

## Inputs

- /media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness@f35d38c61a2bc04b5a81d7678581a161a01f2600
- /media/ak/10E1026C4FA6006E/GitRepos/research_assistant@be4179e5940ac5b61b378a73a6d3810688aea128

## Required skill and version

broad_literature; resolve from registry.yaml before assignment.

## Dependencies

Author review of this plan and a separate programme start decision.

## Output

/media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness/runs/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015; plan: /media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness/runs/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/literature-plan.yaml (SHA-256 9c1d4ffa4a5ef386cfb5b854fe9bf127769f9d4522a066443c8d02452918065f).

## Completion condition

Approved scope and source versions recorded; route and tools checked; worker outputs retain searches, selection, contrary evidence and limitations.

## Scope

Coordinator may prepare this plan. No retrieval or worker dispatch is authorised.

## State, result and next action

Pending author review of the literature plan. The experiment protocol remains unapproved.

## Author literature-plan decision

Approved plan SHA-256 9c1d4ffa4a5ef386cfb5b854fe9bf127769f9d4522a066443c8d02452918065f as ordinary Decision D-3debf5bb28. Next action: record a separate programme start decision before assigning retrieval. No experiment protocol was approved.

## Programme start decision

Author start event at 2026-10-15T22:09:06.006725Z. Next action: check source pins, required tools, dependencies and assignment readiness before assigning one worker.


## Coordinator close, 2026-10-15

Execution `claude-memboundary-20261015` finished with a `finished` event. Coordinator checks: nine output hashes match; plan hash `9c1d4ffa…` unchanged; `sanity_variance_scaling.py` reproduces byte for byte. Budget overrun: 79 records screened against ≤60 (worker crossed the cap at query 15 and ran 6 more; marked in `search-log.md`); 19/20 extracted; 3/12 full-text reads. SSRN and Google Scholar not effectively searched.

Result ([handoff](JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/handoff.md) `a254ca09…`): no paper found mapping size or coverage across memory levels for kernel tests, DML or Hawkes Granger tests. The H = 3/4 line is real only for quadratic-type statistics (dCov, Betken–Dehling Thm 2.2, F; sample covariance, Zhai–Zhong–Wu, F); for two series the condition is H_X + H_Y vs 3/2; for linear statistics (MMD two-sample mean-embedding term, DML population ATE with monotone heterogeneity) any long memory breaks the normal rate (boundary H = 1/2), which corrects C1's boundary to d = 0. HSIC extension is a worker inference. Hawkes Granger depends on branching ratio instead. Financial series mostly sit above 3/4 (durations, counts, volume, spreads 0.75–0.96; order signs straddle). C1 likely a corollary of Robinson–Hidalgo plus Taqqu. Proposed simulation design (not run): 2,000 reps, ≤48 CPU-hours per family, six kill criteria. Evidence class: bounded literature scan plus derivation.

Next action: author decisions — approve/change the design; accept "Hermite-rank map" framing; arm ownership (one job vs E, K/F-CAUSAL, B2-COMP); whether C1 continues; whether to fund an SSRN/Scholar follow-up.
