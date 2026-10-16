# TrackMB3_KernelTestsLongMemory — kernel tests under long memory

**Question.** Where do kernel two-sample, independence and conditional-independence tests (MMD, HSIC, distance covariance) lose size under long memory? Does the break follow the Hermite-rank map (linear statistics at any memory; quadratic statistics at H = 3/4, or H_X + H_Y = 3/2 for pairs) on financial-type series?

**Origin.** Split on 2026-10-16 from proposed track NEW-MEMORY-BOUNDARY. Harness track ID: NEW-MB-KERNEL. Overlaps Track E (kernel tests) and F-KERNEL.

## Layout

- `FINDINGS.md`: what is known so far, with evidence class and source paths.
- `experiments/`: copied simulation outputs, byte for byte, with the original job record.
- `sources/`: copied literature and design jobs this track grew from.
- `hypotheses/`: hypotheses and claims (Agent 1).
- `literature/`: deep literature review, gap analysis and summary (Agent 2).
- `scoping/`: first exploratory experiments (Agent 3).
- `MANIFEST.sha256`: hashes of every copied file at copy time.

## Logs

ResearchHarness (`/media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness`) keeps only the logs: the programme decisions, the job records in `runs/*.md`, and the original job outputs that those records hash. All research content for this track lives here. Copies were made on 2026-10-16. The originals stay in the harness because its records reference them by SHA-256.

## Rules

Exploratory evidence only so far. No accepted claim, locked protocol or Gate decision. Smoke and pilot runs are not evidence.
