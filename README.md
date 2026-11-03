# TrackMB3_KernelTestsLongMemory — kernel tests under long memory

**Question.** Where do kernel two-sample, independence and conditional-independence tests (MMD, HSIC, distance covariance) lose size under long memory? Does the break follow the Hermite-rank map (linear statistics at any memory; quadratic statistics at H = 3/4, or H_X + H_Y = 3/2 for pairs) on financial-type series?

**Origin.** Split on 2026-10-16 from proposed track NEW-MEMORY-BOUNDARY. Harness track ID: NEW-MB-KERNEL. Overlaps Track E (kernel tests) and F-KERNEL.

## Layout

- `paper/long-memory-kernel-tests/`: working paper, PDF, figures, and Overleaf package. The paper separates known-null simulations from descriptive Nasdaq stock-duration placebos.
- `experiments/STOCK-20261023/`: stock-duration extraction, pair results, and request manifests. Raw licensed trade files are excluded.
- `experiments/EXP1-20261016/`: H5, E2, and E3 simulation extensions used in the figures and rule summary table.
- `FINDINGS.md`: what is known so far, with evidence class and source paths.
- `experiments/`: copied simulation outputs, byte for byte, with the original job record.
- `sources/`: copied literature and design jobs this track grew from.
- `hypotheses/`: hypotheses and claims (Agent 1).
- `literature/`: deep literature review, gap analysis and summary (Agent 2).
- `scoping/`: first exploratory experiments (Agent 3).
- `MANIFEST.sha256`: hashes of every copied file at copy time.

## Code: the `kmem` package

`kmem/` holds the shared simulation code: generators, memory kernels, calibration rules, the HSIC, dCov and
split-MMD bundles, and the analysis helpers. `kmem/legacy` is the validated EXP2 core and never changes.

Run everything from the repo root. `kmem` is not installed in the shared anaconda environment.

    python -m kmem env                                   # check numpy, scipy and Python against the pins
    python -I -m pytest tests -q                         # unit tests and the fast reproduction tier
    KMEM_REPRO_TIER=slow python -I -m pytest tests/repro # every remaining saved cell
    python -I tools/analysis_gate.py --check             # rerun every saved analysis, compare hashes

`pyproject.toml` pins the exact versions the saved EXP2-EXP7 rows were made with (Python 3.12.13,
numpy 2.5.1, scipy 1.17.1). Other versions may change the last bits of saved numbers.
`experiments/REPRODUCIBILITY.md` explains what "reproduce exactly" means here.
`kmem/EXTENDING.md` shows how to add a kernel, a generator, a calibration rule or an experiment grid.

To install `kmem` editable, use a separate environment, never the shared base one:

    conda create -n kmem python=3.12.13
    conda activate kmem
    pip install -e ".[analysis,dev]"

## Logs

ResearchHarness (`/media/ak/10E1026C4FA6006E/GitRepos/ResearchHarness`) keeps only the logs: the programme decisions, the job records in `runs/*.md`, and the original job outputs that those records hash. All research content for this track lives here. Copies were made on 2026-10-16. The originals stay in the harness because its records reference them by SHA-256.

## Rules

Exploratory evidence only so far. No accepted claim, locked protocol or Gate decision. Smoke and pilot runs are not evidence.
