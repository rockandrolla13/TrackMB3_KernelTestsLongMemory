# Hermite rank and memory sums shape kernel-test size under long memory

Edit `paper.tex` and `refs.bib`; build with:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper.tex
```

The [public repository](https://github.com/rockandrolla13/TrackMB3_KernelTestsLongMemory) contains this working paper, simulation code and outputs, and derived stock placebo results. The PDF is `build/paper.pdf`; `REPRODUCE.md` gives rebuild commands. Run `python make_figures.py` to regenerate all three figures from saved summaries. The stock-duration section is descriptive because process independence is unverified. The Hermite-rank calculation is a mechanism, not a proved null law for these exact tests.

## Evidence map

| Manuscript section | Primary repo evidence | Qualification |
| --- | --- | --- |
| Theory, Section 2 | `../../sources/JOB-NEW-MEMORY-BOUNDARY-LIT-20261015/boundary-theorems.md`; `../../literature/review.md`; `../../literature/bibliography.bib` | Betken--Dehling's distance-covariance theorem has narrower conditions than the HSIC and mixed-pair predictions. |
| Initial MMD design and results, Sections 3.1 and 4.1 | `../../experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_mmd-dml/results.md`; `../../experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_mmd-dml/run.py`; `../../FINDINGS.md` | Initial wild-bootstrap baseline is oversized at H=0.5. |
| Initial HSIC/dCov/Pearson map, Sections 3.1 and 4.2 | `../../experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear/results.md`; `../../experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear/run_sim.py` | Kernel n=4000 was dropped by author decision. |
| EXP1 H5, E2, E3, Sections 3.2 and 4 | `../../experiments/EXP1-20261016/` scripts, summaries, manifest, and cell CSVs | The slope plot describes adjacent grid changes; it is not a breakpoint estimate. |
| Rule summary table, Section 4 | Initial independence results and `../../experiments/EXP1-20261016/e3_summary.csv` | The table separates baseline error, memory effects, and sample-size changes. |
| Paired multiplier control, Sections 3.2 and 4.2 | `../../investigations/hsic-dcov-calibration-2026-10-23/report.md`; `control_pairs.csv`, `control_summary.json`, `DESIGN.md` in that directory | Uncentred multipliers failed the predeclared non-worsening rule; the stock arm retains the original centred rule. |
| Nasdaq stock duration placebo, Sections 3.2, 4.3, and 5 | `../../experiments/STOCK-20261023/README.md`; source and extraction manifests, `stock_placebo_summary.json`, `stock_placebo_results.md`, and pair CSVs in that directory | The saved-pair audit passed. The market processes' independence remains unverified and days are reused. |
| Figure generation, Sections 4.1 and 4.2 | `make_figures.py`; `figures/manifest.json` | Reads saved E2/E3 summaries. No new simulation or theorem. |
| Scope decisions, Section 4.2 | `../../FINDINGS.md` | The 3/4 description is retained for wild bootstrap, shift, EWC, and Pearson--NW, and dropped for subsampling. These are author decisions, not proofs. |
| Further proof questions, Section 5 | `CLAIM-LEDGER.md`, `LITERATURE-AUDIT.md` | These remain conjectures and research directions. |

The bibliography copies only identities and bibliographic details already recorded in `../../literature/bibliography.bib`. It makes no new DOI or URL claims.

## Items to resolve before external circulation

- The candidate HSIC/distance-covariance multiplier control failed its stated selection rule. Any replacement critical-value rule needs separate validation. The stock day gap does not establish the placebo null.
- Keep the headline as a finite-sample simulation result unless a statistic-specific theorem is proved under explicit assumptions.
- Before journal submission, freeze an archival release and complete verified author, funding, competing-interest, AI-use, and data-access statements.
