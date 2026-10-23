# Reproducing the revised draft

This guide uses saved experiment outputs in the [public repository](https://github.com/rockandrolla13/TrackMB3_KernelTestsLongMemory). Raw licensed market trades are excluded, and the placebo null is not established.

## Inputs and environment

- Manuscript source: `paper.tex`, `refs.bib`.
- Saved simulation summaries: `../../experiments/EXP1-20261016/e2_summary.csv` and `e3_summary.csv`. H5, E2, and E3 cell CSVs, scripts, and per-arm manifests are under the same EXP1 directory.
- Stock-placebo source: `../../experiments/STOCK-20261023/`. Its README lists the Databento request, extraction, saved-pair audit, and local raw-data location.
- Figure input and output SHA-256 hashes: `figures/manifest.json`.
- The source tree is public on `main`; the precise revision is given by its Git commit. This is a working draft, not a frozen archival release.
- The saved simulation arm manifests record Python 3.12.13, NumPy 2.5.2, seeds, and cell counts. The stock arm records its own request and run manifests.
- The local figure/PDF rebuild used Python 3.12.13, NumPy 2.5.1, Matplotlib 3.10.8, TeX Live 2023, and latexmk 4.83. It reads saved simulation summaries rather than regenerating them. LaTeX needs `pdflatex`, `natbib`, `booktabs`, `graphicx`, `float`, and the packages loaded by `paper.tex`.

## Rebuild checks from the repository root

```sh
python experiments/EXP1-20261016/analyze.py h5 e2 e3
python paper/long-memory-kernel-tests/make_figures.py
python experiments/STOCK-20261023/summarize_stock.py
```

The figure generator checks the expected grids before writing all three figures. The stock summary script audits saved pair IDs, day gaps, and within-day window boundaries; it does not query Databento.

The candidate HSIC/dCov multiplier control is separate from the original EXP1 outputs:

```sh
python investigations/hsic-dcov-calibration-2026-10-23/run_check.py --n 500,1000,2000 --reps 1000 --workers 4
```

It compares centred and uncentred Markov multipliers on the same independently generated `H=0.5` pairs. The script writes its own CSV and summary JSON. The predeclared gate failed; see `../../investigations/hsic-dcov-calibration-2026-10-23/report.md`. The stock analysis uses the original centred rule.

## Build the PDF

```sh
cd paper/long-memory-kernel-tests
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper.tex
pdfinfo build/paper.pdf
```

The current build is `build/paper.pdf`. The Overleaf upload package is `overleaf-upload-2026-10-23.zip`.

Before journal submission, the authors must freeze an archival release and complete verified author and disclosure fields. Raw market-data access remains subject to the source licence.
