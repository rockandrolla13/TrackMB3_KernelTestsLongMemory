# One command builds everything: `make` (or `make -B` to force a full rebuild).
# Shared machine: every step runs at nice 19 with at most 4 threads.
RUN = OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 nice -n 19
EXP = ../../../experiments/EXP7-memory-kernel-20261101
FIGS = figures/fig1_paths.png figures/fig2_growth.png figures/fig3_size.png \
       figures/fig4_power.png figures/fig5_seasonal.png figures/fig6_hsplit.png \
       figures/fig7_mmd_size.png figures/fig8_harris.png \
       figures/fig9_talk_size.png

DOC = EXP7-2026-11-03

all: $(DOC).pdf

data/report_numbers.json: build_data.py $(EXP)/summary_main.csv $(EXP)/results.json $(wildcard $(EXP)/main/*.csv)
	$(RUN) python build_data.py

numbers.tex: make_numbers.py data/report_numbers.json
	$(RUN) python make_numbers.py

$(FIGS) &: make_figures.py data/report_numbers.json
	$(RUN) python make_figures.py

$(DOC).pdf: $(DOC).tex numbers.tex $(FIGS)
	$(RUN) latexmk -pdf -interaction=nonstopmode -halt-on-error $(DOC).tex

clean:
	latexmk -C $(DOC).tex

.PHONY: all clean
