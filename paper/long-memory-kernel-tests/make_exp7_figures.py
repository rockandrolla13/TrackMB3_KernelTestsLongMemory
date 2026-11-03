"""Draw the two EXP7 figures for paper.tex.

Reuses fig3 (HSIC size against H) and fig7 (split-MMD size) from the checked EXP7 examples
report, investigations/memory-kernel-2026-11-01/examples-report/make_figures.py. That script
reads only data/report_numbers.json, which build_data.py derives from the EXP7 main-run files.
No simulation is rerun. Writes PDF to figures/ and PNG to the Overleaf copy's figures/.

Run from this folder: OMP_NUM_THREADS=1 nice -n 19 python make_exp7_figures.py
(not -I: the matplotlib that works here lives in the user site)
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "investigations" / "memory-kernel-2026-11-01" / "examples-report" / "make_figures.py"
PDF_DIR = HERE / "figures"
PNG_DIR = HERE / "overleaf-upload-2026-10-23" / "figures"
NAMES = {"fig3_size.png": "exp7_hsic_size", "fig7_mmd_size.png": "exp7_mmd_size"}

spec = importlib.util.spec_from_file_location("exp7_report_figures", SRC)
mf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mf)


def save(fig, name):
    stem = NAMES[name]
    fig.savefig(PDF_DIR / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(PNG_DIR / f"{stem}.png", bbox_inches="tight", facecolor="white", dpi=200)
    mf.plt.close(fig)
    print("wrote", PDF_DIR / f"{stem}.pdf", "and", PNG_DIR / f"{stem}.png")


mf.save = save
mf.fig3()
mf.fig7()
