"""`python -m kmem reproduce` (plan step 4.2): the sample check, the full-run plan, and that it catches changes."""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys

import pytest

from kmem import reproduce as R
from kmem import spec as S

SAVED = {"EXP2": 90000, "EXP3": 42000, "EXP4": 64000, "EXP5": 20788, "EXP6": 3000, "EXP7": 166400}


def _cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "kmem", "reproduce", *args], cwd=S.REPO, capture_output=True, text=True)


def test_full_plan_lists_exactly_the_saved_rows() -> None:
    assert {e: sum(R.planned_rows(e).values()) for e in SAVED} == SAVED
    p = _cli("EXP4", "--full", "--dry-run")
    assert p.returncode == 0 and "20 cells, 64000 rows" in p.stdout


def test_sample_check_passes_through_the_cli() -> None:
    p = _cli("EXP5", "--tier", "fast", "--workers", "4")
    assert p.returncode == 0, p.stdout + p.stderr
    assert "27 of 27 saved rows reproduce exactly" in p.stdout


def test_comparison_catches_a_changed_value(tmp_path) -> None:
    sp = S.load("EXP2")
    name = "mmd__G0__H0.5__n500"
    src = os.path.join(S.REPO, sp.folder, "main", name + ".csv")
    if not os.path.exists(src):
        pytest.skip("saved main outputs are not on this machine")
    shutil.copy(src, tmp_path / f"{name}.csv")
    assert R.passed(R.compare_run("EXP2", str(tmp_path)))
    with open(tmp_path / f"{name}.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    rows[5]["mmd2"] = repr(float(rows[5]["mmd2"]) * (1 + 1e-15))
    with open(tmp_path / f"{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    report = R.compare_run("EXP2", str(tmp_path))
    assert not R.passed(report) and report[name]["diffs"][0][1] == "mmd2"
