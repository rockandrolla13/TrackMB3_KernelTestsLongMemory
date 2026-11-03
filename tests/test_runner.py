"""The shared runner reproduces saved rows through both entry points, and resumes correctly (plan step 4.1)."""
from __future__ import annotations

import csv
import glob
import json
import os
import subprocess
import sys

import pytest

from kmem import reproduce as R
from kmem import spec as S

TIER = os.environ.get("KMEM_REPRO_TIER", "fast")
# (experiment, cells, replications): cheap cells that cover every driver feature
FAST = [("EXP2", "ind__G0__H0.5__n500,mmd__G0__H0.5__n500", 3),
        ("EXP4", "ind__SEAS-V__H0.8__n500", 2),  # four adjustments per replication, reuse of A2 by A3
        ("EXP7", "growth__G0__null0__H0.8__n500,mmd__G0__null0__H0.5__n500,ind__G0__null0__H0.5__n500", 2)]
SLOW = [("EXP3", "ind__G0__lin0.13__H0.8__n500,mmd__G0__shift0.2__H0.5__n500", 2),
        ("EXP5", "mmd__G0__shiftc1.5__H0.8__n500", 2),
        ("EXP6", "ind__G0__null0__H0.5__n2000", 1)]
CASES = FAST + (SLOW if TIER in ("slow", "all") else [])


def _have_saved(exp: str) -> bool:
    return bool(glob.glob(os.path.join(S.REPO, S.EXPERIMENTS[exp], "main", "*.csv")))


def _kmem_run(exp: str, out: str, only: str, reps: int) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "kmem", "run", exp, "--out", out, "--workers", "2", "--only", only,
                           "--reps", str(reps)], cwd=S.REPO, capture_output=True, text=True)


@pytest.mark.parametrize("exp,only,reps", CASES, ids=[c[0] for c in CASES])
def test_kmem_run_reproduces_saved_rows(exp: str, only: str, reps: int, tmp_path) -> None:
    if not _have_saved(exp):
        pytest.skip("saved main outputs are not on this machine")
    p = _kmem_run(exp, str(tmp_path), only, reps)
    assert p.returncode == 0, p.stdout + p.stderr
    report = R.compare_run(exp, str(tmp_path))
    assert set(report) == set(only.split(","))
    assert R.passed(report), {k: v["diffs"][:3] for k, v in report.items() if v["diffs"] or not v["header_equal"]}


def _wrapper(exp: str, args: list[str]) -> subprocess.CompletedProcess:
    sp = S.load(exp)
    return subprocess.run([sys.executable, "-I", sp.driver, *args], cwd=os.path.join(S.REPO, sp.folder),
                          capture_output=True, text=True)


def test_old_driver_flags_reproduce_saved_rows(tmp_path) -> None:
    if not _have_saved("EXP2"):
        pytest.skip("saved main outputs are not on this machine")
    p = _wrapper("EXP2", ["--mode", "main", "--reps", "2", "--out", str(tmp_path), "--workers", "2",
                          "--only", "ind__HT-S__H0.7__n500"])
    assert p.returncode == 0, p.stderr
    assert R.passed(R.compare_run("EXP2", str(tmp_path)))
    new = json.load(open(glob.glob(str(tmp_path / "manifest_*.json"))[0]))
    old = json.load(open(os.path.join(S.REPO, S.load("EXP2").source_manifests[0])))
    assert set(old) <= set(new), set(old) - set(new)
    assert set(old["cells"][0]) <= set(new["cells"][0])


def test_resume_recomputes_truncated_and_error_rows(tmp_path) -> None:
    if not _have_saved("EXP2"):
        pytest.skip("saved main outputs are not on this machine")
    name = "mmd__G0__H0.8__n500"
    assert _kmem_run("EXP2", str(tmp_path), name, 2).returncode == 0
    path = tmp_path / f"{name}.csv"
    with open(path, "a", newline="") as f:  # a crash mid-write, then an error row
        f.write("2,ok,,0.1,0.5\n")
        csv.writer(f).writerow(["3", "error", "RuntimeError()", "0"])
    assert _kmem_run("EXP2", str(tmp_path), name, 4).returncode == 0
    with open(path, newline="") as f:
        ok = [r for r in csv.DictReader(f) if r["status"] == "ok" and None not in r.values()]
    assert sorted(int(r["rep"]) for r in ok) == [0, 1, 2, 3]
    assert R.passed(R.compare_run("EXP2", str(tmp_path)))


def test_a_different_header_is_refused(tmp_path) -> None:
    name = "mmd__G0__H0.8__n500"
    (tmp_path / f"{name}.csv").write_text("rep,status\n0,ok\n")
    p = _kmem_run("EXP2", str(tmp_path), name, 1)
    assert p.returncode != 0 and "header differs" in p.stdout + p.stderr
