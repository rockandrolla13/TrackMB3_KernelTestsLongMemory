"""G-rows: saved replications of EXP2 to EXP7 must reproduce exactly (value-exact, decision D3).

Each experiment runs in its own subprocess (tests/repro/harness.py), all six at once, on at most 8 cores.
The tier comes from KMEM_REPRO_TIER: "fast" (default, about 25 s wall), "slow" or "all".
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
HARNESS = os.path.join(HERE, "harness.py")
EXPERIMENTS = ("EXP2", "EXP3", "EXP4", "EXP5", "EXP6", "EXP7")
TIER = os.environ.get("KMEM_REPRO_TIER", "fast")
WORKERS = {"EXP7": 5 if TIER == "slow" else 3}  # never more than 8 processes in total


def _expected(exp: str) -> set[str]:
    with open(os.path.join(HERE, "fixtures", f"{exp}.json")) as f:
        rows = json.load(f)["rows"]
    return {r["id"] for r in rows if TIER == "all" or r["tier"] == TIER}


def _cmd(exp: str, *extra: str) -> list[str]:
    return [sys.executable, "-I", HARNESS, exp, TIER, "--workers", str(WORKERS.get(exp, 1)), *extra]


@pytest.fixture(scope="session")
def results() -> dict:
    procs = {e: subprocess.Popen(_cmd(e), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for e in EXPERIMENTS if _expected(e)}
    out = {}
    for e, p in procs.items():
        so, se = p.communicate()
        out[e] = (p.returncode, so, se)
    return out


@pytest.mark.parametrize("exp", EXPERIMENTS)
def test_saved_rows_reproduce(exp: str, results: dict) -> None:
    expected = _expected(exp)
    if not expected:
        pytest.skip(f"{exp} has no {TIER}-tier rows")
    code, so, se = results[exp]
    assert code == 0, f"{exp} harness crashed:\n{se[-2000:]}"
    got = json.loads(so)
    assert set(got) == expected, f"{exp}: rows checked differ from the fixture"
    bad = {k: v for k, v in got.items() if v["diff"] or v["error"]}
    shown = {k: bad[k] for k in sorted(bad)[:10]}
    assert not bad, f"{exp}: {len(bad)} of {len(got)} rows differ, e.g. {json.dumps(shown)[:3000]}"


def test_gate_detects_a_seed_change() -> None:
    """Self-test: with every master shifted by one, a cheap EXP2 row must fail."""
    row = "ind__G0__H0.5__n500-r0"
    p = subprocess.run([sys.executable, "-I", HARNESS, "EXP2", "fast", "--only", row, "--perturb-seed", "1"],
                       capture_output=True, text=True, check=True)
    got = json.loads(p.stdout)
    assert set(got) == {row}
    assert got[row]["diff"], "a shifted seed reproduced the saved row: the gate cannot detect changes"
