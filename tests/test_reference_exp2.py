"""Reference test: kmem must reproduce saved EXP2 replication rows exactly.

Each case rebuilds the replication's generator from the EXP2 seed rule,
SeedSequence(entropy=master, spawn_key=(cell_idx, rep)), runs it through kmem, and
compares every output column with experiments/EXP2-20261023/main/<cell>.csv.
Cell indices and master seeds come from the EXP2 main-run manifest, not from code.

Run: PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -I -m pytest tests/
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os

import numpy as np
import pytest

from kmem.legacy import gens
from kmem.legacy import mbcore as mb
from kmem.legacy import mmdcal as M

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXP2 = os.path.join(REPO, "experiments", "EXP2-20261023")
B_WB = 500  # run_exp2.B_WB
SKIP = {"rep", "status", "error", "cpu_s"}

# About 20 replications per family: mostly n = 500 for speed, two n = 2000 cells for coverage.
CASES = {
    "ind": {"ind__G0__H0.8__n500": range(6), "ind__SEAS-V__H0.5__n500": range(6),
            "ind__NL-SV__H0.7__n500": range(6), "ind__HT-L__H0.8__n2000": range(2)},
    "mmd": {"mmd__G0__H0.8__n500": range(6), "mmd__NL-TAR__H0.5__n500": range(6),
            "mmd__SEAS-M__H0.7__n500": range(6), "mmd__HT-S__H0.7__n2000": range(2)},
}


def _manifest() -> dict:
    path = sorted(glob.glob(os.path.join(EXP2, "main", "manifest_*.json")))[0]
    with open(path) as f:
        return json.load(f)


MANIFEST = _manifest()
CELLS = {c["name"]: c for c in MANIFEST["cells"]}


def _saved_rows(name: str) -> dict[int, dict]:
    with open(os.path.join(EXP2, "main", name + ".csv")) as f:
        return {int(r["rep"]): r for r in csv.DictReader(f)}


def _run(cell: dict, rep: int) -> dict:
    master = MANIFEST["masters"][f"{cell['family']}_main"]
    rng = np.random.default_rng(np.random.SeedSequence(entropy=master, spawn_key=(cell["idx"], rep)))
    if cell["family"] == "mmd":
        x = gens.generate(cell["gen"], cell["H"], cell["n"], rng)
        return M.rep_e2(x, rng)
    x, y = gens.pair(cell["gen"], cell["H"], cell["n"], rng)
    out = mb.kind_tests(x, y, rng, B=B_WB)
    out.update(mb.lind_tests(x, y))
    return out


def _same(saved: str, new) -> bool:
    a, b = float(saved), float(new)
    return (math.isnan(a) and math.isnan(b)) or a == b


PARAMS = [(name, rep) for fam in CASES.values() for name, reps in fam.items() for rep in reps]


@pytest.mark.parametrize("name,rep", PARAMS, ids=[f"{n}-r{r}" for n, r in PARAMS])
def test_reproduces_saved_row(name: str, rep: int) -> None:
    saved = _saved_rows(name)[rep]
    assert saved["status"] == "ok"
    new = _run(CELLS[name], rep)
    cols = [k for k in saved if k not in SKIP]
    missing = [k for k in cols if k not in new]
    assert not missing, f"kmem did not return columns {missing}"
    diff = {k: (saved[k], new[k]) for k in cols if not _same(saved[k], new[k])}
    assert not diff, f"{name} rep {rep} differs: {diff}"


def test_manifest_matches_environment() -> None:
    """Exact reproduction is only promised under the numpy version EXP2 used."""
    assert np.__version__ == MANIFEST["numpy"]
