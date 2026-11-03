"""Experiment specs agree with the saved runs and with each driver's grid code (plan step 3.1)."""
from __future__ import annotations

import glob
import importlib.util
import json
import os

import pytest
from pydantic import ValidationError

from kmem import seeds
from kmem import spec as S

EXPS = tuple(S.EXPERIMENTS)


def _driver(exp: str):
    sp = S.load(exp)
    path = os.path.join(S.REPO, sp.folder, sp.driver)
    m = importlib.util.spec_from_file_location(f"testspec_{exp}", path)
    mod = importlib.util.module_from_spec(m)
    m.loader.exec_module(mod)
    return mod


def _driver_grid(exp: str) -> dict[str, dict]:
    """The cells each driver's own code builds for its main run, keyed by name."""
    d = _driver(exp)
    if exp == "EXP2":
        cells = d.cell_list("mmd") + d.cell_list("ind")
    elif exp == "EXP3":
        cells = d.main_cells()
        for i, c in enumerate(cells):
            c["idx"], c["name"] = i, d.name_of(c)
    elif exp in ("EXP4", "EXP5"):
        cells = d.cell_list()
    elif exp == "EXP6":
        cells = d.cells_for("main", (2000,))
    else:
        cells = d.main_cells()
    return {c["name"]: c for c in cells}


@pytest.mark.parametrize("exp", EXPS)
def test_spec_cells_are_the_manifest_cells(exp: str) -> None:
    sp = S.load(exp)
    recorded = {}
    for p in sp.source_manifests:
        with open(os.path.join(S.REPO, p)) as f:
            for c in json.load(f)["cells"]:
                recorded[c["name"]] = c
    assert {c.name for c in sp.cells} == set(recorded)
    for c in sp.cells:
        r = recorded[c.name]
        assert c.seed_index == r["idx"]
        assert all(c.fields[k] == r[k] for k in c.fields if k != "adjustments")


@pytest.mark.parametrize("exp", EXPS)
def test_spec_names_are_the_saved_files(exp: str) -> None:
    sp = S.load(exp)
    saved = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(S.REPO, sp.folder, "main", "*.csv"))}
    if not saved:
        pytest.skip("saved main outputs are not on this machine")
    assert saved == {c.name for c in sp.cells}


@pytest.mark.parametrize("exp", EXPS)
def test_driver_grid_still_gives_the_saved_seed_indices(exp: str) -> None:
    grid = _driver_grid(exp)
    for c in S.load(exp).cells:
        g = grid[c.name]
        assert g["idx"] == c.seed_index, c.name
        shared = set(g) & set(c.fields) - {"adjustments", "target"}
        assert all(g[k] == c.fields[k] for k in shared), c.name


def test_saved_row_totals() -> None:
    want = {"EXP2": 90000, "EXP3": 42000, "EXP4": 64000, "EXP5": 20788, "EXP6": 3000, "EXP7": 166400}
    assert {e: S.load(e).total_rows for e in EXPS} == want


def test_master_registry_has_no_unexplained_reuse() -> None:
    seeds.check()
    assert set(seeds.duplicates()) == set(seeds.ALLOWED_REUSE)


def test_validators_reject_bad_specs() -> None:
    raw = S.load("EXP6").model_dump(mode="json")
    bad = json.loads(json.dumps(raw))
    bad["cells"][1]["seed_index"] = bad["cells"][0]["seed_index"]
    bad["cells"][1]["fields"]["idx"] = bad["cells"][0]["seed_index"]
    with pytest.raises(ValidationError, match="share master and seed index"):
        S.ExperimentSpec.model_validate(bad)
    bad = json.loads(json.dumps(raw))
    bad["cells"][0]["seed_index"] += 1000
    with pytest.raises(ValidationError, match="disagree with fields"):
        S.ExperimentSpec.model_validate(bad)
    bad = json.loads(json.dumps(raw))
    bad["masters"]["main"] = 1
    with pytest.raises(ValidationError, match="not in kmem.seeds.MASTERS"):
        S.ExperimentSpec.model_validate(bad)


def test_prerequisites_are_intact() -> None:
    for e in EXPS:
        for p in S.load(e).prerequisites:
            assert p.verify(), p.path


def test_declared_inputs_exist_and_match() -> None:
    for e in EXPS:
        for i in S.load(e).inputs:
            if not os.path.exists(os.path.join(S.REPO, i.path)):
                pytest.skip(f"{i.path} is not on this machine")
            assert i.verify(), f"{e}: input {i.path} changed since the spec was built"


def test_a_missing_declared_input_is_a_clear_error(tmp_path) -> None:
    """EXP4's analysis in a scratch copy without EXP2's summary must stop and name the file (CR-009)."""
    import shutil
    import subprocess
    import sys

    sp = S.load("EXP4")
    work = tmp_path / "experiments" / os.path.basename(sp.folder)
    work.mkdir(parents=True)
    for n in os.listdir(os.path.join(S.REPO, sp.folder)):
        if n.endswith(".py"):
            shutil.copy(os.path.join(S.REPO, sp.folder, n), work / n)
    os.symlink(os.path.join(S.REPO, sp.folder, "main"), work / "main")
    shutil.copytree(os.path.join(S.REPO, "kmem"), tmp_path / "kmem", ignore=shutil.ignore_patterns("__pycache__"))
    p = subprocess.run([sys.executable, "-I", "analyze_exp4.py", "--in", "main", "--tag", "scratch"], cwd=work,
                       capture_output=True, text=True)
    assert p.returncode != 0
    assert "declared input missing" in p.stderr and "EXP2-20261023/summary_main.csv" in p.stderr
