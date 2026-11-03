"""Recompute fixture rows through an experiment's own row function and compare with the saved values (G-rows).

Run in a fresh process per experiment, so the drivers' bare-name imports (gens, mbcore, ...) cannot collide:
    python -I tests/repro/harness.py EXP7 fast [--workers 3]
It prints one JSON object: {row_id: {"diff": {column: [saved, new]}, "error": str}} for the requested tier
("fast", "slow" or "all"). An empty "diff" and no "error" means the row reproduced.

Comparison (decision D3, value-exact): every saved column except the machine measurements in SKIP must
match. A value matches when float(saved) == float(new), or both are NaN, or, for text, the strings are
equal. A column the row function did not return counts as the empty string, as csv.DictWriter writes it.
EXP7 values first pass through run_exp7.fmt, the formatter EXP7 wrote its CSVs with. EXP4's `reused` flag
depends on how a run was chunked (code review CR-025), so it is skipped for EXP4 only.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIXTURES = os.path.join(REPO, "tests", "repro", "fixtures")
SKIP = {"cpu_s", "wall_s", "maxrss_mb", "error"}
SKIP_BY_EXP = {"EXP4": {"reused"}}

_DRIVER = None
_FIX = None


def load_fixture(exp: str) -> dict:
    with open(os.path.join(FIXTURES, f"{exp}.json")) as f:
        return json.load(f)


def load_driver(fx: dict):
    path = os.path.join(REPO, fx["folder"], fx["driver"])
    spec = importlib.util.spec_from_file_location(f"repro_driver_{fx['experiment']}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_PERTURB = 0  # --perturb-seed: add this to every master, so the self-test can prove the gate fails


def chunk_args(fx: dict, row: dict) -> tuple:
    """The argument tuple the experiment's main() passed to run_chunk for this cell and one replication."""
    exp, cell, reps = fx["experiment"], row["cell"], [row["rep"]]
    opts = dict(fx["options"])
    if _PERTURB:
        if "masters" in opts:
            opts["masters"] = {k: v + _PERTURB for k, v in opts["masters"].items()}
        else:
            opts["master"] = opts["master"] + _PERTURB
    if exp == "EXP2":
        return (cell, reps, opts["masters"][f"{cell['family']}_main"])
    if exp in ("EXP3", "EXP5"):
        return (cell, reps, opts["master"])
    if exp == "EXP4":
        adjs = tuple(s["adj"] for s in row["saved"])
        return (cell, reps, opts["master"], adjs)
    if exp == "EXP6":
        return (cell, reps, opts["master"], tuple(opts["kappas"]), tuple(opts["oracle_kappas"]))
    if exp == "EXP7":
        exact = {int(k): int(v) for k, v in opts["exact_reps"].items()}
        return (cell, reps, opts["master"], False, exact, False)
    raise ValueError(exp)


def same(saved: str, new) -> bool:
    if new is None:
        new = ""
    try:
        a = float(saved)
    except ValueError:
        return saved == str(new)
    try:
        b = float(new)
    except (TypeError, ValueError):
        return False
    return (math.isnan(a) and math.isnan(b)) or a == b


def check_row(row: dict) -> tuple[str, dict]:
    fx, drv = _FIX, _DRIVER
    skip = SKIP | SKIP_BY_EXP.get(fx["experiment"], set())
    fmt = getattr(drv, "fmt", None) if fx["experiment"] == "EXP7" else None
    try:
        out = drv.run_chunk(chunk_args(fx, row))
    except Exception as e:  # report, do not hide
        return row["id"], {"diff": {}, "error": repr(e)[:500]}
    new_rows = out[-1]
    key = "adj" if fx["experiment"] == "EXP4" else None
    new_by = {(r.get(key) if key else None): r for r in new_rows}
    diff = {}
    for saved in row["saved"]:
        new = new_by.get(saved.get(key) if key else None)
        if new is None:
            diff[f"{saved.get(key)}:<row>"] = ["present", "missing"]
            continue
        if new.get("status") != "ok":
            return row["id"], {"diff": {}, "error": f"status {new.get('status')}: {new.get('error')}"}
        for col, sval in saved.items():
            if col in skip:
                continue
            nval = new.get(col, "")
            if fmt is not None and col in new:
                nval = fmt(col, nval)
            if not same(sval, nval):
                tag = f"{saved.get(key)}:{col}" if key else col
                diff[tag] = [sval, str(nval)]
    return row["id"], {"diff": diff, "error": ""}


def _init(exp: str) -> None:
    global _DRIVER, _FIX
    _FIX = load_fixture(exp)
    _DRIVER = load_driver(_FIX)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("tier", choices=["fast", "slow", "all"])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--only", default=None, help="comma list of row ids")
    ap.add_argument("--perturb-seed", type=int, default=0, help="self-test only: shift every master seed")
    a = ap.parse_args()
    global _PERTURB
    _PERTURB = a.perturb_seed
    _init(a.experiment)
    rows = [r for r in _FIX["rows"] if a.tier == "all" or r["tier"] == a.tier]
    if a.only:
        keep = set(a.only.split(","))
        rows = [r for r in rows if r["id"] in keep]
    if a.workers > 1:
        with mp.get_context("fork").Pool(a.workers) as pool:
            results = dict(pool.imap_unordered(check_row, rows))
    else:
        results = dict(map(check_row, rows))
    json.dump(results, sys.stdout, sort_keys=True)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
