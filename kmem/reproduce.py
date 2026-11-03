"""Compare a new run with the saved rows of a past experiment (plan step 4.2, decision D3).

Value-exact: every saved column except machine measurements must match, as float(saved) == float(new), or both
NaN, or equal text. Rows are matched by replication (and EXP4's adjustment), so row order does not matter.
"""
from __future__ import annotations

import csv
import glob
import math
import os

from kmem import spec as S

SKIP = {"cpu_s", "wall_s", "maxrss_mb", "error"}  # machine measurements and error text
SKIP_BY_EXP = {"EXP4": {"reused"}}  # depends on how a run was chunked (code review CR-025)


def same(a: str, b: str) -> bool:
    try:
        x, y = float(a), float(b)
    except ValueError:
        return a == b
    return (math.isnan(x) and math.isnan(y)) or x == y


def _rows(path: str) -> tuple[list[str], dict]:
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        rows = list(rd)
        key = "adj" if "adj" in (rd.fieldnames or []) else None
        return rd.fieldnames, {(r["rep"], r[key] if key else None): r for r in rows}


def compare_cell(exp: str, name: str, run_dir: str) -> dict:
    sp = S.load(exp)
    saved_head, saved = _rows(os.path.join(S.REPO, sp.folder, "main", name + ".csv"))
    new_head, new = _rows(os.path.join(run_dir, name + ".csv"))
    skip = SKIP | SKIP_BY_EXP.get(exp, set())
    diffs = []
    for k, row in new.items():
        if k not in saved:
            diffs.append((k, "<row>", "absent", "present"))
            continue
        for col, sval in saved[k].items():
            if col not in skip and not same(sval, row.get(col, "")):
                diffs.append((k, col, sval, row.get(col)))
    return {"rows": len(new), "saved_rows": len(saved), "header_equal": saved_head == new_head, "diffs": diffs}


def compare_run(exp: str, run_dir: str, full: bool = False) -> dict:
    """Per-cell comparison of every CSV in run_dir. full: also require every saved row to be present."""
    out = {}
    for p in sorted(glob.glob(os.path.join(run_dir, "*.csv"))):
        name = os.path.basename(p)[:-4]
        out[name] = compare_cell(exp, name, run_dir)
    if full:
        for c in S.load(exp).cells:
            if c.name not in out:
                out[c.name] = {"rows": 0, "saved_rows": -1, "header_equal": False, "diffs": [("all", "<file>", "", "")]}
    return out


def passed(report: dict, full: bool = False) -> bool:
    return all(r["header_equal"] and not r["diffs"] and (not full or r["rows"] == r["saved_rows"])
               for r in report.values())


# ---------------------------------------------------------------- the committed sample (tests/repro/fixtures)

FIXTURES = os.path.join(S.REPO, "tests", "repro", "fixtures")
_STATE: dict = {}


def _sample_row(row: dict) -> tuple[str, list]:
    """Recompute one fixture row through the spec and the driver's job_for; list the columns that differ."""
    from kmem import runner

    sp, drv = _STATE["spec"], _STATE["driver"]
    cell = sp.cell(row["cell"]["name"])
    job = drv.job_for(cell.driver_cell(), cell.target, sp.master_for(cell), sp.options)
    subs = (tuple(s["adj"] for s in row["saved"]),) if job.subkeys is not None else ()
    new_rows = drv.run_chunk((job.fields, [row["rep"]], job.master) + tuple(job.extra) + subs)[-1]
    fmt = runner.FORMATTERS[sp.writer] or (lambda k, v: v)
    skip = SKIP | SKIP_BY_EXP.get(sp.id, set())
    by = {r.get("adj"): r for r in new_rows}
    diffs = []
    for saved in row["saved"]:
        new = by.get(saved.get("adj"))
        if new is None or new.get("status") != "ok":
            diffs.append((saved.get("adj"), "<row>", "ok", None if new is None else new.get("error")))
            continue
        for col, sval in saved.items():
            v = new.get(col)
            text = "" if v is None else str(fmt(col, v))
            if col not in skip and not same(sval, text):
                diffs.append((saved.get("adj"), col, sval, text))
    return row["id"], diffs


def sample(exp: str, tier: str = "all", workers: int = 1) -> dict[str, list]:
    """row id -> differing columns, for the committed fixture rows of one experiment (empty lists = reproduced)."""
    import json
    import multiprocessing as mp

    from kmem import drivers

    with open(os.path.join(FIXTURES, f"{exp}.json")) as f:
        rows = [r for r in json.load(f)["rows"] if tier == "all" or r["tier"] == tier]
    _STATE["spec"] = S.load(exp)
    _STATE["driver"] = drivers.load(_STATE["spec"])
    if workers > 1:
        with mp.get_context("fork").Pool(workers) as pool:
            return dict(pool.imap_unordered(_sample_row, rows))
    return dict(map(_sample_row, rows))


def planned_rows(exp: str) -> dict[str, int]:
    """Rows a full rerun writes, per cell: target times the number of adjustments (EXP4) or 1."""
    return {c.name: c.target * len(c.fields.get("adjustments", [None])) for c in S.load(exp).cells}
