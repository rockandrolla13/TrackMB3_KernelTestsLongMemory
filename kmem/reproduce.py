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
