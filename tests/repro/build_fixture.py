"""Build the committed row fixtures for the reproduction gate (G-rows). Run once; the output is committed.

For each of EXP2 to EXP7 this reads the saved main-run manifests (seeds, cell indices, run options) and the
saved main CSVs, picks a deterministic sample of replications, and writes tests/repro/fixtures/<EXP>.json.
Each fixture row holds the cell as the manifest recorded it, the replication number, the tier, and the
saved row values as text, exactly as they appear in the CSV.

Sample rule (every cell is covered once across the two tiers):
  n <= 500  : fast tier, replications 0 and the last saved one.
  n == 2000 : one replication (0). Every k-th such cell is fast, the rest slow; k keeps the fast tier near
              FAST_BUDGET CPU-seconds per experiment, measured from the saved cpu_s column.
  n >= 4000 : slow tier, replication 0.
  EXP7 adds one non-exact-path replication (rep = exact_reps[n]) for the first null cell at each n.

Usage:  python -I tests/repro/build_fixture.py            (writes all six fixtures)
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "tests", "repro", "fixtures")
FAST_BUDGET = 45.0  # CPU-seconds per experiment for the n = 2000 part of the fast tier

EXPERIMENTS = {
    "EXP2": ("experiments/EXP2-20261023", "run_exp2.py"),
    "EXP3": ("experiments/EXP3-power-20261025", "run_exp3.py"),
    "EXP4": ("experiments/EXP4-deseason-20261025", "run_exp4.py"),
    "EXP5": ("experiments/EXP5-scaled-n-20261025", "run_exp5.py"),
    "EXP6": ("experiments/EXP6-bandwidth-20261025", "run_exp6.py"),
    "EXP7": ("experiments/EXP7-memory-kernel-20261101", "run_exp7.py"),
}
DROP_CELL_KEYS = {"already_done", "already_done_pairs", "reps", "seed", "target_done"}


def manifests(folder: str) -> list[dict]:
    paths = sorted(glob.glob(os.path.join(REPO, folder, "main", "manifest_*.json")))
    if not paths:
        sys.exit(f"{folder}: no main manifest")
    return [dict(json.load(open(p)), _file=os.path.relpath(p, REPO)) for p in paths]


def merged_cells(mans: list[dict]) -> dict[str, dict]:
    """name -> cell; every manifest that lists a cell must agree on its seed index."""
    cells: dict[str, dict] = {}
    for m in mans:
        for c in m["cells"]:
            c = {k: v for k, v in c.items() if k not in DROP_CELL_KEYS}
            old = cells.get(c["name"])
            if old is None:
                cells[c["name"]] = c
            elif old["idx"] != c["idx"]:
                sys.exit(f"{c['name']}: manifests disagree on idx ({old['idx']} vs {c['idx']})")
            elif "adjustments" in c:  # EXP4: keep the union of adjustment lists
                old["adjustments"] = sorted(set(old["adjustments"]) | set(c["adjustments"]))
    return cells


def saved_rows(folder: str, name: str) -> list[dict]:
    with open(os.path.join(REPO, folder, "main", name + ".csv"), newline="") as f:
        return list(csv.DictReader(f))


def run_options(exp: str, mans: list[dict]) -> dict:
    m = mans[0]
    if exp == "EXP2":
        opts = {"masters": m["masters"]}
    else:
        opts = {"master": m["master"]}
    if exp == "EXP6":
        opts["kappas"] = m["kappas"]
        opts["oracle_kappas"] = m["oracle_kappas"]
    if exp == "EXP7":
        opts["exact_reps"] = m["exact_reps"]
    for other in mans[1:]:  # options must not differ between resumed runs
        for k, v in opts.items():
            key = {"masters": "masters", "master": "master"}.get(k, k)
            if other.get(key) != v:
                sys.exit(f"{exp}: manifests disagree on {k}")
    opts["manifests"] = [x["_file"] for x in mans]
    return opts


def build(exp: str) -> dict:
    folder, driver = EXPERIMENTS[exp]
    mans = manifests(folder)
    cells = merged_cells(mans)
    opts = run_options(exp, mans)
    by_rep: dict[str, dict[int, list[dict]]] = {}
    for name in cells:
        groups: dict[int, list[dict]] = {}
        for r in saved_rows(folder, name):
            groups.setdefault(int(r["rep"]), []).append(r)
        for rep, rows in groups.items():
            keys = [r.get("adj", "") for r in rows]
            if len(keys) != len(set(keys)):
                sys.exit(f"{exp} {name} rep {rep}: duplicate saved rows")
            if any(r["status"] != "ok" for r in rows):
                groups[rep] = []
        by_rep[name] = {k: v for k, v in groups.items() if v}

    picks: list[tuple[str, str, int]] = []  # (tier, name, rep)
    order = sorted(cells, key=lambda nm: cells[nm]["idx"])
    for nm in order:
        n = cells[nm]["n"]
        reps = sorted(by_rep[nm])
        if n <= 500:
            for rep in sorted({reps[0], reps[-1]}):
                picks.append(("fast", nm, rep))
        elif n >= 4000:
            picks.append(("slow", nm, reps[0]))
    mid = [nm for nm in order if cells[nm]["n"] == 2000]
    cost = sum(float(by_rep[nm][min(by_rep[nm])][0]["cpu_s"]) * len(by_rep[nm][min(by_rep[nm])]) for nm in mid)
    k = max(1, math.ceil(cost / FAST_BUDGET)) if mid else 1
    for i, nm in enumerate(mid):
        picks.append(("fast" if i % k == 0 else "slow", nm, min(by_rep[nm])))
    if exp == "EXP7":
        seen = set()
        for nm in order:
            c = cells[nm]
            key = (c["family"], c["n"])
            if c["family"] == "growth" or c["alt"] != "null" or key in seen:
                continue
            seen.add(key)
            rep = int(opts["exact_reps"][str(c["n"])])
            if rep in by_rep[nm]:
                picks.append(("fast" if c["n"] <= 2000 else "slow", nm, rep))

    rows = []
    seen_ids = set()
    for tier, nm, rep in picks:
        rid = f"{nm}-r{rep}"
        if rid in seen_ids:
            continue
        seen_ids.add(rid)
        rows.append({"id": rid, "tier": tier, "cell": cells[nm], "rep": rep, "saved": by_rep[nm][rep]})
    cpu = {t: round(sum(float(s["cpu_s"] or 0) for r in rows if r["tier"] == t for s in r["saved"]), 1)
           for t in ("fast", "slow")}
    return {"experiment": exp, "folder": folder, "driver": driver, "options": opts,
            "saved_cpu_s": cpu, "rows": rows}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    for exp in (sys.argv[1:] or EXPERIMENTS):
        fx = build(exp)
        path = os.path.join(OUT, f"{exp}.json")
        with open(path, "w") as f:
            json.dump(fx, f, indent=0, sort_keys=True)
            f.write("\n")
        tiers = {t: sum(r["tier"] == t for r in fx["rows"]) for t in ("fast", "slow")}
        print(f"{exp}: {tiers['fast']} fast + {tiers['slow']} slow rows, saved CPU-s {fx['saved_cpu_s']}")


if __name__ == "__main__":
    main()
