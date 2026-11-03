"""Build experiments/EXPn/spec.json for EXP2-EXP7 from the saved main-run manifests (plan step 3.1).

Seed indices, cell names and cell fields come from the manifests, never from code. Targets are the number of
distinct replications in the saved main CSVs. Chunk sizes and worker caps come from each driver's constants.
Run once from the repo root: python -I tools/build_specs.py [EXP2 ...]
"""
from __future__ import annotations

import csv
import glob
import importlib.util
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from kmem.spec import EXPERIMENTS, ExperimentSpec, digest  # noqa: E402

DRIVERS = {"EXP2": "run_exp2.py", "EXP3": "run_exp3.py", "EXP4": "run_exp4.py", "EXP5": "run_exp5.py",
           "EXP6": "run_exp6.py", "EXP7": "run_exp7.py"}
BOOKKEEPING = {"already_done", "already_done_pairs", "reps", "seed"}
INPUTS = {  # (path, used by): files of other experiments that an analysis reads
    "EXP3": [("experiments/EXP2-20261023/main", "analyze.py: null regression against EXP2")],
    "EXP4": [("experiments/EXP2-20261023/summary_main.csv", "analyze_exp4.py: unadjusted A0-EXP2 rows")],
    "EXP5": [("experiments/EXP3-power-20261025/main", "analyze_exp5.py: null regression against EXP3"),
             ("experiments/EXP1-20261016/e3_summary.csv", "analyze_exp5.py: E3 baseline (Part 2)"),
             ("experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear/cells/cell_summary.csv",
              "analyze_exp5.py: KL baseline")],
}
CONSUMERS = {
    "EXP7": [("investigations/memory-kernel-2026-11-01/examples-report/build_data.py",
              ["summary_main.csv", "results.json (hand-made, see RESULTS.md)", "main/*.csv"])],
}
OPTION_KEYS = ("B_wild", "kappas", "oracle_kappas", "exact_reps", "gamma", "scale_c", "shape_a", "tol_chol",
               "tol_pair", "period_true", "numpy", "python")


def load_driver(exp: str):
    path = os.path.join(REPO, EXPERIMENTS[exp], DRIVERS[exp])
    spec = importlib.util.spec_from_file_location(f"specdriver_{exp}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build(exp: str) -> ExperimentSpec:
    folder = EXPERIMENTS[exp]
    paths = sorted(glob.glob(os.path.join(REPO, folder, "main", "manifest_*.json")))
    mans = [json.load(open(p)) for p in paths]
    cells: dict[str, dict] = {}
    for m in mans:
        for c in m["cells"]:
            c = {k: v for k, v in c.items() if k not in BOOKKEEPING}
            old = cells.setdefault(c["name"], c)
            if old["idx"] != c["idx"]:
                raise SystemExit(f"{exp} {c['name']}: manifests disagree on idx")
            if "adjustments" in c:
                old["adjustments"] = sorted(set(old["adjustments"]) | set(c["adjustments"]))
    m0 = mans[0]
    if exp == "EXP2":
        masters = {k: v for k, v in m0["masters"].items() if k.endswith("_main")}
        master_key = "{family}_main"
    else:
        masters, master_key = {"main": m0["master"]}, "main"
    options = {k: m0[k] for k in OPTION_KEYS if k in m0}
    drv = load_driver(exp)
    chunk_src = drv.CHUNK
    families = sorted({c.get("family", "ind") for c in cells.values()})
    if all(isinstance(v, dict) for v in chunk_src.values()):
        chunk = {f: {str(n): int(k) for n, k in chunk_src[f].items()} for f in families}
    else:
        chunk = {f: {str(n): int(k) for n, k in chunk_src.items()} for f in families}
    out_cells = []
    for name, c in sorted(cells.items(), key=lambda kv: kv[1]["idx"]):
        with open(os.path.join(REPO, folder, "main", name + ".csv"), newline="") as f:
            reps = {int(r["rep"]) for r in csv.DictReader(f)}
        if reps != set(range(len(reps))):
            raise SystemExit(f"{exp} {name}: saved reps are not 0..k-1")
        out_cells.append({"name": name, "seed_index": c["idx"], "family": c.get("family", "ind"), "n": c["n"],
                          "target": len(reps), "fields": c})
    prereq = []
    if exp == "EXP7":
        rel = os.path.join(folder, "sn_table.json")
        prereq.append({"path": rel, "sha256": digest(rel),
                       "made_by": "python -I sn_tables.py in the EXP7 folder (master 2026110174, about 1 CPU-hour)"})
    return ExperimentSpec.model_validate({
        "id": exp, "folder": folder, "driver": DRIVERS[exp], "masters": masters, "master_key": master_key,
        "seed_rule": m0["seed_rule"], "options": options, "writer": "fmt12g" if exp == "EXP7" else "str",
        "chunk": chunk, "max_workers": int(getattr(drv, "MAX_WORKERS", 16)), "cells": out_cells,
        "prerequisites": prereq,
        "inputs": [{"path": p, "sha256": digest(p), "used_by": u} for p, u in INPUTS.get(exp, [])],
        "consumers": [{"path": p, "reads": r} for p, r in CONSUMERS.get(exp, [])],
        "source_manifests": [os.path.relpath(p, REPO) for p in paths]})


def main() -> None:
    for exp in sys.argv[1:] or EXPERIMENTS:
        spec = build(exp)
        path = os.path.join(REPO, spec.folder, "spec.json")
        with open(path, "w") as f:
            json.dump(spec.model_dump(mode="json"), f, indent=1)
            f.write("\n")
        print(f"{exp}: {len(spec.cells)} cells, {spec.total_rows} saved rows -> {os.path.relpath(path, REPO)}")


if __name__ == "__main__":
    main()
