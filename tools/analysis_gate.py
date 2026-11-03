"""G-analysis: rerun every saved analysis of EXP2 to EXP7 on the saved raw CSVs and hash what it writes.

The analyses run in a scratch copy of the current working tree (experiments EXP2 to EXP7, kmem, and the two
outside inputs EXP5 reads), so nothing in the repo is ever written. Each job's outputs are deleted in the
copy first, so every hashed file was really regenerated. BLAS is pinned to one thread.

  python -I tools/analysis_gate.py --write-baseline   run twice; record regenerated hashes, which outputs
                                                      also match the saved file, and which are unstable
  python -I tools/analysis_gate.py --check            run once; fail if any stable output's hash moved

The baseline lives in tests/analysis_baseline.json. Unstable outputs (two baseline runs disagree) are
listed there and only checked for existence.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASELINE = os.path.join(REPO, "tests", "analysis_baseline.json")
EXPS = {
    "EXP2": "experiments/EXP2-20261023",
    "EXP3": "experiments/EXP3-power-20261025",
    "EXP4": "experiments/EXP4-deseason-20261025",
    "EXP5": "experiments/EXP5-scaled-n-20261025",
    "EXP6": "experiments/EXP6-bandwidth-20261025",
    "EXP7": "experiments/EXP7-memory-kernel-20261101",
}
OUTSIDE_INPUTS = ["experiments/EXP1-20261016/e3_summary.csv",
                  "experiments/JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear/cells/cell_summary.csv"]


def _tag(t: str, names: list[str]) -> list[str]:
    return [n.format(t=t) for n in names]


# (experiment, command after the interpreter, use -I, outputs relative to the experiment folder)
A2 = ["summary_{t}.csv", "tables_{t}.md", "verdicts_{t}.json"]
A3 = ["null_regression_{t}.csv", "summary_{t}.csv", "results_{t}.json", "tables_{t}.md"]
A5 = ["summary_{t}.csv", "summary_ind_{t}.csv", "null_regression_{t}.csv", "results_{t}.json", "tables_{t}.md"]
A6 = ["summary_{t}.csv", "results_{t}.json", "tables_{t}.md"]
A7 = ["summary_{t}.csv", "results_{t}.json", "tables_{t}.md"]
CHAINS = [
    [("EXP2", ["analyze.py", "--in", "main", "--tag", "main"], True, _tag("main", A2)),
     ("EXP2", ["analyze.py", "--in", "pilot", "--tag", "pilot"], True, _tag("pilot", A2)),
     ("EXP2", ["multiple_testing/mt_analysis.py"], True,
      ["multiple_testing/mt_results.csv", "multiple_testing/mt_claims.csv", "multiple_testing/mt_summary.json"])],
    [("EXP3", ["analyze.py", "--in", "main", "--tag", "main"], True, _tag("main", A3))],
    [("EXP4", ["analyze_exp4.py", "--in", "main", "--tag", "main"], True, _tag("main", A2)),
     ("EXP4", ["analyze_exp4.py", "--in", "pilot", "--tag", "pilot"], True, _tag("pilot", A2))],
    [("EXP5", ["analyze_exp5.py", "--in", "main", "--tag", "main"], True, _tag("main", A5)),
     ("EXP5", ["analyze_exp5.py", "--in", "pilot", "--tag", "pilot"], True, _tag("pilot", A5))],
    [("EXP6", ["analyze_exp6.py", "--in", "main", "--tag", "main"], True, _tag("main", A6)),
     ("EXP6", ["analyze_exp6.py", "--in", "pilot", "--tag", "pilot"], True, _tag("pilot", A6))],
    [("EXP7", ["analyze_exp7.py", "--in", "main", "--tag", "main"], True, _tag("main", A7)),
     ("EXP7", ["analyze_exp7.py", "--in", "pilot", "--extra", "pilot_extras", "--tag", "pilot"], True, _tag("pilot", A7)),
     ("EXP7", ["analyze_exp7.py", "--tune"], True, ["tune_report.txt"]),
     ("EXP7", ["results_tables.py", "main"], True, ["results_tables_main.md"]),
     ("EXP7", ["results_tables.py", "pilot"], True, ["results_tables_pilot.md"]),
     ("EXP7", ["plot_growth.py", "main"], False, ["growth_q95_main.png"]),
     ("EXP7", ["plot_growth.py", "pilot"], False, ["growth_q95_pilot.png"]),
     ("EXP7", ["review_checks.py"], True,
      ["binomial_valid.csv", "trend_main.csv", "growth_bandwidth.csv", "bandwidth_check.csv"])],
]
FIRST = CHAINS[0]  # EXP2 runs before the others: EXP4 reads its summary_main.csv


def sha(path: str) -> str | None:
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def build_tree(work: str) -> None:
    ign = shutil.ignore_patterns("__pycache__", "*.pyc")
    for rel in EXPS.values():
        shutil.copytree(os.path.join(REPO, rel), os.path.join(work, rel), ignore=ign, symlinks=False)
    shutil.copytree(os.path.join(REPO, "kmem"), os.path.join(work, "kmem"), ignore=ign)
    for rel in OUTSIDE_INPUTS:
        os.makedirs(os.path.dirname(os.path.join(work, rel)), exist_ok=True)
        shutil.copy2(os.path.join(REPO, rel), os.path.join(work, rel))


def run_chain(work: str, chain: list) -> list[dict]:
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", MPLBACKEND="Agg")
    log = []
    for exp, cmd, isolated, outs in chain:
        here = os.path.join(work, EXPS[exp])
        for o in outs:
            if os.path.exists(os.path.join(here, o)):
                os.remove(os.path.join(here, o))
        t0 = time.time()
        argv = [sys.executable] + (["-I"] if isolated else []) + cmd
        p = subprocess.run(argv, cwd=here, env=env, capture_output=True, text=True)
        log.append({"exp": exp, "cmd": " ".join(cmd), "code": p.returncode, "seconds": round(time.time() - t0, 1),
                    "stderr": p.stderr[-1500:] if p.returncode else ""})
    return log


def run_all(keep: bool = False) -> tuple[dict, list]:
    work = tempfile.mkdtemp(prefix="kmem-analysis-")
    try:
        build_tree(work)
        log = run_chain(work, FIRST)
        with cf.ThreadPoolExecutor(len(CHAINS) - 1) as ex:
            for part in ex.map(lambda ch: run_chain(work, ch), CHAINS[1:]):
                log += part
        hashes = {}
        for chain in CHAINS:
            for exp, _, _, outs in chain:
                for o in outs:
                    rel = os.path.join(EXPS[exp], o)
                    hashes[rel] = {"regen": sha(os.path.join(work, rel)), "saved": sha(os.path.join(REPO, rel))}
        return hashes, log
    finally:
        if keep:
            print(f"kept scratch tree {work}", file=sys.stderr)
        else:
            shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write-baseline", action="store_true")
    g.add_argument("--check", action="store_true")
    ap.add_argument("--keep", action="store_true", help="keep the scratch tree for inspection")
    a = ap.parse_args()
    first, log = run_all(a.keep)
    failed = [x for x in log if x["code"]]
    for x in log:
        print(f"{x['exp']:5} {x['cmd']:60} exit {x['code']}  {x['seconds']:7.1f}s", flush=True)
    if failed:
        for x in failed:
            print(f"--- {x['exp']} {x['cmd']}\n{x['stderr']}", file=sys.stderr)
        sys.exit("analysis job failed")
    if a.write_baseline:
        second, _ = run_all()
        out = {}
        for rel, h in sorted(first.items()):
            stable = h["regen"] == second[rel]["regen"]
            out[rel] = {"regen": h["regen"], "saved": h["saved"], "stable": stable,
                        "matches_saved": stable and h["regen"] == h["saved"]}
        with open(BASELINE, "w") as f:
            json.dump(out, f, indent=1, sort_keys=True)
            f.write("\n")
        n = len(out)
        print(f"baseline: {n} outputs, {sum(v['matches_saved'] for v in out.values())} match the saved file, "
              f"{sum(not v['stable'] for v in out.values())} unstable")
        return
    base = json.load(open(BASELINE))
    bad = []
    for rel, b in base.items():
        now = first.get(rel, {}).get("regen")
        if now is None:
            bad.append(f"{rel}: not written")
        elif b["stable"] and now != b["regen"]:
            bad.append(f"{rel}: hash changed")
    extra = sorted(set(first) - set(base))
    if extra:
        bad.append(f"outputs not in the baseline: {extra}")
    if bad:
        print("\n".join(bad))
        sys.exit(f"G-analysis FAILED: {len(bad)} problems")
    print(f"G-analysis passed: {len(base)} outputs, {sum(v['stable'] for v in base.values())} checked by hash")


if __name__ == "__main__":
    main()
