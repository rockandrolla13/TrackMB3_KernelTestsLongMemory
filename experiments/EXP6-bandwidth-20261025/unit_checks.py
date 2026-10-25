"""EXP6 unit checks. Usage: python -I unit_checks.py [--regress]  -> prints, exit 1 on failure.

1. bw_core.ind_pair equals EXP3 power_core.ind_pair on the same seed (all gens, alts).
2. kind_bw at kappa = 1 equals EXP3 power_core.kind_ext field by field on the same seed.
3. kappa scaling: sx_k2 = 2 sx_med exactly; kappa = 1 bandwidth equals the median heuristic.
4. --regress: regress/ rows (EXP3 master, EXP3 cell indices) equal the EXP3 main rows at kappa = 1.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import csv  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
EXP3 = os.path.join(ROOT, "experiments", "EXP3-power-20261025")
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.append(EXP3)
import numpy as np  # noqa: E402

import bw_core as K  # noqa: E402
import power_core as P3  # noqa: E402  (EXP3, imports its own gens/mbcore copies)
import run_exp6 as R  # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {name} {detail}")
    if not ok:
        fails.append(name)


def rng(s):
    return np.random.default_rng(np.random.SeedSequence(entropy=777, spawn_key=(s,)))


# 1. generators
for gen in ("G0", "SEAS-V"):
    for alt, rho in (("null", 0.0), ("lin", 0.065), ("sq", 0.075)):
        if gen == "SEAS-V" and alt == "sq":
            continue
        a = K.ind_pair(gen, alt, 0.5, 400, rho, rng(1))
        b = P3.ind_pair(gen, alt, 0.5, 400, rho, rng(1))
        d = max(np.max(np.abs(a[0] - b[0])), np.max(np.abs(a[1] - b[1])))
        check(f"ind_pair {gen} {alt}", d == 0.0, f"maxdiff={d}")

# 2. kappa = 1 equals kind_ext
for n, gen, alt, rho in ((300, "G0", "lin", 0.1), (300, "SEAS-V", "null", 0.0), (2000, "SEAS-V", "lin", 0.065)):
    r1, r2 = rng(2), rng(2)
    x, y = K.ind_pair(gen, alt, 0.5, n, rho, r1)
    P3.ind_pair(gen, alt, 0.5, n, rho, r2)
    a = K.kind_bw(x, y, r1, kappas=(0.5, 1.0, 2.0))
    b = P3.kind_ext(x, y, r2)
    worst = 0.0
    for key, v in b.items():
        ka = key if key in a else f"{key}_k1"
        worst = max(worst, abs(float(a[ka]) - float(v)) if np.isfinite(v) else 0.0)
    check(f"kappa=1 vs kind_ext n={n} {gen} {alt}", worst == 0.0, f"maxdiff={worst}")
    check(f"bandwidth scaling n={n}", a["sx_k2"] == 2 * a["sx_med"] and a["sy_k0p5"] == 0.5 * a["sy_med"]
          and a["sx_k1"] == a["sx_med"])

# 3b. oracle-only path equals the full path at the same kappa
x, y = K.ind_pair("SEAS-V", "lin", 0.5, 500, 0.1, rng(3))
full = K.kind_bw(x, y, rng(4), kappas=(2.0,))
cheap = K.kind_bw(x, y, rng(4), kappas=(1.0,), oracle_kappas=(2.0,))
d = abs(full["hsic_stat_k2"] - cheap["hsic_stat_k2"]) / abs(full["hsic_stat_k2"])
check("oracle-only kappa=2 equals full path", d < 1e-10, f"rel diff={d:.1e}")

# 4. regression against saved EXP3 main rows
if "--regress" in sys.argv:
    for c in R.cells_for("regress", (2000,)):
        mine = {int(r["rep"]): r for r in csv.DictReader(open(os.path.join(HERE, "regress", c["name"] + ".csv")))}
        e3 = [r for r in csv.DictReader(open(os.path.join(EXP3, "main", c["name"] + ".csv")))]
        e3 = {int(r["rep"]): r for r in e3 if int(r["rep"]) in mine}
        worst, nk = 0.0, 0
        for rep, r3 in e3.items():
            for key, v in r3.items():
                if key in ("rep", "status", "error", "cpu_s"):
                    continue
                ka = key if key in mine[rep] else f"{key}_k1"
                a, b = float(mine[rep][ka]), float(v)
                if np.isfinite(b):
                    worst = max(worst, abs(a - b) / max(1.0, abs(b)))
                nk += 1
        check(f"regress {c['name']} ({len(e3)} reps, {nk} values)", len(e3) == len(mine) and worst < 1e-12,
              f"max rel diff={worst:.2e} (CSV round-trip)")

print("ALL PASS" if not fails else f"{len(fails)} FAIL: {fails}")
sys.exit(1 if fails else 0)
