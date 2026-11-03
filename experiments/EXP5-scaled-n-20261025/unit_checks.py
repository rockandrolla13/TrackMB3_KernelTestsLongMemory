"""EXP5 unit checks. Output: unit_checks.txt.

1. sd_half closed form equals the exact fGn-autocovariance value.
2. exp5_core.mmd_ext equals EXP3 power_core.mmd_ext on the same seed for every shared field.
3. exp5_core.ind_rep equals the EXP1 E3 path (mbcore.gen_series G4 + kind_tests + lind_tests).
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import importlib.util  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))  # kmem (D7: not installed)
import numpy as np  # noqa: E402

import exp5_core as P  # noqa: E402

EXP3 = os.path.join(HERE, "..", "EXP3-power-20261025")
EXP1 = os.path.join(HERE, "..", "EXP1-20261016")


def load_from(folder, mod):
    sys.path.insert(0, folder)
    try:
        for m in ("mbcore", "mmdcal", "gens", "common", mod):
            sys.modules.pop(m, None)
        spec = importlib.util.spec_from_file_location(mod, os.path.join(folder, mod + ".py"))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m
    finally:
        sys.path.remove(folder)


def main():
    L = []
    for h in (0.5, 0.7, 0.8):
        for n in (500, 2000):
            a, b = P.sd_half(h, n), P.sd_half_exact(h, n)
            L.append(f"sd_half H={h} n={n}: closed {a:.8f} exact {b:.8f} diff {a - b:.2e}")
            assert abs(a - b) < 1e-9
    pc = load_from(EXP3, "power_core")
    worst = 0.0
    for h, n in ((0.5, 500), (0.8, 500), (0.7, 2000)):
        for r in range(3):
            ss = np.random.SeedSequence(entropy=999, spawn_key=(1, r))
            x1 = P.mmd_series(h, n, 0.1, np.random.default_rng(ss))
            rng1 = np.random.default_rng(ss)
            P.mmd_series(h, n, 0.1, rng1)
            o1 = P.mmd_ext(x1, rng1)
            rng2 = np.random.default_rng(ss)
            x2 = pc.mmd_series("G0", "shift", h, n, 0.1, rng2)
            o2 = pc.mmd_ext(x2, rng2)
            assert np.array_equal(x1, x2)
            for k, v in o2.items():
                d = abs(float(o1[k]) - float(v))
                worst = max(worst, d)
    L.append(f"mmd_ext vs EXP3 power_core.mmd_ext: 9 reps, worst abs diff over shared fields {worst:.3e}")
    assert worst == 0.0
    mb1 = load_from(EXP1, "mbcore")
    worst = 0.0
    for h in (0.5, 0.8):
        ss = np.random.SeedSequence(entropy=999, spawn_key=(2, 0))
        o1 = P.ind_rep(h, 400, np.random.default_rng(ss))
        rng = np.random.default_rng(ss)
        x, y = mb1.gen_series({"proc": "G4", "HX": h, "HY": h, "lam": None}, 400, rng)
        o2 = mb1.kind_tests(x, y, rng, B=500)
        o2.update(mb1.lind_tests(x, y))
        for k, v in o2.items():
            d = abs(float(o1[k]) - float(v))
            if np.isfinite(d):
                worst = max(worst, d)
    L.append(f"ind_rep vs EXP1 mbcore E3 path: 2 reps n=400, worst abs diff {worst:.3e}")
    assert worst == 0.0
    L.append("ALL PASS")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "unit_checks.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
