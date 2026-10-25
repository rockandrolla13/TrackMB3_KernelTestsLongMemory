"""EXP3 unit checks. Run: python -I unit_checks.py > unit_checks.txt"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import math  # noqa: E402
import sys  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import gens  # noqa: E402
import mbcore as mb  # noqa: E402
import power_core as P  # noqa: E402


def rng(k):
    return np.random.default_rng(np.random.SeedSequence(entropy=42, spawn_key=(k,)))


def main():
    ok = True
    # 1. kind_ext reproduces mbcore.kind_tests and lind_tests on the same seed
    for n in (500, 2000):
        x, y = gens.pair("G0", 0.7, n, rng(1))
        a = mb.kind_tests(x, y, rng(2), B=500)
        a.update(mb.lind_tests(x, y))
        b = P.kind_ext(x, y, rng(2), B=500)
        diff = max(abs(float(a[k]) - float(b[k])) for k in a)
        print(f"[1] n={n} max |kind_tests - kind_ext| over {len(a)} shared fields = {diff:.2e}")
        ok &= diff < 1e-12
        # raw p-rule and ratio rule agree up to the counting convention
        print(f"    p_hsic_wb={b['p_hsic_wb']:.4f} r_hsic_wb={b['r_hsic_wb']:.4f}  "
              f"rej_hsic_sub={b['rej_hsic_sub']} r_hsic_sub={b['r_hsic_sub']:.4f}")
        ok &= (b["rej_hsic_sub"] == 1) == (b["r_hsic_sub"] > 1)
    # 2. alt generators: rho = 0 'lin' reproduces the null pair; seasonal null same draw order
    for g in ("G0", "SEAS-V"):
        x0, y0 = P.ind_pair(g, "null", 0.7, 500, 0.0, rng(3))
        x1, y1 = P.ind_pair(g, "lin", 0.7, 500, 0.0, rng(3))
        d = max(np.abs(x0 - x1).max(), np.abs(y0 - y1).max())
        print(f"[2] {g}: max |null - lin(rho=0)| = {d:.2e}")
        ok &= d < 1e-12
    # 3. population checks of the alternatives (large n, H = 0.5)
    r = rng(4)
    xs, ys, zs = [], [], []
    for _ in range(50):
        x, y = P.ind_pair("G0", "lin", 0.5, 4000, 0.3, r)
        xs.append(np.corrcoef(x, y)[0, 1])
        x, y = P.ind_pair("G0", "sq", 0.5, 4000, 0.3, r)
        ys.append(np.corrcoef(x, y)[0, 1])
        zs.append(np.corrcoef(x * x, y)[0, 1])
    print(f"[3] lin rho=0.3: mean corr(X,Y) = {np.mean(xs):.4f} (target 0.3)")
    print(f"    sq  rho=0.3: mean corr(X,Y) = {np.mean(ys):.4f} (target 0), corr(X^2,Y) = {np.mean(zs):.4f} (target 0.3)")
    ok &= abs(np.mean(xs) - 0.3) < 0.01 and abs(np.mean(ys)) < 0.01 and abs(np.mean(zs) - 0.3) < 0.01
    x, y = P.ind_pair("SEAS-V", "lin", 0.5, 4000, 0.3, rng(5))
    print(f"[3] SEAS-V lin rho=0.3 one path: corr = {np.corrcoef(x, y)[0, 1]:.4f}")
    # 4. MMD shift alternative adds delta to the second half only
    x0 = P.mmd_series("G0", "null", 0.7, 500, 0.0, rng(6))
    x1 = P.mmd_series("G0", "shift", 0.7, 500, 0.5, rng(6))
    dd = x1 - x0
    print(f"[4] shift: first half max |d| = {np.abs(dd[:250]).max():.2e}, second half d in [{dd[250:].min():.3f}, {dd[250:].max():.3f}]")
    ok &= np.abs(dd[:250]).max() == 0 and np.allclose(dd[250:], 0.5)
    m = P.mmd_ext(x1, rng(7))
    print(f"[4] mmd_ext on shifted series: " + ", ".join(f"{k}={v:.4g}" for k, v in m.items()))
    ok &= (m["t_n"] > m["q95_sqrt"]) == (m["r_sub_sqrt"] > 1)
    print("ALL OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
