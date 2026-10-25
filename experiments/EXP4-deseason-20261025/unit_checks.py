"""EXP4 unit checks: after the correct known-period adjustment, seasonal series equal the adjusted fGn core.

SEAS-M, SEAS-M-shared: A1(xi + s) = A1(xi) because s is constant within every bin.
SEAS-V: A2((1 + b s) xi) = A2(xi) because 1 + b s > 0 is constant within every bin.
So those adjusted rows equal adjusted G0 in law. Output: unit_checks.txt.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import deseason as D  # noqa: E402
import gens  # noqa: E402
import mbcore as mb  # noqa: E402

rng = np.random.default_rng(2026102501)
out = []
for n in (500, 2000, 1999):
    for H in (0.5, 0.8):
        worst = {"M-A1": 0.0, "M-A2": 0.0, "V-A2": 0.0, "V-A1": 0.0}
        for _ in range(50):
            xi = mb.fgn(H, n, rng)
            s = gens.season(n, int(rng.integers(gens.PERIOD)))
            xm, xv = xi + gens.SEAS_A * s, (1 + gens.SEAS_B * s) * xi
            worst["M-A1"] = max(worst["M-A1"], np.abs(D.adjust(xm, "A1")[0] - D.adjust(xi, "A1")[0]).max())
            worst["M-A2"] = max(worst["M-A2"], np.abs(D.adjust(xm, "A2")[0] - D.adjust(xi, "A2")[0]).max())
            worst["V-A2"] = max(worst["V-A2"], np.abs(D.adjust(xv, "A2")[0] - D.adjust(xi, "A2")[0]).max())
            worst["V-A1"] = max(worst["V-A1"], np.abs(D.adjust(xv, "A1")[0] - D.adjust(xi, "A1")[0]).max())
        out.append(f"n={n} H={H} max|adjusted seasonal - adjusted core| over 50 series: "
                   + ", ".join(f"{k} {v:.2e}" for k, v in worst.items()))
out.append("Expected: M-A1, M-A2, V-A2 at rounding level (1e-15). V-A1 order 1 (A1 leaves the scale).")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "unit_checks.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
