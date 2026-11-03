"""Memory kernels by name. get("herm1") returns the Hermite projection of order 1, and so on.

The ten EXP7 kernels are registered here under their EXP7 names, so saved CSV columns keep their names.
A new kernel is one new module plus one register(...) call; see kmem/kernels/base.py.
"""
from __future__ import annotations

from kmem.kernels.base import Kernel, get, names, register  # noqa: F401
from kmem.kernels.distance import DistanceIncrement
from kmem.kernels.feature_maps import FracOracle, HighPass, Hermite, Increment, Raw
from kmem.kernels.own_series import DATADIFF, SYM
from kmem.legacy import gens

INC_LONG = gens.PERIOD  # inc50: the lag equals the seasonal period of the EXP2 generators
HP_FRAC = 0.1  # hp10: the lowest 10% of DCT frequencies

for _k in (Raw(), Increment(1), Increment(INC_LONG), Hermite(1), Hermite(2), SYM, DATADIFF, HighPass(HP_FRAC, "hp10"),
           DistanceIncrement(1), FracOracle()):
    register(_k)
