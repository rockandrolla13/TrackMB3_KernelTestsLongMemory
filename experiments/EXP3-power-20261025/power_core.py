"""EXP3 rule bundles that return, for every rule, its raw decision input AND a continuous ratio.

The ratio of a rule is R = T / c_ref, where T is the rule's statistic and c_ref is the rule's own
95% reference quantile (wild draws, circular shifts, or sampling windows). The raw rule rejects
iff p <= 0.05 (wild, shift) or T > q95 (subsampling), exactly as in EXP2. Size adjustment
replaces the cut-off 1 on R by the 95% quantile of R over null replications (DESIGN.md).

The statistics and references are the EXP2 code (mbcore.kind_tests, mmdcal.rep_e2), with the
extra quantiles read off the same reference draws. unit_checks.py verifies the p-values match.
Exploratory simulation; not evidence for any accepted claim.
"""
from __future__ import annotations

# The pieces now live in kmem (refactor steps 2.2 and 2.3, 2026-11-03), with the arithmetic and draw order
# unchanged; the saved EXP3 rows still reproduce exactly (tests/repro). The EXP3 names stay as thin wrappers.
import numpy as np

from kmem.calibration.common import ALPHA, pboot, q95  # noqa: F401  (re-exported)
from kmem.generators.independence import pair
from kmem.generators.seasonal import seas_scale  # noqa: F401  (re-exported)
from kmem.generators.split import split_series
from kmem.hsic import hsic_dcov_ext as kind_ext  # noqa: F401  (EXP3 name)
from kmem.legacy import gens
from kmem.mmd import split_mmd_ext as mmd_ext  # noqa: F401  (EXP3 name)

SEAS_PERIOD = gens.PERIOD


def ind_pair(gen: str, alt: str, H: float, n: int, rho: float, rng: np.random.Generator):
    """Independence-arm pair. alt in {null, lin, sq}; gen in {G0, SEAS-V} (kmem.generators.independence.pair).

    null: X, Y independent (gens.pair, the EXP2 code path).
    lin:  latent eta = rho xi + sqrt(1 - rho^2) w, w independent fGn(H).
    sq:   eta = rho (xi^2 - 1)/sqrt(2) + sqrt(1 - rho^2) w  (uncorrelated, rank-two dependence).
    SEAS-V multiplies each latent series by its own seasonal scale with an independent phase.
    """
    return pair(gen, alt, H, H, n, rho, rng)


def mmd_series(gen: str, alt: str, H: float, n: int, delta: float, rng: np.random.Generator):
    """Split-MMD series. alt 'null' is the EXP2 generator; 'shift' adds delta to the second half."""
    if alt not in ("null", "shift"):
        raise ValueError(alt)
    return split_series(gen, alt, H, n, delta, rng)
