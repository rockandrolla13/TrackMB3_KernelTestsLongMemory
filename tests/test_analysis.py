"""Analysis helpers (plan steps 3.3 and 5.1): per-cell bootstrap streams make each cell's intervals independent."""
from __future__ import annotations

import numpy as np

from kmem import analysis as A


def _intervals(cells: dict, rng_for) -> dict:
    out = {}
    for idx, (nv, av) in cells.items():
        rng = rng_for(idx)
        bp = A.boot_power(nv, av, A.boot_idx(nv.size, rng, 400), A.boot_idx(av.size, rng, 400))
        out[idx] = (float(np.quantile(bp, 0.025)), float(np.quantile(bp, 0.975)))
    return out


def _cells() -> dict:
    g = np.random.default_rng(7)
    return {i: (g.standard_normal(300), g.standard_normal(300) + 0.2 * i) for i in range(4)}


def test_per_cell_streams_do_not_depend_on_other_cells() -> None:
    cells = _cells()
    full = _intervals(cells, A.cell_rng)
    del cells[1]
    assert _intervals(cells, A.cell_rng) == {k: v for k, v in full.items() if k != 1}


def test_the_shared_stream_does_depend_on_other_cells() -> None:
    """Why D12 adds cell_rng: with one shared stream, dropping a cell moves the later intervals."""
    cells = _cells()
    shared = np.random.default_rng(A.SEED)
    full = _intervals(cells, lambda i: shared)
    del cells[1]
    shared = np.random.default_rng(A.SEED)
    assert _intervals(cells, lambda i: shared)[2] != full[2]


def test_wilson_and_size_class() -> None:
    lo, hi = A.wilson(50, 1000)
    assert lo < 0.05 < hi and A.size_class(50, 1000) == "valid"
    assert A.size_class(200, 1000) == "over" and A.wilson(0, 0) != A.wilson(0, 0)  # nan pair
