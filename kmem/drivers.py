"""Load an experiment's driver module, which supplies run_chunk (its row function) and job_for (its jobs)."""
from __future__ import annotations

import importlib.util
import os
import sys

from kmem.spec import REPO, ExperimentSpec


def load(sp: ExperimentSpec):
    """Import experiments/<folder>/<driver> under a unique name, registered before it runs, so the fork pool
    can find run_chunk by name in its workers."""
    name = f"kmem_driver_{sp.id}"
    if name in sys.modules:
        return sys.modules[name]
    path = os.path.join(REPO, sp.folder, sp.driver)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        del sys.modules[name]
        raise
    for attr in ("run_chunk", "job_for"):
        if not hasattr(mod, attr):
            raise AttributeError(f"{path} has no {attr}()")
    return mod


def jobs(sp: ExperimentSpec, only: str | None = None, reps: int | None = None) -> list:
    """The driver's jobs for the spec's cells. only: comma list of name substrings. reps: override every target."""
    drv = load(sp)
    cells = [c for c in sp.cells if not only or any(p in c.name for p in only.split(","))]
    return [drv.job_for(c.driver_cell(), c.target if reps is None else reps, sp.master_for(c), sp.options)
            for c in cells]
