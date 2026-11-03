"""Generator equality cases for tests/test_generators.py (plan step 2.2).

Every generator copy in EXP3 to EXP7 is called on the domain its experiment used. The sha256 of each output
array was frozen in tests/repro/generators.json from the pre-move code (git e40180d), so moving the
generators into kmem.generators must give the same bytes.
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SOURCES = {
    "EXP3": "experiments/EXP3-power-20261025/power_core.py",
    "EXP4": "experiments/EXP4-deseason-20261025/gens_exp4.py",
    "EXP5": "experiments/EXP5-scaled-n-20261025/exp5_core.py",
    "EXP6": "experiments/EXP6-bandwidth-20261025/bw_core.py",
    "EXP7": "experiments/EXP7-memory-kernel-20261101/exp7_core.py",
}
NS, HS, SEEDS = (64, 501), (0.5, 0.8), (0, 1)
PROCESS = ("HT-S", "HT-L", "NL-SV", "NL-TAR", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow")


def load(exp: str):
    """Load one experiment module by path; its folder goes first on sys.path for its bare-name imports."""
    rel = SOURCES[exp]
    folder = os.path.join(REPO, os.path.dirname(rel))
    sys.path.insert(0, folder)
    if REPO not in sys.path:
        sys.path.insert(1, REPO)
    try:
        spec = importlib.util.spec_from_file_location(f"gencase_{exp}", os.path.join(REPO, rel))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(folder)
    return mod


def cases() -> list[tuple[str, str, tuple]]:
    """(experiment, function name, arguments without the rng)."""
    out = []
    for n in NS:
        for h in HS:
            for gen in ("G0", "SEAS-V"):
                for alt, par in (("null", 0.0), ("lin", 0.065), ("sq", 0.075)):
                    out += [("EXP3", "ind_pair", (gen, alt, h, n, par)), ("EXP6", "ind_pair", (gen, alt, h, n, par))]
            for gen in ("G0",) + PROCESS:
                out.append(("EXP7", "ind_pair", (gen, "null", h, h, n, 0.0)))
            out += [("EXP7", "ind_pair", ("FACTOR", "null", h, h, n, 0.0)), ("EXP7", "ind_pair", ("SCALE", "null", h, h, n, 0.1))]
            for alt, par in (("lin", 0.065), ("sq", 0.075), ("lag", 0.13), ("flat", 0.13), ("factor", 0.0), ("scale", 0.1)):
                out.append(("EXP7", "ind_pair", ("G0", alt, h, h, n, par)))
            out += [("EXP3", "mmd_series", ("G0", "null", h, n, 0.0)), ("EXP3", "mmd_series", ("G0", "shift", h, n, 0.1)),
                    ("EXP5", "mmd_series", (h, n, 0.0)), ("EXP5", "mmd_series", (h, n, 0.1))]
            for gen, alt, par in (("G0", "null", 0.0), ("HT-S", "null", 0.0), ("G0", "mean", 0.1), ("G0", "mean", 0.5),
                                  ("G0", "scale", 1.1), ("G0", "shape", 0.1)):
                out.append(("EXP7", "mmd_series", (gen, alt, h, n, par)))
            for gen in ("G0", "SEAS-M", "SEAS-V", "SEAS-M-shared", "SEAS-V-slow"):
                out.append(("EXP4", "pair", (gen, h, n)))
        out.append(("EXP7", "ind_pair", ("G0", "null", 0.6, 0.9, n, 0.0)))
    return out


def key(exp: str, fn: str, args: tuple, seed: int) -> str:
    return f"{exp}.{fn}{args!r}@{seed}"


def digest(result) -> str:
    arrays = result if isinstance(result, tuple) else (result,)
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(f"{a.dtype.str}{a.shape}".encode())
        h.update(a.tobytes())
    return h.hexdigest()


def compute(mods: dict) -> dict[str, str]:
    out = {}
    for exp, fn, args in cases():
        for seed in SEEDS:
            rng = np.random.default_rng(np.random.SeedSequence(entropy=2026110302, spawn_key=(seed,)))
            out[key(exp, fn, args, seed)] = digest(getattr(mods[exp], fn)(*args, rng))
    return out
