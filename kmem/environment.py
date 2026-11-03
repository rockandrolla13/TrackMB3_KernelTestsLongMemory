"""The exact dependency pins of pyproject.toml, and a check of the running environment against them."""
from __future__ import annotations

import importlib.metadata
import os
import platform
import re
import tomllib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORE = ("numpy", "scipy")  # the packages whose version can change saved numbers
REFERENCE_PYTHON = "3.12.13"  # the interpreter of every saved EXP2-EXP7 run


def pins() -> dict[str, str]:
    with open(os.path.join(REPO, "pyproject.toml"), "rb") as f:
        deps = tomllib.load(f)["project"]["dependencies"]
    out = {}
    for d in deps:
        m = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^;\s]+)", d.strip())
        if not m:
            raise ValueError(f"pyproject dependency is not an exact pin: {d!r}")
        out[m.group(1).lower()] = m.group(2)
    return out


def check() -> dict[str, tuple[str, str]]:
    """name -> (pinned, installed) for Python and every pinned dependency."""
    out = {"python": (REFERENCE_PYTHON, platform.python_version())}
    for name, want in pins().items():
        try:
            have = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            have = "missing"
        out[name] = (want, have)
    return out
