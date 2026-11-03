"""G-legacy and the pinned raw outputs.

kmem/legacy is the validated EXP2 core and must never change: its files must match tests/legacy.sha256.
experiments/OUTPUTS.sha256 pins every saved output of EXP2 to EXP7 (decision D2). This test checks the
reproduction inputs among them (main-run CSVs and manifests) that exist on this machine. A fresh clone
has none of them, because raw outputs are not committed; the row fixtures in tests/repro carry the sample.
"""
from __future__ import annotations

import hashlib
import os

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _listing(path: str) -> list[tuple[str, str]]:
    with open(path) as f:
        return [(line[:64], line[66:].rstrip("\n")) for line in f if line.strip()]


def test_kmem_legacy_is_unchanged() -> None:
    pinned = _listing(os.path.join(REPO, "tests", "legacy.sha256"))
    on_disk = sorted(f"kmem/legacy/{n}" for n in os.listdir(os.path.join(REPO, "kmem", "legacy"))
                     if n.endswith(".py"))
    assert sorted(p for _, p in pinned) == on_disk, "kmem/legacy gained or lost a .py file"
    changed = [p for digest, p in pinned if _sha256(os.path.join(REPO, p)) != digest]
    assert not changed, f"kmem/legacy changed: {changed}"


PROVENANCE = {  # kmem/legacy/PROVENANCE.md: the EXP2 sha256 of each file, and the one import line kmem changed
    "mbcore.py": ("9c1a3d4af589f5cd2d3cef3433d2988f3085752c1e28eab9317977b61b5e6643", None),
    "common.py": ("0ed009c55c58438a58d7e793082a73b85e33cc6702cfa4f13bac5bef1852e2c1", None),
    "mmdcal.py": ("963f9cd4f45d40baf1df38070c47a0fc5e651a1c3fd38d913d036e3bef5ee43e",
                  (38, "from . import common as C", "import common as C")),
    "gens.py": ("605e2cad5363a1eca2393887f40a724bf754261ea99738142bfd8cacdadd3ff6",
                (13, "from . import mbcore as mb", "import mbcore as mb")),
}


@pytest.mark.parametrize("name", sorted(PROVENANCE))
def test_kmem_legacy_is_the_exp2_core(name: str) -> None:
    """The EXP1-EXP4 copies of the core were deleted (decision D6). kmem/legacy is the same code: undoing the
    documented import line gives the EXP2 file's recorded hash."""
    digest, undo = PROVENANCE[name]
    with open(os.path.join(REPO, "kmem", "legacy", name), "rb") as f:
        data = f.read()
    if undo is not None:
        line, new, old = undo
        lines = data.split(b"\n")
        assert lines[line - 1].startswith(new.encode())
        lines[line - 1] = lines[line - 1].replace(new.encode(), old.encode(), 1)
        data = b"\n".join(lines)
    assert hashlib.sha256(data).hexdigest() == digest


def test_saved_main_outputs_are_unchanged() -> None:
    pinned = [(d, p) for d, p in _listing(os.path.join(REPO, "experiments", "OUTPUTS.sha256")) if "/main/" in p]
    present = [(d, p) for d, p in pinned if os.path.exists(os.path.join(REPO, p))]
    if not present:
        pytest.skip("no saved main outputs on this machine")
    changed = [p for d, p in present if _sha256(os.path.join(REPO, p)) != d]
    assert not changed, f"{len(changed)} saved main outputs changed, e.g. {changed[:5]}"
