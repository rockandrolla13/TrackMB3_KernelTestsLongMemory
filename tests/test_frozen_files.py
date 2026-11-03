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


def test_saved_main_outputs_are_unchanged() -> None:
    pinned = [(d, p) for d, p in _listing(os.path.join(REPO, "experiments", "OUTPUTS.sha256")) if "/main/" in p]
    present = [(d, p) for d, p in pinned if os.path.exists(os.path.join(REPO, p))]
    if not present:
        pytest.skip("no saved main outputs on this machine")
    changed = [p for d, p in present if _sha256(os.path.join(REPO, p)) != d]
    assert not changed, f"{len(changed)} saved main outputs changed, e.g. {changed[:5]}"
