"""Every master seed of the track, in one place, with a uniqueness check (plan step 3.1, AR-EXT-003, CR-016).

Each replication's stream is SeedSequence(entropy=master, spawn_key=(cell seed_index, rep, ...)). Two runs
that share a master must differ in their spawn keys, or they draw the same numbers. A reuse is allowed only
with a written reason in ALLOWED_REUSE. A new experiment adds its masters to MASTERS and runs the tests.
"""
from __future__ import annotations

from collections import Counter

MASTERS: dict[tuple[str, str], int] = {
    ("EXP1-E1", "pilot"): 2026101611, ("EXP1-E1", "main"): 2026101612,
    ("EXP1-E2", "pilot"): 2026101621, ("EXP1-E2", "main"): 2026101622,
    ("EXP1-E3", "pilot"): 2026101631, ("EXP1-E3", "main"): 2026101632,
    ("EXP1-H5", "pilot"): 2026101651, ("EXP1-H5", "main"): 2026101652,
    ("EXP2", "mmd_pilot"): 2026102321, ("EXP2", "mmd_main"): 2026102322,
    ("EXP2", "ind_pilot"): 2026102331, ("EXP2", "ind_main"): 2026102332,
    ("EXP3", "pilot"): 2026102541, ("EXP3", "main"): 2026102542,
    ("EXP4", "pilot"): 2026102541, ("EXP4", "main"): 2026102542,
    ("EXP5", "pilot"): 2026102551, ("EXP5", "main"): 2026102552,
    ("EXP6", "pilot"): 2026102561, ("EXP6", "main"): 2026102562, ("EXP6", "regress"): 2026102542,
    ("EXP7", "pilot"): 2026110171, ("EXP7", "main"): 2026110172, ("EXP7", "tune"): 2026110173,
    ("EXP7", "sn_tables"): 2026110174, ("EXP7", "hcheck_repro"): 2026110175,
    ("EXP7", "layer_memory_check"): 2026110176, ("EXP7", "review_checks_bandwidth"): 2026110177,
}

ALLOWED_REUSE: dict[int, str] = {
    2026102541: "EXP4 shares EXP3's masters; EXP4 seeds data with spawn key (cell, rep, 0) and tests with "
                "(cell, rep, 1), three elements against EXP3's two, so the streams differ (EXP4 DESIGN.md).",
    2026102542: "As 2026102541 for EXP3/EXP4 main. EXP6 regress mode reuses EXP3's main master on purpose: it "
                "re-runs EXP3 cells with their EXP3 seed indices to check bit-for-bit agreement.",
}


def duplicates() -> dict[int, list[tuple[str, str]]]:
    count = Counter(MASTERS.values())
    return {m: [k for k, v in MASTERS.items() if v == m] for m, c in count.items() if c > 1}


def check() -> None:
    """Raise if a master is reused without a reason in ALLOWED_REUSE."""
    bad = {m: users for m, users in duplicates().items() if m not in ALLOWED_REUSE}
    if bad:
        raise ValueError(f"master seeds reused without a recorded reason: {bad}")


def master(experiment: str, key: str) -> int:
    return MASTERS[(experiment, key)]
