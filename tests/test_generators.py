"""Every experiment's generator gives the same arrays as before the move into kmem.generators (plan step 2.2).

The expected hashes in tests/repro/generators.json were computed once from the pre-move code (git e40180d)
with `python -I tests/test_generators.py --freeze`.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):  # also when run as a script
    os.environ[_v] = "1"

import json  # noqa: E402
import sys  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "repro"))
import generator_cases as G  # noqa: E402

FROZEN = os.path.join(HERE, "repro", "generators.json")


def test_generators_match_the_pre_move_arrays() -> None:
    with open(FROZEN) as f:
        expected = json.load(f)
    got = G.compute({exp: G.load(exp) for exp in G.SOURCES})
    assert set(got) == set(expected), "the case list changed"
    bad = sorted(k for k in expected if got[k] != expected[k])
    assert not bad, f"{len(bad)} of {len(expected)} generator outputs changed, e.g. {bad[:5]}"


if __name__ == "__main__" and sys.argv[1:] == ["--freeze"]:
    if os.path.exists(FROZEN):
        sys.exit(f"{FROZEN} exists; the hashes are frozen once, from the pre-move code")
    hashes = G.compute({exp: G.load(exp) for exp in G.SOURCES})
    with open(FROZEN, "w") as f:
        json.dump(hashes, f, indent=0, sort_keys=True)
        f.write("\n")
    print(f"froze {len(hashes)} generator outputs")
