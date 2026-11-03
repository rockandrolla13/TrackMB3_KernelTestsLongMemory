"""The pins in pyproject.toml are exact and match the environment the saved rows were made with (decision D7)."""
from __future__ import annotations

import os
import subprocess
import sys

from kmem import environment

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_core_pins_match_the_running_environment() -> None:
    report = environment.check()
    bad = {k: v for k, v in report.items() if k in environment.CORE + ("python",) and v[0] != v[1]}
    assert not bad, f"environment differs from the pins: {bad}"


def test_python_m_kmem_env_runs_from_the_repo_root() -> None:
    p = subprocess.run([sys.executable, "-m", "kmem", "env"], cwd=REPO, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "numpy" in p.stdout and "scipy" in p.stdout
