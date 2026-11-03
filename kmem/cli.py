"""Command line: `python -m kmem <command>` from the repo root (kmem is not installed in the shared env, decision D7).

BLAS is pinned to one thread before numpy loads: threaded BLAS changes the last bits of Gram-matrix sums,
and exact reproduction of the saved rows needs one thread.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys  # noqa: E402

import typer  # noqa: E402

from kmem import environment  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True, help="Kernel tests under long memory.")


@app.callback()
def _root() -> None:
    """Kernel tests under long memory. Run from the repo root."""


@app.command()
def env() -> None:
    """Check the environment against pyproject.toml. Fails only if Python, numpy or scipy differ, the versions
    that can change saved numbers; other differences are notes."""
    report = environment.check()
    numeric = ("python",) + environment.CORE
    for name, (want, have) in report.items():
        mark = "ok  " if want == have else ("BAD " if name in numeric else "note")
        typer.echo(f"{mark} {name:12} pinned {want:10} installed {have}")
    if any(report[k][0] != report[k][1] for k in numeric):
        typer.echo("Python, numpy or scipy differs from the pins: saved rows may not reproduce exactly")
        raise typer.Exit(1)


def main() -> None:
    app(prog_name=f"{os.path.basename(sys.executable)} -m kmem")
