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


@app.command()
def run(experiment: str = typer.Argument(..., help="EXP2 ... EXP7"),
        out: str = typer.Option(..., help="output folder; an existing one is resumed"),
        workers: int = typer.Option(4, help="worker processes (at most the spec's max_workers)"),
        only: str = typer.Option(None, help="comma list of substrings of cell names"),
        reps: int = typer.Option(None, help="replications per cell instead of the saved targets")) -> None:
    """Run a past experiment's main grid from experiments/EXPn/spec.json, with its saved seeds."""
    import hashlib

    from kmem import drivers, runner, spec

    sp = spec.load(experiment)
    if workers > sp.max_workers:
        raise typer.BadParameter(f"{experiment} allows at most {sp.max_workers} workers")
    bad = [p.path for p in sp.prerequisites if not p.verify()]
    if bad:
        raise typer.BadParameter(f"prerequisite files differ from the spec: {bad}")
    drv = drivers.load(sp)
    jobs = drivers.jobs(sp, only, reps)
    with open(spec.spec_path(experiment), "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    manifest = {"mode": sp.mode, "spec": sp.id, "spec_sha256": digest, "masters": sp.masters,
                "seed_rule": sp.seed_rule, **sp.options}
    runner.run_jobs(jobs, drv.run_chunk, os.path.abspath(out), workers, manifest, writer=sp.writer,
                    max_tasks_per_child=getattr(drv, "MAX_TASKS_PER_CHILD", 50))


@app.command()
def reproduce(experiment: str = typer.Argument(..., help="EXP2 ... EXP7"),
              full: bool = typer.Option(False, "--full", help="rerun every saved cell (resumable) instead of the sample"),
              out: str = typer.Option(None, help="--full: output folder; an existing one is resumed"),
              workers: int = typer.Option(4, help="worker processes"),
              tier: str = typer.Option("all", help="sample tier: fast, slow or all"),
              dry_run: bool = typer.Option(False, "--dry-run", help="--full: list the rows a rerun would write")) -> None:
    """Check that a past experiment still reproduces its saved rows exactly (value-exact, timing excluded)."""
    from kmem import reproduce as R
    from kmem import spec

    sp = spec.load(experiment)
    bad = [p.path for p in sp.prerequisites if not p.verify()]
    if bad:
        typer.echo(f"prerequisite files differ from the spec: {bad}. Rebuild: {[p.made_by for p in sp.prerequisites]}")
        raise typer.Exit(2)
    if not full:
        res = R.sample(experiment, tier, workers)
        bad_rows = {k: v for k, v in res.items() if v}
        for k in sorted(bad_rows)[:10]:
            typer.echo(f"DIFF {k}: {bad_rows[k][:3]}")
        typer.echo(f"{experiment} sample ({tier}): {len(res) - len(bad_rows)} of {len(res)} saved rows reproduce exactly")
        raise typer.Exit(1 if bad_rows else 0)
    plan = R.planned_rows(experiment)
    if dry_run:
        for name, k in plan.items():
            typer.echo(f"{k:7d}  {name}")
        typer.echo(f"{experiment}: {len(plan)} cells, {sum(plan.values())} rows")
        return
    if out is None:
        raise typer.BadParameter("--full needs --out")
    run(experiment, out=out, workers=workers, only=None, reps=None)
    report = R.compare_run(experiment, os.path.abspath(out), full=True)
    ok = R.passed(report, full=True)
    for name, r in report.items():
        if r["diffs"] or not r["header_equal"] or r["rows"] != r["saved_rows"]:
            typer.echo(f"DIFF {name}: rows {r['rows']}/{r['saved_rows']}, header {r['header_equal']}, {r['diffs'][:3]}")
    typer.echo(f"{experiment} full: {sum(r['rows'] for r in report.values())} rows compared, "
               f"{'all reproduce exactly' if ok else 'DIFFERENCES FOUND'}")
    raise typer.Exit(0 if ok else 1)


@app.command()
def lineage() -> None:
    """Which experiment reads which other experiment's outputs (from the specs' inputs and consumers)."""
    from kmem import spec

    for exp in spec.EXPERIMENTS:
        sp = spec.load(exp)
        for i in sp.inputs:
            typer.echo(f"{i.path} -> {exp}   [{i.used_by}]{'' if i.verify() else '   MISSING OR CHANGED'}")
        for c in sp.consumers:
            typer.echo(f"{exp} -> {c.path}   [reads {', '.join(c.reads)}]")


def main() -> None:
    app(prog_name=f"{os.path.basename(sys.executable)} -m kmem")
