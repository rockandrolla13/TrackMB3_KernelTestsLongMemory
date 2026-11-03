# Extending kmem

Each recipe below needs one new file and one registration line. No existing code changes. Each has a
worked toy example that runs as a test. Run everything from the repo root (`python -m kmem ...`).

The one rule that never bends: **saved results must keep reproducing.** After any change, run
`python -I -m pytest tests -q`. Before merging, also run `KMEM_REPRO_TIER=slow python -I -m pytest tests/repro`
and `python -I tools/analysis_gate.py --check`. Never edit `kmem/legacy`.

## Add a kernel

1. Write a small class in a new module under `kmem/kernels/`, with these parts:
   - `name`, `kind` and `lag` (the number of points the kernel loses).
   - A `kind="feature"` kernel maps Gaussian features: `features(f, x, H)`, plus `gram(x, sigma, H)` for the
     exact path.
   - A `kind="own"` kernel builds a Gaussian kernel on its own series: `series(x)`.
   - A `kind="distance"` kernel maps the distance kernel: `from_distances(d)`.
2. Register it: `kernels.register(MyKernel(...))`, in `kmem/kernels/__init__.py` or in your experiment.
3. Pass its name in the kernel list, for example `ind_bundle(..., kernels=("raw", "mine"))`.

Example: `tests/test_kernel_registry.py`. It adds a toy kernel, and EXP7's bundles and exact path run it unchanged.

## Add a data generator

- **A process** (one series at memory H): `def draw(H, n, rng) -> np.ndarray`, then
  `registry.register_process("NAME", draw)`. It gets a null pair (two independent draws) and split-MMD series
  for free.
- **A dependent alternative**: `def alt(gen, H, n, par, rng) -> (x, y)`, then `registry.register_alternative(...)`.
- **A split-MMD change** of the second half: `def change(t, par) -> t_new`, then `registry.register_change(...)`.

Draw in a fixed order from the one `rng` you are given. Never use global random state.
Example: `tests/test_generator_registry.py`.

## Add a calibration rule

Build a `kmem.calibration.registry.CalibrationRule` and `register` it. A rule has these parts:
- the null draws it reads: `"wild"`, `"shift"` or `"subsample"`;
- `decide(stat, draws) -> {column: value}`, which the runner writes;
- `decision` and `reject_value`, which the analysis reads back with `rule.reject(data, column)`.

Then ask for it by name: `rules_fast(fx, fy, w, rules=FAST_RULES + ("mine",))`. The runner and the analysis use
the same object, so they cannot disagree about a rejection.
Example: `tests/test_rule_registry.py`.

## Add an experiment grid

1. Create `experiments/EXPn-<yyyymmdd>/` with a `DESIGN.md` that fixes the grid, the rules and the seeds before
   anything runs.
2. Add the masters to `kmem/seeds.py` (`MASTERS`). A reused master needs a written reason, and the tests enforce it.
3. Write the driver `run_expn.py` with:
   - `cell_list()`, which gives each cell an explicit `idx`, its seed index. Keep the list append-only: never
     insert a cell before existing ones.
   - `run_chunk((cell, reps, master, ...))`, one row per replication, seeded by
     `SeedSequence(entropy=master, spawn_key=(cell["idx"], rep))`. Pin BLAS to one thread before numpy loads.
   - `job_for(cell, target, master, opts) -> kmem.runner.Job`: the CSV columns, chunk size and extra arguments.
   - `main()`, which parses flags and calls `kmem.runner.run_jobs`.
   Compose the statistics from `kmem.generators`, `kmem.kernels`, `kmem.features`, `kmem.rules`, `kmem.hsic`
   and `kmem.mmd`. Do not copy them.
4. In the analysis, use `kmem.analysis` (Wilson, c*, bootstrap power, McNemar) and `cell_rng(seed_index)` for
   bootstrap streams. Read other experiments' outputs only through `require_input`.
5. After the run, add the experiment to `kmem.spec.EXPERIMENTS` and to `DRIVERS` in `tools/build_specs.py`, and
   declare its inputs. Run `python -I tools/build_specs.py EXPn`, then build its reproduction sample with
   `python -I tests/repro/build_fixture.py EXPn` (add it to `EXPERIMENTS` there).
   From then on, `python -m kmem reproduce EXPn` checks it.

Example: `tests/test_grid_recipe.py` runs a two-cell toy experiment, built only from kmem pieces, through the
runner. It checks that a rerun reproduces every value.
