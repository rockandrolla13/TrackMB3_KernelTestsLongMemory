# Code Review Report

**Project:** TrackMB3_KernelTestsLongMemory
**Date:** 2026-11-03
**Scope:** `kmem/` (including `kmem/legacy`, cited but never to be edited), `tests/`, all Python under `experiments/EXP1-20261016` to `experiments/EXP7-memory-kernel-20261101`, `investigations/memory-kernel-2026-11-01/theory_checks.py`, and `investigations/memory-kernel-2026-11-01/examples-report/*.py` as consumers of results.
**Findings:** 31 (Standards 20, Spec 11)

**Method and adaptations.**
- The review follows `~/.claude/skills/code-review/SKILL.md`, with four adaptations the team lead set.
- There is no diff against a fixed point. Much of EXP1 and all of EXP7 is untracked, so a diff would miss it. Each in-scope file was reviewed directly.
- The skill's `docs/agents/issue-tracker.md` requirement was skipped. The spec source is each experiment's `DESIGN.md`. EXP1 has no `DESIGN.md`, so its spec came from `checkpoint.md:18-32`, `results.md:9-15` and `manifest.json`.
- The two axes ran as parallel sub-agents, as the skill prescribes. Their reports are reproduced below, lightly cleaned and given `CR-<NNN>` IDs. The axes are not merged or reranked.
- Word caps were raised from 400 to about 1,300 per axis, because the scope is a whole codebase rather than one diff.
- The aggregator re-read about 30 of the cited line ranges. One citation was wrong and is corrected here: `--oracle-kappas` is defined at `run_exp6.py:81`, not line 231.
- Focus: extensibility and exact reruns. Style nits are left out.
- IDs use `CR-<NNN>`, as the team lead asked. `/home/ak/blueprint/tools/contract_lint.py --skill code-review` therefore fails its ID-format check, which expects `CR-<PILLAR>-<NNN>`. Its other four checks pass. The architecture review and the refactoring plan pass the linter in full.

---

## Executive Summary

The statistical code is careful, and every runner follows its own DESIGN seed rule. The problems are structural. On the **Standards** axis, the worst issue is that each new experiment re-types the HSIC and MMD test bundles instead of composing them (CR-001). On the **Spec** axis, every experiment's row values can in principle be reproduced, but nothing reproduces any of them from one command. Three concrete blockers stand out: EXP5's staged scripts depend on the wall clock (CR-026), EXP7's producing code is not in git (CR-030), and EXP1's driver changed after its run (CR-022).

One hard rule break exists. Two scripts that write saved outputs load numpy without pinning BLAS to one thread (CR-015).

## Standards Sources

The repo has no CLAUDE.md, CONTRIBUTING, `pyproject.toml`, `setup.cfg` or ruff config. These documented rules were applied:

1. `pytest.ini`: `testpaths = tests`, `pythonpath = .`.
2. `tests/conftest.py:1-9`: pin BLAS to one thread before numpy loads, "as every EXP driver does", because threaded BLAS changes last bits.
3. `kmem/__init__.py:3-4`: new code must pass `tests/test_reference_exp2.py`.
4. `kmem/legacy/__init__.py`, `kmem/legacy/PROVENANCE.md`: do not edit legacy.
5. `README.md` "Rules": exploratory evidence only.
6. Any `DESIGN.md` rule on how code must run: seeds, resume, output names. Examples are EXP7 `DESIGN.md:195-198` and EXP5 `DESIGN.md:104`.

The skill's Fowler smell baseline was applied on top. Smells are judgement calls, and a repo rule overrides them.

---

## Standards

Result: 1 hard violation (CR-015) and 19 judgement calls.

### CR-001: Test bundles re-typed, not composed
- **Location:** `experiments/EXP3-power-20261025/power_core.py:83-132` copies `kmem/legacy/mbcore.py:262-316` and adds quantiles. Further copies sit at `experiments/EXP6-bandwidth-20261025/bw_core.py:66-138` and `experiments/EXP7-memory-kernel-20261101/exp7_core.py:697-780`. The MMD bundle follows the same chain: `power_core.py:137-176`, then `exp5_core.py:65-105`, then `exp7_core.py:804-814`.
- **Type:** judgement. Duplicated Code, Data Clumps.
- **Snippet:** `ln = math.floor(n ** 0.5)` / `dn = math.floor(0.1 * n)` / `A = math.floor(0.1 * n)` (`power_core.py:86-88`, `bw_core.py:96-98`, `exp7_core.py:721-722`)
- **Why it matters:** A new calibration rule means editing four bundles, and each one keeps the random draw order by hand. The `(n, ln, dn, A)` clump wants its own type. So does `hsic_fields(x, y, dx, dy, sx, sy, w, n, ln, dn, A)` (`bw_core.py:66`).

### CR-002: Driver skeleton copied six times
- **Location:** `run_exp2.py:81-136`, `run_exp3.py:111-166`, `run_exp4.py:106-163`, `run_exp5.py:102-154`, `run_exp6.py:74-138`, `run_exp7.py:173-246`.
- **Type:** judgement. Duplicated Code.
- **Snippet:** `with mp.get_context("fork").Pool(a.workers, maxtasksperchild=50) as pool:`
- **Why it matters:** Each driver grew its own flags (`--family`, `--adj2000`, `--kappas`, `--extras-only`, `--reps4000`). A one-command rerun needs one driver that reads a spec.

### CR-003: The resume rule differs by runner, and the fix reached only two
- **Location:** `run_exp2.py:102-105`, `run_exp3.py:133-136`, `run_exp5.py:121-124` and `run_exp6.py:102-108` all count any row as done. `run_exp4.py:51-57` and `run_exp7.py:166-170` count only complete `status == "ok"` rows. The rule is stated in EXP7 `DESIGN.md:197`, and the fix is logged in EXP4 `DESIGN.md:154-156`.
- **Type:** judgement. Shotgun Surgery.
- **Snippet:** `done = {int(row["rep"]) for row in csv.DictReader(f)}`
- **Why it matters:** In four drivers, a resumed run never retries error rows or rows truncated by a crash. Only EXP6 checks the header before appending.

### CR-004: Generators copied, with names drifting
- **Location:**
  - `ind_pair`: `power_core.py:41-67`, `bw_core.py:41-61`, `exp7_core.py:109-154`.
  - `seas_scale`: `power_core.py:36`, `bw_core.py:36`, `exp7_core.py:74`.
  - SEAS-V-slow: `gens_exp4.py:15-20` and `exp7_core.py:79-83`.
  - `mmd_series`: `power_core.py:70-78`, `exp5_core.py:56-62`, `exp7_core.py:157-176`.
- **Type:** judgement. Duplicated Code.
- **Snippet:** `if alt == "shift":` (`power_core.py:73`) versus `if alt == "mean":` (`exp7_core.py:165`)
- **Why it matters:** One alternative has two names. The copies keep the same draw order only because cross-folder unit checks compare them.

### CR-005: Two double-centring formulas with different subtraction order
- **Location:** `kmem/memory_kernels.py:17-21` re-types legacy `double_center` (row means first, `mbcore.py:183-186`). Inline column-first versions sit at `power_core.py:156`, `exp5_core.py:85` and `exp7_core.py:811`.
- **Type:** judgement. Duplicated Code, reproduction hazard.
- **Snippet:** `kc = k - k.mean(axis=0, keepdims=True) - k.mean(axis=1, keepdims=True) + k.mean()`
- **Why it matters:** Floating-point subtraction depends on order. Merging the two forms may change the last bits of saved EXP3, EXP5 and EXP7 MMD rows. The refactor must keep both as named variants.

### CR-006: Repeated switches on the kernel name
- **Location:** `exp7_core.py:370-389`, `399-420`, `670-680`, `762-775`, `788-794`, `846-861`. The Hermite order is parsed from the name at `379`, `415`, `850`.
- **Type:** judgement. Repeated Switches, Primitive Obsession.
- **Snippet:** `u = mk.hermite_basis(x, int(kernel[-1]))`
- **Why it matters:** A new kernel needs six edits plus updates to the kernel tuples (`exp7_core.py:41-43`).

### CR-007: Rejection rules and size helpers re-declared per analysis
- **Location:**
  - Rule tables: `EXP2 analyze.py:24-35`, `EXP3 analyze.py:28-43`, `analyze_exp4.py:27`, `analyze_exp5.py:33-41`, `analyze_exp6.py:29-42`, `analyze_exp7.py:146-190`, `examples-report/build_data.py:263-268`.
  - `cstar` and `boot_power`: four copies.
  - `wilson`: `kmem/legacy/common.py:179` plus five copies.
  - Long mains that fuse loading, regression, bootstrap and writing: `EXP3 analyze.py:95-271` (177 lines), `analyze_exp6.py:116-253`, `analyze_exp4.py:97-214`, `analyze_exp7.py:321-407` and `544-655`, `EXP7 results_tables.py:32-168`.
- **Type:** judgement. Duplicated Code.
- **Snippet:** `"hsic_wild": lambda r: float(r["p_hsic_wb"]) <= ALPHA,`
- **Why it matters:** The names already drift (`hsic_wb` at `analyze_exp5.py:41`). No step can be shared between experiments.

### CR-008: A cell is identified by list position plus an f-string filename
- **Location:**
  - Index rules: `run_exp5.py:63-67` (`k = 100`) and `run_exp6.py:43`, `50` (`EXP3_IDX = [0, 2, 6, 18, 24, 30]`).
  - Names rebuilt elsewhere: `EXP2 analyze.py:49`, `EXP3 analyze.py:119`, `analyze_exp4.py:46`, `analyze_exp5.py:95`, `analyze_exp7.py:420`, and `build_data.py:257-260`, which copies `run_exp7.py:113-115`.
- **Type:** judgement. Primitive Obsession.
- **Snippet:** `f"mmd__G0__null0__H{h}__n{n}.csv"`
- **Why it matters:** Seeds depend on list position. Names embed float formatting (`H{h}`, `{param:g}`), so a float change renames files.

### CR-009: Missing input files are skipped silently
- **Location:** `EXP2 analyze.py:50-51`, `EXP3 analyze.py:120`, `analyze_exp5.py:96`, `analyze_exp7.py:78-79`.
- **Type:** judgement. Silent fallback.
- **Snippet:** `old = [r for r in csv.DictReader(open(p2)) if r["status"] == "ok"] if os.path.exists(p2) else []`
- **Why it matters:** After a rename or a move, the cross-experiment regression check comes back empty instead of failing.

### CR-010: One bootstrap stream shared across cells
- **Location:** `analyze_exp7.py:703-711` (draws at `333-334` and `425`), `EXP3 analyze.py:109`, `150-151`, `analyze_exp6.py:148`, `186-187`.
- **Type:** judgement. Determinism hazard.
- **Snippet:** `growth = growth_block(data, summary, rng)` / `power = power_block(data, summary, rng)`
- **Why it matters:** One added or missing cell shifts every later interval. A saved `results_main.json` reproduces only with the identical cell set.

### CR-011: The legacy core is copied into four folders, and EXP1 to EXP4 import the copies
- **Location:** `{common,mbcore,mmdcal,gens}.py` in EXP1 to EXP4 duplicate `kmem/legacy/*`. Import sites include `run_exp2.py:25-30`.
- **Type:** judgement. Duplicated Code.
- **Snippet:** `import mbcore as mb  # noqa: E402`
- **Why it matters:** The package that the reference test protects is not the code EXP2 to EXP4 actually load.

### CR-012: sys.path edits and imports across folders
- **Location:** `exp5_core.py:19-21` (inside a library module), `tests/test_memory_kernels.py:18-21`, `EXP7 unit_checks.py:29-36`, `EXP5 unit_checks.py:23-37`.
- **Type:** judgement. Coupling.
- **Snippet:** `sys.path.insert(0, os.path.join(REPO, "experiments", "EXP7-memory-kernel-20261101"))`
- **Why it matters:** `import gens` resolves by path order, so EXP4's `gens_exp4` gets EXP3's `gens` inside EXP7's unit checks. The package tests cannot run without `experiments/`.

### CR-013: A private cache read from other modules
- **Location:** `power_core.py:131`, `bw_core.py:137`, `exp7_core.py:757`. A second global, `_SN_TABLE`, sits at `exp7_core.py:51-52` and `606-611`.
- **Type:** judgement. Feature Envy, global state.
- **Snippet:** `out["r_ewc"] = abs(lt["t_ewc"]) / mb._TQ[nu]`
- **Why it matters:** Only `lind_tests` fills `_TQ`, as a side effect (`mbcore.py:366-367`). So the column depends on call order. `_SN_TABLE` loads a JSON file from the EXP7 folder.

### CR-014: EXP7 rounds CSV floats to 12 digits, and the manifest does not say so
- **Location:** `run_exp7.py:157-163`; manifest fields at `run_exp7.py:210-215`.
- **Type:** judgement. Determinism hazard.
- **Snippet:** `return f"{float(v):.12g}"`
- **Why it matters:** The other experiments write full precision. One shared writer would change someone's bytes, unless the format becomes part of each experiment's spec.

### CR-015: BLAS not pinned in two scripts that write saved outputs
- **Location:** `experiments/EXP4-deseason-20261025/unit_checks.py:10-11`; `investigations/memory-kernel-2026-11-01/theory_checks.py:3`, `19`.
- **Type:** **hard.** It breaks rule 2 (`tests/conftest.py:1-9`).
- **Snippet:** `Run:  OMP_NUM_THREADS=4 python -I theory_checks.py > theory_checks.out`
- **Why it matters:** `theory_exact.csv` is built with 4 threads, so its last bits may not reproduce. Every EXP2 to EXP7 runner pins BLAS before numpy, which passes. Only `run_exp7.py:18` also pins NumExpr.

### CR-016: Seed masters are hard-coded with no registry, and EXP3 and EXP4 share them
- **Location:** `run_exp3.py:32`, `run_exp4.py:38`; `run_exp6.py:34` reuses EXP3's main master for its `regress` mode by design.
- **Type:** judgement. Primitive Obsession.
- **Snippet:** `MASTER = {"pilot": 2026102541, "main": 2026102542}`
- **Why it matters:** The two streams differ only because EXP4 adds a third spawn-key element (`run_exp4.py:77`, `93`). EXP5 `DESIGN.md:104` and EXP7 `DESIGN.md:196` treat masters as never reused. Every runner follows its own DESIGN seed rule.

### CR-017: A flag-driven bundle, with its column list built in two places
- **Location:** `exp7_core.py:710-712` takes `exact, full, hsplit, dist, strict`, and no caller passes `dist`. Columns are built at `run_exp7.py:118-124` and again at `220-228`.
- **Type:** judgement. Speculative Generality, Divergent Change.
- **Snippet:** `kernels = ("raw", "inc1") if tune else ("raw", "inc1", "herm1", "herm2") if extras_only else K.IND_KERNELS` (`run_exp7.py:147`)
- **Why it matters:** Mode flags change which columns exist. A second column builder can drift from the first.

### CR-018: The reference test pins only EXP2, and it re-types driver constants
- **Location:** `tests/test_reference_exp2.py:27`, `54-63`.
- **Type:** judgement, read against rule 3 (`kmem/__init__.py:3-4`).
- **Snippet:** `B_WB = 500  # run_exp2.B_WB`
- **Why it matters:** The manifest already stores `"B_wild": 500`, and `_run` re-types `run_exp2.py:68-74`. EXP3 to EXP7 have no saved-row pin.

### CR-019: The report reads raw CSVs without the analysis filter
- **Location:** `examples-report/build_data.py:275` and `324`, against `analyze_exp7.py:40-58`. Its `wilson` (`build_data.py:219`) lacks the `m == 0` guard at `analyze_exp7.py:106`.
- **Type:** judgement. Duplicated Code.
- **Snippet:** `pd.read_csv(EXP / "main" / f"{_gname(h, n)}.csv")`
- **Why it matters:** A retried error row would count in the report but not in the analysis. No error rows exist today (EXP7 `DESIGN.md:387`).

### CR-020: Same name, different algorithm
- **Location:** `markov_rademacher` exists at legacy `mbcore.py:189` and at `mmdcal.py:112`. Three median bandwidths exist: `mbcore.py:159`, `mmdcal.py:192`, `exp7_core.py:179`.
- **Type:** judgement. Mysterious Name.
- **Why it matters:** "markov_rademacher" gives a different multiplier law depending on the module it is imported from. A package API should name each by its algorithm and leave legacy untouched.

---

## Spec

Spec sources: each experiment's `DESIGN.md`. For EXP1, the sources are `checkpoint.md`, `results.md` and `manifest.json`.

### Applies to every experiment

#### CR-021: Byte-identical reruns are impossible as the code stands — (c)
- **Spec:** Each DESIGN fixes a per-replication seed rule, so that a run can be repeated (for example EXP7 `DESIGN.md:195-196`). The team lead's goal asks for saved outputs to reproduce exactly.
- **Code:** Every row stores `cpu_s`. Several summaries carry CPU time: EXP2 `cell_cpu_s`, EXP3 `alt_cpu_s`, EXP5 `cpu_s_rep`, and every `tables_*.md`. Rows are appended in `imap_unordered` completion order (`run_exp2.py:123` and the five other drivers).
- **Gap:** Any reproduction check must compare values sorted by `rep`, with timing columns dropped. Values do reproduce. Each replication's seed depends only on (master, cell, rep), so resumes and top-ups do not change values.

### EXP1 (no DESIGN.md)

#### CR-022: The E1 driver changed after the run — (c)
- **Spec:** `manifest.json:71` records `run_e1.py` with sha256 `ef5f7e5c…`.
- **Code:** `run_e1.py` now hashes to `08308ecb…`. It was modified on 2026-10-23, six days after E1 ran. `run_e3.py`, `run_mmd.py`, `mbcore.py` and `mmdcal.py` still match their recorded hashes.
- **Gap:** E1 cannot be rerun from the repo. It also reads market data from `e1/data_cache/*.npz` (`run_e1.py:40`).

**Rerun:** E2, E3 and H5 probably reproduce. E1 does not.

### EXP2

#### CR-023: Design checks not coded as written — (a), (b), (c)
- **(a)+(c)** `DESIGN.md:153-155`: "if a feature cell rejects *less* than nominal by more than 2 pp, or the two H = 0.5 heavy-tail cells (same law) differ … it is flagged". `analyze.py:80-83` compares against G0, not against 5%. No code compares HT-S with HT-L.
- **(a)** `DESIGN.md:162`: "lag-1, 10, 100 ACF of X and |X|". `unit_checks.py:64` computes the |X| ACF at lags 1 and 100 only.
- **(b)** `analyze.py:96-97` adds a verdict label ("breaks alone; dampens memory") that the design does not define. `multiple_testing/mt_analysis.py:1-8` adds an unplanned multiple-testing re-analysis.
- **(c)** `DESIGN.md:125-126` says `manifest.json`. `run_exp2.py:113` writes `manifest_<stamp>.json`.

**Rerun:** Yes, for values. It was one uninterrupted run with 90 CSVs of 1,000 reps. The code has not changed since 97040bd. `tests/test_reference_exp2.py` already checks saved rows exactly.

### EXP3

#### CR-024: Undocumented stop and resume, and the pilot analysis has no code — (a), (c)
- **(a)** `DESIGN.md:97-98`: "The pilot uses the EXP2 G0 null p-values only to pick ρ, δ." `analyze.py:100` loads only `main_cells()`. The A1 table (`DESIGN.md:134-141`) cannot be rebuilt.
- **(c)** `DESIGN.md:65-66` asks the bootstrap to "recompute c* and power". `analyze.py:165-170` drops NaN null ratios for the point c* but keeps them as −inf in the bootstrap.
- **(c)** The main run stopped and resumed (`main.log:1` and `:73`, two manifests). The design does not record this. All 42 cells still hold reps 0 to 999 exactly once.

**Rerun:** Yes, for values. The analysis reads `../EXP2-20261023/main`.

### EXP4

#### CR-025: The `reused` flag depends on resume history — (c)
- **Spec:** `DESIGN.md:62-63`: "When P̂_X = P̂_Y = 50 … The test then reuses the A2 output (flag `reused = 1`)."
- **Code:** `run_exp4.py:90` reuses only when A2 ran in the same chunk. The resume logic at `run_exp4.py:125-140` can run A3 alone, and then it writes `reused = 0`. The values are identical either way.
- **Gap:** An exact column comparison must treat `reused` as schedule-dependent. `main.log` also shows a second stop, with a manifest at 14:10Z and no DONE before the 17:16Z resume. Amendment 4 does not mention it.

**Rerun:** Yes. Commit 4c81f9d changed only task scheduling. Every cell has reps 0 to 999 for each adjustment, with no duplicates and no errors. `summary_main.csv` has no timing column, so it should reproduce byte for byte. The analysis needs `../EXP2-20261023/summary_main.csv`.

### EXP5

#### CR-026: No single script reproduces the run, and Part 3 depends on the wall clock — (c)
- **Spec:** `DESIGN.md` amendments A2 and A3 define the staged run.
- **Code:** `run_main.sh` was stopped (A2). `run_part1b.sh` and `run_part2.sh` wait on `pgrep`. `run_part3.sh:5` waits for `PART2DONE` in `part2.log`. `run_part3.sh:9-10` runs n = 8000 only before 21:25 or 21:35 UTC.
- **Gap:** A rerun at another hour writes no n = 8000 rows and reports nothing. The final targets per cell are MMD 1,000; n = 4000 1,000 at H 0.7 and 0.8 and 500 at H 0.5; n = 8000 144 at H 0.7 and 0.8.

#### CR-027: Three design verdicts are not coded, and one baseline has a different source — (a), (c)
- **(a)** `DESIGN.md:96-97` defines when an excess at sum 1.4 is "flat". `analyze_exp5.py:212-220` stores z values but assigns no label.
- **(a)** `DESIGN.md:71-72`: "A rule that beats the linear oracle by a wide margin … is flagged." Nothing in `analyze_exp5.py:139-163` does this.
- **(a)** `DESIGN.md:49-50`: raw null rates "and the c* values must agree with EXP3". `analyze_exp5.py:99-106` reports both c* values but tests only raw size.
- **(c)** `DESIGN.md:84-85` names "the KL reference row in `e3_tables.md`". `analyze_exp5.py:27` reads `JOB-NEW-MEMORY-BOUNDARY-SIM-20261015_kernel-linear/cells/cell_summary.csv` instead.

**Rerun:** Yes for values, but no script replays the run. The analysis needs EXP3 main, EXP1 `e3_summary.csv` and the kernel-linear folder.

### EXP6

#### CR-028: The documented command is not the one that ran, and the stability verdict is not coded — (a), (c)
- **(c)** `DESIGN.md:146-150` adds oracle-only κ ∈ {0.25, 8} and sets 500 reps. `run_exp6.py:5` documents `--reps 1000`. `run_exp6.py:81` defaults `--oracle-kappas` to empty. The saved manifest argv is `--reps 500 --kappas 0.5,1,2,4 --oracle-kappas 0.25,8`.
- **(a)** `DESIGN.md:25-27`: "the size picture is **stable** if every rule keeps its EXP2 verdict class across κ." `analyze_exp6.py:87-92` assigns a class per κ. Nothing compares classes across κ or with EXP2.

**Rerun:** Yes for values, but only with the argv from `main/manifest_*.json`. All six cells hold reps 0 to 499.

### EXP6 and EXP7

#### CR-029: Saved outputs with no script that writes them — (b)
- **Code:** No file in EXP6 writes `results.json` or `regression_vs_exp3.txt`. No file in EXP7 writes `results.json`, yet `examples-report/build_data.py:337` reads it. `analyze_exp7.py` writes `results_<tag>.json` only.
- **Gap:** The report depends on a file that cannot be regenerated. `results_tables.py` and `plot_growth.py` are also missing from the EXP7 code list (`DESIGN.md:219-222`).

### EXP7

#### CR-030: The producing code is untracked, and a prerequisite is implicit — (c)
- **Spec:** `DESIGN.md:219-222` lists the code. Seeds, cell counts (142), rows (166,400), exact-path rep counts, Harris N (19, 99, 199), Bonferroni levels, the extra chain for inc50 and the SN table (40,000 reps, master 2026110174) all match.
- **Code:** The whole folder, `kmem/memory_kernels.py` and `tests/test_memory_kernels.py` are untracked in git. `sn_tables.py` must run first and takes about one hour. The run itself costs about 44 CPU-hours.
- **Gap:** Nothing pins the code that made 166,400 rows. `run_exp7.py:6` documents `--out pilot_hsplit`, while `DESIGN.md:309` says the extras went to `pilot_extras/`.

#### CR-031: The P3 test is looser than written, and a unit check was relaxed — (c)
- **(c)** `DESIGN.md:180-181`: "P3 holds if herm1 helps in G0 and HT-S at some H ≥ 0.8 where raw is not valid". `analyze_exp7.py:501-502` pools G0 and HT-S with `any`. The saved verdict does not change, because HT-S at H = 0.8 under the wild rule also helps.
- **(c)** `unit_checks.py:327` uses a bound relaxed after the pre-registered one failed. `DESIGN.md:351` and `:369` disclose this.

**Rerun:** Probably yes, for values. It was one invocation (manifest `already_done` is 0), and code timestamps predate the 23:57Z start. Amendment A6 recomputed 16 post-stall rows exactly.

---

## Summary

- **Standards:** 20 findings (1 hard, 19 judgement). Worst within the axis: CR-001. Every new kernel, rule or generator lands in hand-synced copies of the test bundles.
- **Spec:** 11 findings. Worst within the axis: CR-026. EXP5's staged scripts depend on the wall clock, so a scripted rerun can silently drop the n = 8000 rows.

## Handoff

### Summary Table

| Finding ID | Axis | Severity | Pillar | Location | Summary |
|---|---|---|---|---|---|
| CR-001 | Standards | 🔴 | Duplicated Code, Data Clumps | `power_core.py:83-176`, `bw_core.py:66-138`, `exp5_core.py:65-105`, `exp7_core.py:697-814` | HSIC/dCov and MMD bundles re-typed per experiment with hand-kept draw order. |
| CR-002 | Standards | 🟠 | Duplicated Code | `run_exp2.py:81-136` and five other drivers | Driver skeleton copied six times, each with its own flags. |
| CR-003 | Standards | 🟠 | Shotgun Surgery | `run_exp2.py:102-105`, `run_exp3.py:133-136`, `run_exp5.py:121-124`, `run_exp6.py:102-108` | Four drivers count error and truncated rows as done. Only EXP4 and EXP7 have the strict rule. |
| CR-004 | Standards | 🟠 | Duplicated Code | `power_core.py:36-78`, `bw_core.py:36-61`, `exp5_core.py:56-62`, `exp7_core.py:74-176`, `gens_exp4.py:15-20` | Generators copied three times. One alternative is named `shift` and `mean`. |
| CR-005 | Standards | 🟠 | Duplicated Code | `memory_kernels.py:17-21`, `power_core.py:156`, `exp5_core.py:85`, `exp7_core.py:811` | Two double-centring subtraction orders. Merging them can change saved last bits. |
| CR-006 | Standards | 🔴 | Repeated Switches | `exp7_core.py:370-420`, `670-680`, `762-775`, `788-794`, `846-861` | Six switches on the kernel name. Hermite order parsed from the name string. |
| CR-007 | Standards | 🟠 | Duplicated Code | 7 rule tables, 6 Wilson copies, 4 `cstar` copies, long analysis mains | Rejection rules and stats helpers re-declared per analysis, already drifting. |
| CR-008 | Standards | 🟠 | Primitive Obsession | `run_exp5.py:63-67`, `run_exp6.py:43-50`, `build_data.py:257-260` | Seed index is list position. Names rebuilt with f-strings in five places. |
| CR-009 | Standards | 🟡 | Silent fallback | `EXP2 analyze.py:50-51`, `EXP3 analyze.py:120`, `analyze_exp5.py:96`, `analyze_exp7.py:78-79` | Missing upstream files are skipped, not reported. |
| CR-010 | Standards | 🟡 | Determinism hazard | `analyze_exp7.py:703-711`, `EXP3 analyze.py:150-151`, `analyze_exp6.py:186-187` | One bootstrap stream across cells. One cell change shifts every later interval. |
| CR-011 | Standards | 🟠 | Duplicated Code | `experiments/EXP1..EXP4/{common,mbcore,mmdcal,gens}.py`, `run_exp2.py:25-30` | EXP1 to EXP4 load local copies of the frozen core, not `kmem.legacy`. |
| CR-012 | Standards | 🟠 | Coupling | `exp5_core.py:19-21`, `tests/test_memory_kernels.py:18-21`, EXP5 and EXP7 `unit_checks.py` | `sys.path` edits in a library module and in package tests. Module names collide. |
| CR-013 | Standards | 🟡 | Feature Envy | `power_core.py:131`, `bw_core.py:137`, `exp7_core.py:757`, `606-611` | Private `_TQ` read across modules. SN table loaded from a fixed folder. |
| CR-014 | Standards | 🟡 | Determinism hazard | `run_exp7.py:157-163`, `210-215` | EXP7 writes floats at 12 digits. The manifest does not record the format. |
| CR-015 | Standards | 🟠 | Hard: conftest BLAS rule | `EXP4 unit_checks.py:10-11`, `theory_checks.py:3` | Two scripts that write saved outputs load numpy without pinning BLAS. |
| CR-016 | Standards | 🟡 | Primitive Obsession | `run_exp3.py:32`, `run_exp4.py:38`, `run_exp6.py:34` | Master seeds hard-coded with no registry. EXP3 and EXP4 share masters. |
| CR-017 | Standards | 🟡 | Speculative Generality | `exp7_core.py:710-712`, `run_exp7.py:118-124`, `220-228` | Five-flag bundle with one unused flag. Column list built in two places. |
| CR-018 | Standards | 🟠 | Test coverage | `tests/test_reference_exp2.py:27`, `54-63` | Only EXP2 rows are pinned. The test re-types driver constants. |
| CR-019 | Standards | 🟡 | Duplicated Code | `build_data.py:219`, `275`, `324` | Report reads raw CSVs without the analysis filter. Its Wilson lacks a guard. |
| CR-020 | Standards | 🟡 | Mysterious Name | `mbcore.py:159`, `189`; `mmdcal.py:112`, `192`; `exp7_core.py:179` | Same function name, different algorithm, depending on the module. |
| CR-021 | Spec | 🟠 | Spec (c) | all drivers; summaries with CPU columns | Timing columns and completion-order rows make byte-identical reruns impossible. |
| CR-022 | Spec | 🔴 | Spec (c) | `EXP1 manifest.json:71`, `run_e1.py:40` | E1 driver hash no longer matches its manifest. E1 also needs market data. |
| CR-023 | Spec | 🟡 | Spec (a)(b)(c) | `EXP2 analyze.py:80-83`, `96-97`, `unit_checks.py:64`, `run_exp2.py:113` | Suspicion check and one ACF lag not coded as written. Extra label and re-analysis. |
| CR-024 | Spec | 🟡 | Spec (a)(c) | `EXP3 analyze.py:100`, `165-170`, `main.log:1`, `:73` | Pilot effect-size analysis has no code. NaN handling differs. Stop and resume undocumented. |
| CR-025 | Spec | 🟠 | Spec (c) | `run_exp4.py:90`, `125-140` | `reused` flag depends on resume history. A second stop is undocumented. |
| CR-026 | Spec | 🔴 | Spec (c) | `EXP5 run_part3.sh:5-10`, `run_part1b.sh`, `run_part2.sh` | Staged scripts wait on processes and the wall clock. A rerun can silently drop n = 8000. |
| CR-027 | Spec | 🟡 | Spec (a)(c) | `analyze_exp5.py:27`, `99-106`, `139-163`, `212-220` | Three DESIGN verdicts not coded. KL baseline read from a different source. |
| CR-028 | Spec | 🟠 | Spec (a)(c) | `run_exp6.py:5`, `81`; `analyze_exp6.py:87-92` | Documented command differs from the saved run. Stability verdict across κ not coded. |
| CR-029 | Spec | 🟠 | Spec (b) | EXP6 `results.json`, `regression_vs_exp3.txt`; EXP7 `results.json`; `build_data.py:337` | Saved outputs with no writing script. The report reads one of them. |
| CR-030 | Spec | 🔴 | Spec (c) | `experiments/EXP7-memory-kernel-20261101/`, `kmem/memory_kernels.py`, `exp7_core.py:51` | EXP7 code is untracked in git. The SN table prerequisite is implicit. |
| CR-031 | Spec | 🟡 | Spec (c) | `analyze_exp7.py:501-502`, `unit_checks.py:327` | P3 pooled with `any`, looser than written. Unit check relaxed after failing (disclosed). |

**Totals:** Standards 20 (worst CR-001). Spec 11 (worst CR-026).
