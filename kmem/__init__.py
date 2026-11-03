"""kmem: kernel tests under long memory, shared by the experiments of this track.

kmem.legacy is the validated EXP2 core, frozen (tests/legacy.sha256). The rest is built on it:
generators, kernels, calibration rules, the HSIC/dCov/MMD bundles and the analysis helpers.
Every change must keep the saved EXP2-EXP7 rows reproducing exactly (tests/repro, experiments/REPRODUCIBILITY.md).
Run from the repo root: `python -m kmem --help`.
"""
__version__ = "0.2.0"
