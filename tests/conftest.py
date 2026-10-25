"""Pin BLAS to one thread before numpy loads, as every EXP driver does.

Multi-threaded BLAS sums in a different order and changes the last bits of
large Gram-matrix reductions, so exact reproduction needs one thread.
"""
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
