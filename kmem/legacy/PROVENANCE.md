# Provenance of kmem.legacy

Copied on 2026-10-25 from `experiments/EXP2-20261023/`. The EXP2 copies are byte-identical
to the copies in EXP1, EXP3 and EXP4 (checked by md5 on 2026-10-25).

The only change is one import line per file, so the modules work as a package:

| file | source sha256 (EXP2) | change |
|---|---|---|
| mbcore.py | 9c1a3d4af589f5cd2d3cef3433d2988f3085752c1e28eab9317977b61b5e6643 | none |
| common.py | 0ed009c55c58438a58d7e793082a73b85e33cc6702cfa4f13bac5bef1852e2c1 | none |
| mmdcal.py | 963f9cd4f45d40baf1df38070c47a0fc5e651a1c3fd38d913d036e3bef5ee43e | line 38: `import common as C` -> `from . import common as C` |
| gens.py | 605e2cad5363a1eca2393887f40a724bf754261ea99738142bfd8cacdadd3ff6 | line 13: `import mbcore as mb` -> `from . import mbcore as mb` |

`tests/test_reference_exp2.py` re-runs saved EXP2 replications through this package and
requires every output column to match the saved CSV exactly.
