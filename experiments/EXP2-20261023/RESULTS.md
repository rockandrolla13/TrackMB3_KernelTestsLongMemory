# EXP2 results (main run, 2026-10-25)

Run: 1,000 reps per cell, 90 cells, 90,000 replications, 0 errors, 7,453 s wall, 26.4 CPU-h.
Tables: `tables_main.md`. Rates: `summary_main.csv`. Verdicts: `verdicts_main.json`.
Exploratory simulation. Monte Carlo SE is about 0.7 pp at 5% and 1.6 pp at 50%.

## 1. Shared seasonality breaks every single-path test, including Pearson

SEAS-M (periodic mean, independent random phase per series), n = 2000, H = 0.5:
HSIC wild 91.2%, HSIC subsampling 90.9%, dCov wild 91.7%, dCov subsampling 91.3%,
Pearson-NW 91.4%, Pearson-EWC 91.7%. HSIC circular shift: 5.8%.

Reading: the two series are independent, but each path carries a fixed sinusoid with the same
period. The sample correlation of two such sinusoids tends to cos(phase gap), not zero. The
random phase acts like a component with infinite memory, so the paths are not ergodic. This is
Yule-type nonsense correlation. It is not kernel-specific, because Pearson fails equally.
Circular shift survives because a random shift re-randomises the phase, which matches the null.
With a common fixed phase (SEAS-M-shared) every rule rejects 100% and circular shift 28–34%.

The "dampens memory" verdicts for SEAS cells are a ceiling effect: rates near 90% cannot rise.

## 2. Seasonal variance breaks the kernel tests but not Pearson

SEAS-V, n = 2000, H = 0.5: HSIC wild 73.3%, HSIC subsampling 71.5%, dCov wild 51.7%,
dCov subsampling 45.6%. Pearson-NW 4.7%, Pearson-EWC 4.5%. Circular shift 5.0%.

Reading: a shared seasonal scale makes |X| and |Y| move together along each path, while the
signs stay unrelated. Pearson sees only the linear layer and stays at 5%. Kernel tests see the
second (squared) layer and reject. This mirrors the rank-two Hermite story in the paper: the
kernel test's extra sensitivity is also an extra exposure.

## 3. Heavy tails: instantaneous tails do nothing; linear-filter tails add a little

HT-S (t3 transform of fGn, Hermite rank one) tracks G0 within noise in every rule.
HT-L (ARFIMA with t3 innovations, not subordinated) adds about +5 to +9 pp at H = 0.8 under
subsampling (MMD 57.3% vs 49.1%; HSIC 25.9% vs 17.5%). It adds nothing at H = 0.5.

## 4. Nonlinear dynamics: SV hides memory from HSIC; TAR adds to it

NL-SV puts the memory in the volatility, so in the squared (rank-two) layer. HSIC and dCov
subsampling on independent SV pairs stay at or below 6% at every H (dCov 1.9% at H = 0.8).
For split MMD, SV breaks subsampling even with no memory (13.3% at H = 0.5). The AR(0.9)
volatility gives strong short-range persistence that a sqrt(n) block does not cover.
NL-TAR raises HSIC and dCov subsampling by about +11 pp at H = 0.8. Per Amendment (1), TAR adds
short-range persistence, so this excess is not a pure nonlinearity effect.

## 5. Self-normalised block MMD stays closest to nominal

sn_block stays in 3–13% across all process families at n = 2000. The worst cell is NL-SV at
H = 0.8 (13.0%). Uncentred wild MMD reaches 26–31% at H = 0.8 for all families.

## What this means for the paper

Memory is still the main driver for G0, HT and NL processes. Process features change the size
of the excess but not its direction, except SV, which moves memory into a layer the test sees
less. Seasonality is a separate failure mode. It breaks every single-path test, not only kernel
tests, and circular shift is the one rule that survives it.
