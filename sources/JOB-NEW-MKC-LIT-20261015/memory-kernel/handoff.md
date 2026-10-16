# Handoff — memory-kernel workstream

To: coordinator. Execution `claude-mkc-kernel-20261015`, 2026-10-15. Bounded scan. No novelty claim.

## What the field looks like

Kernel methods work well under short memory. Long memory breaks them in a set order.

- Consistency mostly survives. Kernel estimates still converge, just more slowly.
- Valid tests do not survive. The wild-bootstrap theorems behind HSIC and MMD tests need dependence
  to die off fast (about lag^−6). Long memory, and even slow power-law dependence, fails that.
- There is a clean boundary at Hurst H = 3/4. Below it, kernel-type statistics still have a normal
  limit. Above it, the limit changes and its law is unknown. Two independent papers (2021–26 and
  Sep 2026) land on this same line.
- Finance straddles that line. Order signs sit just below it (H about 0.7). Realised volatility sits
  well above it (H about 0.9).

On memory kernels (Hawkes, impact propagators), recent RKHS estimators exist (2025, 2026). Both cut
the kernel off at a fixed lag. Neither gives a convergence rate. So power-law tails near criticality
are not covered. That regime is exactly what links Hawkes models to rough volatility.

## The 2–3 most promising gaps

1. ⚠ **Signature-kernel causal discovery from one path with fractional noise (K6, K7).** This is
   the triple intersection. Manten et al. (ICLR 2025) show their signature-kernel CI test keeps its
   power on fBm-driven systems. But their discovery algorithm stops being valid without independent
   noise increments. All their guarantees also need many independent paths, and markets give one.
   - Cheapest check: rerun their fBm example on blocks of one long path.
   - Overlap: partly Track E; otherwise none.
2. **RKHS estimation of power-law Hawkes and propagator kernels near criticality, with rates (K4,
   K5).** This is the best match to "RKHS kernels meet memory kernels".
   - Cheapest check: simulate near-critical power-law Hawkes. Measure how much the support cut-off
     biases the branching ratio.
   - Overlap: Track A (Hawkes power-law item, Paper 2 propagators) and Track G (mixing lemma).
3. **Size map for kernel tests across memory levels (K1, K2, contested C1).** This is a
   simulation with no data needed. It decides whether long memory is a scope note (Track E's current
   view) or the main obstacle for volatility-channel kernel tests. Run it first, whichever direction
   is chosen.
   - Overlap: Track E (N3 audit, G5).

## Overlap with existing tracks

- E owns short-memory kernel CI (CL-E4). Its G5 "kernel CI under long memory" is ranked low there.
  This scan suggests it matters for volatility series.
- F-KERNEL owns kernel conformal scores. Nothing found covers long memory (K9).
- J owns kernel proxy methods under β-mixing (CL-J1). The long-memory version is unclaimed (K10).
- A and G already hold the Hawkes mixing question. Who owns the Hawkes lemma is still an open
  author decision.

## For the coordinator (triple intersection)

- K3: no nonlinear kernel CI or Granger test under long memory was found. The only one is a linear
  Gaussian test (Zhai, Zhong & Wu, Sep 2026), valid only for H < 3/4. memory-causal may hold the
  causal side.
- K12: one unrefereed 2026 paper uses signatures to detect the direction of Hawkes excitation. It
  covers exponential kernels only. Power-law and single-path versions are open.
- Two unverified leads on near-critical Hawkes inference are passed on in search-log.md.

Details: field-map.md, gap-map.md, directions.md, search-log.md, bibliography.bib.
