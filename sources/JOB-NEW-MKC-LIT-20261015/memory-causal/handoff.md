# Handoff — memory-causal workstream

To the coordinator, for the author. Execution `claude-mkc-causal-20261015`, 2026-10-15.
Plan hash verified. No novelty claims. No decisions taken.

## What the field looks like

Causal inference on time series assumes short memory almost everywhere. Every single-series
DML paper found assumes fast (geometric or summable) mixing. One of them, Ballinari and Wehrli,
says outright that its condition is never met by a long-memory process.

Some things do survive memory. Randomised designs survive, as in Bojinov and Shephard, because
their validity comes from the randomisation, not from the data. Linear Granger tests have also
been made robust to long memory (Bauer and Maynard, 2012). So that question is mostly settled.

Self-exciting processes are a separate world. Granger structure for Hawkes is well defined. The
one test with a single-path guarantee (Wang, Kolar and Shojaie) needs a fixed margin from
criticality and bounded intensity. Markets may sit near criticality. That claim is itself
disputed: with a flexible background rate, criticality is rejected.

The newest finance work treats market impact as a counterfactual, but only through a model.
One 2026 paper splits order flow into "core" orders, which carry the memory, and "reaction"
orders. The authors note that the two cannot be told apart in anonymous data. That is an
identification problem, though they do not call it one.

## The three most promising gaps

1. **Which part of a causal estimator carries memory.** This is my own simple derivation, not a
   result from the literature. The DML score splits into a residual part and a state part. The
   residual part has no memory. The state part does, but only when the effect varies with a
   long-memory state such as volatility. If this holds, population-average effects lose the
   usual √T rate, while effects averaged over the observed path keep it. A one-hour simulation
   would settle it. It may be folklore, and it overlaps tracks **K** and **F-CAUSAL**.
2. **Hawkes Granger tests near criticality.** Nothing found shows how test size behaves as the
   branching ratio nears 1. Nothing found covers a slowly moving background rate that fakes
   cross-excitation either. A no-data simulation would settle it. Overlaps **B2-COMP** and the
   Track A ρ-dependence check.
3. **Core vs reaction flow as a causal question.** This fits the author's interests best:
   memory, rough dynamics, near-critical Hawkes, impact. It also carries the most data risk,
   because it likely needs account-level labels or an instrument. Overlaps **B2-HTE** and
   **B2-QUEUE**.

## Flags for the intersection synthesis (memory × kernels × causal)

- Track **J** (kernel proxy causal) assumes β-mixing. One 2026 paper (Hao, PLOS One) claims
  proxy-based identification under long-range dependence. Its proofs are unread, and it uses
  time itself as the proxy. This is the clearest case where all three themes meet.
- Local-independence tests (Christgau et al.) use flexible ML nuisances but need i.i.d. copies
  of the process. One long self-exciting path is open.
- "Kernel" means two things here: the Hawkes excitation kernel and an RKHS kernel. Keep them
  apart.

## Overlap with existing tracks

K, F-CAUSAL, B2-COMP, B2-HTE, B2-QUEUE, J and E all touch this workstream. No direction here
needs a new track. C1 belongs with K or F-CAUSAL, C2 with B2-COMP, C3 with B2, and C4 with J.

## Limits

- I screened 74 records. That is about 4 over the 70 limit; I found the overrun when I recounted
  at write-up. I extracted 24 records and did 9 new full-text reads.
- No SSRN or Google Scholar direct search. 2026 point-process causal papers are appearing
  monthly, so the pre-emption risk is real.
- No simulations or data were run. Every "check" is a proposal.
