## Coordinator assessment and coverage

Both workers loaded the registered `research-gap-finder` SKILL.md (SHA-256 `e81cb1de…1d25919`) at matching hashes; the coordinator checked every `finished`-event output hash for both jobs on 2026-10-23. Retrieval was live WebSearch plus WebFetch primary-source opens. This is a bounded scan, not a systematic review; no-paper-found is a search result, not a novelty claim. Where the two workers disagree, the reviewer verdict is adopted only where the reviewer opened the primary text the discovery worker did not; a verdict resting on either worker's own algebra stays marked unresolved. No track state, hold, claim, protocol or experiment changed.

**Adjudicated.**
1. MB-G1-HSIC-MIXED-NULL: reviewer's **narrow** adopted with a qualification. Coordinator opened arXiv 1704.08609 (Düker): it proves functional limit theorems for sample means and cross-autocovariances of multivariate linear processes with mixed short/long-range components (limits include operator fBm and matrix Rosenblatt). It does not mention HSIC. The "linear case is done" reading therefore holds through the reduction of linear-kernel HSIC to squared cross-covariance, which the coordinator accepts but did not check against Theorem 4.1's exact conditions. The non-linear bounded-kernel mixed case remains the open part. Lévy-Leduc & Taqqu (Metrika 2014) is univariate by abstract; full text unread.
2. MB-G2-LOWFREQ-KERNEL: **narrow** adopted; the reviewer sourced discovery's "may fail above the line" (Müller–Watson need Gaussian averages; Düker gives non-Gaussian limits there). Priority lowered.
3. MB-UNIFORM-SIZE-LIMITS: **provisional, unranked** adopted; it is a reading task, and discovery's wording on Pötscher 2002 overreached per Preinerstorfer–Pötscher (read in full by the reviewer).
4. MB-G3-LMSV-RANK1: **reject** (both agree). Qualification from the reviewer kept: rank 2 is possible without leverage when the volatility function is symmetric, so G3's premise fails in the standard LMSV model but not in every model.
5. MB-HAWKES-NO-H-LINE: **reject** (both agree); near-critical Hawkes gives H below 1/2, so no 3/4 line arises.
6. MB-CROSS-CLASS-MAP: **provisional** (both agree); forward-citation sweeps found no cross-class map.

**Best next investigation.** G1 restricted to non-linear bounded kernels in the mixed-memory case, after confirming Düker Thm 4.1's conditions. Blocked: Pötscher 2002 and Lévy-Leduc–Taqqu full texts; Semantic Scholar rate limits.
