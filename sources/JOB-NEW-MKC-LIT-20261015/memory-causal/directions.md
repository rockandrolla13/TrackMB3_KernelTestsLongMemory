# Candidate directions — memory-causal workstream

Execution `claude-mkc-causal-20261015`, 2026-10-15. Four candidates, ranked by fit to the
author's interests and by how cheaply each can be killed. None is a hypothesis or protocol; none
is approved. "Cheapest decisive check" means the smallest step whose result would change the
decision to continue. Checks marked *simulation* are proposals only; no experiment was run.

---

## C1. Where memory enters a causal estimator: population vs path-conditional effects

**Idea.** Split the orthogonal (AIPW/DML) score into a state term τ(X_t) − θ and a residual
term. Under sequential ignorability the residual term is a martingale difference (an MA(h) at
horizon h), so it carries no long memory. Only the state term does, and only when the effect
is heterogeneous in a long-memory state (volatility, volume, spread). Prediction: population-ATE
intervals under-cover and converge at T^{d−1/2}; path-conditional ATE intervals keep √T and
nominal coverage. (field-map §8, D.)

**Assumptions.** Sequential ignorability given the past; positivity; stationary long-memory state
with Hermite-rank-1 effect function; nuisances learned at rates that may themselves be slower
under long memory.

**Closest prior work.** Ballinari–Wehrli arXiv:2411.10009 (F; condition fails for persistent
processes, but stated on the score); Semenova et al. arXiv:1712.09988 (F, reused); Yu et al.
arXiv:2609.27531 (F); Balkus et al. arXiv:2601.10899 (F; as-i.i.d. folds under m-dependence);
Bojinov–Shephard arXiv:1706.07840 (design-based, survives); i.i.d. sample-vs-population ATE
literature (not searched here — pre-emption risk that this is considered obvious).

**Cheapest decisive check.** *Simulation*, under an hour of compute, no data: ARFIMA(0,d,0)
state X_t with d ∈ {0, 0.2, 0.4}; binary treatment with known logistic propensity in X_t;
outcomes with (a) constant effect, (b) effect = X_t. Report HAC-Wald coverage for the population
ATE and for the path-conditional ATE, with as-i.i.d. and buffered cross-fitting. Kill if
coverage of the population ATE stays nominal at d = 0.4 under (b). Second kill criterion: a
literature hit stating the decomposition for long-memory series.

**Fit.** High: memory × causal × finance (effects that scale with volatility are the norm). It
also answers G3 (do bespoke folds matter beyond m-dependence) in the same run.

**Overlap.** Track **K** (DML reliability — nuisance alignment; this adds a dependence axis),
**F-CAUSAL** (DR-ACI horizon-h pivot uses mixing), JOB-NEW-CAUSAL-RISK-EXCEED (its B1 baseline
assumes geometric β-mixing). Should be owned by K or F-CAUSAL, not a new track.

---

## C2. Granger causality for Hawkes near criticality, and the endo–exo confounder

**Idea.** Existing single-path Hawkes Granger inference (Wang–Kolar–Shojaie) holds for a fixed
stability margin γ < 1 and bounded intensity. Markets are fitted at branching ratios near 1
(contested). Two questions: (i) how does test size degrade as ρ → 1 and with power-law
kernels; (ii) does a slowly varying exogenous baseline create spurious cross-excitation
(endo–exo confounding), as it creates spurious criticality in the univariate case?

**Assumptions.** Stationary multivariate linear Hawkes; for (ii), a smooth or piecewise-constant
common baseline shared by two streams with no true cross-excitation.

**Closest prior work.** Wang–Kolar–Shojaie arXiv:2007.07448 (F); Eichler–Dahlhaus–Dueck
arXiv:1605.06759 (F, consistency only); Christgau et al. arXiv:2203.13559 (A, i.i.d. replicates);
Filimonov–Sornette arXiv:1308.6756 (S); Wheatley–Wehrli–Sornette QF 2019 (A); Cordi et al.
arXiv:1709.08516 (A); JOB-A-HAWKES-CHECK heuristic mixing time ≈ 2/(1−ρ)² events.

**Cheapest decisive check.** *Simulation*, no data: bivariate exponential-kernel Hawkes with zero
cross-excitation, ρ ∈ {0.5, 0.8, 0.95, 0.99}, T fixed; plus a variant with a shared slow
baseline. Run a likelihood-ratio or de-correlated-score Granger test. Kill if size stays within
±1 pp of nominal across ρ and baseline variants. Before simulating, read the ICLR 2026 Hawkes
identification paper and Rong & Nair (ICASSP 2026) to check pre-emption.

**Fit.** High for memory (self-excitation, criticality) and causal (Granger); direct
microstructure use (bid/ask/cancel streams).

**Overlap.** **B2-COMP** (coverage vs branching ratio — same simulator family), Track A/G Hawkes
mixing items, **E** (graph methods). Natural home: B2-COMP or the Track A ρ-dependence check
already proposed in JOB-A-HAWKES-CHECK.

---

## C3. Core vs reaction flow as a causal decomposition of impact with memory

**Idea.** Muhle-Karbe et al. (2026) write order flow as autonomous "core" orders (metaorders —
the source of memory, H₀ ≈ 3/4) plus endogenous "reaction" flow, and note the two cannot be
told apart in anonymous data. In causal terms, core flow is the treatment and reaction flow plus
price are responses; memory in the treatment is caused by splitting (Sato–Kanazawa). The
question: what in this decomposition is identified from aggregate data, and what needs
account-level labels or an instrument?

**Assumptions.** The structural Hawkes-on-Hawkes model; either labelled core flow (account IDs)
or an instrument for core flow (Track B2-QUEUE's queue position is one candidate).

**Closest prior work.** arXiv:2601.23172 (F); Leclère et al. arXiv:2607.03239 (F; model-based
counterfactual via shared noise); Naviglio et al. arXiv:2501.17096 (S); Maitrier et al.
arXiv:2502.16246 (A); Sato–Kanazawa arXiv:2301.13505, 2308.01112 (A); Bojinov–Shephard (only
randomised trading evidence).

**Cheapest decisive check.** Paper-only identification check: in the Muhle-Karbe et al. model,
is the reaction kernel φ₁ identified from the second-order structure of aggregate signed flow
alone? Write the spectral density of S_t in terms of (φ₀, φ₁, φ₂) and count. If identified, the
causal content is model-borne and the direction reduces to estimation; if not, the direction
needs labelled data — then check whether any account-level data are accessible (marketdb holds
LOB, TRACE, CRSP; no account IDs known). Kill if neither identification nor data is available.

**Fit.** Very high on interests (memory, rough dynamics, near-critical Hawkes, causal,
microstructure). Highest data risk.

**Overlap.** **B2-HTE** ("when does order flow move prices") and **B2-QUEUE** (instrument). This
should be offered to B2 rather than started separately.

---

## C4. Proxy identification under long memory (bridge to Track J) — triple intersection

**Idea.** Track J's kernel-proxy claims assume β-mixing. Check whether the proxy bridge-function
estimator's rate and inference survive a long-memory confounder, and whether the one paper
claiming long-memory proximal identification (Hao 2026) holds up.

**Assumptions.** Proximal completeness; stationary confounder with memory d; RKHS bridge
function (kernel route) or sieve (Hao's Bernstein route).

**Closest prior work.** Hao, doi:10.1371/journal.pone.0347847 (S; single author; time used as the
proxy; fBm limit claimed; proofs unread); Shi et al. arXiv:2108.13935 (A); Ying et al.
arXiv:2109.07030 (A); Track J README (catalog only).

**Cheapest decisive check.** Read Hao's supplementary proofs for Thm 4–5 (one session) and
Track J's CL-J1 mixing assumption side by side. Kill C4 if (a) Hao's result depends on time
being a valid proxy in a way that does not generalise, and (b) CL-J1's proof uses mixing only
through a concentration step that a functional-dependence bound already replaces.

**Fit.** Medium-high: the only candidate that touches all three themes. Weaker evidence base.

**Overlap.** **J** directly. Hand to the coordinator for the intersection synthesis and to
Track J's assumption audit (already its next action).

---

## Not recommended as directions

- Linear Granger testing under long memory (G5): largely handled (Bauer–Maynard).
- Neural causal discovery on event sequences (G13): no error control; low fit to a
  theory-first programme.
- Rough volatility as a causal variable (G12): premise contested; revisit only if C1 shows
  memory in the score matters.
