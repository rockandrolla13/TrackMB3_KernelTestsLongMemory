# What a kernel built to remove long memory can and cannot do

Theory agent, 1 November 2026. Exact numerics are in `theory_checks.py`, `theory_exact.csv` and `theory_exponents.csv`. Monte Carlo appears only in `theory_path_check.py`, and it is labelled there.

**Marks.** Every result carries one mark.
- **[Proved]**: the proof is here, or follows in a few lines from a source I opened.
- **[Proved, preprint]**: my proof is complete, but it rests on a theorem in an unreviewed preprint.
- **[Sketch]**: the argument is here, but a step is not written out.
- **[Conjecture]**: I believe it, and I have no proof.
- **[Refuted]**: the claim as stated is false, and I give the counterexample or the number.

## 1. Bottom line

1. The increment kernel does what the brief hoped, and more. Under independence, `E T_n^D` converges to a finite limit at every `H_X, H_Y < 1`. `T_n^D` has a limit law. The circular-shift rule is asymptotically exact for it, at a fixed bandwidth. The limit law and the shift result rest on one preprint theorem and on my own contraction bound.
2. The price is what the increment kernel cannot see. It is useless for split MMD. It is exactly blind to flat three-lag dependence. It is weak against a slow common factor.
3. Deflation helps at finite `n`, but it costs rate. It needs a deflation order that grows without bound as `H` approaches 1. Under stochastic volatility (SV) it reduces the memory but never removes it.
4. The Hermite-split design recovers the first-layer power that plain deflation loses. It recovers the raw rate, not the `sqrt(n)` rate. A Bonferroni combination is valid only if each part is valid. No hsplit part has a proof, so this result does not cover hsplit1 or hsplit2. The proved version, Harris on each layer plus Bonferroni, is not an EXP7 arm.
5. The cap on location power is a property of the problem. It is specific to a level shift in the layer that carries the memory.

| Claim | Verdict | Main change from the brief | EXP7 arms affected |
|---|---|---|---|
| C1 rank arithmetic | **Proved** | The pair covariance is exact, not asymptotic. Equality gives `log n`. | all |
| C2 PSD, rank, Gram form | **Proved** | `Q K Q` needs regression coefficients. The sharp rank condition is new. | herm1, herm2 |
| C2 SV removes nothing | **Weakened** | True for odd `r` only. `r = 2` removes 97% of the rank-one weight `c_1`. The rate stays. | herm1, herm2 |
| C2 blind set | **Weakened** | The blind set exists for a fixed reference. It is empty for pooled scores. | herm1, herm2 |
| C2 estimated `F` spoils deflation | **Refuted by path check** | Estimating `F` lowers the mean statistic. It also changes the limit law. | herm1, herm2 |
| C3 `E T_n^D = O(1)` | **Proved** | Also `E T_n^D = Lambda - c/n`. One differenced side is enough. | inc1, inc50, dcov_inc1 |
| C3 "any process" | **Refuted outside fGn-type spectra** | A memory pole away from frequency zero defeats lag 1. | inc1, inc50 |
| C3 limit law | **Proved, preprint** | Gaussian limit in Hilbert-Schmidt norm, fixed bandwidth. | inc1, inc50 |
| C4 target and blind set | **Proved** | Exact blind condition `r(0) = r(1) = r(-1)`. | inc1 |
| C4 oracle gives sum 1 | **Refuted** | The fractional filter gives 1.0026 to 1.0371 on fGn. | frac_oracle |
| C5 telescoping and cap | **Proved** | Exact order `I_n ~ n^(2-2H)`. The cap does not extend to SV level shifts or to scale shifts. | split-MMD arms |
| C6 lag-P and lag-1 | **Proved, with two corrections** | `T_n` attenuates by `(2 pi/P)^4`. Lag P does not remove the phase-dependent scale effect. | inc1, inc50 |
| C7 shift rule | **Proved, preprint** | Asymptotically exact for increments. Raw is covered only for `H_X + H_Y < 3/2`. | inc1, inc50, raw |
| C8 Hermite split | **Mixed** | See Section 9. | herm1, herm2 |

Section 10 lists every claim I could not verify and every claim the checks refuted.

## 2. Setting and the one exact identity

`X_t = g(xi_t, e_t)` and `Y_t = f(eta_t, u_t)` are independent. `xi` and `eta` are unit-variance fGn with `H_X, H_Y` in `[1/2, 1)`. `e_t` and `u_t` are iid and independent of the drivers. `D = 2 - 2H`. `rho(j)` is the fGn autocovariance, and `R = rho(1) = 2^(2H-1) - 1`.

Expand the feature in the driver: `E[phi(X_t) | xi_t] - mu = sum_{p>=1} a_p He_p(xi_t) / p!`. Set `c_p = ||a_p||^2 / p!`. Then `sum c_p <= E||phi(X) - mu||^2 =: kappa(0)`. Extras add to `kappa(0)` only.

**Identity [Proved].** Mehler's formula and conditional independence of the extras give, for `j != 0`,
`kappa(j) = E k(X_0, X_j) - E k(X, X') = sum_p c_p rho(j)^p`.
Independence of the two series then gives
`E n||C_n||^2 = sum_{|j|<n} (1 - |j|/n) kappa_X(j) kappa_Y(j)`.
Sample centring adds two terms that the script also computes exactly. Every number in `theory_exact.csv` comes from this identity. The script matches it to the literal double sum at `n = 1000`.

For the Gaussian kernel with `g = id`, `kappa` has a closed form. With `q = 2/(sig^2+2)`, `c_p = sig (sig^2+2)^(-1/2) [(1/2)_p / p!] q^p`.

## 3. C1: rank arithmetic. [Proved]

The pair `(p, q)` contributes `c_p^X c_q^Y rho_X(j)^p rho_Y(j)^q` to the summand. This holds exactly at every lag, not only in the tail. The pair is summable if and only if `p D_X + q D_Y > 1`. That is the same as `p(1-H_X) + q(1-H_Y) > 1/2`. Above the line it adds a bounded amount. On the line it grows like `log n`. Below the line it grows like `n^(1 - p D_X - q D_Y)`.

- The pair `(1,1)` gives `H_X + H_Y = 3/2`, as in the paper.
- For split MMD, rank `q` is summable if and only if `q(1-H) > 1/2`. Rank 2 gives `H = 3/4`. Rank 1 is never summable for `H > 1/2`. Both match the paper.
- Deflating `r` layers on both sides, at equal memory, needs `2(r+1)D > 1`. So the smallest order is `r*(H) = floor(1/(4(1-H)))`. That gives 0 below `H = 3/4`, 1 up to `7/8`, 2 up to `11/12`, and 5 at `H = 0.95`. **`r*` grows without bound as `H` approaches 1.**

"Short memory" here means a summable covariance. It gives the rate of `E T_n`. It does not give a limit law.

**Numbers.** Each pure pair matches its exponent with an error below `5e-5` at `n` up to about `2e6`. Centring keeps the exponent, with an error of 0.001. Whole-kernel growth exponents match within 0.02 uncentred and 0.04 centred for the Gaussian kernel, and within 0.05 and 0.10 for a kernel on `g = 2 sinh(xi/2)`. Where the prediction is a finite limit, the check asks only that the differences shrink. One qualification applies. At `H = 0.95`, neighbouring powers differ by only 0.1. The whole-kernel fit there brackets the prediction and does not pin it. The pure-pair check is what pins that case.

## 4. C2: the Hermite-deflated kernel

### What holds

**PSD [Proved].** `k_r(x,y) = <phi_r(x), phi_r(y)>` is an inner product, so it is positive semidefinite. It is not bounded, because `He_q(T(x))` is not. It is not translation invariant either.

**Gram form [Proved].** With regression coefficients on `B = [1, He_1(z), ..., He_r(z)]`, the residual Gram matrix is `Q K Q`, and `Q` is a projection. With sample-moment coefficients it is `(I-A) K (I-A)`, and `A` is not a projection unless the scores are exactly orthogonal. The EXP7 herm arms use `Q K Q`, which is the clean form.

**Sharp rank condition [Proved].** Let `b_p = E[He_p(xi) phi_r(X)]`. Write `Pi_r` for the `L^2(P_X)` projection onto `He_1(T), ..., He_r(T)`. Then `b_p` is the embedding of the signed measure `(I - Pi_r) E[He_p(xi) | X] dP_X`. For a characteristic kernel, `b_p = 0` if and only if `E[He_p(xi) | X]` lies in that span. So:
- For strictly monotone `g`, `T(X) = xi` or `-xi`, and the feature rank is at least `r+1`. This matches the brief.
- For anything else, the deflated rank is at least `r+1` only if `E[He_p(xi) | X]` is a polynomial of degree at most `r` in `T` for every `p <= r`.

**Non-monotone case [Proved].** For `X = xi^2`, `He_2(xi)` tends to `-1` as `T` tends to `-infinity`. A nonconstant polynomial in `T` is unbounded there. So `c_2` stays positive for every `r`, and the rank in `xi` stays 2. The size falls quickly: `c_2` keeps 10.8% at `r = 1`, 0.76% at `r = 2`, and 0.011% at `r = 3`.

### Stochastic volatility: "removes nothing" holds only for odd `r`

Take `X_t = exp(xi_t/2) e_t` with symmetric `e_t`.

- **Odd layers [Proved].** The law of `X` given `xi` is symmetric, and `T(-x) = -T(x)`. So `E[He_q(T(X)) | xi] = 0` for odd `q`. Each odd layer has zero correlation with every function of `xi`. Removing it changes `kappa(0)` and leaves `kappa(j)` unchanged for `j != 0`. The script confirms this with a largest change of about `1e-19` in `c_p`.
- **Even layers. The brief's claim is [Refuted] as a general statement.** Removing `He_2(T)` removes most of the memory. `He_2(T)` is a noisy proxy for `xi`, with correlation 0.488. With a Gaussian kernel on `X`, `c_1` falls to 2.87% at `r = 2`, 0.32% at `r = 4`, and 0.06% at `r = 6`.
- **The rate stays [Sketch].** The rank-one layer vanishes only if `E[xi | X]` is a polynomial in `T`. It is not, and `c_1 > 0` at each `r` I computed, up to `r = 6`. I have no complete proof for every `r`. The rank-one layer alone grows at exactly `n^(1-2D)`.

**What this means.** At `H = 0.95` and `n = 64000`, `E T_n` is 8.04 raw and 0.079 for population `r = 2`. Its local slope is still 0.69, and the limit slope is 0.80. So herm2 can look cured on SV at practical `n` without being cured. herm1 should not help at all. This agrees with EXP7 prediction P3 for herm1. It does not extend to herm2.

### Blind set and local rates

**Fixed reference [Proved].** Let `P` be a fixed law with continuous `F`. Let `w = dP_A/dP - dP_B/dP` be in `L^2(P)`. For a characteristic kernel, the deflated MMD vanishes if and only if `w` lies in `span{He_1(T), ..., He_r(T)}`. Blind pairs exist. Take `f = c T_+ + 1 - c/sqrt(2 pi)` and `f - c T` as the two densities, with `c = 1`. The raw MMD squared is 0.19 and the deflated one is below `1e-27`.

**Pooled reference. The claim "the kernel loses r directions" is [Refuted] for this case, and the blind set is empty [Proved].** With `P = lam P_A + (1-lam) P_B`, the density ratios are bounded, so `|w|` is bounded. `T` is unbounded under `P`, so a nonconstant polynomial in `T` is unbounded. Therefore `w` lies in the span only if `w = 0`. The empirical scores use the pooled law. So **deflated MMD on pooled scores stays consistent against every fixed alternative with continuous `F`.** What it loses is local sensitivity, not consistency.

**Local rates for a Gaussian mean shift, balanced split.**
- Fixed reference [Proved]. The generating function gives `E_A phi_r - E_B phi_r = -2 sum_{q>r, q odd} (delta/2)^q a_q / q!`. Even terms vanish exactly. The first surviving power is `delta^(r+1)` for even `r` and `delta^(r+2)` for odd `r`.
- Pooled reference [Sketch, confirmed numerically]. `w = -2 tanh(delta x / 2)`, and `T = x/s + O(delta^4)`. So `(I - Pi_1) w = (delta^3/12) He_3(T) + O(delta^5)`.
- At `r = 1`, the MMD squared is `delta^6 ||a_3||^2 / 144` pooled and `delta^6 ||a_3||^2 / 576` fixed. The script gets ratios of 0.995 and 0.999 to these constants. The pooled reference therefore gives a value four times larger, with the opposite sign in the embedding.
- An unbalanced split brings back the even term. At `lam = 0.25` and `r = 1`, the MMD squared has slope 4 in `delta`.

**Consequence for the review's heuristic.** The review says a rank-two feature detects a shift at `n^(H-1)` above `H = 3/4`. That needs a `delta^2` signal, so it holds only for an unbalanced split. The paper's balanced split has first surviving signal `delta^3`. Detection there needs `delta ~ n^(-1/6)` below `H = 3/4` and `n^(-(2/3)(1-H))` above it. Both are slower than raw `n^(H-1)` **[Sketch]**. Signal-equals-null-scale values are in Section 9.

**HSIC under a Gaussian copula [Proved].** Mehler gives `C^r = sum_{p>r} rho^p a_p tensor b_p / p!`. So `||C^r||^2 = rho^(2(r+1)) c_(r+1)^X c_(r+1)^Y (1 + O(rho))`. Series and closed form agree within `1e-15`. Once `r >= r*`, signal equals null scale at `rho ~ n^(-1/(2(r+1)))` **[Sketch]**. That is `n^(-1/4)` for `r = 1`. It is always slower than the increment rate `n^(-1/2)`.

### Estimated scores do not spoil the deflation

**Reduction [Sketch].** Dehling-Taqqu reduction gives `z_t` close to `(xi_t - xibar)/s`. I did not open that paper, so its use is **not checked**. In the path check, the largest gap between `z_t` and `(xi_t - xibar)/s` is about 0.35 to 0.47 at `n` from 1000 to 4000. It does not shrink over this range. I expect it sits at the extreme ranks, but I did not check where it sits.

**Rank-one leak [Path check, not a proof].** After `Q`, the norm of the coefficient on the true driver falls to about 1% of its raw size. It keeps falling with `n`, from about 1.2% at `n = 1000` to about 0.5% at `n = 4000` for `r = 2`.

**Mean statistic [Path check, not a proof].** Compared with `T` known, the ratio of means is 0.89 to 1.03 at `H = 0.7`, 0.84 to 1.01 at `H = 0.8`, and 0.34 to 0.50 at `H = 0.9`. These use 60 replications per cell. The job hit its time limit before `H = 0.9, n = 4000`, so those two cells are missing.

**Why estimation helps [Proved for the exact formula, numeric for the limits].** The second-layer contrast with a sample centre has variance `Var(Delta_2 - 2 xibar Delta_1)`, which has a closed form. Relative to a known centre, this ratio moves toward 1 below `H = 3/4`. Above `H = 3/4` it settles below 1. At `n` about `1e6` it is 0.54 at `H = 0.8`, 0.146 at `H = 0.9`, and 0.055 at `H = 0.95`. The sample centre removes the level component of the second layer.

**The cost is the limit law.** Bai-Taqqu 2019 Theorem 3.4 (opened, polynomial `G`) says that subtracting the sample mean can change the limit law without changing the order. So estimated scores help the mean, but a null law taken from known `T` would be wrong. No theorem I found covers an RKHS feature with estimated scores. Duker-Zoubouloglou Theorem 3.1 does not cover it.

**Bai-Taqqu 2019 Proposition 4.1 [checked].** It says Hermite rank equals power rank, for scalar `G` in `L^2` of a standard Gaussian. It applies to a coordinate of `phi_r` when `T` is fixed and the marginal is Gaussian. It does not apply to pooled scores, because `T` moves with the shift. It does not cover the balanced-split cancellation of even terms. It does not apply to SV, because the shift acts on `X` and not on `xi`. It does not give the HSIC trade-off. Mehler gives that.

## 5. C3: the increment kernel

Write `g_p(j) = 2 rho(j)^p - rho(j+1)^p - rho(j-1)^p`.

### Lemma 1: uniform in `p` [Proved]

For `j >= 2` and every `p >= 1`, with `R = rho(1) < 1`,
`|g_p(j)| <= A |Delta^2 rho(j)| + (B/2)(u^2 + v^2) <= C(H) j^(2H-4)`.
Here `u = rho(j+1) - rho(j)`, `v = rho(j-1) - rho(j)`, `A = sup_p p R^(p-1)` and `B = sup_p p(p-1) R^(p-2)`.

*Proof.* Taylor's theorem at `x = rho(j)` gives `(x+u)^p = x^p + p x^(p-1) u + (1/2) p(p-1) z^(p-2) u^2`, with `z` between `x` and `x+u`. For `j >= 2`, every argument lies in `[0, rho(1)]`, because fGn with `H > 1/2` has a positive, decreasing `rho`. Adding the two expansions gives the first bound. `Delta^2 rho(j)` is half a fourth difference of `|j|^(2H)`, so it is of order `j^(2H-4)`. `u^2 + v^2` is of order `j^(4H-6)`, which is no larger because `H < 1`. QED.

`A` and `B` blow up as `H` approaches 1, but `C(H)` does not. The bound `C(H)` is 0.25 to 0.83 for `H` from 0.55 to 0.99, because the differences shrink as fast. The largest value seen over `p < 500`, `2 <= j < 2000` is 0.23 to 0.54. The script confirms the inequality at each of those cells.

### Theorem A: the expectation [Proved]

For any measurable `g, f`, any kernels with finite `E k(X,X)`, and any iid extras,
`E T_n^D = sum_{|j|<n} (1 - |j|/n) kappa_X^D(j) kappa_Y^D(j)`, which tends to `Lambda`, with `|E T_n^D - Lambda| <= C/n`.

*Proof.* For `|j| >= 2`, `kappa^D(j) = sum_p c_p g_p(j)`, so `|kappa^D(j)| <= kappa(0) C(H) |j|^(2H-4)`. Extras enter only at `|j| <= 1`, where `kappa^D` is bounded by `4 kappa(0)`. The product of the two sides is of order `|j|^(2H_X + 2H_Y - 8)`, which is summable even with a factor `|j|`. QED.

So `T_n^D` is tight, by Markov's inequality. Boundedness of `k` enters only through `sum c_p < infinity`. So the expectation result should cover `dcov_inc1` once `E|X| < infinity` **[Sketch]**. The limit law below assumes bounded kernels, so it does not cover that arm as written.

**Numbers.** The slope of `E_(2n) - E_n` equals -1 within 0.02 at every `H` from 0.6 to 0.95, for lag 1 and lag 50, and for the mixed pairs. Between `j = 1e8` and `1e9`, `kappa^D` has slope -2.800 at `H = 0.6` and -2.115 at `H = 0.95`, against theory values -2.8 and -2.1. At `H = 0.95` the next power, -2.2, still pulls the slope. `Lambda` is 1.166, 0.981, 0.762, 0.638, 0.502, 0.357, 0.208, 0.0714, at `H` = 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95. It falls to 0.00004 at `H = 0.999`, because the features themselves stop moving. So a fixed cutoff from another `H` would be wrong. The rule must adapt.

**One side is enough [Proved].** With only `X` differenced, the summand is of order `|j|^(2H_X + 2H_Y - 6)`, which is summable for all `H < 1`. The target is `C_0 - C_(-1)`. It costs efficiency: 0.50 at `H = 0.5` and 0.36 at `H = 0.95`, against 0.67 and 0.86 for two sides.

### Scope limit [Refuted beyond this scope]

"Any process with long memory" is false. The proof uses the smoothness of `rho`. Take `rho(k) = cos(w k) rho_H(k)`, which is a valid covariance. Second differences do not gain decay away from frequency zero. At `H = 0.9`, `w = pi/3`, the exact slope for lag 1 is 0.600, which equals the rate `n^(1-2D)`. Lag 6 restores -1.000. The theorem holds for memory at frequency zero only. It needs a spectral singularity at zero or at a harmonic of `1/L`.

### Theorem B: the limit law [Proved, preprint]

**Review check.** I opened Duker-Zoubouloglou arXiv 2405.11452v3 and confirmed the review's claims.
- Assumption 3.1 asks `sum_v ||R(v)||_op^q < infinity`, with `q` the Hermite rank. Assumption 3.2 asks the same of `Gamma(v)` in trace norm.
- For `k_D`, `R(v)` is of order `v^(-D)`, and the HSIC product has rank `q = 2`. So Theorem 3.1 gives nothing past `H = 3/4`. The condition cannot see the cancellation.
- For deflated kernels it applies to the population form, with `T` known and `X` a function of the Gaussian alone, and with `q D_min > 1`.
- The review's coverage claim needs two more limits. It fails with extras. Encoding `e_t` as a Gaussian makes the rank 2 whatever `r` is. It also fails with estimated scores.

**Statement.** Take bounded continuous kernels on separable spaces, a fixed bandwidth, and `H_X, H_Y` in `[1/2, 1)`. Then `F_n = sqrt(n) C_n^D` converges in law in the Hilbert-Schmidt space to a centred Gaussian `Z`. Its covariance is `S = sum_j Gamma_X^D(j) tensor Gamma_Y^D(j)`, which converges in trace norm, with trace `Lambda`. So `T_n^D` converges in law to `||Z||^2 = sum_i lam_i N_i^2`. The same holds at any fixed lag `L`.

**Proof.** I use Theorem 3.4 of arXiv 2509.20342v3, which I opened. It has four conditions.

*Encoding.* Write `e_t = G(eps_t)` with `eps_t` iid standard normal. Then `D phi_t` expands in `He_p(xi_t) He_l(eps_t)` minus the same term at `t-1`. Layers are indexed by `(p, l, q, l')`. Chaos order is `r = p + l + q + l'`. Fixed `r` holds finitely many layers. Different layers are orthogonal.

*Conditions (i) and (ii).* The covariance of chaos `r` is a finite sum of rank-one terms. Each coefficient is `sum_{|j|<n} (1-|j|/n) g^X(j) g^Y(j)`, which converges by Lemma 1 and dominated convergence. Lemma 1 gives `sum_j sup_(p,q) |g_p^X(j) g_q^Y(j)| =: K < infinity`. So total trace is at most `K kappa_X(0) kappa_Y(0)`.

*Condition (iv).* The tail above chaos `N` is at most `K` times the sum of `c_p c_q` over `p + q > N`. This is uniform in `n` and tends to 0.

*Condition (iii), contractions.* The identity `||f tensor_m g||^2 = <f tensor_(r-m) f, g tensor_(r-m) g>` reduces cross terms to self-contractions of one layer. A self-contraction has squared norm `n^-2 sum_{s,t,u,v} P`, with `P` a product of edge terms around the cycle `s-t-v-u-s`. Each edge carries multiplicities that sum to `m` on `s-t` and `u-v`, and to `r-m` on `s-u` and `t-v`. Since `0 < m < r`, each edge has positive total multiplicity.

- *Case 1: no extras in the layer.* The discrete Leibniz rule moves each difference to one adjacent edge. Each vertex has two differences, one for `X` and one for `Y`. An edge with `d` differences decays like `<j>^(-a)` with `a >= d + D_*`, where `D_* = min(D_X, D_Y) > 0`. An edge with `d >= 1` is therefore summable. The total is 8 differences, at most 4 per edge. Two adjacent edges cannot both have `d = 0`, because their shared vertex must place its differences on one of them. So at most two edges have `d = 0`, and they are opposite. Then the sum over `j_1, j_2, j_3` is at most `C n^((1-2 D_*)_+) log n`. Hence the squared norm is at most `C n^(-min(1, 2 D_*)) log n`, which tends to 0 for every `H < 1`.
- *Case 2: extras in the layer.* An extra argument pairs only with an extra argument at an equal time, up to shift 1. By symmetry of the self-contraction, this makes two opposite edges compactly supported. The same bound follows, with no cancellation needed.

The theorem then gives `F_n -> Z`. The continuous mapping theorem gives the law of `T_n^D`. QED.

**Exact support for condition (iii).** For the `(1,1)` layer, `n^-2 tr((Gamma_X Gamma_Y)^2)` has slope exactly -1.000 for increments at `(0.95, 0.95)`, `(0.9, 0.9)` and `(0.95, 0.6)`. For the `(2,2)` layer with `m = (1,1)` at `H = 0.95`, the raw contraction grows from 20.8 to 78.4 for `n` from 12 to 36, and the increment contraction falls from 0.0068 to 0.0041.

**What rests on what.** Theorem 3.4 is a preprint result. The contraction bound is mine, and a referee should check it. An RKHS-valued Breuer-Major theorem as stated in 2405.11452v3 does not give this limit. The fourth-moment route does.

**Median bandwidth [Sketch].** All EXP7 arms use a data-driven bandwidth. Gaussian features of all bandwidths embed in one `L^2` space. The derivative in `sig` has the same Hermite structure, with a finite `sum c_p`. So `sup` of `T_n^D` over a compact bandwidth range is tight, and `F_n(sig-hat) - F_n(sig*)` tends to 0. I have not written out this step. Theorems A and B and C7 below are proved for a **fixed** bandwidth only.

### Theorem C: raw HSIC below the line [Proved, preprint]

The same proof, with no differences, covers raw `T_n` for `H_X + H_Y < 3/2`. Conditions (i), (ii) and (iv) use `|rho_X rho_Y|` summable. Condition (iii) uses Young's inequality, which gives a bound of order `n^((3 - 2(a_1 + a_2))_+)`, with `a_1 + a_2 = p D_X + q D_Y > 1`. That is `o(n)`. A short case check handles `a_i > 1` and layers with extras. The centring term has squared norm of order `n^(1 - D_X - D_Y)`, which tends to 0.

This covers **mixed pairs with one `H` above 3/4**, which Betken-Dehling Remark 2.1 leaves out, and transformed features with extras. It addresses the first open proof problem at `paper.tex:438`, for a fixed bandwidth and below the line. The exact contraction has slope -0.200 at `(0.9, 0.55)` and `(0.85, 0.6)`, which equals the theory value. So convergence there is slow, and a finite-sample law may look far from the limit.

**datadiff [Proved for Gaussian data, Sketch for SV].** For Gaussian data, raw increments have a summable `R(v)`, so Duker-Zoubouloglou Theorem 3.1 applies at every `H < 1`. This agrees with Bardet-Surgailis as the review reports it. I did not open Bardet-Surgailis. Under SV, `phi(X_t - X_(t-1))` has scale memory through `xi_t + xi_(t-1)`, so its rank-one layer keeps `n^(1-2D)` growth. This supports EXP7 prediction P2.

## 6. C4: target and power

**Target [Proved].** `C_n^D` has population target `2 C_0 - C_1 - C_(-1)`. For Gaussian pairs, the `p`-th term is `(2 r(0)^p - r(1)^p - r(-1)^p) a_p tensor b_p / p!`. The series matches a closed form built from 16 Gaussian integrals within `4e-16`.

**Consistency at equal memory [Proved].** For `Y = rho X + sqrt(1-rho^2) Z`, the target equals `2 Phi (T_rho - T_(rho R)) Psi^*`, where `T_s` is the Ornstein-Uhlenbeck operator. That operator difference has eigenvalue `rho^p (1 - R^p) != 0` on every layer. If both RKHSs are dense in `L^2(gamma)` modulo constants, as for Gaussian kernels, the target vanishes only at `rho = 0`.

**Exact blind set [Proved].** For Gaussian pairs and a dense RKHS, the target vanishes if and only if `2 r(0)^p = r(1)^p + r(-1)^p` for all `p`. The cases `p = 1, 2` and strict convexity give `r(0) = r(1) = r(-1)`. Flat dependence `r = 1/sqrt(3)` gives exactly 0.

**Near-blind set.** Take a common factor `F` with share 0.3. The ratio of target to raw HSIC is 1.91 at `H_F = 0.7`, 0.29 at 0.9, and 0.0034 at 0.99.

**Lag dependence, the reverse case [Proved].** If `r(0) = 0` and `r(1) = rho`, the raw HSIC target is 0 for a Gaussian pair, and the increment target is `-rho^p` per layer. So the increment kernel sees it and raw does not.

**First-layer efficiency.** The numbers in the brief hold.

| H | R | `sum r_D^2` (exact) | increment | one side | raw layer | fractional filter sum |
|---|---|---|---|---|---|---|
| 0.5 | 0 | 1.50000 | 0.667 | 0.500 | 1.000 | 1.00000 |
| 0.6 | 0.149 | 1.41542 | 0.707 | 0.495 | 0.924 | 1.00262 |
| 0.7 | 0.320 | 1.33520 | 0.749 | 0.479 | 0.519 | 1.00933 |
| 0.75 | 0.414 | 1.29716 | 0.771 | 0.465 | 0 | 1.01382 |
| 0.8 | 0.516 | 1.26075 | 0.793 | 0.447 | 0 | 1.01892 |
| 0.85 | 0.625 | 1.22616 | 0.816 | 0.423 | 0 | 1.02455 |
| 0.9 | 0.741 | 1.19360 | 0.838 | 0.394 | 0 | 1.03063 |
| 0.95 | 0.866 | 1.16327 | 0.860 | 0.358 | 0 | 1.03712 |

Lag sums and spectral integrals agree. The "should give sum 1" claim is **[Refuted]**. Fractional differencing with `d = H - 1/2` whitens ARFIMA exactly, but not fGn. The gap is small: up to 3.7% in variance. Only the exact whitening filter gives 1.

**Bigger cost than the table shows.** `z = sqrt(n) rho / sqrt(sum r_D^2)` looks only at the first layer. In the full kernel, the first-layer signal shrinks by `(1-R)^2`, while higher-layer noise shrinks by `(1-R^p)^2`, which is much less. The next table gives the `rho` at which the first-order signal equals the null mean. It is a proxy, not power.

| H, n | raw | herm r=1 | herm r=2 | inc1 | inc50 | stationary fractional |
|---|---|---|---|---|---|---|
| 0.5, 2000 | 0.051 | 0.237 | 0.389 | 0.063 | 0.063 | 0.051 |
| 0.8, 2000 | 0.065 | 0.241 | 0.390 | 0.085 | 0.076 | 0.065 |
| 0.95, 2000 | 0.118 | 0.348 | 0.474 | 0.116 | 0.160 | 0.097 |
| 0.8, 64000 | 0.014 | 0.101 | 0.219 | 0.015 | 0.013 | 0.011 |
| 0.95, 64000 | 0.064 | 0.219 | 0.315 | 0.021 | 0.028 | 0.017 |

**What this means.** On this alternative, increments beat deflation by a factor of 3 to 20. At `n = 2000`, raw and increments are close, and raw is slightly ahead up to `H = 0.9`. Raw has no proved cutoff above the line, and its null mean keeps growing. At `n = 64000` and `H = 0.95`, increments lead raw by a factor of 3. The fractional column uses a stationary filter with the whole past. EXP7's frac_oracle truncates at the sample start, so this column does not describe that arm.

## 7. C5: telescoping and the location cap

**Telescoping [Proved].** `sum_{t<=m} D phi_t = phi_m - phi_0`. So `E T_n^D = (n/m^2) [6 kappa(0) - 8 kappa(m) + 2 kappa(2m)]`, which is `O(1/n)`. Under a level shift the expected contrast is one spike of `D phi`, also of order `1/n`. The statistic depends on three observations, so it is not consistent against any fixed alternative. **Use no increment arm for split MMD.**

**The cap has exact order [Proved].** Paper Proposition 2(iii) caps every location-invariant test through `I_n`. I add the order.
- Lower bound: the half-mean contrast gives `I_n >= m^(2-2H) / (4 - 2^(2H))`.
- Upper bound: for `w` with `w'1 = 0`, replace `c` by `c - 1/2` and apply Cauchy-Schwarz in frequency. This gives `I_n <= (2 pi)^-2 integral sin^4(m l/2) / (sin^2(l/2) f(l)) dl`. With `1/f <= C |l|^(2H-1)` and `l = x/m`, this is at most `K m^(2-2H) integral_0^inf sin^4(x/2) x^(2H-3) dx`. The integral converges for `H < 1`.

So `I_n ~ n^(2-2H)`. **No location-invariant test has nontrivial power below `delta ~ n^(H-1)`.** Both bounds hold numerically. For example, at `H = 0.8`, `n = 16000`, the three values are 37.59, 48.37 and 55.15.

**Filtering cannot change this [Proved].** Any invertible filter of the invariant data leaves `I_n` unchanged. GLS on first differences gives the same information to within `1e-7`. A plain contrast on the differences is the telescoped one. At `H = 0.7, n = 600`, it has `I = 0.73` against 24.64.

**One-line reason at the first layer [Proved].** By Stein's identity, the signal of any feature is `delta a_1`, and its rank-one noise is `a_1 Delta_m`. Their ratio is `R_m` for every feature. Kernel choice cannot change it.

**The half-mean test does not reach the cap at high `H`.** The efficiency against `I_n` is 0.98, 0.91, 0.85, 0.77, 0.67, 0.54, 0.38 at `H` = 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95. The rate is the same. The constant is not. A contrast with weights near the change point gains up to a factor of 1.6 in signal-to-noise at `H = 0.95`, if `H` and the location are known.

**Where the impossibility stops.** The brief's general wording is **[Refuted]** in two ways.
1. *SV level shift [Proved].* `E[X_t | xi] = 0`, so the level of `X` is a martingale difference. A half-mean test detects at `delta ~ n^(-1/2)`. The memory is not in that layer.
2. *Scale shift [Proved].* The Fisher information for `(1 + delta 1{t>m}) xi_t` is `m + tr(Gamma^-1 E Gamma E) >= m`. The exact value at `n = 800` is `I/n = 1.00` to `1.01`. So information allows `n^(-1/2)` at every `H`. I do not claim a kernel achieves it.

The exact statement is this. **The cap holds for a level shift in the variable whose rank-one layer carries the memory.**

## 8. C6 and C7: seasonality and the shift rule

### C6. Lag `P` under a season of period `P`

**Lag `P` removes the seasonal mean embedding exactly [Proved].** If `P` divides `L`, the conditional mean of `phi(X_t) - phi(X_(t-L))` is 0 for any phase, including a shared fixed phase. This holds for a seasonal mean and for a seasonal scale. Hermite coefficients then depend on `t mod P`, and they are equal at `t` and `t-L`. So Lemma 1 applies, with Cauchy-Schwarz on the periodic coefficients, and `E T_n` stays bounded. It also removes the pole at frequency zero. At `H = 0.8`, `E T_n` moves from 0.934 to 0.938 under SEAS-M, and from 1.310 to 1.315 under SEAS-V, between `n = 2000` and 4000.

**The brief's "removes exactly" needs one correction [Refuted for the null level].** Lag `P` removes the mean, not the pairing of two scale patterns. Under SEAS-V with lag 50 and `n = 2000`, `E T_n` is 1.148 at gap 0 and 1.038 at gap 12 for `H = 0.5`. It is 1.310 and 1.195 for `H = 0.8`. This difference is bounded and does not grow. A shift reference that averages over gaps sees a different mean when the true phase is shared and fixed.

### C6. Lag 1 under a season of period 50

**Attenuation [Proved].** Per series, lag 1 multiplies harmonic `k` by `2 sin(pi k / P)` in amplitude. The brief's "(2 pi/P)^2 per series" is correct in power. In `T_n` the factor is `(2 sin(pi/P))^4 = 2.487e-4`. The exact ratio at `a = 1` is 2.643e-4. It sits 6% above, because lag 1 favours the higher harmonics. For seasonal scale it is 2.61e-4.

**Size tends to 1 for fixed `P` [Sketch].** `h_1(d) > 0` at every gap, with minimum 7.06e-8 at `d = 38`. The ergodic lemma in the paper's proofs (`lem:ergodic`) applies with a window of length 1, so `T_n^D / n -> h_1(D) > 0`. Any rule with a cutoff of order `o(n)` then has size tending to 1. The exact slope of `E T_n` in `n` matches `h_1(0) = 1.551e-6` with ratios 1.02 and 1.05.

**It matters only at huge `n`.** At `n = 2000`, gap 0 adds 0.003 to a null mean of 1.17. The seasonal term reaches the null mean near `n = 7.5e5` for `H = 0.5` and `n = 3.2e5` for `H = 0.8`, at gap 0. At the smallest gap it needs about 20 times more.

### C6. Slow drift, SEAS-V-slow, period `4n`

Here I expected lag 1 to cure the problem. It removes most of it. It does not remove all of it.

| n = 2000, 16 phase pairs | E T_n(0) / mean over offsets of E T_n(s): min | median | max |
|---|---|---|---|
| H = 0.5, raw | 1.032 | 1.173 | 1.737 |
| H = 0.5, lag 1 | 0.927 | 1.001 | 1.078 |
| H = 0.8, raw | 1.056 | 1.167 | 1.449 |
| H = 0.8, lag 1 | 0.887 | 1.001 | 1.127 |

**Raw [Exact].** The ratio is above 1 at every phase pair, with a median near 1.17. This is a one-sided bias, and it fits the paper's 10.1%. I did not compute the size.

**Lag 1 [Exact].** The median is 1.001. The one-sided bias has gone. A spread of about plus or minus 10% remains, from the pairing of two scale patterns. I computed the ratio only at `n = 2000`. `E T_n(0)` at lag 1 and gap 0 is the same at `n = 2000` and 4000 (1.7677 and 1.7678 at `H = 0.5`), so I expect the spread to persist. So **the paper's "no correction" for drift may ease on inc1, but I cannot say by how much.** EXP7 measures it.

### C7. The shift rule

**Review check.** I opened Harris and Liang et al.

**Harris 2020 [Proved there, checked].** `Y` is stationary and independent of `X`. For any measurable `V` on trimmed segments, and `s` from `-N` to `N`, Harris's own count `m` includes `s = 0`, and it satisfies `P(m <= M) <= M/(N+1)`. The review prints `P(m < M)` without saying whether `m` includes `s = 0`. With Harris's own `m`, that form would be vacuous at `M = 1`, and Harris's example rejects at `m = 1` with `N = 19`. This correction applies to Harris's own `m` only. EXP7 amendment A3 counts `s != 0`, so its `m < M` is the same event and is correct.

What it gives:
- It holds at finite `n` for raw and for increments, at every `H`, for any `g`. It needs no limit theory and works with a data-driven bandwidth. So **a conservative shift rule already has proved size for every arm that is a function of the segments.**
- Theorem 3 of Harris gives the rule with twice the power, but only as `N` tends to infinity with `D` fixed. That leaves the test with no growing sample. With `D` proportional to `T`, it does not apply.
- **New corollary [Proved].** If `Y` is `P`-periodic in law and independent of `X`, then `V_(kP)` is stationary in `k`. So Harris's bound holds for shifts on multiples of `P`, including a shared fixed phase. If the phase of `Y` is uniform, as in SEAS-M, SEAS-V and SEAS-V-slow, `Y` is stationary and the plain bound applies.
- Increments help power here, not size. Under `Y = rho X + ...`, the raw target at shift 1 keeps `R^2` of its height. That is 0.27 at `H = 0.8` and 0.75 at `H = 0.95`. The increment target has a sharp peak, with `r_D(1)^2` equal to 0.121 at `H = 0.8` and 0.064 at `H = 0.95`. Harris says the same in words. I did not compute the power.

**Liang et al. 2025 [checked].** Theorem 1 covers a local window statistic with inflation near 4 at mid-sample, under stationarity alone. For whole-series wrap-around, Section 4.2 says a guarantee is not possible without further assumptions. Theorem 3 bounds the size by `alpha + min{2 beta_X(tau) + 2 beta_Y(tau)}`. fGn is not beta-mixing (audit M1). The review says differencing does not restore mixing. I did not recheck that calculation. So Theorem 3 gives nothing for raw or for increments.

**Neither source covers the circular rule the paper uses.** The next two results do, for fixed bandwidth.

**Theorem D: the circular rule on increments is asymptotically exact [Proved, preprint].** Let `S_n` be the integer offsets in `[eps n, (1-eps) n]`. Under the assumptions of Theorem B, the `p`-value tends to uniform on `(0,1)`.

*Proof.*
1. *Same law at each offset.* The covariance of chaos `r` changes by at most `(1/n) sum_j min(j, n-j) |g^X(j)| |g^Y(j) - g^Y(n-j)|`, which is `O(1/n)`. Exact check: `E T_n(s) / E T_n(0) - 1` is about `-1e-4` at `n = 2000` and `-4e-6` at `n = 64000`.
2. *Decoupling.* For offsets `s, s'` with circular distance `m`, the cross-covariance trace is at most `C <m>^(-min(a,b))`, with `a = 4-2H_X`, `b = 4-2H_Y`. Exact check: the ratio is `-5e-12` at `n = 2000` and `-4e-17` at `n = 64000`, at `s = n/2`, `H = 0.8`.
3. *Joint limit.* Apply Theorem 3.4 of 2509.20342 on `H x H`. Condition (iii) holds uniformly, because X and Y factors are both bounded by functions of circular distance. So the Young bound of Theorem B holds on the group `Z_n`. Then `(T_n(s_n), T_n(s'_n))` converges to two iid copies of `T`.
4. *Self-averaging.* So `(1/|S_n|) sum_s 1{T_n(s) <= x}` has mean tending to `F(x)` and variance tending to 0. `F` is continuous, because `Lambda > 0` gives a nonatomic weighted chi-square law. Polya's theorem gives uniform convergence. Then `p_n = 1 - F(T_n(0)) + o_p(1)`. QED.

Both missing steps in Remark `rem:shift` close here. Uniform convergence comes from steps 1 to 3. Ties do not arise, because there is no order-`n` term with discrete values when no season is present.

**Theorem D applies without a season only.** With a season, the limit depends on the gap. I did not prove the rule with seasonal terms for lag `P` or lag 1.

**Raw below the line [Proved, preprint].** For `H_X + H_Y < 3/2`, the same steps work, because `D_X + D_Y > 1` gives decay `<m>^(1 - D_X - D_Y)` and a wrap effect of order `n^(1 - D_X - D_Y)`. So the paper's rule is asymptotically exact below the line, at a fixed bandwidth.

**Raw above the line [Sketch, numbers exact].** The argument fails, and exact second moments show why.

| H_X, H_Y | E T_n(n/2)/E T_n(0), n=2000 | n=64000 | centred cross ratio at s=n/2, n=2000 | n=64000 |
|---|---|---|---|---|
| 0.7, 0.7 | 0.9991 | 0.9996 | -0.016 | -0.006 |
| 0.8, 0.8 | 0.9910 | 0.9894 | -0.089 | -0.096 |
| 0.9, 0.9 | 0.9626 | 0.9525 | -0.263 | -0.310 |
| 0.95, 0.95 | 0.9395 | 0.9308 | -0.364 | -0.394 |
| 0.7, 0.9 | 0.9921 | 0.9911 | -0.091 | -0.092 |

Below the line (0.7, 0.7), both columns approach their no-memory values. Above the line they do not. **My own first guess was that wrap-around bias drives the failure. That is [Refuted] as the main cause.** The mean bias is small: 1% at `H = 0.8`. The larger effect is that shifted statistics stay correlated across distant offsets. So the reference distribution does not self-average, and its limit is random **[Sketch]**. This fits the paper's slow rise in rate at memory sum 1.6. I have not computed the limit size.

**What is still open for the paper's rule.** The rule with median bandwidth. Any size statement above the line. The rule with a seasonal term. Finite-sample accuracy: below the line, contraction falls only like `n^(-0.2)` for mixed pairs near the line.

## 9. C8: the Hermite-split design

**Design.** Map each series to normal scores. Remove low Hermite layers from the Gram matrix. Test those low layers separately with self-normalised long-memory tests. Combine the results.

### The prior rough check

The statistic is `n` times split MMD squared on scores, with `m = r + 1` removed terms. Inflation is the mean at `H` over the mean at `H = 0.5`, at `n = 800`. "Earlier repeat" is an earlier, noisier rerun of mine with 600 replications, against the note's 150. EXP7's 2,000-replication repeat is the more precise one. "Exact" takes `T` known and uses no simulation.

| Data, m | Note, H=0.7, 0.8 | Earlier repeat, H=0.7 (se) | Earlier repeat, H=0.8 (se) | Exact, H=0.7, 0.8 | Exact, H=0.8, large n |
|---|---|---|---|---|---|
| fGn, 1 | 3.0, 5.9 | 3.23 (0.17) | 7.87 (0.47) | 3.12, 6.91 | 78.2 at n=64000 |
| fGn, 2 | 1.1, 1.3 | 1.08 (0.05) | 1.42 (0.07) | 1.18, 2.14 | 4.26 |
| fGn, 3 | 1.1, 1.1 | 1.03 (0.05) | 1.13 (0.06) | 1.03, 1.28 | 1.35 |
| fGn, 4 | 1.0, 1.1 | 0.93 (0.05) | 0.93 (0.05) | 1.01, 1.09 | 1.10 |
| SV, 1 | 1.4, 2.4 | 1.37 (0.06) | 2.00 (0.10) | 1.41, 2.08 | 83.1 at n=1,024,000 |
| SV, 2 | 1.7, 3.1 | 1.70 (0.09) | 2.81 (0.15) | 1.60, 2.55 | 119.5 |
| SV, 3 | 1.2, 1.2 | 1.08 (0.06) | 1.16 (0.06) | 1.05, 1.16 | 7.90 |
| SV, 4 | 1.0, 1.2 | 1.14 (0.05) | 1.19 (0.06) | 1.06, 1.22 | 10.22 |

**Verdict on "x5.9 falls to x1.1".** The direction is right and the top number is low.
- At `m = 1`, the exact value with `T` known is x6.91. EXP7's 2,000-replication repeat gives x6.65 [6.22, 7.09]. The note's x5.9 had 150 replications. An earlier, noisier run of mine gave x7.9 (600 replications, se 0.5).
- The fall to about x1.1 holds at `m = 3` and `m = 4`, at `n = 800`.
- At `m = 2`, the note's x1.3 and my x1.42 reflect estimated scores. With `T` known the value is x2.14. The second-layer ratio in Section 4 explains the direction. That ratio is my explanation and I did not separate it from other estimation effects.
- For `m = 3` at `H = 0.8` the inflation is bounded, because `3D = 1.2 > 1` **[Proved]**. For `m = 1, 2` it grows without bound **[Proved, T known]**.
- As the note says, this is not a size result.

### Which layer holds the SV memory [Proved by exact quadrature]

Correlation of `He_p(xi)` with `He_q(T(X))`, both normalised:

| p | q=1 | q=2 | q=3 | q=4 | q=5 | q=6 | q=7 | q=8 |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | +0.488 | 0 | -0.113 | 0 | +0.051 | 0 | -0.029 |
| 2 | 0 | +0.199 | 0 | +0.214 | 0 | -0.107 | 0 | +0.067 |
| 3 | 0 | +0.049 | 0 | +0.222 | 0 | +0.032 | 0 | -0.040 |

- Odd score layers carry no memory, exactly.
- **The second score layer holds 93.6% of the rank-one memory.** The fourth holds 5.0%, the sixth 1.0%, and the eighth 0.3%. So the note's reading is right about where most of it sits.
- It is not all there. With the kernel on scores, removing `He_2(T)` leaves 4.0% of the rank-one weight. So the rate stays `n^(1-2D)`, and inflation at `m = 3` reaches x7.9 at `n` about `1e6` for `H = 0.8`.
- Scores recover only part of the driver: `R^2 = 0.254` across layers up to 8.
- Removing an odd layer raises inflation. It lowers the lag-0 variance and leaves the memory alone. So inflation goes up from `m = 1` to `m = 2`, and from `m = 3` to `m = 4`. Both exact values and the repeat show this.

### Does separate testing recover the power that deflation loses?

**Yes for first-layer signals, at the raw rate [Sketch at the level of a test].** The numbers below use an oracle scale. They exclude the cost of self-normalisation, which I did not compute.

*Location.* The rank contrast keeps the signal: `E_A T - E_B T` over `delta` is 0.999, 0.970, 0.893, 0.692 at `delta` = 0.1, 0.5, 1, 2. By Stein's identity its signal-to-noise ratio is the paper's `R_m` to first order. The table gives the shift at which signal equals noise scale at `n = 2000`. The first column uses the standard deviation of a one-degree contrast. The others use the null mean of `T`, so they are not on the same footing.

| H | first layer | raw | deflated r=1 | deflated r=2 |
|---|---|---|---|---|
| 0.5 | 0.045 | 0.081 | 0.799 | 0.758 |
| 0.8 | 0.247 | 0.269 | 0.928 | 0.792 |
| 0.9 | 0.361 | 0.435 | 1.227 | 1.014 |
| 0.95 | 0.366 | 0.544 | 1.452 | 1.294 |

The local regime is `delta` of order `n^(H-1)`. There the first layer has power, up to the cap of Section 7, and the deflated remainder has none. Plain deflation needs `delta ~ n^(-1/6)` or `n^(-(2/3)(1-H))`. No layer test can pass the cap.

*Independence.* The first-layer Pearson on scores has `z = sqrt(n) rho / sqrt(S)`. For `H_X + H_Y > 3/2`, `S` grows like `n^(2(H_X+H_Y)-3)`. So detection is at `rho ~ n^(-(2-H_X-H_Y))`, which equals the raw rate. Plain deflation needs `n^(-1/(2(r+1)))`. Between those two scales, separate testing recovers the signal and plain deflation does not. The same `S` for the increment kernel stays bounded, so the increment rate is `n^(-1/2)`.

| H | S first layer, n=2000 | S first layer, n=64000 | S increment |
|---|---|---|---|
| 0.6 | 1.07 | 1.08 | 1.42 |
| 0.7 | 1.49 | 1.71 | 1.34 |
| 0.8 | 3.21 | 7.39 | 1.26 |
| 0.9 | 7.33 | 59.3 | 1.19 |
| 0.95 | 6.40 | 102.6 | 1.16 |

### Is the combination valid when the parts are dependent?

**Bonferroni is valid only if each part is valid [Proved].** If each part has valid size, a Bonferroni combination is valid under any dependence. It costs only the larger quantile. No hsplit part has a proof, so this result does not give hsplit1 or hsplit2 valid size.

**Fisher or a chi-square sum is not justified above `H = 3/4` [Proved for known centre].** Exact cumulants give `Corr(Delta_1^2, Delta_2^2)` for the half contrasts of `He_1` and `He_2`. Scores behave like the centred row.

| H | known centre, n=500 | n=2000 | sample centre, n=500 | n=2000 |
|---|---|---|---|---|
| 0.6 | 0.033 | 0.014 | 0.0003 | 0.0001 |
| 0.7 | 0.179 | 0.133 | 0.0031 | 0.0021 |
| 0.8 | 0.397 | 0.386 | 0.0142 | 0.0128 |
| 0.9 | 0.480 | 0.479 | 0.0363 | 0.0357 |
| 0.95 | 0.493 | 0.493 | 0.0507 | 0.0504 |

With a known centre, the correlation is large above `H = 3/4` and does not fall with `n`. A sample centre shrinks it to about 1% to 5%, and it still barely moves with `n` above `H = 3/4`. So the two parts stay dependent there, but weakly on this measure. Zero correlation of squares would not prove independence. An independent-combination rule on rank scores may work in practice **[Conjecture]**.

**Remainder against low layers [Sketch].** If `r >= r*`, the remainder has a Gaussian limit and lives in higher chaos. Peccati-Tudor reasoning gives asymptotic independence from a Gaussian first layer. If `r < r*`, both limits are non-central functionals of one random measure, so they are dependent. Bai-Taqqu 2013 reportedly covers the scalar case. I did not open it.

**The marginals are the real gap. I have no proof that any self-normalised low-layer test holds size [Open].**
- The first-layer limit is a functional of fractional Brownian motion, so its law depends on `H`.
- The second-layer limit is a Rosenblatt-type law above `H = 3/4`. Bai-Taqqu Theorem 3.4 says the sample mean changes it. The sample scale in rank scores is not covered.
- The paper already reports 6.0% to 8.5% for self-normalised blocks on split MMD.

**A valid version exists for independence [Proved].** Apply Harris's bound to each layer statistic and to the remainder, and combine with Bonferroni. No self-normalisation and no value of `H` is needed. It is conservative. Split MMD has one series, so this route does not apply there. This version is not an EXP7 arm, because A3 applies Harris to raw and inc1 only.

### Head to head with the increment kernel

| Case | Hermite split | Increment lag 1 | Status |
|---|---|---|---|
| Level shift in a split | First layer sees it at `n^(H-1)`, up to the cap | Blind, telescopes | Proved for increment, sketch for split |
| Scale shift in a split | Second layer sees it | Blind, telescopes | Same |
| Equal-memory correlation, `H = 0.8`, `n = 64000` | `S = 7.39`, rate `n^(-0.4)` | `S = 1.26`, rate `n^(-1/2)` | Exact S, oracle scale |
| Common factor, white idiosyncratic noise, `s = 0.3`, `H_F = 0.9`, `n = 2000` | `z = 11.0` | `z = 3.7`, correlation 0.10 | Exact, oracle scale |
| Same with `H_F = 0.99` | `z = 17.1` | `z = 0.43`, correlation 0.012 | Exact, oracle scale |
| Common SV factor, `H = 0.9`, `n = 2000`, layer (2,2) | `z = 10.8` | `z = 3.7` | Exact, oracle scale |
| Flat three-lag dependence | Sees it | Exactly blind | Proved |
| Lag-one dependence, `r(0) = 0` | Blind for a Gaussian pair | Sees it | Proved |
| Seasonal mean | Size tends to 1, as in paper Proposition 1 | Lag `P` removes it. Lag 1 shrinks it by `2.6e-4` | Sketch for split, proved for lag `P` |
| SV null | Memory falls by 96% to 97% at `r = 2`. Rate stays | Bounded `E T_n` | Proved |
| Non-monotone `g` | Partial removal. Rate stays | Bounded `E T_n` | Proved |
| Memory pole away from frequency zero | Not studied | Lag 1 fails. Matching lag works | Exact |
| Null limit | Depends on `H`. Non-Gaussian above `H = 3/4`. Changes with estimated scores | Gaussian in Hilbert-Schmidt norm | Proved, preprint, for increment |
| Size proof | Only Harris on each layer plus Bonferroni, for independence, which is not an EXP7 arm. None for hsplit1 or hsplit2 | Circular shift, fixed bandwidth, plus Harris | Proved, preprint |

The dividing line is where the dependence lives. The increment correlation equals the share of the dependence at high frequency. If idiosyncratic noise is as persistent as the shared part, increments lose little. If the shared part is slow and the rest is fast, increments lose nearly everything, and the first layer keeps it. **They answer different questions. Neither replaces the other. Report both.** The Granger-Newbold advice on levels and changes, as the review quotes it, says the same.

## 10. Refuted and unverified claims, by arm

### Refuted or weakened

1. **"Deflation removes nothing under SV."** True for odd `r` only. `r = 2` removes 96% to 97%. Arms: herm1 keeps full memory, and herm2 looks cured at finite `n` while the rate stays.
2. **"The kernel loses `r` directions."** False for pooled scores, where the blind set is empty. It holds for a fixed reference. Arms: herm1, herm2.
3. **"Estimating `F` may spoil deflation."** The path check shows lower means, with ratios 0.34 to 1.03. The limit law does change. Arms: herm1, herm2.
4. **"`E T_n^D = O(1)` for any process."** False for memory poles away from frequency zero. True for fGn-type spectra. Arms: inc1, inc50, dcov_inc1.
5. **"The oracle filter gives sum 1."** It gives 1.0026 to 1.0371 on fGn. Arm: frac_oracle benchmark.
6. **"Memory cannot be filtered out of split MMD."** Holds for a level shift in the memory layer. False as a general statement: SV level shifts and scale shifts have no such cap. Arms: split-MMD arms.
7. **"(2 pi/P)^2 per series" read as the `T_n` factor.** `T_n` scales with the fourth power. Arm: inc1 under seasonal cells.
8. **"Lag `P` removes the season exactly."** True for the mean embedding. False for the phase-dependent null level under seasonal scale. Arm: inc50 under SEAS-V with a shared phase.
9. **My guess that lag 1 cures slow drift.** It removes the one-sided bias, and a spread of about 10% remains. Arm: inc1 under SEAS-V-slow.
10. **My guess that wrap bias drives the raw shift failure.** The mean bias is 1% at `H = 0.8`. The non-vanishing correlation across offsets is larger. Arm: raw with circular shifts.
11. **The note's "x5.9".** The exact value with `T` known is x6.91. EXP7's 2,000-replication repeat gives x6.65 [6.22, 7.09]. An earlier, noisier run of mine gave x7.9 (600 replications, se 0.5). The note's x1.3 at `m = 2` needs estimated scores. Arms: herm1, herm2.
12. **Review items.** With Harris's own `m`, which includes `s = 0`, the inequality needs `<=`. A3 counts `s != 0`, so its `m < M` is correct. The rank-two rate holds only for an unbalanced split. Theorem 3.1 does not cover extras or estimated scores. Arms: herm1, herm2 and split-MMD arms.

### Could not verify

1. **Median bandwidth.** All limit and shift results are for a fixed bandwidth. The extension is a sketch. All arms use the median. Harris's bound is unaffected if the bandwidth is computed from the segments.
2. **Theorem 3.4 of arXiv 2509.20342v3.** I read the statement and part of its proof. It is an unreviewed preprint. Theorems B, C and D rest on it. Arms: inc1, inc50, raw.
3. **My contraction bound.** It has no independent check beyond one exact second-chaos family and one fourth-chaos case at `n <= 36`. A referee should check it. Same arms.
4. **Finite-sample size of any rule.** I give means and second moments only. All arms.
5. **Wild multipliers on increments.** No theory here. Increments have lag-one correlation near -1/2 at `H = 0.5`, so a positive-only block scheme may misbehave. Arms: inc1, inc50 with wild rules.
6. **`c_1 > 0` for every `r` under SV.** Checked up to `r = 6`. Arm: herm2.
7. **Deflation with estimated scores.** Path check only, with 60 replications. The `H = 0.9, n = 4000` cells are missing because the job hit its time limit. Arms: herm1, herm2.
8. **frac_oracle as implemented.** My exact values use a stationary filter with the whole past. EXP7 truncates at the sample start. Arm: frac_oracle.
9. **sym, hp10 and dcov_inc1.** I did no exact calculation for sym or hp10. For dcov_inc1 the expectation result is a sketch, and the limit law does not cover it. Arms: sym, hp10, dcov_inc1.
10. **Size-adjusted power.** All power statements are oracle `z` or signal-to-null ratios. They exclude self-normalisation. Arms: all power cells.
11. **Self-normalised low-layer tests.** No size proof. Their cost in power is not computed. EXP7 amendment A2 measures this design as arms hsplit1 and hsplit2. It has a plug-in table rule, block-sampling rules and a pre-registered Bonferroni split. A2 itself says no cited theorem covers its low-layer statistics under long memory. So those arms give evidence on size, not a proof. A2 compares with this document as follows.
    - *Bonferroni form: no mismatch.* A2 gives the remainder `alpha/2` and splits `alpha/2` equally over the low-layer statistics. The levels sum to 5%, and the code uses those levels. For split MMD, the remainder uses the wild rule and no shift rule, which fits Section 9.
    - *Bonferroni validity: mismatch.* Section 9 proves validity only if each part is valid. No hsplit part has that proof. Theorem D does not cover the deflated remainder under the shift rule, and no theorem covers the low-layer table.
    - *The proved version is not an arm: mismatch.* Section 9's proved version applies Harris to each layer and to the remainder. A3 applies Harris to raw and inc1 only.
    - *Correlation table: mismatch.* Section 9's table uses squared half contrasts. A2 uses self-normalised `G^2` and, for HSIC, products. Bonferroni needs no such table, so validity is unaffected. The table does not describe the dependence among A2's statistics.
    - *SV generator: mismatch.* My exact work uses `exp(xi/2) e` with fGn `xi`, which is EXP7's SV-hcheck. P6 reads NL-SV, whose log-volatility is AR(0.9)-filtered fGn rescaled to unit variance. The marginal law is the same. So the layer correlations, the 93.6% share, the zero odd layers and the 2.87% carry over. `E T_n`, inflation values and slopes do not, because they use the fGn `rho`. Lemma 1 and Theorem A also assume fGn. On EXP7's own repeat at `m = 3, H = 0.8`, inflation is x1.23 to x1.26 on NL-SV and x1.05 to x1.06 on SV-hcheck.
    - *"hsplit2 removes the SV memory": mismatch.* P6 says it does. With A2's remainder, which keeps K on raw `x`, `r = 2` leaves 2.87% of `c_1`, and the rate stays. The 4.0% in Section 9 is for K on scores. So P6 can hold only at finite `n`. The memory also moves into the second low layer, and it does not leave the test.
    - *Table law under SV: A2 flags this as an approximation, and I add the direction.* Under SV the memory in `He_2(T)` is rank one in the driver, so its limit type differs from the rank-two law in `split2` **[Sketch]**. About 72% of the variance of `He_2(T)` is white noise, so local Whittle may read `H-hat` too low **[Conjecture]**. A low reading would pick a low cut-off and could raise the low-layer rate on NL-SV.
12. **The review's mixing calculation for differenced fGn.** I did not recheck it. It matters only for the claim that Liang et al. Theorem 3 gives nothing for increments.

## 11. Sources and how to reproduce

**Opened by me.**
- Duker-Zoubouloglou arXiv 2405.11452v3: Assumptions 3.1, 3.2, Theorem 3.1.
- Duker-Zoubouloglou arXiv 2509.20342v3: Theorems 3.3, 3.4 and part of the proof of 3.4.
- Harris arXiv 2012.06862: full text.
- Liang et al. arXiv 2512.17372v2: Theorems 1, 3, Section 4.2.
- Bai-Taqqu arXiv 1710.01612: Theorem 3.4, Proposition 4.1, the definition of power rank.
- Betken-Dehling arXiv 2107.03041v3: Remark 2.1.

**Not checked.** Dehling-Taqqu 1989. Breuer-Major 1983. Arcones 1994. Bai-Taqqu 2013. Peccati-Tudor. Bardet-Surgailis 2011. Nourdin-Peccati's book. Granger-Newbold 1974, which I cite only through the review. The product-moment contraction identity and Mehler's formula are standard, and I did not open a source for them.

**Run.** All from this folder.

```
OMP_NUM_THREADS=8 python -I theory_checks.py > theory_checks.out
OMP_NUM_THREADS=2 python -I theory_path_check.py > theory_path_check.out
OMP_NUM_THREADS=2 python -I theory_path_check.py --note > theory_note_repro.out
```

`theory_checks.py` runs 74 checks. All 74 pass, in 222 seconds on the final run. An earlier run with other jobs on the machine took 630 seconds, so run time depends on load. `theory_exact.csv` has 2988 rows. `theory_exponents.csv` has 211 fits, and all pass. Tolerances are in the script. I chose some of them after seeing a first run, and the bracket rule at `H = 0.95` is one of those.
