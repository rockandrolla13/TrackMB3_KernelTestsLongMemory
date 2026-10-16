"""DML: AIPW estimate of the ATE with a long-memory state.

Design: X_t from G1 (fGn, H grid); D_t ~ Bernoulli(logit^-1(0.5 X_t)); Y_t = tau(X_t) D_t + X_t + eps_t,
eps i.i.d. N(0,1). Effects: (i) tau = 1; (ii) tau = 1 + X; (iii) tau = 1 + (X^2 - 1).
Population ATE = 1 in all three. Path-conditional ATE = mean_t tau(X_t). T in {500, 2000, 8000}.

Arms (all computed on the same simulated path):
  oracle        known propensity e(x) = logit^-1(0.5 x) and true outcome regressions
                mu0(x) = x, mu1(x) = x + tau(x). AIPW score = (tau(X_t) - theta) + martingale
                difference exactly, so this arm isolates the state term.
  kp_iid        known propensity; outcome regressions by gradient boosting, 5-fold cross-fitting,
                folds assigned at random (as-i.i.d.).
  kp_blk        known propensity; gradient boosting; 5 contiguous folds; training set for fold k
                excludes observations within b = floor(T^0.5) of fold k.
  lrn_iid       learned propensity (logistic regression, numpy IRLS, features [1, x]) and gradient
                boosting outcome regressions; as-i.i.d. folds.
  lrn_blk       as lrn_iid with blocked folds and buffer b.
Gradient boosting (numpy, histogram trees on the scalar state; squared loss): quantile bins of
the training-fold state (255 bins, as in histogram GB defaults), depth-2 regression trees with splits on bin edges, min 20 observations
per leaf, learning rate 0.1, 100 rounds, initial value = training mean. mu1 is fit on treated and
mu0 on control units of the training fold. Learned propensities clipped to [0.01, 0.99]
(count recorded). IRLS: up to 50 Newton steps, converged if max |step| < 1e-8 (failures recorded,
final iterate used).
Variance of the mean score: Newey-West Bartlett with Andrews AR(1) bandwidth (normal critical
value) and EWC with nu = floor(0.4 T^{2/3}) (t_nu critical value). The same interval is scored
against both estimands.
"""
import math

import numpy as np

import common as C

T_GRID = (500, 2000, 8000)
TAUS = ("const", "rank1", "rank2")
K_FOLDS = 5
GB_BINS = 255
GB_ROUNDS = 100
GB_LR = 0.1
GB_MIN_LEAF = 20
CLIP = 0.01
ARMS = ("oracle", "kp_iid", "kp_blk", "lrn_iid", "lrn_blk")
COLUMNS = (("theta_path", "treated_share", "irls_nonconv_iid", "n_clip_iid", "irls_nonconv_blk", "n_clip_blk")
           + tuple(f"{a}_{k}" for a in ARMS for k in ("est", "se_nw", "se_ewc", "bw_nw")))
SIZE_KEY = "T"


def cell_list():
    cells = []
    for tau in TAUS:
        for h in C.HURST_GRID:
            for t in T_GRID:
                cells.append({"tau": tau, "H": h, "T": t})
    for i, c in enumerate(cells):
        c["cell_idx"] = i
        c["cell_id"] = f"dml_{c['tau']}_H{c['H']}_T{c['T']}"
    return cells


def tau_fn(kind, x):
    if kind == "const":
        return np.ones_like(x)
    if kind == "rank1":
        return 1.0 + x
    if kind == "rank2":
        return 1.0 + (x * x - 1.0)
    raise ValueError(kind)


def expit(z):
    return 1.0 / (1.0 + np.exp(-z))


def logistic_irls(x, d, max_iter=50, tol=1e-8):
    z = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2)
    for _ in range(max_iter):
        p = expit(z @ beta)
        w = p * (1.0 - p)
        g = z.T @ (d - p)
        hss = (z * w[:, None]).T @ z
        step = np.linalg.solve(hss, g)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            return beta, True
    return beta, False


def _best_split(cnt, sm, lo, hi):
    """Best split of bins [lo, hi) into [lo, j) and [j, hi); returns (gain, j) or (0, None)."""
    c = np.cumsum(cnt[lo:hi])
    s = np.cumsum(sm[lo:hi])
    ct, st = c[-1], s[-1]
    cl, sl = c[:-1], s[:-1]
    cr, sr = ct - cl, st - sl
    ok = (cl >= GB_MIN_LEAF) & (cr >= GB_MIN_LEAF)
    if not ok.any():
        return 0.0, None
    with np.errstate(divide="ignore", invalid="ignore"):
        gain = np.where(ok, sl * sl / cl + sr * sr / cr - st * st / ct, -np.inf)
    j = int(np.argmax(gain))
    if gain[j] <= 0:
        return 0.0, None
    return float(gain[j]), lo + j + 1


def gb_fit_predict(xtr, ytr, xte):
    """Histogram gradient boosting, squared loss, scalar feature. Exact on binned data: every tree
    splits on bin edges, so the ensemble is a function of the bin and only per-bin counts and sums
    are needed."""
    edges = np.unique(np.quantile(xtr, np.arange(1, GB_BINS) / GB_BINS))
    btr = np.searchsorted(edges, xtr, side="right")
    nb = len(edges) + 1
    cnt = np.bincount(btr, minlength=nb).astype(float)
    ysum = np.bincount(btr, weights=ytr, minlength=nb)
    g = np.full(nb, ytr.mean())
    for _ in range(GB_ROUNDS):
        res = ysum - cnt * g
        upd = np.zeros(nb)
        _, j = _best_split(cnt, res, 0, nb)
        leaves = [(0, nb)] if j is None else []
        if j is not None:
            for lo, hi in ((0, j), (j, nb)):
                _, j2 = _best_split(cnt, res, lo, hi)
                leaves.extend([(lo, hi)] if j2 is None else [(lo, j2), (j2, hi)])
        for lo, hi in leaves:
            cc = cnt[lo:hi].sum()
            if cc > 0:
                upd[lo:hi] = res[lo:hi].sum() / cc
        g = g + GB_LR * upd
    return g[np.searchsorted(edges, xte, side="right")]


def folds(t, scheme, rng):
    if scheme == "iid":
        f = np.empty(t, dtype=int)
        f[rng.permutation(t)] = np.arange(t) % K_FOLDS
        return f, 0
    f = (np.arange(t) * K_FOLDS) // t
    return f, int(math.floor(math.sqrt(t)))


def cross_fit(x, d, y, scheme, rng):
    t = len(x)
    f, buf = folds(t, scheme, rng)
    idx = np.arange(t)
    e_hat = np.empty(t)
    m0 = np.empty(t)
    m1 = np.empty(t)
    nonconv = 0
    for k in range(K_FOLDS):
        te = f == k
        tr = ~te
        if buf:
            lo, hi = idx[te].min(), idx[te].max()
            tr &= (idx < lo - buf) | (idx > hi + buf)
        beta, ok = logistic_irls(x[tr], d[tr])
        nonconv += int(not ok)
        e_hat[te] = expit(beta[0] + beta[1] * x[te])
        t1 = tr & (d == 1)
        t0 = tr & (d == 0)
        m1[te] = gb_fit_predict(x[t1], y[t1], x[te])
        m0[te] = gb_fit_predict(x[t0], y[t0], x[te])
    n_clip = int(np.sum((e_hat < CLIP) | (e_hat > 1 - CLIP)))
    return np.clip(e_hat, CLIP, 1 - CLIP), m0, m1, nonconv, n_clip


def aipw(d, y, e, m0, m1):
    return m1 - m0 + d * (y - m1) / e - (1 - d) * (y - m0) / (1 - e)


def run_rep(cell, rng):
    t, h = cell["T"], cell["H"]
    x = C.fgn(h, t, rng)
    e = expit(0.5 * x)
    d = (rng.random(t) < e).astype(float)
    eps = rng.standard_normal(t)
    tau = tau_fn(cell["tau"], x)
    y = tau * d + x + eps
    out = {"theta_path": float(tau.mean()), "treated_share": float(d.mean())}
    scores = {"oracle": aipw(d, y, e, x, x + tau)}
    for scheme in ("iid", "blk"):
        e_hat, m0, m1, nonconv, n_clip = cross_fit(x, d, y, "iid" if scheme == "iid" else "blk", rng)
        scores[f"kp_{scheme}"] = aipw(d, y, e, m0, m1)
        scores[f"lrn_{scheme}"] = aipw(d, y, e_hat, m0, m1)
        out[f"irls_nonconv_{scheme}"] = nonconv
        out[f"n_clip_{scheme}"] = n_clip
    for arm in ARMS:
        psi = scores[arm]
        est = float(psi.mean())
        u = psi - est
        lrv_nw, s_nw = C.nw_lrv(u)
        lrv_ewc, _ = C.ewc_lrv(u)
        out[f"{arm}_est"] = est
        out[f"{arm}_se_nw"] = math.sqrt(max(lrv_nw, 0.0) / t)
        out[f"{arm}_se_ewc"] = math.sqrt(max(lrv_ewc, 0.0) / t)
        out[f"{arm}_bw_nw"] = s_nw
    return out
