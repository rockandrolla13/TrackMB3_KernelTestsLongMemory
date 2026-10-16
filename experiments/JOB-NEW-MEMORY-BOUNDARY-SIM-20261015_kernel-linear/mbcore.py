"""Core generators and tests for the kernel-linear arm of the memory breakdown map.

Job JOB-NEW-MEMORY-BOUNDARY-SIM-20261015, worker claude-mb-kernel-20261015.
numpy only. Every test is implemented here from its published definition.
"""
from __future__ import annotations

import math

import numpy as np

# ---------------------------------------------------------------- normal / t helpers

_SQRT2 = math.sqrt(2.0)
_erf = np.frompyfunc(math.erf, 1, 1)


def norm_cdf(x: np.ndarray) -> np.ndarray:
    """Standard normal CDF via math.erf (vectorised with frompyfunc)."""
    x = np.asarray(x, dtype=float)
    return 0.5 * (1.0 + np.asarray(_erf(x / _SQRT2), dtype=float))


def t_cdf(x: float, nu: int) -> float:
    """Student-t CDF by numerical integration of the density (Simpson, 20k panels)."""
    if x == 0.0:
        return 0.5
    c = math.exp(math.lgamma((nu + 1) / 2) - math.lgamma(nu / 2)) / math.sqrt(nu * math.pi)
    a = abs(x)
    m = 20000
    grid = np.linspace(0.0, a, m + 1)
    f = c * (1.0 + grid**2 / nu) ** (-(nu + 1) / 2)
    h = a / m
    integral = h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
    return 0.5 + integral if x > 0 else 0.5 - integral


def t_quantile(p: float, nu: int) -> float:
    """Student-t quantile by bisection on t_cdf."""
    lo, hi = 0.0, 50.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, nu) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- fGn (Davies-Harte)

def fgn_acov(H: float, n: int) -> np.ndarray:
    """Autocovariance of unit-variance fGn at lags 0..n-1."""
    k = np.arange(n, dtype=float)
    return 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))


_DH_CACHE: dict[tuple[float, int], np.ndarray] = {}


def dh_eigs(H: float, n: int) -> np.ndarray:
    """Eigenvalues of the minimal circulant embedding (size m = 2(n-1)) of the fGn Toeplitz matrix.

    Raises if any eigenvalue is negative beyond rounding (embedding would not be exact).
    """
    key = (H, n)
    if key not in _DH_CACHE:
        g = fgn_acov(H, n)
        c = np.concatenate([g, g[-2:0:-1]])  # length 2n-2
        lam = np.fft.fft(c).real
        if lam.min() < -1e-10 * lam.max():
            raise ValueError(f"Davies-Harte embedding not nonnegative: H={H}, n={n}, min={lam.min()}")
        _DH_CACHE[key] = np.clip(lam, 0.0, None)
    return _DH_CACHE[key]


def fgn(H: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Exact unit-variance fGn path of length n by Davies-Harte. Uses the real part only."""
    if n < 2:
        raise ValueError("n >= 2")
    lam = dh_eigs(H, n)
    m = lam.size
    z = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    v = np.fft.fft(np.sqrt(lam / m) * z)
    return v.real[:n].copy()


def fou(H: float, lam: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Euler fOU: x_t = (1 - lam) x_{t-1} + fGn_t(H), started at 0, burn-in ceil(10/lam)."""
    burn = int(math.ceil(10.0 / lam))
    e = fgn(H, n + burn, rng)
    x = np.empty(n + burn)
    prev = 0.0
    a = 1.0 - lam
    for t in range(n + burn):
        prev = a * prev + e[t]
        x[t] = prev
    return x[burn:]


def phi_transform(xi: np.ndarray) -> np.ndarray:
    """G2: X = 2 Phi(xi) - 1."""
    return 2.0 * norm_cdf(xi) - 1.0


# ---------------------------------------------------------------- process settings

def make_settings() -> list[dict]:
    """Fixed process settings (design table)."""
    s: list[dict] = []
    for proc in ("G1", "G2", "G3"):
        for H in (0.5, 0.6, 0.7, 0.75, 0.8, 0.9):
            s.append({"proc": proc, "HX": H, "HY": H, "lam": None})
    for hx, hy in ((0.6, 0.85), (0.7, 0.8), (0.7, 0.9), (0.6, 0.9)):
        s.append({"proc": "G4", "HX": hx, "HY": hy, "lam": None})
    for lam in (0.01, 0.1):
        s.append({"proc": "R1", "HX": 0.1, "HY": 0.1, "lam": lam})
    for st in s:
        st["sid"] = setting_id(st)
    return s


def setting_id(st: dict) -> str:
    if st["proc"] == "R1":
        return f"R1_lam{st['lam']}"
    if st["proc"] == "G4":
        return f"G4_HX{st['HX']}_HY{st['HY']}"
    return f"{st['proc']}_H{st['HX']}"


def gen_series(st: dict, n: int, rng: np.random.Generator, r: float = 0.0):
    """Generate (X, Y). Independent under the null (r = 0).

    Power check (r != 0): Betken-Dehling 'linear' design, Y-latent = r*xi + sqrt(1-r^2)*eta,
    which gives Cov(X_i, Y_j) = r Cov(X_i, X_j) on the latent scale.
    """
    proc = st["proc"]
    if proc == "R1":
        x = fou(st["HX"], st["lam"], n, rng)
        y = fou(st["HY"], st["lam"], n, rng)
        return x, y
    xi = fgn(st["HX"], n, rng)
    eta = fgn(st["HY"], n, rng)
    if r != 0.0:
        if st["HX"] != st["HY"]:
            raise ValueError("power design needs equal H")
        eta = r * xi + math.sqrt(1.0 - r * r) * eta
    if proc in ("G1", "G4"):
        return xi, eta
    if proc == "G2":
        return phi_transform(xi), phi_transform(eta)
    if proc == "G3":
        return np.sign(xi), np.sign(eta)
    raise ValueError(proc)


# ---------------------------------------------------------------- kernel helpers

def median_bandwidth(d: np.ndarray) -> float:
    """Median of off-diagonal pairwise distances from a full n x n distance matrix.

    The n diagonal zeros are the n smallest entries, so the off-diagonal median is read off
    the full sorted array after skipping them. If that median is 0 (ties, e.g. sign data),
    use the median of the strictly positive distances.
    """
    n = d.shape[0]
    flat = d.ravel()
    N = n * (n - 1)
    k1 = n + N // 2 - 1
    k2 = n + N // 2
    part = np.partition(flat, (k1, k2))
    med = 0.5 * (part[k1] + part[k2])
    if med <= 0.0:
        pos = flat[flat > 0]
        med = float(np.median(pos)) if pos.size else 1.0
    return float(med)


def gauss_gram(d: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-(d * d) / (2.0 * sigma * sigma))


def double_center(a: np.ndarray) -> np.ndarray:
    rm = a.mean(axis=1, keepdims=True)
    cm = a.mean(axis=0, keepdims=True)
    return a - rm - cm + a.mean()


def markov_rademacher(n: int, B: int, ln: float, rng: np.random.Generator) -> np.ndarray:
    """Rademacher block-Markov multipliers: W_1 = +-1, flip sign w.p. p = (1 - e^{-1/ln})/2.

    Then Corr(W_s, W_t) = e^{-|s-t|/ln}, the autocorrelation of the Gaussian AR wild-bootstrap
    process in Chwialkowski-Sejdinovic-Gretton (2014). Returned empirically centred (W tilde).
    """
    p = 0.5 * (1.0 - math.exp(-1.0 / ln))
    first = np.where(rng.random(B) < 0.5, -1.0, 1.0)
    flips = (rng.random((n - 1, B)) < p).astype(np.int8)
    parity = np.concatenate([np.zeros((1, B), np.int64), np.cumsum(flips, axis=0, dtype=np.int64)], axis=0)
    w = first[None, :] * (1.0 - 2.0 * (parity & 1))
    return w - w.mean(axis=0, keepdims=True)


def wild_boot_stats(m: np.ndarray, w: np.ndarray) -> np.ndarray:
    """(1/n) sum_ij W_i W_j M_ij for each column of w."""
    n = m.shape[0]
    return np.einsum("ib,ib->b", w, m @ w) / n


def shift_stats(kc: np.ndarray, lc: np.ndarray, chunk: int = 500) -> np.ndarray:
    """S(c) = (1/n) sum_ij Kc_ij Lc_{i+c, j+c} (indices mod n) for every circular shift c.

    Computed by FFT along circulant diagonals: S(c) = sum_d xcorr(Kd_d, Ld_d)(c), with
    Kd[d, i] = Kc[i, (i+d) mod n]. Centred Gram matrices commute with the shift permutation,
    so this equals n * HSIC between X and Y shifted by c.
    """
    n = kc.shape[0]
    i = np.arange(n)
    acc = np.zeros(n // 2 + 1, dtype=complex)
    for d0 in range(0, n, chunk):
        d = np.arange(d0, min(n, d0 + chunk))
        cols = (i[None, :] + d[:, None]) % n
        kd = kc[i[None, :], cols]
        ld = lc[i[None, :], cols]
        acc += (np.conj(np.fft.rfft(kd, axis=1)) * np.fft.rfft(ld, axis=1)).sum(axis=0)
    return np.fft.irfft(acc, n=n) / n


def _windows(x: np.ndarray, start: np.ndarray, l: int) -> np.ndarray:
    return x[start[:, None] + np.arange(l)[None, :]]


def subsample_stats(x: np.ndarray, y: np.ndarray, l: int, dn: int, kind: str,
                    sx: float | None = None, sy: float | None = None, chunk: int = 400) -> np.ndarray:
    """Betken-Dehling sampling-window statistics T_{l,k}, k = 0..m-1, m = n - l - dn.

    Block of X starting at k, block of Y starting at k + dn. kind 'hsic' uses the full-sample
    bandwidths sx, sy; kind 'dcov' uses distance covariance. T_l = l * V-statistic.
    """
    n = x.size
    m = n - l - dn
    out = np.empty(m)
    for k0 in range(0, m, chunk):
        ks = np.arange(k0, min(m, k0 + chunk))
        xw = _windows(x, ks, l)
        yw = _windows(y, ks + dn, l)
        dx = np.abs(xw[:, :, None] - xw[:, None, :])
        dy = np.abs(yw[:, :, None] - yw[:, None, :])
        if kind == "hsic":
            a = np.exp(-(dx * dx) / (2 * sx * sx))
            b = np.exp(-(dy * dy) / (2 * sy * sy))
        else:
            a, b = dx, dy
        bc = b - b.mean(axis=2, keepdims=True) - b.mean(axis=1, keepdims=True) + b.mean(axis=(1, 2), keepdims=True)
        if kind == "dcov":
            a = a - a.mean(axis=2, keepdims=True) - a.mean(axis=1, keepdims=True) + a.mean(axis=(1, 2), keepdims=True)
        out[ks] = np.einsum("kij,kij->k", a, bc) / l
    return out


# ---------------------------------------------------------------- K-ind test bundle

def kind_tests(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, B: int = 500) -> dict:
    """All K-ind tests on one (x, y) pair. Returns p-values and statistics.

    HSIC: wild bootstrap, shift test, subsampling. dCov: wild bootstrap, subsampling.
    """
    n = x.size
    ln = math.floor(n ** 0.5)  # wild-bootstrap ln and subsampling block length
    dn = math.floor(0.1 * n)
    A = math.floor(0.1 * n)  # shift range [A, n - A]

    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    sx = median_bandwidth(dx)
    sy = median_bandwidth(dy)

    kc = double_center(gauss_gram(dx, sx))
    lc = double_center(gauss_gram(dy, sy))
    mh = kc * lc
    t_h = mh.sum() / n

    ac = double_center(dx)
    bc = double_center(dy)
    del dx, dy
    md = ac * bc
    del ac, bc
    t_d = md.sum() / n

    w = markov_rademacher(n, B, ln, rng)
    bh = wild_boot_stats(mh, w)
    del mh
    bd = wild_boot_stats(md, w)
    del md, w

    s = shift_stats(kc, lc)
    del kc, lc
    shifts = s[A:n - A + 1]

    sub_h = subsample_stats(x, y, ln, dn, "hsic", sx, sy)
    sub_d = subsample_stats(x, y, ln, dn, "dcov")

    def pboot(t, b):
        return (1.0 + np.sum(b >= t)) / (1.0 + b.size)

    return {
        "hsic_stat": t_h,
        "dcov_stat": t_d,
        "p_hsic_wb": pboot(t_h, bh),
        "p_hsic_shift": pboot(s[0], shifts),
        "p_hsic_sub": float(np.mean(sub_h >= t_h)),
        "rej_hsic_sub": int(t_h > np.quantile(sub_h, 0.95)),
        "p_dcov_wb": pboot(t_d, bd),
        "p_dcov_sub": float(np.mean(sub_d >= t_d)),
        "rej_dcov_sub": int(t_d > np.quantile(sub_d, 0.95)),
        "shift_check": float(s[0] - t_h),
    }


# ---------------------------------------------------------------- L-ind (Pearson)

def andrews_bartlett_bw(u: np.ndarray) -> float:
    """Andrews (1991) AR(1) plug-in bandwidth for the Bartlett kernel."""
    n = u.size
    rho = float(np.dot(u[1:], u[:-1]) / np.dot(u[:-1], u[:-1]))
    alpha1 = 4 * rho * rho / ((1 - rho) ** 2 * (1 + rho) ** 2)
    return 1.1447 * (alpha1 * n) ** (1.0 / 3.0)


def lrv_nw(u: np.ndarray, S: float) -> float:
    """Newey-West (Bartlett) long-run variance with bandwidth S: weights 1 - j/S, j < S."""
    n = u.size
    v = float(np.dot(u, u) / n)
    jmax = min(int(math.ceil(S)) - 1, n - 1)
    for j in range(1, jmax + 1):
        wgt = 1.0 - j / S
        if wgt <= 0:
            break
        v += 2 * wgt * float(np.dot(u[j:], u[:-j]) / n)
    return v


def lrv_ewc(u: np.ndarray, nu: int) -> float:
    """Equal-weighted cosine LRV (Lazarus-Lewis-Stock-Watson 2018), nu cosine projections."""
    n = u.size
    t = np.arange(1, n + 1)
    j = np.arange(1, nu + 1)
    basis = math.sqrt(2.0 / n) * np.cos(np.pi * j[:, None] * (t[None, :] - 0.5) / n)
    lam = basis @ u
    return float(np.mean(lam * lam))


_TQ: dict[int, float] = {}


def lind_tests(x: np.ndarray, y: np.ndarray) -> dict:
    n = x.size
    u = (x - x.mean()) * (y - y.mean())
    ubar = u.mean()
    uc = u - ubar
    S = min(andrews_bartlett_bw(uc), float(n - 1))
    v_nw = lrv_nw(uc, S)
    nu = math.floor(0.4 * n ** (2.0 / 3.0))
    v_ewc = lrv_ewc(uc, nu)
    t_nw = math.sqrt(n) * ubar / math.sqrt(v_nw) if v_nw > 0 else float("nan")
    t_ewc = math.sqrt(n) * ubar / math.sqrt(v_ewc)
    if nu not in _TQ:
        _TQ[nu] = t_quantile(0.975, nu)
    p_nw = float(2 * (1 - norm_cdf(abs(t_nw)))) if np.isfinite(t_nw) else float("nan")
    return {
        "pearson_r": float(np.corrcoef(x, y)[0, 1]),
        "nw_bw": S,
        "t_nw": t_nw,
        "p_nw": p_nw,
        "rej_nw": int(abs(t_nw) > 1.959963984540054) if np.isfinite(t_nw) else -1,
        "ewc_nu": nu,
        "t_ewc": t_ewc,
        "rej_ewc": int(abs(t_ewc) > _TQ[nu]),
        "nw_lrv_nonpos": int(not v_nw > 0),
    }
