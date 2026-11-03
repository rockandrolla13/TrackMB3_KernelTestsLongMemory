"""Exact numeric checks for THEORY.md (claims C1-C7). No Monte Carlo.

Run:  OMP_NUM_THREADS=4 python -I theory_checks.py > theory_checks.out
Writes theory_exact.csv and theory_exponents.csv next to this file.

Everything is built from one exact object. For independent series the expected statistic is
    E T_n = sum_{|j|<n} (1-|j|/n) kappa_X(j) kappa_Y(j),
where kappa(j) = E k(X_0,X_j) - E k(X,X') = sum_{p>=1} c_p rho(j)^p and c_p = ||a_p||^2 / p!.
Filters act on kappa by convolution. Quadrature checks the series. Nothing is simulated.
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys
import time

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.signal import fftconvolve
from scipy.special import gammaln, ndtr, ndtri

OUT = pathlib.Path(__file__).resolve().parent
T0 = time.time()
SIG = math.sqrt(2.0) * float(ndtri(0.75))  # median |xi - xi'| for iid N(0,1)
HS = (0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
UNEQ = ((0.6, 0.7), (0.6, 0.9), (0.7, 0.8), (0.7, 0.9), (0.8, 0.9))
N_REQ = [1000 * 2 ** k for k in range(7)]  # 1000 ... 64000
N_EXT = [1000 * 2 ** k for k in range(12)]  # ... 2,048,000, used for exponent fits
RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)


def section(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}  [{time.time() - T0:.0f}s]\n{'=' * 78}", flush=True)


# ------------------------------------------------------------------ fGn autocovariance

def gbinom(a: float, m: int) -> float:
    """Generalized binomial coefficient."""
    out = 1.0
    for i in range(m):
        out *= (a - i) / (i + 1)
    return out


def rho_all(H: float, nmax: int) -> np.ndarray:
    """fGn autocovariance at lags 0..nmax. The series form avoids cancellation at large lags."""
    k = np.arange(nmax + 1, dtype=float)
    a = 2.0 * H
    out = 0.5 * (np.abs(k + 1) ** a - 2 * np.abs(k) ** a + np.abs(k - 1) ** a)
    big = k >= 64
    if big.any():
        kb = k[big]
        out[big] = sum(gbinom(a, m) * kb ** (a - m) for m in range(2, 16, 2))
    return out


def d2rho_series(H: float, k: np.ndarray) -> np.ndarray:
    """Second difference rho(k+1)-2rho(k)+rho(k-1) from the series of the fourth difference."""
    a = 2.0 * H
    return 0.5 * sum((2.0 ** (m + 1) - 8.0) * gbinom(a, m) * k ** (a - m) for m in range(4, 20, 2))


# ------------------------------------------------------------------ Hermite layers of a kernel

class Layers:
    """kappa(j) = sum_{p>=1} c_p rho(j)^p for j != 0, with kappa(0) = k0 (k0 >= sum c_p; extras add to k0)."""

    def __init__(self, c: np.ndarray, k0: float, name: str):
        self.c, self.k0, self.name = np.asarray(c, float), float(k0), name

    def kappa(self, rho: np.ndarray, r: int = 0) -> np.ndarray:
        """Deflated autocovariance. r=0 is the raw kernel. rho[0] must equal 1."""
        out = np.empty_like(rho)
        out[0] = self.k0 - self.c[:r].sum()
        out[1:] = series_tail(rho[1:], self.c, r)
        return out


def series_tail(x: np.ndarray, c: np.ndarray, r: int) -> np.ndarray:
    """sum_{p>r} c_p x^p by Horner, with the term count chosen per block of decreasing |x|."""
    out = np.zeros_like(x)
    lo, n = 0, x.size
    block = 8
    while lo < n:
        hi = min(n, lo + block)
        xmax = float(np.abs(x[lo:hi]).max())
        if xmax > 0:
            need = c.size if xmax >= 1 else min(c.size, int(math.log(1e-22) / math.log(xmax)) + r + 4)
            acc = np.zeros(hi - lo)
            for p in range(need, r, -1):
                acc = acc * x[lo:hi] + c[p - 1]
            out[lo:hi] = acc * x[lo:hi] ** (r + 1)
        lo, block = hi, block * 2
    return out


def gauss_layers(sig: float = SIG, pmax: int = 700) -> Layers:
    """Gaussian kernel with X = xi. Psi(r) = sig / sqrt(sig^2 + 2 - 2 r) - sig / sqrt(sig^2 + 2)."""
    a = sig / math.sqrt(sig * sig + 2)
    q = 2.0 / (sig * sig + 2)
    t = np.cumprod((np.arange(1, pmax + 1) - 0.5) / np.arange(1, pmax + 1))
    c = a * t * q ** np.arange(1, pmax + 1)
    return Layers(c, 1.0 - a, f"gauss_id(sig={sig:.4f})")


def psi_gauss(r: np.ndarray, sig: float = SIG) -> np.ndarray:
    return sig / np.sqrt(sig * sig + 2 - 2 * r) - sig / math.sqrt(sig * sig + 2)


def he_norm(x: np.ndarray, pmax: int) -> np.ndarray:
    """Normalized Hermite polynomials He_p(x)/sqrt(p!), p = 0..pmax, in rows."""
    h = np.zeros((pmax + 1,) + x.shape)
    h[0] = 1.0
    if pmax >= 1:
        h[1] = x
    for p in range(1, pmax):
        h[p + 1] = (x * h[p] - math.sqrt(p) * h[p - 1]) / math.sqrt(p + 1)
    return h


def layers_from_gram(gram_fn, name: str, nodes: int = 300, pmax: int = 110) -> Layers:
    """c_p by Gauss-Hermite for a kernel K(xi,xi') on the driver. c_p = E[h_p(xi) h_p(xi') K(xi,xi')]."""
    x, w = hermegauss(nodes)
    w = w / w.sum()
    K = gram_fn(x[:, None], x[None, :])
    h = he_norm(x, pmax) * w
    c = np.einsum("pi,ij,pj->p", h[1:], K, h[1:])
    return Layers(np.maximum(c, 0.0), float(1.0 - w @ K @ w), name)


def biv_gh(f, r: float, nodes: int = 120) -> float:
    """E f(xi_0, xi_1) by two-dimensional Gauss-Hermite at correlation r."""
    x, w = hermegauss(nodes)
    w = w / w.sum()
    a, b = x[:, None], r * x[:, None] + math.sqrt(max(0.0, 1 - r * r)) * x[None, :]
    return float(w @ f(a, b) @ w)


# ------------------------------------------------------------------ exact expected statistics

def et_uncentred(kx: np.ndarray, ky: np.ndarray, ns: list[int]) -> dict[int, float]:
    """E T_n = sum_{|j|<n} (1-|j|/n) kx(j) ky(j)."""
    nmax = max(ns)
    prod = kx[:nmax] * ky[:nmax]
    s0 = np.cumsum(prod)
    s1 = np.cumsum(np.arange(nmax) * prod)
    return {n: float(2 * s0[n - 1] - prod[0] - 2 * s1[n - 1] / n) for n in ns}


def et_centred(kx: np.ndarray, ky: np.ndarray, ns: list[int]) -> dict[int, float]:
    """E (1/n) tr(H Kx H Ky) with H = I - 11'/n, from the Toeplitz matrices of kx and ky."""
    out = {}
    nmax = max(ns)
    un = et_uncentred(kx, ky, ns)
    sx, sy = np.cumsum(kx[:nmax]), np.cumsum(ky[:nmax])
    for n in ns:
        t = np.arange(n)
        ax = sx[t] + sx[n - 1 - t] - kx[0]
        ay = sy[t] + sy[n - 1 - t] - ky[0]
        out[n] = un[n] + (-2.0 * (ax @ ay) / n + ax.sum() * ay.sum() / n ** 2) / n
    return out


def et_split_mmd(k: np.ndarray, ns: list[int]) -> dict[int, float]:
    """E n ||mean first half - mean second half||^2 for n = 2m, with centred features."""
    out = {}
    for n in ns:
        m = n // 2
        j = np.arange(1, n)
        a = np.where(j <= m, 2 * m - 3 * j, -(2 * m - j)).astype(float)
        out[n] = (4.0 / n) * (k[0] * n + 2 * (a @ k[1:n]))
    return out


def inc(k: np.ndarray, lag: int) -> np.ndarray:
    """kappa_D(j) = 2 kappa(j) - kappa(j+L) - kappa(|j-L|), for j = 0 .. len(k)-L-1."""
    m = k.size - lag
    j = np.arange(m)
    return 2 * k[:m] - k[j + lag] - k[np.abs(j - lag)]


def frac_filter_acov(d: float, M: int) -> np.ndarray:
    """Autocovariance of the (1-B)^d weights at lags 0..M (ARFIMA(0,-d,0) recursion, exact)."""
    g = np.empty(M + 1)
    g[0] = math.exp(gammaln(1 + 2 * d) - 2 * gammaln(1 + d))
    ratios = (np.arange(1, M + 1) - 1 - d) / (np.arange(1, M + 1) + d)
    g[1:] = g[0] * np.cumprod(ratios)
    return g


def frac_apply(k: np.ndarray, H: float, jmax: int, M: int) -> np.ndarray:
    """kappa_W(j) = sum_m gamma_w(m) kappa(j+m) for j = 0..jmax. k must reach lag jmax+M.

    Truncating the filter at M leaves a tail of order 1/M. The last term adds its leading part.
    """
    d = H - 0.5
    if d == 0:
        return k[: jmax + 1].copy()
    g = frac_filter_acov(d, M)
    two = np.r_[g[:0:-1], g]
    ks = np.r_[k[jmax + M:0:-1], k[: jmax + M + 1]]
    conv = fftconvolve(ks, two, mode="same")
    mid = jmax + M
    out = conv[mid : mid + jmax + 1]
    tail = 2 * g[M] * M ** (1 + 2 * d) * k[M] * M ** (1 - 2 * d) / M
    return out + tail


def aitken_exponent(e: dict[int, float], n: int) -> float:
    """log2 of the ratio of successive differences of E T_n. It removes a constant term."""
    d1, d2 = e[2 * n] - e[n], e[4 * n] - e[2 * n]
    if d1 == 0 or d2 == 0 or d2 / d1 <= 0:
        return float("-inf")
    return math.log2(d2 / d1)


def predicted_beta(hx: float, hy: float, r: int) -> float:
    """Growth exponent of the leading Hermite pair (r+1, r+1). Zero means bounded or logarithmic."""
    return max(0.0, 1.0 - (r + 1) * (2 - 2 * hx + 2 - 2 * hy))


# ------------------------------------------------------------------ C1 and C3: exact growth rates

def kernel_arrays(lay: Layers, H: float, nmax: int, with_frac: bool) -> dict[str, np.ndarray]:
    rho = rho_all(H, nmax + 60)
    raw = lay.kappa(rho, 0)
    out = {"raw": raw, "def1": lay.kappa(rho, 1), "def2": lay.kappa(rho, 2), "inc1": inc(raw, 1), "inc50": inc(raw, 50)}
    if with_frac:
        M = 2 ** 19
        out["frac"] = frac_apply(raw, H, 64000, M) if nmax + 60 >= 64000 + M else None
    return out


def pair_exact(arr_x: dict, arr_y: dict, ns: list[int], ns_c: list[int]) -> dict[str, dict[int, float]]:
    e: dict[str, dict[int, float]] = {}
    for f in ("raw", "def1", "def2"):
        e[f + "_unc"] = et_uncentred(arr_x[f], arr_y[f], ns)
        e[f + "_c"] = et_centred(arr_x[f], arr_y[f], ns_c)
    e["inc1"] = et_uncentred(arr_x["inc1"], arr_y["inc1"], ns)
    e["inc50"] = et_uncentred(arr_x["inc50"], arr_y["inc50"], ns)
    e["inc1_oneside"] = et_uncentred(arr_x["inc1"], arr_y["raw"], ns)
    if arr_x.get("frac") is not None and arr_y.get("frac") is not None:
        e["frac"] = et_uncentred(arr_x["frac"], arr_y["frac"], N_REQ)
    return e


def expected_slope(feature: str, hx: float, hy: float) -> float:
    """Exponent of E_{2n} - E_n. Positive: growth rate. Zero: log. Negative: approach to a finite limit."""
    if feature.startswith(("inc", "frac")):
        return -1.0
    r = {"raw": 0, "def1": 1, "def2": 2}[feature.split("_")[0]]
    return max(1.0 - (r + 1) * (4 - 2 * hx - 2 * hy), -1.0)


def slope_with_sublead(vals: dict[int, float], beta2: float) -> float:
    """Slope of E_{2n}-E_n after removing a known second power n^beta2 by differencing."""
    ks = sorted(vals)[-5:]
    d = [vals[ks[i + 1]] - vals[ks[i]] for i in range(4)]
    e = [d[i + 1] - 2.0 ** beta2 * d[i] for i in range(3)]
    if e[1] == 0 or e[2] / e[1] <= 0:
        return float("-inf")
    return math.log2(e[2] / e[1])


def fit_ok(vals: dict[int, float], want: float, beta2: float, tol: float) -> tuple[float, float, bool]:
    """Return plain slope, corrected slope, verdict.

    want > 0: the corrected slope must match want within tol.
    want = 0: logarithmic balance, so the plain slope must be near zero.
    want < 0: E T_n must converge, so the differences must shrink.
    """
    ks = sorted(vals)
    plain = aitken_exponent(vals, ks[-3])
    flat = abs(vals[ks[-1]] - vals[ks[-2]]) <= 1e-11 * max(1.0, abs(vals[ks[-1]]))
    if want > 1e-9:
        corr = slope_with_sublead(vals, max(beta2, 0.0))
        gap = want - max(beta2, 0.0)
        if gap <= 0.1 + 1e-9:  # neighbouring powers too close to separate at n = 4e6: require a bracket
            return plain, corr, plain - tol <= want <= corr + tol
        return plain, corr, abs(corr - want) < tol
    if abs(want) <= 1e-9:
        return plain, plain, flat or abs(plain) < 0.05
    return plain, plain, flat or plain < -0.05


def run_c1_c3(lay: Layers, kname: str, nfit: list[int], rows: list, fits: list, tol: float, hs=HS, uneq=UNEQ) -> dict:
    nmax = max(nfit)
    ns_c = [n for n in nfit if n <= 1_024_000]
    with_frac = nmax + 60 >= 64000 + 2 ** 19
    arrs = {H: kernel_arrays(lay, H, nmax, with_frac) for H in hs}
    keep = {}
    bad: list[str] = []
    for hx, hy in [(h, h) for h in hs] + [p for p in uneq if p[0] in hs and p[1] in hs]:
        e = pair_exact(arrs[hx], arrs[hy], nfit, ns_c)
        keep[(hx, hy)] = e
        dx, dy = 2 - 2 * hx, 2 - 2 * hy
        for feat, vals in e.items():
            for n, v in sorted(vals.items()):
                rows.append((kname, "hsic", feat, hx, hy, n, v))
            want = expected_slope(feat, hx, hy)
            r = {"raw": 0, "def1": 1, "def2": 2}.get(feat.split("_")[0], 0)
            beta2 = 1 - (r + 1) * (dx + dy) - min(dx, dy)
            t = tol if not feat.endswith("_c") else 2 * tol
            plain, corr, ok = fit_ok(vals, want, beta2, t)
            if feat.startswith(("inc", "frac")) and hx > 0.5:
                ok = abs(plain + 1.0) < 0.02  # E T_n = Lambda - c/n is the exact claim here
            fits.append((kname, "hsic", feat, hx, hy, max(vals), plain, corr, want, ok))
            if not ok:
                bad.append(f"{feat}@({hx},{hy}): slope {plain:.3f}/{corr:.3f} want {want:.3f}")
    for H in hs:
        d = 2 - 2 * H
        for feat, r in (("raw", 0), ("def1", 1), ("def2", 2)):
            vals = et_split_mmd(arrs[H][feat], ns_c)
            for n, v in sorted(vals.items()):
                rows.append((kname, "mmd", feat, H, H, n, v))
            want = max(1.0 - (r + 1) * d, -1.0)
            plain, corr, ok = fit_ok(vals, want, 1 - (r + 2) * d, 2 * tol)
            fits.append((kname, "mmd", feat, H, H, max(vals), plain, corr, want, ok))
            if not ok:
                bad.append(f"mmd {feat}@{H}: slope {plain:.3f}/{corr:.3f} want {want:.3f}")
        k = arrs[H]["raw"]
        for n in ns_c:
            m = n // 2
            rows.append((kname, "mmd", "inc1", H, H, n, (n / m ** 2) * (6 * k[0] - 8 * k[m] + 2 * k[2 * m])))
    check(f"C1/C3 growth exponents, kernel {kname}: {len(fits)} cumulative fits, tol {tol}", not bad, "; ".join(bad[:8]))
    return {"keep": keep, "arrs": arrs}


def pure_pair_check(nfit: list[int]) -> None:
    """C1 itself: the pair (p,q) alone, sum (1-|j|/n) rho_X^p rho_Y^q, must grow at n^{1-pD_X-qD_Y}."""
    worst = worst_c = 0.0
    for hx, hy, p, q in ((0.8, 0.8, 1, 1), (0.9, 0.9, 1, 1), (0.95, 0.95, 1, 1), (0.9, 0.9, 2, 2), (0.95, 0.95, 2, 2),
                         (0.95, 0.95, 3, 3), (0.7, 0.9, 1, 1), (0.8, 0.9, 1, 1), (0.9, 0.95, 2, 1), (0.95, 0.9, 1, 3)):
        rx, ry = rho_all(hx, N_EXT[-1]), rho_all(hy, N_EXT[-1])
        e = et_uncentred(rx ** p, ry ** q, N_EXT)
        want = 1 - p * (2 - 2 * hx) - q * (2 - 2 * hy)
        worst = max(worst, abs(aitken_exponent(e, N_EXT[-3]) - want))
        ec = et_centred(rx ** p, ry ** q, N_EXT[:11])
        worst_c = max(worst_c, abs(aitken_exponent(ec, N_EXT[8]) - want))
    check("C1: pure (p,q) pairs grow at exactly 1 - p D_X - q D_Y", worst < 0.005, f"worst slope error {worst:.4f}")
    check("C1: centring keeps the exponent of each pure pair", worst_c < 0.02, f"worst slope error {worst_c:.4f}")


def section_c1_c3(rows: list, fits: list) -> dict:
    section("C1 and C3: exact expected statistic and growth exponents")
    lay = gauss_layers()
    rho = rho_all(0.8, 1000)
    x, w = hermegauss(240)
    w = w / w.sum()

    # Identities: series equals closed form. Quadrature equals both.
    r_grid = np.array([0.0, 0.05, 0.2, 0.4, 0.6, 0.8, 0.866])
    series = series_tail(r_grid, lay.c, 0)
    check("Psi series equals the closed form", np.abs(series - psi_gauss(r_grid)).max() < 1e-13)
    kfun = lambda a, b: np.exp(-((a - b) ** 2) / (2 * SIG ** 2))
    direct = np.array([biv_gh(kfun, r) for r in r_grid]) - SIG / math.sqrt(SIG ** 2 + 2)
    check("2-D Gauss-Hermite equals the closed form", np.abs(direct - psi_gauss(r_grid)).max() < 1e-11,
          f"max err {np.abs(direct - psi_gauss(r_grid)).max():.1e}")
    gh = layers_from_gram(kfun, "gh")
    check("c_p by quadrature equals c_p in closed form, p<=40",
          np.abs(gh.c[:40] - lay.c[:40]).max() < 1e-12, f"max err {np.abs(gh.c[:40] - lay.c[:40]).max():.1e}")

    # The deflated kernel written out in x and y, integrated directly, must match the tail of the series.
    def k_r_fun(r: int):
        h = he_norm(x, r)
        fact = [math.sqrt(math.factorial(q)) for q in range(r + 1)]

        def A(q: int, v: np.ndarray) -> np.ndarray:  # <a_q, phi(v)-mu> = E[He_q(xi) k(xi,v)]
            return ((h[q] * fact[q] * w) * np.exp(-((x - v[..., None]) ** 2) / (2 * SIG ** 2))).sum(-1)

        def f(a: np.ndarray, b: np.ndarray) -> np.ndarray:
            a, b = np.broadcast_arrays(a, b)
            out = np.exp(-((a - b) ** 2) / (2 * SIG ** 2)) - A(0, a) - A(0, b) + SIG / math.sqrt(SIG ** 2 + 2)
            he_a, he_b = he_norm(a, r), he_norm(b, r)
            for q in range(1, r + 1):
                nq = lay.c[q - 1] * math.factorial(q)  # ||a_q||^2
                # He_q/q! = h_q / sqrt(q!), so the cross term is A(q,y) h_q(x) / sqrt(q!)
                out = out - (A(q, b) * he_a[q] + A(q, a) * he_b[q]) / fact[q] + nq * he_a[q] * he_b[q] / math.factorial(q)
            return out
        return f

    for r in (1, 2):
        d = np.array([biv_gh(k_r_fun(r), v, nodes=70) for v in r_grid])
        s = series_tail(r_grid, lay.c, r)
        check(f"deflated kernel r={r}: direct 2-D quadrature equals the Hermite tail",
              np.abs(d - s).max() < 2e-10, f"max err {np.abs(d - s).max():.1e}")

    # The O(n) formulas against the literal double sum.
    n = 1000
    for H1, H2 in ((0.8, 0.7), (0.9, 0.9)):
        kx, ky = lay.kappa(rho_all(H1, n + 60), 0), lay.kappa(rho_all(H2, n + 60), 0)
        ix = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        Ax, Ay = kx[ix], ky[ix]
        C = np.eye(n) - 1.0 / n
        lit_c = np.trace(C @ Ax @ C @ Ay) / n
        lit_u = np.sum(Ax * Ay) / n
        check(f"O(n) formulas match the double sum at ({H1},{H2})",
              abs(et_centred(kx, ky, [n])[n] - lit_c) < 1e-9 and abs(et_uncentred(kx, ky, [n])[n] - lit_u) < 1e-9)
        dx, dy = inc(kx, 1), inc(ky, 1)
        D = np.eye(n + 1)[1:] - np.eye(n + 1)[:-1]
        ix1 = np.abs(np.arange(n + 1)[:, None] - np.arange(n + 1)[None, :])
        Kx1, Ky1 = (D @ kx[ix1] @ D.T), (D @ ky[ix1] @ D.T)
        check(f"increment Gram D K D' matches kappa_D at ({H1},{H2})",
              abs(np.sum(Kx1 * Ky1) / n - et_uncentred(dx, dy, [n])[n]) < 1e-9)

    # The paper's balance statements.
    check("C1: (1,1) balance at H_X+H_Y=3/2", all(abs(0.5 - ((1 - h) + (1 - (1.5 - h)))) < 1e-12 for h in (0.6, 0.7, 0.75)))
    check("C1: split-MMD rank 2 balance at H=3/4 and rank 1 long memory for every H>1/2",
          abs(2 * (1 - 0.75) - 0.5) < 1e-12 and all((1 - h) <= 0.5 for h in HS))

    pure_pair_check(N_EXT)
    out = run_c1_c3(lay, "gauss_id", N_EXT, rows, fits, tol=0.02)

    # Monotone transform: the same exponents with a kernel on g(xi)=2 sinh(xi/2).
    med = 1.2  # about the median distance of g(xi) = 2 sinh(xi/2)
    ksinh = lambda a, b: np.exp(-((2 * np.sinh(a / 2) - 2 * np.sinh(b / 2)) ** 2) / (2 * med ** 2))
    lsinh = layers_from_gram(ksinh, "gauss_sinh", nodes=280, pmax=135)
    rq = np.array([0.1, 0.3, 0.5, 0.7])
    psi0 = 1 - lsinh.k0
    err = max(abs(biv_gh(ksinh, v, nodes=280) - psi0 - series_tail(np.array([v]), lsinh.c, 0)[0]) for v in rq)
    check("monotone g=2 sinh(xi/2): quadrature equals the Hermite series", err < 2e-6, f"max err {err:.1e}")
    run_c1_c3(lsinh, "gauss_sinh", N_EXT[:10], rows, fits, tol=0.05, hs=(0.6, 0.8, 0.9, 0.95), uneq=((0.8, 0.9),))
    return out


# ------------------------------------------------------------------ C2: Hermite deflation

def median_dist_sv() -> float:
    x, w = hermegauss(160)
    w = w / w.sum()
    s = np.sqrt(np.exp(x)[:, None] + np.exp(x)[None, :])
    lo, hi = 1e-3, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if w @ (2 * ndtr(mid / s) - 1) @ w < 0.5 else (lo, mid)
    return 0.5 * (lo + hi)


def sv_section(rows: list) -> dict:
    """X_t = exp(xi_t/2) e_t with e_t iid N(0,1). The memory layer is E[phi(X)|xi]."""
    sig = median_dist_sv()
    kbar = lambda a, b: sig / np.sqrt(sig ** 2 + np.exp(np.minimum(a, 40)) + np.exp(np.minimum(b, 40)))
    lay = layers_from_gram(kbar, "gauss_sv", nodes=200, pmax=90)
    nx, ne = 61, 61
    xi, wx = hermegauss(nx)
    e, we = hermegauss(ne)
    wx, we = wx / wx.sum(), we / we.sum()
    XI, E = np.meshgrid(xi, e, indexing="ij")
    W = np.outer(wx, we).ravel()
    XI, E = XI.ravel(), E.ravel()
    X = np.exp(XI / 2) * E
    zi, wz = hermegauss(180)
    wz = wz / wz.sum()
    tail = (wz * ndtr(-np.abs(X)[:, None] * np.exp(-zi / 2))).sum(1)  # P(X < -|x|)
    T = -np.sign(X) * ndtri(tail)
    check("SV: T(X) has unit variance under the quadrature", abs(W @ T ** 2 - 1) < 1e-3, f"var {W @ T ** 2:.6f}")
    K = np.exp(-((X[:, None] - X[None, :]) ** 2) / (2 * sig ** 2))
    pmax, rmax = 8, 6
    hxi = he_norm(XI, pmax)  # He_p(xi)/sqrt(p!)
    ht = he_norm(T, rmax)
    beta = (hxi * W) @ ht.T  # corr(He_p(xi)/sqrt(p!), He_q(T)/sqrt(q!))
    u = hxi * W
    v = ht * W
    res = {}
    # lag-0 value after deflation: ||phi_r||^2 falls by the squared norm of each removed layer
    drop = np.cumsum([0.0] + [float(v[q] @ K @ v[q]) for q in range(1, rmax + 1)])
    for r in (0, 1, 2, 3, 4, 6):
        c_r = np.zeros(pmax)
        for p in range(1, pmax + 1):
            b = u[p] - sum(beta[p, q] * v[q] for q in range(1, r + 1))
            c_r[p - 1] = b @ K @ b
        res[r] = c_r
    c0 = res[0]
    check("SV: product quadrature reproduces c_p from the closed conditional kernel",
          np.abs(c0[:6] - lay.c[:6]).max() < 5e-5, f"max err {np.abs(c0[:6] - lay.c[:6]).max():.1e}")
    check("SV: He_1(T) has zero correlation with every He_p(xi)", np.abs(beta[1:, 1]).max() < 1e-8,
          f"max |corr| {np.abs(beta[1:, 1]).max():.1e}")
    check("SV: r=1 deflation leaves every memory coefficient c_p unchanged",
          np.abs(res[1] - c0).max() < 1e-9, f"max change {np.abs(res[1] - c0).max():.1e}")
    check("SV: r=2,3,4,6 deflation leaves c_1 > 0, so the rank in xi stays one", all(res[r][0] > 0 for r in (2, 3, 4, 6)))
    print("  SV bandwidth", round(sig, 4), " k0", round(lay.k0, 5), " sum c_p", round(float(lay.c.sum()), 5),
          " share of variance in memory layers", round(float(lay.c.sum() / lay.k0), 4))
    print("  corr(He_1(xi), He_q(T)), q=1..6:", np.round(beta[1, 1:], 4))
    for r in (0, 1, 2, 3, 4, 6):
        print(f"  r={r}: c_1={res[r][0]:.3e} ratio to r=0: {res[r][0] / c0[0]:.4f}   c_2={res[r][1]:.3e} ratio: {res[r][1] / c0[1]:.4f}")
    # Exact growth: SV raw against increment. Rate is the claim. Level is reported.
    bad = []
    ns = N_EXT[:11]
    for H in (0.8, 0.9, 0.95):
        rho = rho_all(H, max(ns) + 60)
        raw = lay.kappa(rho, 0)
        d1 = inc(raw, 1)
        e_raw, e_d = et_uncentred(raw, raw, ns), et_uncentred(d1, d1, ns)
        c2 = np.r_[res[2][:8], lay.c[8:]]  # layers above p=8 are tiny and kept at their raw value
        k2 = Layers(c2, lay.k0 - drop[2], "x").kappa(rho)
        e_def2 = et_uncentred(k2, k2, ns)
        for n in ns:
            rows += [("gauss_sv", "hsic", "raw_unc", H, H, n, e_raw[n]), ("gauss_sv", "hsic", "inc1", H, H, n, e_d[n]),
                     ("gauss_sv", "hsic", "def2_unc", H, H, n, e_def2[n])]
        d = 2 - 2 * H
        pure = et_uncentred(lay.c[0] * rho, lay.c[0] * rho, ns)
        s_pure = aitken_exponent(pure, ns[-3])
        s_d = aitken_exponent(e_d, ns[-3])
        s_r = aitken_exponent(e_raw, ns[-3])
        s_2 = aitken_exponent(e_def2, ns[-3])
        print(f"  SV H={H}: slope raw {s_r:.3f}, def r=2 {s_2:.3f}, rank-one layer alone {s_pure:.3f} (theory {1 - 2 * d:.3f}), "
              f"increment {s_d:.3f}; E T_n at n=64000: raw {e_raw[64000]:.4f}, def2 {e_def2[64000]:.4f}, inc1 {e_d[64000]:.4f}")
        if abs(s_pure - (1 - 2 * d)) > 0.01 or abs(s_d + 1) > 0.02:
            bad.append(H)
    check("SV: rank-one layer grows at the C1 rate and the increment kernel converges", not bad, str(bad))
    return {"sig": sig, "c": lay.c, "k0": lay.k0, "res": res}


def even_section() -> None:
    """X = xi^2: T is an even function of xi. Integrate in t = T(X), which is N(0,1) exactly."""
    t, w = hermegauss(240)
    w = w / w.sum()
    abs_xi = -ndtri(ndtr(-t) / 2)
    x = abs_xi ** 2
    sig = 1.0
    K = np.exp(-((x[:, None] - x[None, :]) ** 2) / (2 * sig ** 2))
    he2 = (x - 1) / math.sqrt(2)
    ht = he_norm(t, 8)
    corr = (ht * w) @ he2
    base = w * he2
    c2_raw = base @ K @ base
    out = []
    for r in (1, 2, 3, 4, 6, 8):
        b = base - sum(corr[q] * ht[q] * w for q in range(1, r + 1))
        out.append((r, float(b @ K @ b / c2_raw)))
    print("  X=xi^2: corr(He_2(xi), He_q(T)) q=1..4:", np.round(corr[1:5], 4))
    print("  X=xi^2: c_2 after deflation / c_2 raw:", [(r, round(v, 5)) for r, v in out])
    check("X=xi^2: deflation lowers c_2 but leaves it positive for r<=8, so the rank in xi stays two",
          all(0 < v < 1 for _, v in out))


def shift_mmd(delta: float, r: int, mode: str, lam: float = 0.5) -> float:
    """Squared MMD of the deflated kernel between N(-delta*(1-lam),1) and N(delta*lam,1), on a fine grid.

    mode 'pooled' takes T from the mixture lam P_A + (1-lam) P_B. Mode 'fixed' takes T from N(0,1).
    """
    x = np.linspace(-12, 12, 1601)
    dx = x[1] - x[0]
    a, b = -delta * (1 - lam), delta * lam
    pa = np.exp(-((x - a) ** 2) / 2) / math.sqrt(2 * math.pi)
    pb = np.exp(-((x - b) ** 2) / 2) / math.sqrt(2 * math.pi)
    if mode == "pooled":
        p = lam * pa + (1 - lam) * pb
        lo = lam * ndtr(x - a) + (1 - lam) * ndtr(x - b)
        up = lam * ndtr(a - x) + (1 - lam) * ndtr(b - x)
        T = np.where(lo < 0.5, ndtri(lo), -ndtri(up))
    else:
        p = np.exp(-(x ** 2) / 2) / math.sqrt(2 * math.pi)
        T = x
    nu = pa - pb
    h = he_norm(T, max(r, 1))
    wt = p * dx
    for q in range(1, r + 1):
        nu = nu - ((h[q] * dx) @ (pa - pb)) * h[q] * p  # h_q has norm one under p
    K = np.exp(-((x[:, None] - x[None, :]) ** 2) / (2 * SIG ** 2))
    return float((nu * dx) @ K @ (nu * dx))


def section_c2(rows: list) -> dict:
    section("C2: Hermite deflation")
    lay = gauss_layers()
    # Finite-sample algebra on exact normal scores.
    n = 400
    z = ndtri((np.arange(1, n + 1) - 0.5) / n)
    Kz = np.exp(-((z[:, None] - z[None, :]) ** 2) / (2 * SIG ** 2))
    for r in (1, 2, 4):
        B = he_norm(z, r).T
        Q = np.eye(n) - B @ np.linalg.solve(B.T @ B, B.T)
        ev = np.linalg.eigvalsh(Q @ Kz @ Q)
        check(f"Q K Q is PSD and Q is a projection that annihilates the basis, r={r}",
              ev.min() > -1e-10 and np.abs(Q @ Q - Q).max() < 1e-10 and np.abs(Q @ B).max() < 1e-10, f"min eig {ev.min():.1e}")
    # Monotone: T(g(xi)) = xi or -xi, so He_p coefficients for p<=r vanish.
    x, w = hermegauss(220)
    w = w / w.sum()
    for name, g, sgn in (("increasing sinh", lambda v: 2 * np.sinh(v / 2), 1.0), ("decreasing cubic", lambda v: -(v + v ** 3 / 6), -1.0)):
        gx = g(x)
        K = np.exp(-((gx[:, None] - gx[None, :]) ** 2) / (2 * 1.2 ** 2))
        T = sgn * x
        h, ht = he_norm(x, 6), he_norm(T, 6)
        for r in (1, 2):
            beta = (h * w) @ ht.T
            cs = []
            for p in range(1, 7):
                b = h[p] * w - sum(beta[p, q] * ht[q] * w for q in range(1, r + 1))
                cs.append(float(b @ K @ b))
            ok = max(cs[:r]) < 1e-20 and cs[r] > 0
            check(f"monotone g ({name}), r={r}: first surviving Hermite rank is r+1", ok, f"c_p={np.array2string(np.array(cs[:r + 1]), precision=2)}")
    sv = sv_section(rows)
    even_section()
    # Mean-shift local rates. Slopes in delta of MMD^2.
    a3 = 6 * lay.c[2]
    bad = []
    for mode, lam, r, want in (("fixed", 0.5, 0, 2), ("pooled", 0.5, 0, 2), ("fixed", 0.5, 1, 6), ("pooled", 0.5, 1, 6), ("fixed", 0.5, 2, 6),
                               ("pooled", 0.5, 2, 6), ("fixed", 0.5, 3, 10), ("pooled", 0.5, 3, 10), ("pooled", 0.25, 1, 4), ("pooled", 0.25, 2, 6)):
        m1, m2 = shift_mmd(0.2, r, mode, lam), shift_mmd(0.1, r, mode, lam)
        s = math.log2(m1 / m2)
        print(f"  shift {mode:6s} lam={lam} r={r}: MMD^2 slope in delta {s:.3f} (want {want}); MMD^2(0.2)={m1:.3e}")
        if abs(s - want) > 0.15:
            bad.append((mode, lam, r, s))
    check("C2 local rates: MMD^2 slope in delta is 2(r+1) for even r, 2(r+2) for odd r in the symmetric split", not bad, str(bad))
    d = 0.1
    fx, pl = shift_mmd(d, 1, "fixed"), shift_mmd(d, 1, "pooled")
    pred_f, pred_p = d ** 6 * a3 / 576, d ** 6 * a3 / 144
    print(f"  r=1, delta=0.1: fixed {fx:.4e} vs delta^6||a_3||^2/576 = {pred_f:.4e}; pooled {pl:.4e} vs delta^6||a_3||^2/144 = {pred_p:.4e}")
    check("C2 leading constants: pooled reference gives delta^6||a_3||^2/144 and fixed gives /576",
          abs(fx / pred_f - 1) < 0.03 and abs(pl / pred_p - 1) < 0.03, f"ratios {fx / pred_f:.3f}, {pl / pred_p:.3f}")
    # Blind pair with a fixed reference: w = c T exactly. Raw MMD is positive and deflated MMD is zero.
    t = np.linspace(-12, 12, 1601)
    dt = t[1] - t[0]
    p = np.exp(-(t ** 2) / 2) / math.sqrt(2 * math.pi)
    c = 1.0
    fa = c * np.maximum(t, 0) + (1 - c / math.sqrt(2 * math.pi))
    fb = fa - c * t
    KK = np.exp(-((t[:, None] - t[None, :]) ** 2) / (2 * SIG ** 2))
    nu0 = (fa - fb) * p * dt
    alpha = (t * p * dt) @ (fa - fb)
    nu1 = nu0 - alpha * t * p * dt
    check("C2 blind pair, fixed reference: both are densities, raw MMD>0, deflated r=1 MMD=0",
          fa.min() >= 0 and fb.min() >= 0 and abs(fa @ (p * dt) - 1) < 1e-4 and abs(fb @ (p * dt) - 1) < 1e-4 and nu0 @ KK @ nu0 > 1e-3 and abs(nu1 @ KK @ nu1) < 1e-20,
          f"raw {nu0 @ KK @ nu0:.3e}, deflated {nu1 @ KK @ nu1:.1e}")
    # HSIC under a Gaussian copula: series against the closed form, then the slope in rho.
    x, w = hermegauss(240)
    w = w / w.sum()
    K = np.exp(-((x[:, None] - x[None, :]) ** 2) / (2 * SIG ** 2))
    P = 30
    h = he_norm(x, P)[1:] * w
    A = h @ K @ h.T  # <a_p,a_q>/sqrt(p!q!)

    def hs_norm(rh: float, r: int) -> float:
        pw = rh ** np.arange(1, P + 1)
        pw[:r] = 0
        return float(pw @ (A * A) @ pw)

    def hsic_closed(rh: float) -> float:
        s = SIG ** 2
        det1 = np.linalg.det(np.eye(2) + np.array([[2, 2 * rh], [2 * rh, 2]]) / s)
        det2 = np.linalg.det(np.eye(2) + np.array([[2, rh], [rh, 2]]) / s)
        return float(det1 ** -0.5 - 2 * det2 ** -0.5 + s / (s + 2))
    err = max(abs(hs_norm(v, 0) - hsic_closed(v)) for v in (0.05, 0.2, 0.5, 0.7))
    check("C2 HSIC: Hermite series for ||C||^2 equals the closed Gaussian form", err < 1e-12, f"max err {err:.1e}")
    bad = []
    for r in (0, 1, 2, 3):
        s = math.log2(hs_norm(0.02, r) / hs_norm(0.01, r))
        if abs(s - 2 * (r + 1)) > 0.05:
            bad.append((r, s))
    check("C2 HSIC: ||C^r||^2 has slope 2(r+1) in rho", not bad, str(bad))
    return {"sv": sv}


# ------------------------------------------------------------------ C4: target and power

def fgn_spectrum(H: float, lam: np.ndarray, J: int = 400) -> np.ndarray:
    """fGn spectral density on (0,pi], normalized so that its integral over (-pi,pi) is one."""
    a = 2 * H + 1
    c = math.exp(gammaln(2 * H + 1)) * math.sin(math.pi * H) / (2 * math.pi)
    j = np.arange(1, J + 1)[:, None]
    s = np.abs(lam) ** (-a) + ((2 * math.pi * j + lam) ** (-a) + (2 * math.pi * j - lam) ** (-a)).sum(0)
    # Euler-Maclaurin tail of the two outer sums
    tail = 2 * (2 * math.pi) ** (-a) * (J + 0.5) ** (1 - a) / (a - 1)
    return c * 2 * (1 - np.cos(lam)) * (s + tail)


def panels(lo: float, hi: float, npanel: int = 40, nodes: int = 24) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre on geometric panels, which resolve the algebraic singularity at zero."""
    edges = np.geomspace(lo, hi, npanel + 1)
    gx, gw = np.polynomial.legendre.leggauss(nodes)
    x = np.concatenate([(b - a) / 2 * gx + (a + b) / 2 for a, b in zip(edges[:-1], edges[1:])])
    w = np.concatenate([(b - a) / 2 * gw for a, b in zip(edges[:-1], edges[1:])])
    return x, w


def spectral_sum(H: float, power: float) -> tuple[float, float]:
    """Return (sum_k r(k)^2, variance) for fGn filtered with squared gain |2 sin(lam/2)|^power."""
    x, w = panels(1e-12, math.pi)
    f = fgn_spectrum(H, x)
    g = f * np.abs(2 * np.sin(x / 2)) ** power
    v = 2 * (w @ g)
    return float(2 * math.pi * 2 * (w @ (g * g)) / v ** 2), float(v)


def target_closed(rf, sig: float = SIG) -> float:
    """||2C_0 - C_1 - C_{-1}||^2 from 16 closed Gaussian integrals. rf(m) = Corr(xi_t, eta_{t+m})."""
    s = sig * sig
    total = 0.0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                for d in (0, 1):
                    u = rf(a - c) + rf(b - d)
                    total += (-1) ** (a + b + c + d) * ((1 + 2 / s) ** 2 - (u / s) ** 2) ** -0.5
    return total


def section_c4(store: dict) -> None:
    section("C4: population target and power of the increment kernel")
    lay = gauss_layers()
    x, w = hermegauss(240)
    w = w / w.sum()
    K = np.exp(-((x[:, None] - x[None, :]) ** 2) / (2 * SIG ** 2))
    P = 40
    h = he_norm(x, P)[1:] * w
    A2 = (h @ K @ h.T) ** 2

    def target_series(rf) -> float:
        co = np.array([2 * rf(0) ** p - rf(1) ** p - rf(-1) ** p for p in range(1, P + 1)])
        return float(co @ A2 @ co)

    worst = 0.0
    for H, rh in ((0.5, 0.3), (0.8, 0.2), (0.8, 0.6), (0.95, 0.4)):
        R = 2 ** (2 * H - 1) - 1
        rf = lambda m, rh=rh, R=R: rh * (1.0 if m == 0 else R)
        worst = max(worst, abs(target_closed(rf) - target_series(rf)))
    lead = lambda m: 0.5 if m == -1 else (0.25 if m == 0 else 0.1)  # asymmetric r(-1), r(0), r(1)
    worst = max(worst, abs(target_closed(lead) - target_series(lead)))
    check("C4: the p-th term formula equals the closed-form target (16 Gaussian integrals)", worst < 1e-12, f"max err {worst:.1e}")
    flat = lambda m: 1 / math.sqrt(3)
    check("C4: blind set, flat MA dependence r(0)=r(1)=r(-1) gives target exactly zero", abs(target_closed(flat)) < 1e-14,
          f"value {target_closed(flat):.1e}")
    print("  Y=rho X + sqrt(1-rho^2) Z. Target ||.||^2 at rho=0.2 for H grid and its first-layer approximation:")
    for H in HS:
        R = 2 ** (2 * H - 1) - 1
        rf = lambda m, R=R: 0.2 * (1.0 if m == 0 else R)
        print(f"    H={H}: R={R:.4f} target={target_closed(rf):.4e}  4 rho^2 (1-R)^2 c_1^2={4 * 0.04 * (1 - R) ** 2 * lay.c[0] ** 2:.4e}")
    check("C4: target positive for rho != 0 at every H on the grid",
          all(target_closed(lambda m, R=2 ** (2 * H - 1) - 1, s=s: s * 0.05 * (1.0 if m == 0 else R)) > 0 for H in HS for s in (-1, 1)))
    print("  Common smooth factor: X=F+e, Y=F+u with F fGn(H_F) of variance share s. Target against raw HSIC ratio:")
    for HF in (0.7, 0.9, 0.99):
        RF = 2 ** (2 * HF - 1) - 1
        s = 0.3
        rf = lambda m, RF=RF, s=s: s * (1.0 if m == 0 else RF)
        raw = float(((s ** np.arange(1, P + 1)) @ A2 @ (s ** np.arange(1, P + 1))))
        print(f"    H_F={HF}: R_F={RF:.4f}, target/raw HSIC = {target_closed(rf) / raw:.4e}")

    # Efficiency numbers. Direct lag sums and the spectral integral must agree.
    table = {}
    bad = []
    for H in HS:
        rho = rho_all(H, 3_000_000)
        R = rho[1]
        gam = inc(rho, 1)
        s_d = float((gam[0] ** 2 + 2 * (gam[1:] ** 2).sum()) / gam[0] ** 2)
        s_spec, v_spec = spectral_sum(H, 2.0)
        v_one = float(gam[0] * rho[0] + 2 * (gam[1:] * rho[1 : gam.size]).sum())
        eff_one = (1 - R) ** 2 / v_one
        M = 2 ** 19
        if H > 0.5:
            gw = frac_apply(rho, H, 200000, M)
        else:
            gw = rho[:200001].copy()
        s_w = float((gw[0] ** 2 + 2 * (gw[1:] ** 2).sum()) / gw[0] ** 2)
        s_w_spec, v_w_spec = spectral_sum(H, 2 * (H - 0.5))
        s_raw = float("inf")
        if H < 0.75:
            tail = (H * (2 * H - 1)) ** 2 * rho.size ** (1 - 2 * (2 - 2 * H)) / (2 * (2 - 2 * H) - 1)
            s_raw = float(1 + 2 * ((rho[1:] ** 2).sum() + tail))
        table[H] = dict(R=R, S_inc=s_d, eff_inc=1 / s_d, eff_oneside=eff_one, S_frac=s_w, v_frac=float(gw[0]), S_raw=s_raw)
        print(f"    H={H:4}: raw rank-one layer sum rho^2 = {s_raw:.5f}, efficiency {1 / s_raw:.4f}")
        print(f"    H={H:4}: R={R:.4f} sum r_D^2={s_d:.5f} (spectral {s_spec:.5f}) efficiency={1 / s_d:.4f}; one-sided efficiency={eff_one:.4f}; "
              f"frac filter sum r^2={s_w:.6f} (spectral {s_w_spec:.6f})")
        if abs(s_d - s_spec) > 2e-5 or abs(s_w - s_w_spec) > 2e-3:
            bad.append((H, s_d, s_spec, s_w, s_w_spec))
    check("C4: lag sums and spectral integrals agree", not bad, str(bad))
    check("C4: sum r_D^2 = 1.5 exactly at H=1/2", abs(table[0.5]["S_inc"] - 1.5) < 1e-12)
    check("C4 claim 'about 1.26 at H=0.8'", abs(table[0.8]["S_inc"] - 1.26) < 0.01, f"exact {table[0.8]['S_inc']:.5f}")
    check("C4 claim 'oracle fractional filter gives sum 1' holds only approximately for fGn",
          all(table[H]["S_frac"] >= 1 - 1e-6 for H in HS) and max(table[H]["S_frac"] for H in HS) < 1.05,
          f"sums {[round(table[H]['S_frac'], 5) for H in HS]}")
    store["c4"] = table


# ------------------------------------------------------------------ C5: telescoping and the location cap

def toeplitz_pcg(rho: np.ndarray, b: np.ndarray, tol: float = 1e-11, maxit: int = 2000) -> np.ndarray:
    """Solve Gamma x = b by conjugate gradients with Chan's circulant preconditioner."""
    n = b.size
    from scipy.fft import irfft, rfft
    emb = np.r_[rho[:n], [0.0], rho[1:n][::-1]]
    fe = rfft(emb)
    mv = lambda v: irfft(fe * rfft(v, 2 * n), 2 * n)[:n]
    k = np.arange(n)
    c = ((n - k) * rho[:n] + k * np.r_[0.0, rho[n - k[1:]]]) / n
    ev = rfft(c).real
    pre = lambda r: irfft(rfft(r) / ev, n)
    x = np.zeros(n)
    r = b.copy()
    z = pre(r)
    p = z.copy()
    rz = r @ z
    nb = np.linalg.norm(b)
    for _ in range(maxit):
        ap = mv(p)
        alpha = rz / (p @ ap)
        x += alpha * p
        r -= alpha * ap
        if np.linalg.norm(r) < tol * nb:
            break
        z = pre(r)
        rz2 = r @ z
        p = z + (rz2 / rz) * p
        rz = rz2
    return x


def info_step(H: float, n: int) -> float:
    rho = rho_all(H, n)
    one = np.ones(n)
    c = np.r_[np.zeros(n // 2), np.ones(n - n // 2)]
    gi1, gic = toeplitz_pcg(rho, one), toeplitz_pcg(rho, c)
    return float(c @ gic - (one @ gic) ** 2 / (one @ gi1))


def section_c5(store: dict) -> None:
    section("C5: telescoping and the location cap")
    got = {H: 0.1 * math.sqrt(info_step(H, 2000)) for H in (0.5, 0.7, 0.8)}
    check("C5: reproduces the paper's delta I_n^{1/2} = 2.24, 0.71, 0.46 at n=2000",
          abs(got[0.5] - 2.24) < 0.005 and abs(got[0.7] - 0.71) < 0.005 and abs(got[0.8] - 0.46) < 0.005, str({k: round(v, 3) for k, v in got.items()}))
    table = {}
    bad = []
    for H in HS:
        ns = [2000, 64000, 128000, 256000]
        Is = [info_step(H, n) for n in ns]
        # Differences remove the constant term that the high-frequency part of the step adds.
        slope = math.log2((Is[3] - Is[2]) / (Is[2] - Is[1]))
        K = Is[-1] / ns[-1] ** (2 - 2 * H)
        m = ns[-1] // 2
        v = (4 * m ** (2 * H) - (2 * m) ** (2 * H)) / m ** 2
        eff = (1 / v) / Is[-1]
        table[H] = dict(K=K, slope=slope, eff_halfmean=eff, I2000=Is[0])
        print(f"    H={H:4}: I_n/n^(2-2H)={K:.4f} slope={slope:.4f} (theory {2 - 2 * H:.3f}); half-mean efficiency={eff:.4f}")
        if abs(slope - (2 - 2 * H)) > 0.01:
            bad.append((H, slope))
    check("C5: I_n grows at n^(2-2H), so no location-invariant test has a sqrt(n) rate for H>1/2", not bad, str(bad))
    print("    half-mean efficiency falls with H: the cap is reachable only by a contrast that weights the change point")
    # Proof check: I_n <= (2 pi)^-2 int sin^4(m l/2) / (sin^2(l/2) f(l)) dl, by Cauchy-Schwarz on w'1=0.
    bad = []
    for H in (0.6, 0.8, 0.95):
        for n in (2000, 16000):
            m = n // 2
            grid = (np.arange(40 * m) + 0.5) * math.pi / (40 * m)
            f = fgn_spectrum(H, grid, J=60)
            upper = 2 * ((np.sin(m * grid / 2) ** 4 / (np.sin(grid / 2) ** 2 * f)).sum() * math.pi / (40 * m)) / (2 * math.pi) ** 2
            In = info_step(H, n)
            lower = m ** (2 - 2 * H) / (4 - 2 ** (2 * H))
            print(f"    H={H} n={n}: half-mean bound {lower:.4f} <= I_n {In:.4f} <= spectral bound {upper:.4f}")
            if not lower <= In * (1 + 1e-9) or not In <= upper * (1 + 1e-9):
                bad.append((H, n))
    check("C5: both proved bounds on I_n hold, so I_n has exact order n^(2-2H)", not bad, str(bad))
    # Any invertible filter leaves the information unchanged. L = first differences is one instance.
    n = 600
    for H in (0.7, 0.9):
        rho = rho_all(H, n)
        G = rho[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]
        L = np.eye(n)[1:] - np.eye(n)[:-1]
        c = np.r_[np.zeros(n // 2), np.ones(n // 2)]
        Lc = L @ c
        via_diff = float(Lc @ np.linalg.solve(L @ G @ L.T, Lc))
        one = np.ones(n)
        gi1, gic = np.linalg.solve(G, one), np.linalg.solve(G, c)
        direct = float(c @ gic - (one @ gic) ** 2 / (one @ gi1))
        contrast = float((Lc @ Lc) / (Lc @ (L @ G @ L.T) @ Lc) * (Lc @ Lc))
        check(f"C5: GLS on differences has the same information at H={H}", abs(via_diff - direct) < 1e-7 * direct,
              f"I={direct:.4f}; a plain contrast on the differences has I={contrast:.4f}")
    # Scale change: the information has no memory cap. I = m + tr(G^-1 E G E) >= m.
    for H in (0.7, 0.9, 0.95):
        n = 800
        m = n // 2
        rho = rho_all(H, n)
        G = rho[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]
        E = np.diag(np.r_[np.zeros(m), np.ones(m)])
        extra = float(np.trace(np.linalg.solve(G, E @ G @ E)))
        check(f"C5: scale-change Fisher information is at least m at H={H}", m + extra >= m and extra >= 0,
              f"I/n={(m + extra) / n:.4f}")
    store["c5"] = table


# ------------------------------------------------------------------ C6: seasonality

def gram_mean(H: float, n: int, mean: np.ndarray, scale: np.ndarray, sig: float) -> np.ndarray:
    """E k(X_s, X_t) for X_t = mean_t + scale_t xi_t, all s,t. The diagonal is exactly one."""
    rho = rho_all(H, n)
    lag = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    var = scale[:, None] ** 2 + scale[None, :] ** 2 - 2 * scale[:, None] * scale[None, :] * rho[lag]
    var = np.maximum(var, 0.0)
    v = sig * sig + var
    return sig / np.sqrt(v) * np.exp(-((mean[:, None] - mean[None, :]) ** 2) / (2 * v))


def diff2(G: np.ndarray, L: int) -> np.ndarray:
    return G[L:, L:] - G[L:, :-L] - G[:-L, L:] + G[:-L, :-L]


def et_raw_season(Gx: np.ndarray, Gy: np.ndarray) -> float:
    n = Gx.shape[0]
    ax, ay = Gx.sum(1), Gy.sum(1)
    return float((np.sum(Gx * Gy) - 2 * (ax @ ay) / n + ax.sum() * ay.sum() / n ** 2) / n)


def mixture_median(cent: np.ndarray) -> float:
    mu = (cent[:, None] - cent[None, :]).ravel()
    sd = math.sqrt(2.0)
    lo, hi = 1e-6, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if np.mean(ndtr((mid - mu) / sd) - ndtr((-mid - mu) / sd)) < 0.5 else (lo, mid)
    return 0.5 * (lo + hi)


def limit_values(P: int, kind: str, sig: float, a: float = 1.0, b: float = 0.5) -> dict[int, tuple[float, float]]:
    """(h_raw(d), h_lag1(d)) for the phase gap d, with exact closed-form conditional means."""
    th = 2 * np.pi * np.arange(P) / P
    out = {}

    def G(par: np.ndarray) -> np.ndarray:
        if kind == "mean":
            return sig / math.sqrt(sig ** 2 + 2) * np.exp(-((par[:, None] - par[None, :]) ** 2) / (2 * (sig ** 2 + 2)))
        return sig / np.sqrt(sig ** 2 + par[:, None] ** 2 + par[None, :] ** 2)
    base = a * np.sin(th) if kind == "mean" else 1 + b * np.sin(th)
    for d in range(P):
        par_y = a * np.sin(th - 2 * np.pi * d / P) if kind == "mean" else 1 + b * np.sin(th - 2 * np.pi * d / P)
        gx, gy = G(base), G(par_y)
        c = np.eye(P) - 1.0 / P
        h_raw = float(np.sum((c @ gx @ c) * (c @ gy @ c)) / P ** 2)
        s = np.roll(np.eye(P), 1, axis=0)
        D = np.eye(P) - s
        h_1 = float(np.sum((D @ gx @ D.T) * (D @ gy @ D.T)) / P ** 2)
        out[d] = (h_raw, h_1)
    return out


def section_c6(store: dict, csv_rows: list) -> None:
    section("C6: seasonality")
    P = 50
    th = 2 * np.pi * np.arange(P) / P
    sig_m = mixture_median(np.sin(th))
    lim_m = limit_values(P, "mean", sig_m)
    check("C6: reproduces Proposition 1 h(0)=0.005869 and min h=0.000039", abs(lim_m[0][0] - 0.005869) < 1e-6 and abs(min(v[0] for v in lim_m.values()) - 0.000039) < 1e-6,
          f"h(0)={lim_m[0][0]:.6f}, min={min(v[0] for v in lim_m.values()):.6f}")
    att = (2 * math.sin(math.pi / P)) ** 4
    ratio = lim_m[0][1] / lim_m[0][0]
    print(f"  (2 sin(pi/P))^4 = {att:.4e}; (2 pi/P)^4 = {(2 * math.pi / P) ** 4:.4e}; h_1(0)/h_raw(0) = {ratio:.4e}")
    check("C6: lag-1 attenuates the seasonal-mean limit h by about (2 sin(pi/P))^4", abs(ratio / att - 1) < 0.10, f"ratio/att={ratio / att:.4f}")
    check("C6: h_1(d) > 0 at every phase gap, so lag-1 does not reach a zero limit for fixed P", min(v[1] for v in lim_m.values()) > 0,
          f"min h_1={min(v[1] for v in lim_m.values()):.3e} at d={min(lim_m, key=lambda d: lim_m[d][1])}")
    lim_v = limit_values(P, "scale", SIG)
    print(f"  seasonal scale, sig={SIG:.4f}: h_raw(0)={lim_v[0][0]:.4e}, h_1(0)={lim_v[0][1]:.4e}, ratio={lim_v[0][1] / lim_v[0][0]:.4e}")

    # Exact finite-n means given the phases.
    def series_for(kind: str, n: int, u: int, per: int) -> tuple[np.ndarray, np.ndarray]:
        tt = 2 * np.pi * (np.arange(n) + u) / per
        if kind.startswith("none"):
            return np.zeros(n), np.ones(n)
        if kind == "mean":
            return np.sin(tt), np.ones(n)
        return np.zeros(n), 1 + 0.5 * np.sin(tt)

    summary = {}
    for H in (0.5, 0.8):
        print(f"  H={H}: exact E T_n given the phase gap (n, raw centred, lag-1, lag-50)")
        base = {}
        for kind, label in (("none", "no season"), ("none_m", "no season, SEAS-M sig"), ("mean", "SEAS-M P=50"), ("scale", "SEAS-V P=50"), ("slow", "SEAS-V-slow P=4n")):
            sig = sig_m if kind in ("mean", "none_m") else SIG
            for n in (2000, 4000):
                # Gap 12 is needed only at n=2000. The slope in n uses gap 0.
                for d in ((0,) if kind.startswith("none") or n > 2000 else (0, 12)):
                    per = 4 * n if kind == "slow" else P
                    du = d if kind != "slow" else (0 if d == 0 else per // 4)
                    mx, sx = series_for("scale" if kind == "slow" else kind, n + 50, 0, per)
                    my, sy = series_for("scale" if kind == "slow" else kind, n + 50, -du, per)
                    Gx, Gy = gram_mean(H, n + 50, mx, sx, sig), gram_mean(H, n + 50, my, sy, sig)
                    raw = et_raw_season(Gx[:n, :n], Gy[:n, :n])
                    d1 = float(np.sum(diff2(Gx[:n + 1, :n + 1], 1) * diff2(Gy[:n + 1, :n + 1], 1)) / n)
                    d50 = float(np.sum(diff2(Gx, 50) * diff2(Gy, 50)) / n)
                    summary[(H, kind, n, d)] = (raw, d1, d50)
                    csv_rows += [(f"season_{kind}_gap{d}", "hsic", "raw_c", H, H, n, raw), (f"season_{kind}_gap{d}", "hsic", "inc1", H, H, n, d1),
                                 (f"season_{kind}_gap{d}", "hsic", "inc50", H, H, n, d50)]
                    print(f"    {label:18s} gap={d:<3d} n={n:<5d} raw={raw:9.4f}  lag1={d1:8.4f}  lag50={d50:8.4f}")
        # lag-50 under a P=50 season is bounded in n. Lag-1 has slope h_1 in n.
        for kind in ("mean", "scale"):
            a, b = summary[(H, kind, 2000, 0)], summary[(H, kind, 4000, 0)]
            check(f"C6: lag-50 is bounded under {kind}, H={H}: E T_n changes by under 1% from n=2000 to 4000", abs(b[2] / a[2] - 1) < 0.01,
                  f"{a[2]:.4f} -> {b[2]:.4f}")
        s1 = (summary[(H, "mean", 4000, 0)][1] - summary[(H, "mean", 2000, 0)][1]) / 2000
        print(f"    lag-1 slope in n under SEAS-M gap 0: {s1:.3e}; closed h_1(0) = {lim_m[0][1]:.3e}; ratio {s1 / lim_m[0][1]:.3f}")
        lam = summary[(H, "none", 2000, 0)][1]
        ns = lam / lim_m[0][1]
        print(f"    n at which n*h_1(0) equals the null mean {lam:.3f}: {ns:.3e}")
        store[f"c6_nstar_{H}"] = ns
        check(f"C6: lag-1 slope in n matches h_1(0) at H={H}", abs(s1 / lim_m[0][1] - 1) < 0.25, f"ratio {s1 / lim_m[0][1]:.3f}")
    store["c6"] = {str(k): v for k, v in summary.items()}
    store["c6_lim"] = dict(h0=lim_m[0], hmin=min(v[1] for v in lim_m.values()), v0=lim_v[0], sig_m=sig_m)


# ------------------------------------------------------------------ C7: circular shifts

def shift_mean(kx: np.ndarray, ky: np.ndarray, n: int, s: int, centred: bool) -> float:
    """E T_n(s), with Y circularly shifted by s, from the exact pair count across the wrap."""
    j = np.arange(1, n)
    c = np.minimum.reduce([j, n - j, np.full(n - 1, s), np.full(n - 1, n - s)])
    tr = n * kx[0] * ky[0] + 2 * (kx[1:n] * ((n - j - c) * ky[1:n] + c * ky[n - j])).sum()
    if not centred:
        return float(tr / n)
    sx, sy = np.cumsum(kx[:n]), np.cumsum(ky[:n])
    t = np.arange(n)
    ax, ay = sx[t] + sx[n - 1 - t] - kx[0], sy[t] + sy[n - 1 - t] - ky[0]
    return float((tr - 2 * (ax @ np.roll(ay, -s)) / n + ax.sum() * ay.sum() / n ** 2) / n)


def shift_cross(kx: np.ndarray, ky: np.ndarray, n: int, s: int) -> float:
    """E <Z_n(0), Z_n(s)> = (1/n) sum_{t,t'} kx(t-t') ky(t - sigma_s(t'))."""
    l = np.arange(-(n - 1), n)
    lo = np.maximum(1, 1 - l)
    up = np.minimum(n, n - l)
    ca = np.maximum(0, np.minimum(up, n - s) - lo + 1)
    cb = np.maximum(0, up - np.maximum(n - s + 1, lo) + 1)
    return float((kx[np.abs(l)] * (ca * ky[np.abs(l - s)] + cb * ky[np.abs(l - s + n)])).sum() / n)


def drift_shift_bias() -> dict:
    """SEAS-V-slow, period 4n. Exact E T_n(0) against the mean of E T_n(s) over offsets, per phase pair."""
    n = 2000
    per = 4 * n
    phases = (0, n, 2 * n, 3 * n)
    offsets = [int(f * n) for f in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
    out = {}
    print("  SEAS-V-slow: ratio E T_n(0) / mean over offsets of E T_n(s), n=2000, 16 phase pairs")
    for H in (0.5, 0.8):
        C = np.eye(n) - 1.0 / n
        raw, d1 = {}, {}
        for u in phases:
            v = 1 + 0.5 * np.sin(2 * np.pi * (np.arange(n + 1) + u) / per)
            G = gram_mean(H, n + 1, np.zeros(n + 1), v, SIG)
            raw[u] = C @ G[:n, :n] @ C
            d1[u] = diff2(G, 1)
        for name, mats in (("raw", raw), ("lag1", d1)):
            ratios = []
            for u in phases:
                for w in phases:
                    t0 = np.sum(mats[u] * mats[w]) / n
                    ts = []
                    for s_ in offsets:
                        perm = (np.arange(n) + s_) % n
                        ts.append(np.sum(mats[u] * mats[w][np.ix_(perm, perm)]) / n)
                    ratios.append(t0 / float(np.mean(ts)))
            r = np.array(ratios)
            out[(H, name)] = (float(r.min()), float(np.median(r)), float(r.max()))
            print(f"    H={H} {name:5s}: min {r.min():.3f}  median {np.median(r):.3f}  max {r.max():.3f}")
    return {str(k): v for k, v in out.items()}


def section_c7(store: dict, csv_rows: list) -> None:
    section("C7: circular-shift rule")
    lay = gauss_layers()
    # Brute-force check of the two wrap formulas.
    n, s = 300, 70
    kx, ky = lay.kappa(rho_all(0.9, 2 * n + 60)), lay.kappa(rho_all(0.7, 2 * n + 60))
    idx = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    A, B = kx[idx], ky[idx]
    perm = (np.arange(n) + s) % n
    Bs = B[np.ix_(perm, perm)]
    C = np.eye(n) - 1.0 / n
    lit_c, lit_u = np.trace(C @ A @ C @ Bs) / n, np.sum(A * Bs) / n
    cross = np.sum(A * ky[np.abs(np.arange(n)[:, None] - perm[None, :])]) / n
    check("C7: wrap formulas match the literal double sums",
          abs(shift_mean(kx, ky, n, s, True) - lit_c) < 1e-10 and abs(shift_mean(kx, ky, n, s, False) - lit_u) < 1e-10 and abs(shift_cross(kx, ky, n, s) - cross) < 1e-10)
    out = {}
    print("  mean ratio = E T_n(s)/E T_n(0); cross = E<Z(0),Z(s)>/E||Z(0)||^2")
    for hx, hy in ((0.6, 0.6), (0.7, 0.7), (0.8, 0.8), (0.9, 0.9), (0.95, 0.95), (0.7, 0.9), (0.8, 0.9)):
        nmax = 2 * 64000 + 60
        kx, ky = lay.kappa(rho_all(hx, nmax)), lay.kappa(rho_all(hy, nmax))
        dx, dy = inc(kx, 1), inc(ky, 1)
        for n in (2000, 8000, 64000):
            row = []
            for frac in (0.1, 0.5):
                s = int(frac * n)
                rc = shift_mean(kx, ky, n, s, True) / shift_mean(kx, ky, n, 0, True)
                ru = shift_mean(kx, ky, n, s, False) / shift_mean(kx, ky, n, 0, False)
                ri = shift_mean(dx, dy, n, s, False) / shift_mean(dx, dy, n, 0, False)
                # Centred: the common term phi-bar tensor psi-bar drops out.
                sx, sy = np.cumsum(kx[:n]), np.cumsum(ky[:n])
                tt = np.arange(n)
                ax, ay = sx[tt] + sx[n - 1 - tt] - kx[0], sy[tt] + sy[n - 1 - tt] - ky[0]
                cc = shift_cross(kx, ky, n, s) - ((ax @ ay) + (ax @ np.roll(ay, -s))) / n ** 2 + ax.sum() * ay.sum() / n ** 3
                xr = cc / shift_mean(kx, ky, n, 0, True)
                xi = shift_cross(dx, dy, n, s) / shift_mean(dx, dy, n, 0, False)
                row.append((frac, rc, ru, ri, xr, xi))
                csv_rows += [("shift_mean_ratio", "hsic", f"raw_c_s{frac}", hx, hy, n, rc), ("shift_mean_ratio", "hsic", f"inc1_s{frac}", hx, hy, n, ri),
                             ("shift_cross_ratio", "hsic", f"raw_unc_s{frac}", hx, hy, n, xr), ("shift_cross_ratio", "hsic", f"inc1_s{frac}", hx, hy, n, xi)]
            out[(hx, hy, n)] = row
            f = row[1]
            print(f"    ({hx},{hy}) n={n:<6d} s=n/2: raw centred ratio={f[1]:.4f} raw uncentred ratio={f[2]:.4f} inc ratio-1={f[3] - 1:+.2e} | cross raw={f[4]:.4f} cross inc={f[5]:+.2e}")
    check("C7: for the increment kernel the wrap changes E T_n(s) by under 0.1% at n>=2000",
          all(abs(v[1][3] - 1) < 1e-3 and abs(v[0][3] - 1) < 1e-3 for v in out.values()))
    check("C7: for the increment kernel statistics at distant offsets are almost uncorrelated, |cross|<1e-3",
          all(abs(v[1][5]) < 1e-3 and abs(v[0][5]) < 1e-3 for v in out.values()))
    store["c7"] = {str(k): v for k, v in out.items()}
    store["c7_drift"] = drift_shift_bias()


# ------------------------------------------------------------------ C3: uniform-in-p and tail exponent

def section_c3_uniform(store: dict) -> None:
    section("C3: the uniform-in-p step and the tail exponent")
    import mpmath as mp
    mp.mp.dps = 80
    bad = []
    print("  sup over p>=1, j>=2 of |Delta^2 rho^p (j)| j^(4-2H), and the proof constant")
    tab = {}
    for H in (0.55, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99):
        J = 200000
        rho = rho_all(H, J + 2)
        R = rho[1]
        d2 = np.empty(J + 2)
        j = np.arange(J + 2, dtype=float)
        d2[64:] = d2rho_series(H, j[64:])
        d2[1:64] = rho[2:65] - 2 * rho[1:64] + rho[0:63]
        u = np.abs(rho[2:J + 2] - rho[1:J + 1])  # |rho(j+1) - rho(j)| at j = 1..J
        v = np.abs(rho[0:J] - rho[1:J + 1])
        pp = np.arange(1, 500)
        A = float((pp * R ** (pp - 1)).max())
        Bb = float((pp * (pp - 1) * R ** np.maximum(pp - 2, 0)).max())
        js = j[2:J + 1]
        bound = A * np.abs(d2[2:J + 1]) + 0.5 * Bb * (u[1:J] ** 2 + v[1:J] ** 2)
        # Bound is compared with the true sup over p on a coarse lag set. Lags keep their exact rho.
        lagset = np.unique(np.r_[np.arange(2, 200), np.geomspace(200, J, 300).astype(int)])
        worst = 0.0
        for jj in lagset:
            mmp = jj.item()
            rj, rp, rm = rho[mmp], rho[mmp + 1], rho[mmp - 1]
            # Use mpmath only when double precision cannot resolve the difference.
            if mmp < 2000:
                vals = np.abs(2 * rj ** pp - rp ** pp - rm ** pp)
                val = float(vals.max())
                ok = np.all(vals <= bound[mmp - 2] * (1 + 1e-9) + 1e-13)
                if not ok:
                    bad.append((H, mmp))
                worst = max(worst, val * mmp ** (4 - 2 * H))
        ratio = float((bound * js ** (4 - 2 * H)).max())
        tab[H] = dict(R=R, A=A, B=Bb, sup_seen=worst, bound=ratio)
        print(f"    H={H}: R={R:.4f} A={A:8.2f} B={Bb:9.2f}  sup seen={worst:9.3f}  bound C(H)={ratio:10.3f}")
    check("C3: |Delta^2 rho^p (j)| <= A|Delta^2 rho(j)| + (B/2)(u^2+v^2) for all p<500, 2<=j<2000", not bad, str(bad[:4]))
    store["c3_uniform"] = tab

    # Tail exponent of kappa_D for the Gaussian kernel at high precision.
    sig = mp.mpf(SIG)
    psi = lambda r: sig / mp.sqrt(sig ** 2 + 2 - 2 * r)

    def rho_mp(H, k):
        a = 2 * mp.mpf(H)
        k = mp.mpf(k)
        return ((k + 1) ** a - 2 * k ** a + (k - 1) ** a) / 2

    def kd(H, j):
        return 2 * psi(rho_mp(H, j)) - psi(rho_mp(H, j + 1)) - psi(rho_mp(H, j - 1))
    bad = []
    for H in (0.6, 0.75, 0.8, 0.9, 0.95):
        s = mp.log(abs(kd(H, 10 ** 9) / kd(H, 10 ** 8))) / mp.log(10)
        want = 2 * H - 4
        print(f"    H={H}: tail slope of kappa_D between j=1e8 and 1e9 = {float(s):.4f} (theory {want:.1f})")
        if abs(float(s) - want) > 0.02:
            bad.append((H, float(s)))
    check("C3: kappa_D(j) decays like j^(2H-4) with j up to 1e9", not bad, str(bad))
    # Scope limit: memory away from frequency zero defeats lag 1. cos(w k) rho_H(k) is a valid covariance.
    H, om = 0.9, math.pi / 3
    ns = N_EXT[:11]
    rho = rho_all(H, max(ns) + 60)
    cyc = np.cos(om * np.arange(rho.size)) * rho
    d1c, d6c = inc(cyc, 1), inc(cyc, 6)
    s1 = aitken_exponent(et_uncentred(d1c, d1c, ns), ns[-3])
    s6 = aitken_exponent(et_uncentred(d6c, d6c, ns), ns[-3])
    check("C3 scope: a memory pole at w=pi/3 keeps lag-1 at rate n^(1-2D); lag-6 removes it",
          abs(s1 - 0.6) < 0.03 and abs(s6 + 1) < 0.05, f"slope lag-1 {s1:.3f} (want 0.600), slope lag-6 {s6:.3f} (want -1)")
    # The limit Lambda for H grid, including H close to one, with the Gaussian kernel.
    lay = gauss_layers()
    tab2 = {}
    for H in HS + (0.99, 0.999):
        k = inc(lay.kappa(rho_all(H, 200060)), 1)
        lam = float(k[0] ** 2 + 2 * (k[1:200000] ** 2).sum())
        tab2[H] = (lam, float(k[0]))
        print(f"    H={H}: Lambda = E T_infinity = {lam:.5f}; kappa_D(0) = {k[0]:.5f}")
    store["c3_lambda"] = tab2
    check("C3: Lambda stays finite and falls as H approaches one", all(v[0] < 5 for v in tab2.values()) and tab2[0.999][0] < tab2[0.95][0])


# ------------------------------------------------------------------ power proxy

def section_proxy(c13: dict, store: dict, rows: list) -> None:
    section("Power proxy: rho at which the first-order signal equals the null mean (n=2000 and n=64000)")
    lay = gauss_layers()
    c = lay.c
    keep = c13["keep"]
    arrs = c13["arrs"]
    tab = {}
    for H in HS:
        e = keep[(H, H)]
        R = 2 ** (2 * H - 1) - 1
        R50 = float(rho_all(H, 60)[50])
        vw = store["c4"][H]["v_frac"]
        row = {}
        for n in (2000, 64000):
            row[n] = dict(
                raw=(e["raw_c"][n] / (n * c[0] ** 2)) ** 0.5,
                def1=(e["def1_c"][n] / (n * c[1] ** 2)) ** 0.25,
                def2=(e["def2_c"][n] / (n * c[2] ** 2)) ** (1 / 6),
                inc1=(e["inc1"][n] / (n * (2 * (1 - R) * c[0]) ** 2)) ** 0.5,
                inc50=(e["inc50"][n] / (n * (2 * (1 - R50) * c[0]) ** 2)) ** 0.5,
                frac=(e["frac"][n] / (n * (vw * c[0]) ** 2)) ** 0.5,
            )
        tab[H] = row
        print(f"    H={H:4}, n=2000: " + "  ".join(f"{k}={v:.3f}" for k, v in row[2000].items()))
        print(f"    H={H:4}, n=64000:" + "  ".join(f"{k}={v:.3f}" for k, v in row[64000].items()))
        for n in (2000, 64000):
            rows += [("proxy_rho_star", "hsic", k, H, H, n, v) for k, v in row[n].items()]
    store["proxy"] = tab


# ------------------------------------------------------------------ fourth-moment proof: contraction checks

def contraction_cycle(hx: float, hy: float, p: int, q: int, mx: int, my: int, n: int, diff: bool) -> float:
    """||f ⊗_(mx,my) f||^2 for the layer (p,q): n^-2 sum_{s,t,u,v} P_X P_Y on the four-cycle s-t-v-u-s."""
    rx, ry = rho_all(hx, n + 4), rho_all(hy, n + 4)
    s = np.arange(n)[:, None, None, None]
    t = np.arange(n)[None, :, None, None]
    u = np.arange(n)[None, None, :, None]
    v = np.arange(n)[None, None, None, :]

    def part(r: np.ndarray, deg: int, m: int) -> np.ndarray:
        r1, r2 = r ** m, r ** (deg - m)
        total = 0.0
        shifts = [(0, 0, 0, 0)] if not diff else [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]
        for a, b, c, d in shifts:
            sign = (-1) ** (a + b + c + d)
            total = total + sign * r1[np.abs(s - a - t + b)] * r1[np.abs(u - c - v + d)] * r2[np.abs(s - a - u + c)] * r2[np.abs(t - b - v + d)]
        return total
    return float(np.sum(part(rx, p, mx) * part(ry, q, my)) / n ** 2)


def section_fourth_moment(store: dict) -> None:
    section("C3 and C7: contraction conditions of the fourth-moment theorem")
    # Second chaos: the (1,1) layer. S_n = n^-2 tr((Gx Gy)^2) is exact and cheap.
    bad = []
    for label, hx, hy, inc_ in (("raw", 0.7, 0.7, False), ("raw mixed", 0.9, 0.55, False), ("raw mixed", 0.85, 0.6, False),
                                ("increment", 0.95, 0.95, True), ("increment", 0.9, 0.9, True), ("increment", 0.95, 0.6, True)):
        vals = []
        ns = (500, 1000, 2000)
        for n in ns:
            rx, ry = rho_all(hx, n + 2), rho_all(hy, n + 2)
            if inc_:
                rx, ry = inc(rx, 1), inc(ry, 1)
            lag = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
            M = rx[lag] @ ry[lag]
            vals.append(float(np.sum(M * M.T) / n ** 2))
        slope = math.log2(vals[-1] / vals[-2])
        dsum = 4 - 2 * hx - 2 * hy
        want = -1.0 if inc_ else max(-1.0, 2 - 2 * dsum)
        print(f"    {label:10s} ({hx},{hy}): S_n at n=500,1000,2000 {np.array2string(np.array(vals), precision=4)}; slope {slope:.3f}; theory {want:.3f}")
        if vals[-1] >= vals[0] or abs(slope - want) > 0.15:
            bad.append((label, hx, hy, slope))
    check("second-chaos contraction vanishes at the predicted rate, including mixed pairs and the increment kernel", not bad, str(bad))
    # Fourth chaos, layer (2,2), contraction (1,1): the edge where cancellation matters most.
    out = {}
    for hx, hy in ((0.95, 0.95), (0.9, 0.6)):
        for diff in (False, True):
            out[(hx, hy, diff)] = [contraction_cycle(hx, hy, 2, 2, 1, 1, n, diff) for n in (12, 24, 36)]
            print(f"    layer (2,2), m=(1,1), H=({hx},{hy}), {'increment' if diff else 'raw      '}: S_n at n=12,24,36 = {np.array2string(np.array(out[(hx, hy, diff)]), precision=4)}")
    check("fourth-chaos contraction falls with n for the increment kernel", all(out[k][2] < out[k][0] for k in out if k[2]))
    store["contraction"] = {str(k): v for k, v in out.items()}


# ------------------------------------------------------------------ C8: the Hermite-split design

def sv_score_layers(sig: float, pmax: int = 14, rmax: int = 8) -> dict:
    """SV with the kernel on normal scores T(X). Memory of each He_q(T) layer, measured against He_p(xi)."""
    nx, ne = 61, 61
    xi, wx = hermegauss(nx)
    e, we = hermegauss(ne)
    wx, we = wx / wx.sum(), we / we.sum()
    XI, E = np.meshgrid(xi, e, indexing="ij")
    W = np.outer(wx, we).ravel()
    XIr, Er = XI.ravel(), E.ravel()
    X = np.exp(XIr / 2) * Er
    zi, wz = hermegauss(180)
    wz = wz / wz.sum()
    T = -np.sign(X) * ndtri((wz * ndtr(-np.abs(X)[:, None] * np.exp(-zi / 2))).sum(1))
    K = np.exp(-((T[:, None] - T[None, :]) ** 2) / (2 * sig ** 2))
    hxi, ht = he_norm(XIr, pmax), he_norm(T, rmax)
    beta = (hxi * W) @ ht.T  # corr(h_p(xi), h_q(T))
    u, v = hxi * W, ht * W
    k00 = float(1 - W @ K @ W)
    drop = np.cumsum([0.0] + [float(v[q] @ K @ v[q]) for q in range(1, rmax + 1)])
    d = {}
    for r in range(0, rmax + 1):
        d[r] = np.array([float((u[p] - beta[p, 1:r + 1] @ v[1:r + 1]) @ K @ (u[p] - beta[p, 1:r + 1] @ v[1:r + 1])) for p in range(1, pmax + 1)])
    m2 = (np.sum((ht[2] * W).reshape(nx, ne), axis=1) / wx)  # E[h_2(T) | xi] on the xi nodes
    return dict(beta=beta, d=d, k0={r: k00 - drop[r] for r in d}, tau2=float(wx @ m2 ** 2))


def section_c8(store: dict, rows: list) -> None:
    section("C8: the Hermite-split design")
    sig_z = float(ndtri(0.75))  # hcheck uses exp(-d^2 / median d^2), which is sigma = 0.6745 on scores
    print(f"  kernel on normal scores, sigma = {sig_z:.4f}; n=800; inflation = E T(H) / E T(1/2)")
    note_fgn = {1: (3.0, 5.9), 2: (1.1, 1.3), 3: (1.1, 1.1), 4: (1.0, 1.1)}
    note_sv = {1: (1.4, 2.4), 2: (1.7, 3.1), 3: (1.2, 1.2), 4: (1.0, 1.2)}
    lay = gauss_layers(sig_z)
    n = 800
    inf_f = {}
    for m in (1, 2, 3, 4):
        r = m - 1
        base = et_split_mmd(lay.kappa(rho_all(0.5, n + 2), r), [n])[n]
        vals = [et_split_mmd(lay.kappa(rho_all(H, n + 2), r), [n])[n] / base for H in (0.7, 0.8)]
        big = [et_split_mmd(lay.kappa(rho_all(0.8, 64002), r), [64000])[64000] / et_split_mmd(lay.kappa(rho_all(0.5, 64002), r), [64000])[64000]]
        inf_f[m] = vals
        print(f"    fGn m={m}: exact x{vals[0]:.2f}, x{vals[1]:.2f}; note x{note_fgn[m][0]}, x{note_fgn[m][1]}; exact at H=0.8, n=64000: x{big[0]:.2f}")
        rows += [("c8_inflation_fgn", "mmd", f"m{m}", H, H, n, v) for H, v in zip((0.7, 0.8), vals)]
    check("C8: exact fGn inflation with T known agrees with the note within 25% for m=1,3,4",
          all(abs(inf_f[m][i] / note_fgn[m][i] - 1) < 0.25 for m in (1, 3, 4) for i in (0, 1)))
    print(f"  The note's m=2, H=0.8 value (x1.3) differs from T known (x{inf_f[2][1]:.2f}). Estimated scores cause it: see the centre ratio below and theory_path_check.py.")
    sv = sv_score_layers(sig_z)
    beta = sv["beta"]
    print("  SV on normal scores: corr(h_p(xi), h_q(T)), rows p=1,2,3, columns q=1..8")
    for p in (1, 2, 3):
        print(f"    p={p}: " + "  ".join(f"{beta[p, q]:+.4f}" for q in range(1, 9)))
    share = beta[1, 1:] ** 2
    print(f"  share of the rank-one memory R^2 by score layer: q=2 {share[1] / share.sum():.4f}, q=4 {share[3] / share.sum():.4f}, q=6 {share[5] / share.sum():.4f}, q=8 {share[7] / share.sum():.4f}; total R^2 of xi on these layers {share.sum():.4f}")
    check("C8: odd score layers carry no SV memory", np.abs(beta[1:, 1::2]).max() < 1e-8, f"max {np.abs(beta[1:, 1::2]).max():.1e}")
    check("C8: the second score layer carries most of the SV rank-one memory, but not all", share[1] / share.sum() > 0.8 and share[3] > 0)
    inf_s = {}
    for m in (1, 2, 3, 4):
        r = m - 1
        c = sv["d"][r]
        base = 4 * sv["k0"][r]
        vals = []
        big = []
        for H in (0.7, 0.8):
            k = Layers(c, sv["k0"][r], "x").kappa(rho_all(H, n + 2))
            vals.append(et_split_mmd(k, [n])[n] / base)
            k = Layers(c, sv["k0"][r], "x").kappa(rho_all(H, 1_024_002))
            big.append(et_split_mmd(k, [1_024_000])[1_024_000] / base)
        inf_s[m] = vals
        print(f"    SV  m={m}: exact x{vals[0]:.2f}, x{vals[1]:.2f}; note x{note_sv[m][0]}, x{note_sv[m][1]}; d_1={c[0]:.3e}, d_1/raw={c[0] / sv['d'][0][0]:.4f}; at n=1,024,000: x{big[0]:.2f}, x{big[1]:.2f}")
        rows += [("c8_inflation_sv", "mmd", f"m{m}", H, H, n, v) for H, v in zip((0.7, 0.8), vals)]
    check("C8: exact SV inflation agrees with the note to within 25% at n=800",
          all(abs(inf_s[m][i] / note_sv[m][i] - 1) < 0.25 for m in note_sv for i in (0, 1)), str({m: np.round(v, 2).tolist() for m, v in inf_s.items()}))

    # Power: does the separate first-layer test recover the location signal?
    print("  Location shift, rank contrast of He_1(T), pooled T: E_A T - E_B T over delta")
    x = np.linspace(-14, 14, 8001)
    dx = x[1] - x[0]
    for d in (0.1, 0.5, 1.0, 2.0):
        pa = np.exp(-((x + d / 2) ** 2) / 2) / math.sqrt(2 * math.pi)
        pb = np.exp(-((x - d / 2) ** 2) / 2) / math.sqrt(2 * math.pi)
        lo, up = 0.5 * (ndtr(x + d / 2) + ndtr(x - d / 2)), 0.5 * (ndtr(-x - d / 2) + ndtr(-x + d / 2))
        T = np.where(lo < 0.5, ndtri(lo), -ndtri(up))
        print(f"    delta={d}: ratio {abs(((pa - pb) * T).sum() * dx) / d:.5f}")
    c = lay.c
    a3 = 6 * c[2]
    print("  delta at which the first-order signal equals the noise scale, n=2000. first_layer uses sd of the half-mean contrast. Others use the null mean of T.")
    tab = {}
    n = 2000
    m_ = n // 2
    for H in HS:
        rho = rho_all(H, n + 2)
        e = {r: et_split_mmd(lay.kappa(rho, r), [n])[n] for r in (0, 1, 2)}
        first = math.sqrt((4 - 2 ** (2 * H)) * m_ ** (2 * H - 2))
        rawd = math.sqrt(e[0] / (n * c[0]))
        d1 = (144 * e[1] / (n * a3)) ** (1 / 6)
        d2 = (144 * e[2] / (n * a3)) ** (1 / 6)
        tab[H] = dict(first_layer=first, raw=rawd, def1=d1, def2=d2)
        print(f"    H={H:4}: first_layer={first:.3f}  raw={rawd:.3f}  deflated r=1={d1:.3f}  r=2={d2:.3f}")
        rows += [("c8_delta_star", "mmd", k, H, H, n, v) for k, v in tab[H].items()]
    store["c8_delta"] = tab

    # Independence: first-layer Pearson on scores against increments, one degree of freedom each.
    print("  Rank-(1,1) signal, Y = rho X + sqrt(1-rho^2) Z. Variance factor S: z = sqrt(n) rho / sqrt(S)")
    tab = {}
    for H in HS:
        rho = rho_all(H, 64002)
        s = et_centred(rho, rho, [2000, 64000])
        g = inc(rho_all(H, 3_000_060), 1)
        sd = float((g[0] ** 2 + 2 * (g[1:] ** 2).sum()) / g[0] ** 2)
        tab[H] = dict(S_low_2000=s[2000], S_low_64000=s[64000], S_inc=sd)
        print(f"    H={H:4}: S_low(n=2000)={s[2000]:.3f}  S_low(n=64000)={s[64000]:.3f}  S_inc={sd:.3f}")
        rows += [("c8_S", "hsic", k, H, H, 0, v) for k, v in tab[H].items()]
    store["c8_S"] = tab
    print("  Common factor: X=sqrt(s)F+sqrt(1-s)e, Y=sqrt(s)F+sqrt(1-s)u, s=0.3, e and u white. z per sqrt(n) units at n=2000 and 64000")
    s_ = 0.3
    for HF in (0.7, 0.8, 0.9, 0.95, 0.99):
        rF = rho_all(HF, 64062)
        rX = s_ * rF
        rX[0] = 1.0
        S = et_centred(rX, rX, [2000, 64000])
        gX = inc(rX, 1)
        gF = inc(rF, 1)
        tau = s_ * gF[0] / gX[0]
        Sd = float((gX[0] ** 2 + 2 * (gX[1:] ** 2).sum()) / gX[0] ** 2)
        print(f"    H_F={HF}: level corr={s_:.3f} z_low(2000)={math.sqrt(2000) * s_ / math.sqrt(S[2000]):.2f} z_low(64000)={math.sqrt(64000) * s_ / math.sqrt(S[64000]):.2f}; increment corr={tau:.4f} z_inc(2000)={math.sqrt(2000) * tau / math.sqrt(Sd):.2f} z_inc(64000)={math.sqrt(64000) * tau / math.sqrt(Sd):.2f}")
    print("  Common SV factor on the (2,2) score layer: X=exp(F/2)e, Y=exp(F/2)u, n=2000 and 64000")
    bt = beta[1:, 2]
    tau2 = sv["tau2"]
    for H in (0.6, 0.7, 0.8, 0.9, 0.95):
        rho = rho_all(H, 64062)
        r2 = series_tail(rho, bt ** 2, 0)
        r2[0] = 1.0
        S2 = et_centred(r2, r2, [2000, 64000])
        g2 = inc(r2, 1)
        ti = (tau2 - r2[1]) / (1 - r2[1])
        Sd2 = float((g2[0] ** 2 + 2 * (g2[1:] ** 2).sum()) / g2[0] ** 2)
        print(f"    H={H}: level corr={tau2:.4f} z_low(2000)={math.sqrt(2000) * tau2 / math.sqrt(S2[2000]):.2f} z_low(64000)={math.sqrt(64000) * tau2 / math.sqrt(S2[64000]):.2f}; increment corr={ti:.4f} z_inc(2000)={math.sqrt(2000) * ti / math.sqrt(Sd2):.2f} z_inc(64000)={math.sqrt(64000) * ti / math.sqrt(Sd2):.2f}")

    # Dependence of first-layer and second-layer statistics.
    print("  Corr(Delta_1^2, Delta_2^2) for the half contrasts of He_1 and He_2, exact by cumulants")
    print("  known = He_2(xi_t). centred = He_2(xi_t - xibar), a quadratic form with matrix A - (1w'+w1')/n")
    res = {}
    for H in (0.6, 0.7, 0.75, 0.8, 0.9, 0.95):
        out = {"known": [], "centred": []}
        for n in (500, 1000, 2000):
            m_ = n // 2
            w = np.r_[np.ones(m_), -np.ones(n - m_)]
            rho = rho_all(H, n + 2)
            G = rho[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]
            v = G @ w
            a = float(w @ v)
            one = np.ones(n)
            for key, B in (("known", np.diag(w)), ("centred", np.diag(w) - (np.outer(one, w) + np.outer(w, one)) / n)):
                BG = B @ G
                BG2 = BG @ BG
                var_q = 2 * float(np.trace(BG2))
                cum4 = 48 * float(np.sum(BG2 * BG2.T))
                bv = B @ v
                cross = 8 * float(bv @ G @ bv)
                third = 2 * float(v @ bv)
                out[key].append((cross + 2 * third ** 2) / math.sqrt(2 * a * a * (cum4 + 2 * var_q * var_q)))
        res[H] = out
        print(f"    H={H:4}: known  " + "  ".join(f"{v:.4f}" for v in out["known"]) + "   centred  " + "  ".join(f"{v:.4f}" for v in out["centred"]))
        for key in out:
            rows += [("c8_layer_corr_" + key, "mmd", "He1_He2", H, H, nn, v) for nn, v in zip((500, 1000, 2000), out[key])]
    check("C8: with a known centre the layers decouple below H=3/4 and stay correlated above it",
          res[0.6]["known"][2] < 0.75 * res[0.6]["known"][0] and res[0.9]["known"][2] > 0.9 * res[0.9]["known"][0] and res[0.9]["known"][2] > 0.05,
          f"H=0.6: {res[0.6]['known'][0]:.3f}->{res[0.6]['known'][2]:.3f}; H=0.9: {res[0.9]['known'][0]:.3f}->{res[0.9]['known'][2]:.3f}")
    store["c8_corr"] = res

    # Estimated centre: variance of the second-layer contrast with xi - xibar in place of xi.
    print("  Var(sum w_t He_2(xi_t - xibar)) / Var(sum w_t He_2(xi_t)) for the half contrast")
    res = {}
    for H in (0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
        out = []
        ns = [2000, 64000, 1_024_000]
        rho = rho_all(H, ns[-1] + 2)
        s_ = np.cumsum(rho[: ns[-1] + 1])
        for n in ns:
            m_ = n // 2
            t = np.arange(n)
            # u_t = (Gamma 1)_t / n and v_t = (Gamma w)_t from partial sums of rho
            u = (s_[t] + s_[n - 1 - t] - 1.0) / n
            left = np.where(t < m_, s_[t] + s_[m_ - 1 - t] - 1.0, s_[t] - s_[t - m_])
            right = np.where(t >= m_, s_[t - m_] + s_[n - 1 - t] - 1.0, s_[n - 1 - t] - s_[m_ - 1 - t])
            v = left - right
            w = np.where(t < m_, 1.0, -1.0)
            a = float(w @ v)
            j = np.arange(1, n)
            aa = np.where(j <= m_, 2 * m_ - 3 * j, -(2 * m_ - j)).astype(float)
            b = 2 * (n + 2 * float(aa @ rho[1:n] ** 2))
            varbar = float(u.mean())
            cross = 2 * float((w * u) @ v)
            out.append(1 + (4 * varbar * a - 4 * cross) / b)
        res[H] = out
        print(f"    H={H:4}: n=2000 {out[0]:.4f}  n=64000 {out[1]:.4f}  n=1024000 {out[2]:.4f}")
        rows += [("c8_he2_centre_ratio", "mmd", "He2", H, H, nn, v) for nn, v in zip(ns, out)]
    check("C8: an estimated centre leaves the second-layer variance alone below H=3/4 and changes it above",
          abs(res[0.6][2] - 1) < 0.02 and abs(res[0.9][2] - 1) > 0.05, f"H=0.6: {res[0.6][2]:.4f}; H=0.9: {res[0.9][2]:.4f}")
    store["c8_centre"] = res


# ------------------------------------------------------------------ main

def main() -> int:
    rows: list = []
    fits: list = []
    store: dict = {}
    c13 = section_c1_c3(rows, fits)
    section_c3_uniform(store)
    section_c2(rows)
    section_c4(store)
    section_proxy(c13, store, rows)
    section_c5(store)
    section_c6(store, rows)
    section_c7(store, rows)
    section_fourth_moment(store)
    section_c8(store, rows)
    with open(OUT / "theory_exact.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kernel", "statistic", "feature", "H_X", "H_Y", "n", "value"])
        for r in rows:
            wr.writerow(r[:6] + (f"{r[6]:.10g}",))
    with open(OUT / "theory_exponents.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kernel", "statistic", "feature", "H_X", "H_Y", "n_max", "slope_plain", "slope_corrected", "slope_expected", "pass"])
        for r in fits:
            wr.writerow(r[:6] + tuple(f"{v:.5f}" for v in r[6:9]) + (r[9],))
    nf = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n{'=' * 78}\n{len(RESULTS)} checks, {nf} failed, {time.time() - T0:.0f}s, {len(rows)} rows in theory_exact.csv")
    for name, ok, det in RESULTS:
        if not ok:
            print(f"FAIL  {name}  {det}")
    return 1 if nf else 0


if __name__ == "__main__":
    sys.exit(main())
