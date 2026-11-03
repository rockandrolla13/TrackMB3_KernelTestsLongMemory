"""Unit tests for kmem.memory_kernels and the fast paths moved from EXP7 into kmem.

Run: PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -I -m pytest tests/test_memory_kernels.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

from kmem import bandwidth, features, hsic, kernels, layers, mmd, rules
from kmem import memory_kernels as mk
from kmem.calibration import harris, selfnorm, shift, subsample, wild
from kmem.generators import independence
from kmem.legacy import mbcore as mb

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "experiments", "EXP7-memory-kernel-20261101"))
import exp7_core as K  # noqa: E402  (kernel switches and the EXP7 bundle; they move in steps 3.2 and 5.2)

N = 240


def rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(entropy=20261101, spawn_key=(seed,)))


def series(seed: int = 1, H: float = 0.8, n: int = N) -> np.ndarray:
    return mb.fgn(H, n, rng(seed))


def gaussian_gram(x: np.ndarray) -> np.ndarray:
    d = np.abs(x[:, None] - x[None, :])
    return mb.gauss_gram(d, mb.median_bandwidth(d))


# ---------------------------------------------------------------- Gram transforms

@pytest.mark.parametrize("kernel", ["raw", "inc1", "inc50", "herm1", "herm2", "sym", "datadiff", "hp10", "dcov_inc1", "frac_oracle"])
@pytest.mark.parametrize("H", [0.5, 0.9])
def test_transformed_gram_is_psd_and_centred(kernel: str, H: float) -> None:
    x = series(2, H)
    d = np.abs(x[:, None] - x[None, :])
    g = kernels.get(kernel).gram(x, mb.median_bandwidth(d), H)
    assert np.allclose(g, g.T, atol=1e-10)
    ev = np.linalg.eigvalsh(g)
    assert ev.min() >= -1e-9 * max(ev.max(), 1.0), f"{kernel}: min eigenvalue {ev.min()}"
    assert ev.max() > 1e-6, f"{kernel}: Gram matrix is numerically zero"
    assert np.abs(g.sum(axis=0)).max() < 1e-8, f"{kernel}: rows do not sum to zero"


def test_double_center_matches_legacy() -> None:
    k = gaussian_gram(series(3))
    assert np.array_equal(mk.double_center(k), mb.double_center(k))


@pytest.mark.parametrize("lag", [1, 2, 50])
def test_increment_equals_feature_difference_for_linear_kernel(lag: int) -> None:
    x = series(4)
    k = np.outer(x, x)
    dx = x[lag:] - x[:-lag]
    assert np.allclose(mk.increment_raw(k, lag), np.outer(dx, dx), atol=1e-12)
    dc = dx - dx.mean()
    assert np.allclose(mk.increment(k, lag), np.outer(dc, dc), atol=1e-12)
    assert mk.increment(k, lag).shape == (N - lag, N - lag)


def test_increment_equals_explicit_d_matrix() -> None:
    k = gaussian_gram(series(5))
    d = np.zeros((N - 1, N))
    d[np.arange(N - 1), np.arange(1, N)] = 1.0
    d[np.arange(N - 1), np.arange(N - 1)] = -1.0
    assert np.allclose(mk.increment_raw(k, 1), d @ k @ d.T, atol=1e-12)


def test_normal_scores_are_symmetric_and_monotone() -> None:
    x = series(6)
    z = mk.normal_scores(x)
    assert np.all(np.diff(z[np.argsort(x)]) > 0)
    assert abs(z.mean()) < 1e-12
    assert abs(np.sort(z)[0] + np.sort(z)[-1]) < 1e-12


def test_hermite_design_recursion() -> None:
    z = np.array([-1.5, 0.0, 0.7, 2.0])
    h = mk.hermite_design(z, 3)
    assert np.allclose(h[:, 1], z)
    assert np.allclose(h[:, 2], z ** 2 - 1)
    assert np.allclose(h[:, 3], z ** 3 - 3 * z)


def test_herm1_leaves_almost_no_linear_signal() -> None:
    x = mb.fgn(0.5, 2000, rng(7))
    k = np.outer(x, x)
    g = mk.hermite(k, x, 1)
    share = np.trace(g) / np.trace(k)
    assert share < 0.02, f"share of the linear kernel left after herm1: {share}"
    assert np.trace(mk.double_center(k)) / np.trace(k) > 0.99


def test_herm2_removes_the_square_of_the_scores() -> None:
    x = mb.fgn(0.5, 2000, rng(8))
    z = mk.normal_scores(x)
    he2 = z * z - 1.0
    assert np.trace(mk.hermite(np.outer(he2, he2), x, 2)) < 1e-18 * np.trace(np.outer(he2, he2)) + 1e-9
    assert np.trace(mk.hermite(np.outer(he2, he2), x, 1)) > 0.99 * np.trace(np.outer(he2, he2))


def test_project_out_equals_explicit_projector() -> None:
    x = series(9)
    k = gaussian_gram(x)
    u = mk.hermite_basis(x, 2)
    q = np.eye(N) - u @ u.T
    assert np.allclose(mk.project_out(k, u), q @ k @ q, atol=1e-11)
    assert np.allclose(mk.hermite(k, x, 2), q @ k @ q, atol=1e-11)


def test_dct_basis_is_orthonormal_and_low_frequency() -> None:
    from scipy import fft
    b = mk.dct_basis(N, mk.n_low(N, 0.1))
    assert b.shape == (N, 24)
    assert np.allclose(b.T @ b, np.eye(24), atol=1e-12)
    e = np.eye(N)
    dct = fft.dct(e, type=2, norm="ortho", axis=0)
    assert np.allclose(b.T, dct[:24], atol=1e-12)
    assert np.ptp(b[:, 0]) < 1e-12


def test_highpass_removes_constant_and_slow_cosine() -> None:
    t = np.arange(N)
    slow = np.cos(np.pi * 3 * (t + 0.5) / N)
    k = np.outer(slow + 1.0, slow + 1.0)
    assert np.abs(mk.highpass(k, 0.1)).max() < 1e-10


def test_fracdiff_weights_and_special_orders() -> None:
    assert np.allclose(mk.fracdiff_weights(1.0, 5), [1, -1, 0, 0, 0])
    assert np.allclose(mk.fracdiff_weights(0.0, 4), [1, 0, 0, 0])
    assert np.allclose(mk.fracdiff_weights(0.3, 3), [1, -0.3, -0.105])


def test_fracdiff_gram_equals_explicit_filter_matrix() -> None:
    k = gaussian_gram(series(10))
    d = 0.3
    pi = mk.fracdiff_weights(d, N)
    f = np.zeros((N, N))
    for t in range(N):
        f[t, :t + 1] = pi[t::-1]
    h = np.eye(N) - 1.0 / N
    assert np.allclose(mk.fracdiff(k, d), h @ f @ (h @ k @ h) @ f.T @ h, atol=1e-10)
    assert np.allclose(mk.fracdiff(k, 0.0), mk.double_center(k), atol=1e-14)


def test_even_map_and_first_difference() -> None:
    x = np.array([1.0, 3.0, 7.0, 2.0])
    assert np.allclose(mk.even_map(x), [1.5, 0.5, 4.5, 0.5])
    assert np.allclose(mk.first_difference(x), [2.0, 4.0, -5.0])


def test_distance_kernel_increment_is_positive_semidefinite_on_differences() -> None:
    x = series(11)
    g = mk.increment(mk.distance_kernel(x), 1)
    assert np.linalg.eigvalsh(g).min() > -1e-9


# ---------------------------------------------------------------- telescoping

@pytest.mark.parametrize("n,h", [(200, 99), (201, 100), (64, 10)])
def test_increment_kernel_telescopes_in_split_mmd(n: int, h: int) -> None:
    x = series(12, 0.8, n)
    k = gaussian_gram(x)
    m = n - 1
    s = mmd.contrast(m, h)
    v = mk.telescoped_contrast(n, h)
    assert np.count_nonzero(v) == 3
    assert abs(s @ mk.increment_raw(k, 1) @ s - v @ k @ v) < 1e-12
    # the statistic can never exceed the squared l1 norm of v, whatever the data: about 64 / n^2 for a mid split
    assert v @ k @ v <= np.abs(v).sum() ** 2 + 1e-12


def test_telescoped_statistic_shrinks_as_one_over_n_squared() -> None:
    for n in (100, 400, 1600):
        v = mk.telescoped_contrast(n, n // 2 - 1)
        assert (np.abs(v).sum() * (n - 1)) ** 2 < 66.0  # ||v||_1 = 2/h + 2/(m - h), about 8/m at a mid split


# ---------------------------------------------------------------- raw reproduces EXP6

@pytest.mark.parametrize("gen,alt,rho,H,n", [("G0", "null", 0.0, 0.8, 300), ("SEAS-V", "lin", 0.1, 0.5, 300)])
def test_raw_exact_reproduces_exp6_hsic_fields(gen: str, alt: str, rho: float, H: float, n: int) -> None:
    r6, r7 = rng(13), rng(13)
    x, y = independence.pair(gen, alt, H, H, n, rho, r6)
    independence.pair(gen, alt, H, H, n, rho, r7)
    ln, dn = math.floor(n ** 0.5), math.floor(0.1 * n)
    dx, dy = np.abs(x[:, None] - x[None, :]), np.abs(y[:, None] - y[None, :])
    w = mb.markov_rademacher(n, 500, ln, r6)  # EXP6 kind_bw draws its multipliers this way
    f = hsic.hsic_full_rules(x, y, dx, dy, mb.median_bandwidth(dx), mb.median_bandwidth(dy), w, n, ln, dn, dn)
    want = {f"{key}_k1": v for key, v in f.items()}
    got = K.ind_bundle(x, y, r7, kernels=("raw",), exact=True, hsplit=False)
    pairs = {"hsic_stat": "t", "sx": "sx", "sy": "sy", "p_hsic_wb": "p_wb", "r_hsic_wb": "r_wb", "p_hsic_shift": "p_shift",
             "r_hsic_shift": "r_shift", "p_hsic_sub": "p_sub", "rej_hsic_sub": "rej_sub", "r_hsic_sub": "r_sub",
             "shift_check": "check"}
    for a, b in pairs.items():
        assert float(want[f"{a}_k1"]) == float(got[f"xraw_{b}"]), a


# ---------------------------------------------------------------- fast paths equal full paths

def test_median_bw_fast_equals_full() -> None:
    for n in (240, 241):
        x = series(14, 0.7, n)
        full = mb.median_bandwidth(np.abs(x[:, None] - x[None, :]))
        assert abs(bandwidth.median_bw_fast(x) - full) <= 1e-11 * full


def test_wild_prefix_equals_mbcore() -> None:
    x, y = series(15), series(16)
    w = mb.markov_rademacher(N, 200, math.floor(N ** 0.5), rng(17))
    m = mk.double_center(np.abs(x[:, None] - x[None, :])) * mk.double_center(np.abs(y[:, None] - y[None, :]))
    a, b = mb.wild_boot_stats(m, w), wild.wild_prefix(m, w)
    assert np.abs(a - b).max() < 1e-10 * max(np.abs(a).max(), 1.0)
    v = w[:N - 7] - w[:N - 7].mean(axis=0, keepdims=True)
    m2 = m[:N - 7, :N - 7]
    assert np.abs(mb.wild_boot_stats(m2, v) - wild.wild_prefix(m2, v)).max() < 1e-10


@pytest.mark.parametrize("n", [240, 241])
def test_shift_stats_sym_equals_mbcore(n: int) -> None:
    a = mk.double_center(np.abs(series(18, 0.8, n)[:, None] - series(18, 0.8, n)[None, :]))
    b = mk.double_center(np.abs(series(19, 0.8, n)[:, None] - series(19, 0.8, n)[None, :]))
    assert np.abs(mb.shift_stats(a, b) - shift.shift_stats_sym(a, b)).max() < 1e-11


@pytest.mark.parametrize("kernel", ["raw", "inc1", "inc50", "herm1", "herm2", "hp10", "frac_oracle", "sym", "datadiff"])
def test_fast_rules_equal_full_rules(kernel: str) -> None:
    n = 600
    H = 0.8
    x, y = mb.fgn(H, n, rng(20)), mb.fgn(H, n, rng(21))
    sx, sy = bandwidth.median_bw_fast(x), bandwidth.median_bw_fast(y)
    mult = wild.Multipliers(n, 200, rng(22))
    kc, lc = kernels.get(kernel).gram(x, sx, H), kernels.get(kernel).gram(y, sy, H)
    if kernel in ("sym", "datadiff"):
        ux, uy = kernels.get(kernel).series(x), kernels.get(kernel).series(y)
        fx = features.gauss_features(ux, bandwidth.median_bw_fast(ux))
        fy = features.gauss_features(uy, bandwidth.median_bw_fast(uy))
    else:
        fx = kernels.get(kernel).features(features.gauss_features(x, sx), x, H)
        fy = kernels.get(kernel).features(features.gauss_features(y, sy), y, H)
    w = mult.get(kc.shape[0])
    full = rules.rules_full(kc, lc, w)
    fast = rules.rules_fast(features.compress(fx), features.compress(fy), w)
    assert abs(fast["t"] / full["t"] - 1.0) < 2e-5
    assert fast["p_wb"] == full["p_wb"]
    assert fast["p_shift"] == full["p_shift"]
    assert abs(fast["r_wb"] / full["r_wb"] - 1.0) < 2e-5
    assert abs(fast["r_shift"] / full["r_shift"] - 1.0) < 2e-5
    # the shift function at zero shift equals the statistic, on the full Gram matrices and in the fast path
    assert abs(full["check"]) < 1e-9
    assert abs(fast["check"]) < 1e-9


def test_rules_distance_equals_rules_full() -> None:
    n = 500
    x, y = mb.fgn(0.8, n, rng(23)), mb.fgn(0.8, n, rng(24))
    ka, la = mk.increment(mk.distance_kernel(x), 1), mk.increment(mk.distance_kernel(y), 1)
    mult = wild.Multipliers(n, 200, rng(25))
    a, b = rules.rules_full(ka, la, mult.get(n - 1)), rules.rules_distance(ka, la, mult.get(n - 1))
    for key in ("p_wb", "p_shift"):
        assert a[key] == b[key]
    for key in ("t", "r_wb", "r_shift"):
        assert abs(a[key] / b[key] - 1.0) < 1e-9


def _window_stats(gx: np.ndarray, gy: np.ndarray, l: int, dn: int) -> np.ndarray:
    """Brute force: double-centre each window's Gram matrices and sum their product."""
    m = gx.shape[0]
    out = []
    h = np.eye(l) - 1.0 / l
    for k in range(m - l - dn):
        a = h @ gx[k:k + l]
        b = h @ gy[k + dn:k + dn + l]
        out.append(((a @ a.T) * (b @ b.T)).sum() / l)
    return np.array(out)


@pytest.mark.parametrize("kernel", ["raw", "inc1"])
def test_sub_fast_equals_brute_force_windows(kernel: str) -> None:
    n = 400
    x, y = mb.fgn(0.8, n, rng(26)), mb.fgn(0.8, n, rng(27))
    fx = kernels.get(kernel).features(features.gauss_features(x, bandwidth.median_bw_fast(x)), x, 0.8)
    fy = kernels.get(kernel).features(features.gauss_features(y, bandwidth.median_bw_fast(y)), y, 0.8)
    (px, lx), (py, ly) = features.compress(fx), features.compress(fy)
    i, j = features.select_pairs(lx, ly)
    m = px.shape[0]
    l, dn = math.floor(m ** 0.5), math.floor(0.1 * m)
    fast = subsample.sub_fast(px, py, i, j, l, dn)
    brute = _window_stats(fx, fy, l, dn)
    assert fast.size == m - l - dn
    assert np.abs(fast - brute).max() < 1e-4 * brute.max()


def test_raw_sub_fast_equals_mbcore_subsample() -> None:
    n = 500
    x, y = mb.fgn(0.8, n, rng(28)), mb.fgn(0.8, n, rng(29))
    sx, sy = bandwidth.median_bw_fast(x), bandwidth.median_bw_fast(y)
    (px, lx), (py, ly) = features.compress(features.gauss_features(x, sx)), features.compress(features.gauss_features(y, sy))
    i, j = features.select_pairs(lx, ly)
    l, dn = math.floor(n ** 0.5), math.floor(0.1 * n)
    legacy = mb.subsample_stats(x, y, l, dn, "hsic", sx, sy)
    assert np.abs(subsample.sub_fast(px, py, i, j, l, dn) - legacy).max() < 1e-4 * legacy.max()


def test_multipliers_share_one_chain() -> None:
    mult = wild.Multipliers(2000, 50, rng(30))
    w = mult.get(1999)
    assert np.allclose(w.mean(axis=0), 0.0)
    assert np.allclose(w[1:] - w[:-1], mult.w[1:1999] - mult.w[:1998], atol=1e-12)
    assert mult.w_long is mult.w
    assert wild.Multipliers(500, 50, rng(31)).ln_long == 21


# ---------------------------------------------------------------- self-normalised statistics of hsplit

def test_sn_split_stat_is_betken_statistic_as_written() -> None:
    v = rng(32).standard_normal(301)
    k = 120
    num = v[:k].sum() - k / v.size * v.sum()
    s1, s2 = np.cumsum(v[:k] - v[:k].mean()), np.cumsum(v[k:] - v[k:].mean())
    g = abs(num) / math.sqrt(((s1 * s1).sum() + (s2 * s2).sum()) / v.size)
    assert abs(math.sqrt(selfnorm.sn_split_stat(v, k)) - g) < 1e-10


def test_sn_mean_stat_is_invariant_to_scale() -> None:
    u = rng(33).standard_normal(300)
    assert abs(selfnorm.sn_mean_stat(3.7 * u) - selfnorm.sn_mean_stat(u)) < 1e-8


def test_sn_windows_equal_direct_computation() -> None:
    v = rng(34).standard_normal(300)
    for l in (20, 21):
        direct = np.array([selfnorm.sn_split_stat(v[k:k + l]) for k in range(300 - l + 1)])
        assert np.abs(selfnorm.sn_split_windows(v, l) / direct - 1.0).max() < 1e-6
        direct = np.array([selfnorm.sn_mean_stat(v[k:k + l]) for k in range(300 - l + 1)])
        assert np.abs(selfnorm.sn_mean_windows(v, l) / direct - 1.0).max() < 1e-6


def test_hermite_layers_are_centred_and_match_the_design() -> None:
    x = series(35)
    v = layers.hermite_layers(x)
    z = mk.normal_scores(x)
    assert np.abs(v.mean(axis=0)).max() < 1e-12
    assert np.allclose(v[:, 0], z - z.mean())
    assert np.allclose(v[:, 1], z * z - (z * z).mean())


# ---------------------------------------------------------------- Harris (2020) shift test

@pytest.mark.parametrize("kernel", ["raw", "inc1"])
@pytest.mark.parametrize("big_n", [19, 49])
def test_harris_shifts_equal_brute_force(kernel: str, big_n: int) -> None:
    n = 400
    x, y = mb.fgn(0.8, n, rng(36)), mb.fgn(0.8, n, rng(37))
    sx, sy = bandwidth.median_bw_fast(x), bandwidth.median_bw_fast(y)
    fx = kernels.get(kernel).features(features.gauss_features(x, sx), x, 0.8)
    fy = kernels.get(kernel).features(features.gauss_features(y, sy), y, 0.8)
    kx, ky = fx @ fx.T, fy @ fy.T
    m = kx.shape[0]
    d = m - 2 * big_n
    a = mk.double_center(kx[big_n:m - big_n, big_n:m - big_n])
    brute = np.empty(2 * big_n + 1)
    for s in range(-big_n, big_n + 1):
        b = mk.double_center(ky[s + big_n:s + m - big_n, s + big_n:s + m - big_n])
        brute[s + big_n] = (a * b).sum() / d
    fast = harris.harris_shifts(features.compress(fx), features.compress(fy), big_n)
    assert np.abs(fast / brute - 1.0).max() < 1e-4


def test_harris_ratio_matches_the_counting_rule() -> None:
    r = rng(38)
    for big_n, limit in ((19, 1.0), (99, 5.0), (99, 9.95), (19, 1.95)):
        for _ in range(200):
            v = r.standard_normal(2 * big_n + 1) ** 2
            m = int(np.sum(np.delete(v, big_n) >= v[big_n]))
            assert (harris.harris_ratio(v, big_n, limit) > 1) == (m < limit)


def test_harris_lemma_counts_local_maxima() -> None:
    """Lemma 1 of the paper: at most M (M, N)-local maxima in any N + 1 consecutive points."""
    x = rng(39).standard_normal(600)
    for big_n, big_m in ((19, 1), (19, 3), (40, 2)):
        local = np.array([np.sum(np.delete(x[t - big_n:t + big_n + 1], big_n) >= x[t]) < big_m for t in range(big_n, x.size - big_n)])
        windows = np.convolve(local.astype(int), np.ones(big_n + 1, int), mode="valid")
        assert windows.max() <= big_m


def test_harris_big_n_gives_whole_number_limit() -> None:
    assert [harris.harris_big_n(n) for n in (500, 2000, 4000)] == [19, 99, 199]
    for n in (500, 2000, 4000):
        assert abs(0.05 * (harris.harris_big_n(n) + 1) - round(0.05 * (harris.harris_big_n(n) + 1))) < 1e-9


def test_harris_strict_form_is_close_to_the_sliding_form() -> None:
    n = 500
    x, y = mb.fgn(0.7, n, rng(40)), mb.fgn(0.7, n, rng(41))
    fx, fy = features.gauss_features(x, bandwidth.median_bw_fast(x)), features.gauss_features(y, bandwidth.median_bw_fast(y))
    strict = harris.harris_strict19(x, y)
    a = harris.harris_rules(features.compress(fx), features.compress(fy), n)
    assert abs(a["rh19"] / strict["xh19_raw"] - 1.0) < 0.1
