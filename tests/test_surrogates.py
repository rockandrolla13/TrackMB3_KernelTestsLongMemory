import numpy as np

from kmem.generators.registry import get_process
from kmem.surrogates import hurst_lw, iaaft, iaaft_drift, variance_trend


def _fgn(h, n, seed):
    return get_process("G0")(h, n, np.random.default_rng(seed))


def _pgram_corr(a, b):
    pa, pb = np.abs(np.fft.rfft(a - a.mean()))[1:], np.abs(np.fft.rfft(b - b.mean()))[1:]
    return float(np.corrcoef(np.log(pa ** 2), np.log(pb ** 2))[0, 1])


def test_iaaft_keeps_values_and_spectrum():
    y = _fgn(0.8, 2000, 1) ** 3  # non-Gaussian marginal
    s = iaaft(y, np.random.default_rng(2))
    assert np.array_equal(np.sort(s), np.sort(y))
    assert not np.array_equal(s, y)
    smooth = lambda v: np.convolve(np.abs(np.fft.rfft(v - v.mean()))[1:] ** 2, np.ones(20) / 20, "valid")
    assert np.corrcoef(np.log(smooth(s)), np.log(smooth(y)))[0, 1] > 0.95


def test_iaaft_keeps_memory():
    y = _fgn(0.8, 2000, 3)
    hs = [hurst_lw(iaaft(y, np.random.default_rng(k))) for k in range(5)]
    assert abs(np.mean(hs) - hurst_lw(y)) < 0.08


def test_iaaft_is_reproducible():
    y = _fgn(0.7, 500, 4)
    assert np.array_equal(iaaft(y, np.random.default_rng(9)), iaaft(y, np.random.default_rng(9)))


def test_drift_surrogate_keeps_trend():
    y = _fgn(0.5, 2000, 5) * np.exp(np.linspace(0, 1, 2000))  # variance rises with position
    s = iaaft_drift(y, np.random.default_rng(6))
    assert abs(variance_trend(s)[0] - variance_trend(y)[0]) < 0.5
    assert abs(variance_trend(iaaft(y, np.random.default_rng(6)))[0]) < 0.5 * abs(variance_trend(y)[0])


def test_variance_trend_matches_diurnal_diagnostic():
    y = _fgn(0.6, 2000, 7)
    v = y.reshape(8, -1).var(axis=1)
    assert np.isclose(variance_trend(y)[0], np.polyfit(np.linspace(0, 1, 8), np.log(v), 1)[0])
