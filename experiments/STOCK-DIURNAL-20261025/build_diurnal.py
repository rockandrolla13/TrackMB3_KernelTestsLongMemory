"""Build D1 (Engle-Russell diurnal factor) and D2 (plus diurnal scale) duration series.

Reads the STOCK-20261023 cache and writes adjusted series to a git-ignored cache. The fitted
curves and scale diagnostics, which are derived statistics only, go to this folder.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SRC = ROOT / "data/raw/stock_duration_cache_2024"
OUT = ROOT / "data/raw/stock_diurnal_cache_2024"
SYMBOLS = ("AAPL", "MSFT", "NVDA")
SESSION = 23400
KNOTS_MIN = (15, 30, 60, 120, 210, 300, 360)
DEGREE = 3


def basis(t: np.ndarray) -> np.ndarray:
    inner = np.asarray(KNOTS_MIN, float) * 60
    knots = np.r_[[0.0] * (DEGREE + 1), inner, [float(SESSION)] * (DEGREE + 1)]
    return BSpline.design_matrix(np.clip(t, 0, SESSION - 1e-6), knots, DEGREE).toarray()


GRID = np.arange(SESSION) + 0.5
B_GRID = basis(GRID)


def smooth_fit(sec: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Least-squares spline coefficients from per-second sums (regressor rounded to 1 s)."""
    s = sec.astype(np.int64)
    counts = np.bincount(s, minlength=SESSION).astype(float)
    sums = np.bincount(s, weights=y, minlength=SESSION)
    xtx = B_GRID.T @ (B_GRID * counts[:, None])
    xty = B_GRID.T @ sums
    return np.linalg.solve(xtx, xty)


def evaluate(coef: np.ndarray, sec: np.ndarray) -> np.ndarray:
    out = np.empty(sec.size)
    for k in range(0, sec.size, 2_000_000):
        out[k:k + 2_000_000] = basis(sec[k:k + 2_000_000]) @ coef
    return out


def scale_range(sec: np.ndarray, bins: np.ndarray, v: np.ndarray) -> dict:
    """Smooth and clock-bin standard-deviation ranges across the session."""
    m = evaluate(smooth_fit(sec, v), sec)
    var_curve = B_GRID @ smooth_fit(sec, (v - m) ** 2)
    sd_curve = np.sqrt(np.maximum(var_curve, 1e-12))
    bin_sd = np.array([v[bins == b].std() for b in range(26)])
    return {
        "smooth_sd_max_over_min": float(sd_curve.max() / sd_curve.min()),
        "smooth_sd_min": float(sd_curve.min()), "smooth_sd_max": float(sd_curve.max()),
        "smooth_sd_argmin_minute": float(GRID[sd_curve.argmin()] / 60),
        "smooth_sd_argmax_minute": float(GRID[sd_curve.argmax()] / 60),
        "bin_sd_max_over_min": float(bin_sd.max() / bin_sd.min()),
        "bin_sd": [round(float(x), 4) for x in bin_sd],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"knots_minutes": list(KNOTS_MIN), "degree": DEGREE, "symbols": {}}
    for symbol in SYMBOLS:
        with np.load(SRC / f"{symbol}.npz") as z:
            d0 = z["dur"].astype(float)
            ts = z["event_ts_ns"]
            raw = z["raw_dur_ns"].astype(float) / 1e9
            bins = z["clock_bin"].astype(int)
            day = z["dur_day"]
            days = z["days"]
        local = pd.to_datetime(ts, utc=True).tz_convert("America/New_York")
        sec = ((local.hour * 3600 + local.minute * 60 + local.second).to_numpy() - 34200
               + local.nanosecond.to_numpy() * 1e-9 + local.microsecond.to_numpy() * 1e-6)
        if sec.min() < 0 or sec.max() >= SESSION:
            raise ValueError(f"{symbol}: time outside session")
        if not np.array_equal((sec // 900).astype(int), bins):
            raise ValueError(f"{symbol}: seconds disagree with cached clock bins")

        phi_coef = smooth_fit(sec, raw)
        phi = evaluate(phi_coef, sec)
        if phi.min() <= 0:
            raise ValueError(f"{symbol}: nonpositive diurnal factor")
        r = np.log(raw / phi)
        m_coef = smooth_fit(sec, r)
        m = evaluate(m_coef, sec)
        v_coef = smooth_fit(sec, (r - m) ** 2)
        v = evaluate(v_coef, sec)
        if v.min() <= 0:
            raise ValueError(f"{symbol}: nonpositive diurnal variance")
        d2 = (r - m) / np.sqrt(v)

        np.savez(OUT / f"{symbol}.npz", D1=r.astype(np.float32), D2=d2.astype(np.float32),
                 dur_day=day, days=days)
        phi_grid = B_GRID @ phi_coef
        report["symbols"][symbol] = {
            "phi_seconds_at_minute": {str(mn): round(float(phi_grid[mn * 60]), 5)
                                      for mn in (0, 5, 15, 30, 60, 120, 195, 300, 360, 385, 389)},
            "phi_max_over_min": float(phi_grid.max() / phi_grid.min()),
            "scale_D0": scale_range(sec, bins, d0),
            "scale_D1": scale_range(sec, bins, r),
            "scale_D2": scale_range(sec, bins, d2),
        }
        print(symbol, json.dumps({k: report["symbols"][symbol][k]["smooth_sd_max_over_min"]
                                  for k in ("scale_D0", "scale_D1", "scale_D2")}), flush=True)
    (HERE / "diurnal_fit.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
