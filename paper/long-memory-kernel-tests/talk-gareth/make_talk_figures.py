"""Charts for the talk. Result charts use numbers from paper.tex and saved experiment summaries.

Panels marked illustrative are single simulated draws. They explain an idea and are not results.
Run: python make_talk_figures.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE / "figures"
INK, MUTED, GRID = "#0b0b0b", "#6b6a66", "#e4e3df"
SHIFT, WILD, SUB, DCOV, PEARSON = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 14, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "legend.frameon": False, "savefig.dpi": 220, "figure.facecolor": "white",
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)


def fgn(n, H, rng):
    """Davies-Harte draw of unit-variance fGn."""
    k = np.arange(n + 1)
    g = 0.5 * (np.abs(k + 1.0) ** (2 * H) - 2 * np.abs(k * 1.0) ** (2 * H) + np.abs(k - 1.0) ** (2 * H))
    lam = np.clip(np.fft.fft(np.r_[g, g[-2:0:-1]]).real, 0, None)
    m = len(lam)
    z = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    return np.fft.fft(np.sqrt(lam / m) * z).real[:n]


def false_alarm_dots():
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))
    for ax, hits, title in ((axes[0], 5, "What a 5% test promises:\n5 alarms per 100 true nulls"),
                            (axes[1], 73, "HSIC with wild multipliers\nunder a seasonal scale: 73"),):
        xy = [(i, 9 - j) for j in range(10) for i in range(10)]
        for idx, (i, j) in enumerate(xy):
            hit = idx < hits
            ax.scatter(i, j, s=120, color=WILD if hit else GRID, marker="X" if hit else "o")
        ax.set_title(title, fontsize=14, color=INK)
        ax.set_aspect("equal"), ax.axis("off")
    save(fig, "false_alarm")


def memory_paths():
    rng = np.random.default_rng(20261101)
    fig, axes = plt.subplots(2, 1, figsize=(9, 4.4), sharex=True, sharey=True)
    for ax, H, label in ((axes[0], 0.5, "No memory (H = 0.5)"), (axes[1], 0.8, "Long memory (H = 0.8)")):
        x = fgn(1000, H, rng)
        ax.plot(x, color=GRID, lw=0.8)
        ax.plot(np.convolve(x, np.ones(50) / 50, mode="same"), color=SHIFT, lw=2.2)
        ax.axhline(0, color=MUTED, lw=0.7)
        ax.text(0.01, 0.98, label, transform=ax.transAxes, va="top", color=INK)
        ax.set_yticks([])
    axes[1].set_xlabel("Time (illustrative draws; blue line is a 50-point average)")
    save(fig, "memory_paths")


def shift_diagram():
    rng = np.random.default_rng(7)
    n, s = 400, 130
    x, y = fgn(n, 0.7, rng), fgn(n, 0.7, rng)
    fig, axes = plt.subplots(3, 1, figsize=(9, 4.0), sharex=True)
    axes[0].plot(x, color=MUTED, lw=1.2)
    axes[1].plot(y[:n - s], color=SHIFT, lw=1.2)
    axes[1].plot(range(n - s, n), y[n - s:], color=WILD, lw=1.2)
    axes[2].plot(range(0, s), y[n - s:], color=WILD, lw=1.2)
    axes[2].plot(range(s, n), y[:n - s], color=SHIFT, lw=1.2)
    axes[2].axvline(s, color=INK, lw=1, ls=":")
    for ax, label in zip(axes, ("Series X stays put", "Series Y", "Series Y after a wrapped shift")):
        ax.text(0.01, 1.0, label, transform=ax.transAxes, va="top", color=INK, fontsize=13)
        ax.set_yticks([]), ax.spines["left"].set_visible(False), ax.margins(y=0.3)
    axes[2].set_xticks([]), axes[2].set_xlabel("Time")
    axes[0].spines["bottom"].set_visible(False), axes[1].spines["bottom"].set_visible(False)
    axes[0].tick_params(bottom=False), axes[1].tick_params(bottom=False)
    save(fig, "shift_diagram")


def yardsticks():
    fig, ax = plt.subplots(figsize=(9, 3.4))
    rows = [(2, "Subsampling", SUB), (1, "Wild multipliers", WILD), (0, "Circular shift", SHIFT)]
    n = 2000
    w = int(np.sqrt(n))
    for y, _, c in rows:
        ax.add_patch(Rectangle((0, y - 0.18), n, 0.36, color=GRID))
    ax.add_patch(Rectangle((800, 2 - 0.18), w, 0.36, color=SUB))
    ax.annotate("one window of about √n points", (800 + w, 2.18), (900, 2.42), color=INK, fontsize=13,
                arrowprops={"arrowstyle": "-", "color": MUTED})
    t = np.linspace(0, 3 * w, 200)
    ax.fill_between(800 + t, 1 - 0.18, 1 - 0.18 + 0.36 * np.exp(-t / w), color=WILD)
    ax.text(800, 1.25, "weights on nearby points fade on a √n scale", color=INK, fontsize=13)
    ax.add_patch(Rectangle((0, -0.18), n, 0.36, color=SHIFT))
    ax.text(0, 0.25, "every shifted copy uses the whole series", color=INK, fontsize=13)
    ax.set_yticks([r[0] for r in rows], [r[1] for r in rows])
    ax.tick_params(left=False, labelcolor=INK)
    ax.set_xlim(0, n), ax.set_ylim(-0.4, 2.7), ax.set_xticks([0, n], ["0", "n = 2000"])
    ax.spines["left"].set_visible(False)
    save(fig, "yardsticks")


def rule_lines():
    x = [0, 1, 2]
    series = [("dCov subsampling", DCOV, [5.3, 10.6, 23.0], 0.0), ("HSIC subsampling", SUB, [5.6, 9.2, 19.7], 0.0),
              ("Pearson NW", PEARSON, [5.8, 8.2, 13.3], 0.0), ("HSIC wild", WILD, [6.5, 7.4, 8.9], 0.5),
              ("HSIC shift", SHIFT, [5.0, 5.5, 5.7], 0.3)]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    ax.axhline(5, color=INK, lw=1, ls="--")
    ax.text(0.02, 3.4, "5% target", color=INK, fontsize=12)
    for label, c, y, dy in series:
        ax.plot(x, y, color=c, lw=2.5, marker="o", ms=8, label=label)
        ax.text(2.05, y[-1] + dy, f"{label}  {y[-1]:.1f}", va="center", color=INK, fontsize=13)
    ax.set_xticks(x, ["No memory\n(H = 0.5 each)", "Moderate memory\n(H = 0.7 each)", "Strong memory\n(H = 0.8 each)"])
    ax.tick_params(axis="x", labelcolor=INK)
    ax.set_xlim(-0.1, 3.0), ax.set_ylim(0, 25), ax.set_ylabel("False alarms (%)")
    ax.yaxis.grid(True, color=GRID), ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=12, ncol=2)
    save(fig, "rule_lines")


def mmd_ranks():
    rows = list(csv.DictReader((REPO / "experiments/EXP1-20261016/e2_summary.csv").open())) \
        if (REPO / "experiments/EXP1-20261016/e2_summary.csv").exists() else None
    if rows is None:
        raise SystemExit("e2 summary not found")
    fig, ax = plt.subplots(figsize=(9, 4.4))
    for proc, label, c in (("G1", "Plain series (rank one)", WILD), ("G5", "Squared series (rank two)", SHIFT)):
        sub = sorted((r for r in rows if r["process"] == proc and r["calibration"] == "wild_unc_05"
                      and r["n"] == "2000" and r["power"] == "False"), key=lambda r: float(r["H"]))
        h, y = [float(r["H"]) for r in sub], [100 * float(r["size"]) for r in sub]
        ax.plot(h, y, color=c, lw=2.5, marker="o", ms=8, label=label)
        ax.text(h[-1] + 0.005, y[-1], f"{label}\n{y[-1]:.1f}", va="center", color=INK, fontsize=13)
    ax.axhline(5, color=INK, lw=1, ls="--")
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9]), ax.set_xlim(right=max(h) + 0.14), ax.set_xlabel("Memory (Hurst index H)"), ax.set_ylabel("False alarms (%)")
    ax.tick_params(axis="x", labelcolor=INK)
    ax.yaxis.grid(True, color=GRID), ax.set_axisbelow(True), ax.legend(loc="upper left", fontsize=12)
    save(fig, "mmd_ranks")


def feature_heatmap():
    gens = ["Gaussian", "Heavy tails, instantaneous", "Heavy tails, filtered", "Stochastic volatility",
            "Threshold AR", "Seasonal mean", "Seasonal scale"]
    cols = ["HSIC\nwild", "HSIC\nshift", "HSIC\nsub", "dCov\nsub", "Pearson\nNW"]
    h5 = [[6.8, 5.3, 5.8, 6.8, 5.0], [6.1, 4.5, 4.3, 3.2, 5.7], [8.9, 6.4, 6.8, 3.2, 4.2], [7.9, 5.8, 4.9, 2.2, 5.2],
          [7.1, 4.7, 7.2, 6.4, 6.6], [91.2, 5.8, 90.9, 91.3, 91.4], [73.3, 5.0, 71.5, 45.6, 4.7]]
    h8 = [[8.1, 5.3, 17.5, 22.8, 14.1], [7.0, 3.8, 15.7, 16.5, 11.5], [10.6, 5.3, 25.9, 25.9, 16.2], [6.7, 3.9, 6.2, 1.9, 4.6],
          [11.9, 7.0, 29.8, 33.5, 12.1], [85.2, 5.9, 85.4, 86.6, 84.5], [68.0, 5.1, 73.0, 53.3, 13.4]]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    for ax, data, title in ((axes[0], h5, "No memory (H = 0.5)"), (axes[1], h8, "Strong memory (H = 0.8)")):
        a = np.array(data)
        ax.imshow(np.sqrt(a), cmap="Blues", vmin=np.sqrt(3), vmax=np.sqrt(120), aspect="auto")
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                ax.text(j, i, f"{a[i, j]:.0f}" if a[i, j] >= 10 else f"{a[i, j]:.1f}", ha="center", va="center",
                        color="white" if a[i, j] > 40 else INK, fontsize=12.5)
        ax.set_xticks(range(5), cols, fontsize=12), ax.set_title(title, color=INK, fontsize=14)
        ax.tick_params(length=0, labelcolor=INK), ax.xaxis.set_ticks_position("top"), ax.xaxis.set_label_position("top")
        ax.set_title(title, color=INK, fontsize=14, pad=44)
        for s in ax.spines.values():
            s.set_visible(False)
    axes[0].set_yticks(range(7), gens)
    save(fig, "feature_heatmap")


def seasonal_scale():
    rng = np.random.default_rng(11)
    n, P = 400, 50
    t = np.arange(n)
    sx = 1 + 0.5 * np.sin(2 * np.pi * (t + 3) / P)
    sy = 1 + 0.5 * np.sin(2 * np.pi * t / P)
    x, y = sx * rng.standard_normal(n), sy * rng.standard_normal(n)
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 4.4), sharex=True)
    axes[0].plot(x, color=SHIFT, lw=0.9, label="X"), axes[0].plot(y, color=WILD, lw=0.9, label="Y", alpha=0.85)
    axes[0].text(0.01, 1.0, "The values: signs look unrelated", transform=axes[0].transAxes, va="top", color=INK, fontsize=13)
    axes[0].margins(y=0.3), axes[0].set_yticks([]), axes[0].legend(loc="lower right", ncol=2, fontsize=12)
    k = np.ones(10) / 10
    axes[1].plot(np.convolve(np.abs(x), k, "same"), color=SHIFT, lw=2)
    axes[1].plot(np.convolve(np.abs(y), k, "same"), color=WILD, lw=2)
    axes[1].text(0.01, 1.0, "Their sizes: swing together", transform=axes[1].transAxes, va="top", color=INK, fontsize=13)
    axes[1].margins(y=0.3), axes[1].set_yticks([]), axes[1].set_xlabel("Time (one illustrative draw, small phase gap)")
    save(fig, "seasonal_scale")


def hbars(ax, labels, values, colors):
    y = range(len(values))[::-1]
    ax.barh(list(y), values, color=colors, height=0.68)
    for yi, v in zip(y, values):
        ax.text(v + 1.2, yi, f"{v:.1f}", va="center", color=INK)
    ax.set_yticks(list(y), labels), ax.tick_params(axis="y", length=0, labelcolor=INK)
    ax.axvline(5, color=INK, lw=1, ls="--")
    ax.spines["left"].set_visible(False)


def seasonal_bars():
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    hbars(ax, ["HSIC wild", "HSIC sub", "dCov sub", "Pearson NW", "HSIC shift"], [73.3, 71.5, 45.6, 4.7, 5.0],
          [WILD, SUB, DCOV, PEARSON, SHIFT])
    ax.set_xlim(0, 90), ax.set_xlabel("False alarms (%). Dashed line: 5%")
    save(fig, "seasonal_bars")


def phase_gap():
    d = np.arange(50)
    c = (0.5 * np.cos(2 * np.pi * d / 50)) / 1.5
    small = np.isin(d, [12, 13, 37, 38])
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.axhline(0, color=MUTED, lw=1)
    ax.vlines(d, 0, c, color=GRID, lw=1.5)
    ax.scatter(d[~small], c[~small], color=SHIFT, s=55, label="46 gaps: limit clearly away from zero")
    ax.scatter(d[small], c[small], color=WILD, marker="D", s=70, label="4 gaps: limit close to zero (|c| = 0.021)")
    ax.set_xlabel("Phase gap between the two series (random, fixed for the whole path)")
    ax.set_ylabel("Where the sample\ncorrelation ends up")
    ax.tick_params(labelcolor=INK)
    ax.legend(loc="upper center", fontsize=12, ncol=1)
    ax.set_ylim(-0.38, 0.62)
    save(fig, "phase_gap")


def adjustment():
    labels = ["No adjustment", "Remove mean only", "Remove mean and scale"]
    wild, shift = [73.3, 73.8, 7.2], [5.0, 4.3, 5.0]
    y = np.arange(3)[::-1]
    fig, ax = plt.subplots(figsize=(9, 3.9))
    ax.barh(y + 0.19, wild, height=0.34, color=WILD, label="HSIC wild")
    ax.barh(y - 0.19, shift, height=0.34, color=SHIFT, label="HSIC shift")
    for yi, a, b in zip(y, wild, shift):
        ax.text(a + 1, yi + 0.19, f"{a:.1f}", va="center", color=INK)
        ax.text(b + 1, yi - 0.19, f"{b:.1f}", va="center", color=INK)
    ax.axvline(5, color=INK, lw=1, ls="--")
    ax.set_yticks(y, labels), ax.tick_params(axis="y", length=0, labelcolor=INK), ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 88), ax.set_xlabel("False alarms (%). Dashed line: 5%")
    ax.legend(loc="lower right", fontsize=13)
    save(fig, "adjustment")


def power_panels():
    data = {"Linear link": ([56.0, 28.7], [57.2, 30.6], [54.4, 27.1]), "Squared link": ([43.0, 19.1], [44.0, 21.2], [40.9, 17.1]),
            "Seasonal scale,\nlinear link": ([13.8, 12.7], [15.2, 17.4], [14.6, 15.6])}
    fig, axes = plt.subplots(1, 3, figsize=(10, 4.2), sharey=True)
    for ax, (title, (w, s, u)) in zip(axes, data.items()):
        x = np.arange(2)
        for off, vals, c, lab in ((-0.27, s, SHIFT, "HSIC shift"), (0, w, WILD, "HSIC wild"), (0.27, u, SUB, "HSIC sub")):
            ax.bar(x + off, vals, width=0.25, color=c, label=lab)
        for xi, v in zip(x - 0.27, s):
            ax.text(xi, v + 1.2, f"{v:.0f}", ha="center", color=INK, fontsize=13)
        ax.set_xticks(x, ["H = 0.5", "H = 0.8"]), ax.tick_params(axis="x", length=0, labelcolor=INK)
        ax.set_title(title, color=INK, fontsize=14)
    axes[0].set_ylabel("Detection at a true 5% (%)")
    axes[2].legend(loc="upper right", fontsize=12)
    save(fig, "power")


def half_mean():
    m = 1000
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), sharex=True)
    for ax, H in zip(axes, (0.5, 0.8)):
        sd = np.sqrt((4 - 2 ** (2 * H)) * m ** (2 * H - 2))
        x = np.linspace(-0.75, 0.75, 800)
        y = norm.pdf(x, 0, sd)
        ax.fill_between(x, y, color=GRID)
        ax.plot(x, y, color=MUTED, lw=1.5)
        ax.axvline(0.1, color=WILD, lw=2.5)
        ax.set_yticks([]), ax.spines["left"].set_visible(False)
        ax.set_title(f"H = {H}: the shift is {0.1 / sd:.2f} null SDs", color=INK, fontsize=14)
        ax.set_xlabel("Difference of the two half-means")
    save(fig, "half_mean")


def limits():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.0))
    ax = axes[0]
    n = ["500", "1000", "2000", "4000", "8000"]
    y, lo, hi = [5.9, 5.9, 5.7, 6.3, 8.3], [5.0, 4.9, 4.8, 5.0, 4.8], [7.1, 7.0, 6.8, 8.0, 14.0]
    ax.errorbar(range(5), y, yerr=[np.subtract(y, lo), np.subtract(hi, y)], color=SHIFT, marker="o", ms=8, lw=2.2, capsize=4)
    ax.axhline(5, color=INK, lw=1, ls="--")
    ax.set_xticks(range(5), n), ax.set_ylim(0, 16), ax.tick_params(axis="x", labelcolor=INK)
    ax.set_xlabel("Series length n"), ax.set_ylabel("False alarms (%)")
    ax.set_title("Strong memory, HSIC shift\n(n = 8000 has 144 replications)", color=INK, fontsize=13)
    ax = axes[1]
    x = np.arange(2)
    a, b = [6.2, 10.1], [8.8, 16.8]
    ax.bar(x - 0.2, a, width=0.38, color=SHIFT, label="HSIC shift"), ax.bar(x + 0.2, b, width=0.38, color=WILD, label="HSIC wild")
    for xi, v in list(zip(x - 0.2, a)) + list(zip(x + 0.2, b)):
        ax.text(xi, v + 0.3, f"{v:.1f}", ha="center", color=INK)
    ax.axhline(5, color=INK, lw=1, ls="--")
    ax.set_xticks(x, ["n = 500", "n = 2000"]), ax.tick_params(axis="x", length=0, labelcolor=INK), ax.set_ylim(0, 20)
    ax.set_title("Slow scale drift, no memory\n(window sees a quarter cycle)", color=INK, fontsize=13)
    ax.legend(loc="upper left", fontsize=12)
    save(fig, "limits")


def stock():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    ax = axes[0]
    pairs = ["AAPL-MSFT", "AAPL-NVDA", "MSFT-NVDA", "AAPL-AAPL"]
    d0, d2 = [8.4, 9.2, 10.2, 7.2], [8.8, 9.0, 9.8, 7.4]
    x = np.arange(4)
    ax.bar(x - 0.2, d0, width=0.38, color=MUTED, label="Original adjustment")
    ax.bar(x + 0.2, d2, width=0.38, color=SUB, label="Time-of-day mean and scale removed")
    ax.axhline(5, color=INK, lw=1, ls="--")
    ax.set_xticks(x, pairs, fontsize=11), ax.tick_params(axis="x", length=0, labelcolor=INK)
    ax.set_ylabel("Rejected pairs (%)"), ax.set_ylim(0, 15), ax.legend(loc="upper left", fontsize=11)
    ax.set_title("HSIC subsampling, n = 2000", color=INK, fontsize=13)
    ax = axes[1]
    x = np.arange(3)
    sub, sh = [7.5, 8.8, 11.0], [3.6, 5.4, 4.3]
    ax.bar(x - 0.2, sub, width=0.38, color=SUB, label="HSIC sub"), ax.bar(x + 0.2, sh, width=0.38, color=SHIFT, label="HSIC shift")
    for xi, v in list(zip(x - 0.2, sub)) + list(zip(x + 0.2, sh)):
        ax.text(xi, v + 0.25, f"{v:.1f}", ha="center", color=INK)
    ax.axhline(5, color=INK, lw=1, ls="--")
    ax.set_xticks(x, ["0 windows\n(689 pairs)", "1 window\n(966 pairs)", "2 windows\n(345 pairs)"], fontsize=11)
    ax.tick_params(axis="x", length=0, labelcolor=INK), ax.legend(loc="upper left", fontsize=11)
    ax.set_title("By number of windows with a variance trend", color=INK, fontsize=13)
    save(fig, "stock")


def bandwidth():
    k = ["0.5", "1", "2", "4"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
    ax = axes[0]
    ax.plot(k, [10.2, 15.6, 42.6, 72.8], color=WILD, lw=2.5, marker="o", ms=8, label="Seasonal scale, linear link")
    ax.plot(k, [22.0, 37.2, 41.0, 7.2], color=SUB, lw=2.5, marker="s", ms=8, label="Squared link")
    ax.set_ylim(0, 95), ax.set_xlabel("Bandwidth, as a multiple of the median"), ax.set_ylabel("HSIC wild detection at a true 5% (%)", fontsize=12)
    ax.legend(loc="upper left", fontsize=12), ax.tick_params(axis="x", labelcolor=INK)
    ax.yaxis.grid(True, color=GRID), ax.set_axisbelow(True)
    ax = axes[1]
    ax.plot(k, [5.0, 5.6, 5.8, 4.8], color=SHIFT, lw=2.5, marker="o", ms=8, label="HSIC shift")
    ax.plot(k, [72.0, 74.4, 56.8, 6.6], color=WILD, lw=2.5, marker="s", ms=8, label="HSIC wild")
    ax.axhline(5, color=INK, lw=1, ls="--"), ax.set_ylim(0, 95)
    ax.set_xlabel("Bandwidth, as a multiple of the median"), ax.set_ylabel("False alarms, seasonal scale (%)", fontsize=12)
    ax.legend(loc="upper right", fontsize=12), ax.tick_params(axis="x", labelcolor=INK)
    save(fig, "bandwidth")


if __name__ == "__main__":
    for f in (false_alarm_dots, memory_paths, shift_diagram, yardsticks, rule_lines, mmd_ranks, feature_heatmap, seasonal_scale,
              seasonal_bars, phase_gap, adjustment, power_panels, half_mean, limits, stock, bandwidth):
        f()
