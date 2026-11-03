"""Draw the six report figures from data/report_numbers.json. No other input.

Palette: the dataviz reference categorical slots, validated with validate_palette.js on a
white surface (adjacent pairs for line charts, all pairs for the 4-series dot plot).
raw is always blue and inc1 always orange. Other series take the next slots in fixed order
within a figure. Each series also has its own marker, because three slots are below 3:1
contrast on white (the validator's relief rule).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
D = json.loads((HERE / "data" / "report_numbers.json").read_text())
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.labelsize": 8.5, "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.fontsize": 8, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.dpi": 200, "lines.linewidth": 1.6, "lines.markersize": 5.5, "text.color": INK,
})

# Entity style: colour, marker, label.
STY = {
    "raw": (BLUE, "o", "raw"),
    "inc1": (ORANGE, "s", "inc1"),
    "herm1": (AQUA, "^", "herm1"),
    "datadiff": (YELLOW, "D", "datadiff"),
    "dcov_inc1": (YELLOW, "D", "dcov_inc1"),
    "hsplit1": (MAGENTA, "v", "hsplit1"),
    "harris": (MAGENTA, "P", "raw, Harris rule"),
    "hsplit2": (YELLOW, "X", "hsplit2"),
}
RING = dict(markeredgecolor="white", markeredgewidth=0.9)


def save(fig, name):
    fig.savefig(FIG / name, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", FIG / name)


def five_line(ax, ymax=None):
    ax.axhline(5, color=INK2, lw=0.9, zorder=1)
    ax.text(1.0, 5, " 5%", transform=ax.get_yaxis_transform(), color=INK2, va="center", fontsize=7.5)


# ---------------------------------------------------------------- 1. one pair of paths

def fig1():
    e = D["ex1"]
    x, dx = np.array(e["x"]), np.array(e["dx"])
    L = np.arange(e["maxlag"] + 1)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    a = axs[0]
    t = np.arange(x.size)
    a.plot(t, x + 4.5, color=BLUE, lw=0.6)
    a.plot(t[1:], dx - 4.5, color=ORANGE, lw=0.6)
    a.set_yticks([])
    a.text(10, 9.3, "raw series $X_t$", color=INK, fontsize=8)
    a.text(10, -0.4, "increments $X_{t+1}-X_t$", color=INK, fontsize=8)
    a.set_xlabel("time t")
    a.set_title("(a) one path, H = 0.9, n = 2000", loc="left")
    a.grid(False)

    b = axs[1]
    k = L[1:201]
    b.plot(k, np.array(e["acf_x"])[1:201], color=BLUE, lw=0.7, alpha=0.55)
    b.plot(k, np.array(e["acf_dx"])[1:201], color=ORANGE, lw=0.7, alpha=0.55)
    b.plot(k, np.array(e["acf_theory_raw"])[1:201], color=BLUE, lw=1.8, label="raw")
    b.plot(k, np.array(e["acf_theory_inc"])[1:201], color=ORANGE, lw=1.8, label="increments")
    b.axhline(0, color=AXIS, lw=0.8)
    b.set_xscale("log")
    b.set_xlabel("lag j")
    b.set_ylabel("autocorrelation")
    b.set_title("(b) autocorrelation", loc="left")
    b.legend(loc="upper right", handlelength=1.5)

    c = axs[2]
    k = L[1:]
    c.plot(k, np.array(e["pairsum_theory_raw"]), color=BLUE, lw=1.8, label="raw")
    c.plot(k, np.array(e["pairsum_theory_inc"]), color=ORANGE, lw=1.8, label="increments")
    c.plot(k, np.array(e["pairsum_sample_raw"]), color=BLUE, lw=0.7, alpha=0.55)
    c.plot(k, np.array(e["pairsum_sample_inc"]), color=ORANGE, lw=0.7, alpha=0.55)
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xlabel("lags summed, |j| up to L")
    c.set_ylabel(r"$\sum_{|j|\leq L} r_X(j)\, r_Y(j)$")
    c.set_title("(c) pair sum", loc="left")
    c.text(1.5, e["pairsum_theory_raw"][-1] * 0.5, "keeps growing", color=INK2, fontsize=7.5)
    c.text(30, 1.6, "levels off", color=INK2, fontsize=7.5)
    fig.text(0.01, -0.04, "Thick lines: exact fGn theory. Thin lines: this simulated pair (seed %d)." % e["seed"],
             color=INK2, fontsize=7.5)
    fig.tight_layout()
    save(fig, "fig1_paths.png")


# ---------------------------------------------------------------- 2. null growth

def fig2():
    g = D["growth"]
    ns = [500, 1000, 2000, 4000, 8000]
    Hs = ["0.5", "0.7", "0.8", "0.9", "0.95"]
    fig, axs = plt.subplots(1, 5, figsize=(7.2, 2.4), sharey=True)
    for ax, H in zip(axs, Hs):
        for k in ["raw", "inc1"]:
            c, mk, lab = STY[k]
            q = g[f"{k}|{H}"]["q95"]
            base = q["500"]["v"]
            v = np.array([q[str(n)]["v"] for n in ns]) / base
            lo = np.array([q[str(n)]["lo"] for n in ns]) / base
            hi = np.array([q[str(n)]["hi"] for n in ns]) / base
            ax.fill_between(ns, lo, hi, color=c, alpha=0.12, lw=0)
            ax.plot(ns, v, color=c, marker=mk, label=lab, **RING)
            dy = 5 if k == "raw" else -5
            ax.annotate(f"x{v[-1]:.2f}", (ns[-1], v[-1]), xytext=(3, dy), textcoords="offset points",
                        fontsize=7, color=INK2, va="center")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(ns)
        ax.set_xticklabels(["500", "", "2k", "", "8k"])
        ax.minorticks_off()
        ax.set_title(f"H = {H}", loc="left")
        ax.set_xlim(400, 16000)
        ax.set_xlabel("n")
    axs[0].set_yticks([0.5, 0.7, 1, 1.5, 2, 3, 4.5])
    axs[0].set_yticklabels(["0.5", "0.7", "1", "1.5", "2", "3", "4.5"])
    axs[0].set_ylabel("95th pct of null T\n(relative to n = 500)")
    axs[0].legend(loc="upper left")
    fig.tight_layout()
    save(fig, "fig2_growth.png")


# ---------------------------------------------------------------- 3. false alarms against H

def fig3():
    s = D["size"]
    Hf = [float(h) for h in ["0.5", "0.6", "0.7", "0.8", "0.9", "0.95"]]
    panels = [("wild", "wild bootstrap"), ("shift", "circular shift"), ("sub", "subsampling"), ("harris", "Harris shift")]
    kern = {"wild": ["raw", "inc1", "herm1", "datadiff", "hsplit1"], "shift": ["raw", "inc1", "herm1", "datadiff", "hsplit1"],
            "sub": ["raw", "inc1"], "harris": ["raw", "inc1"]}
    fig, axs = plt.subplots(1, 4, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1.5, 1.5, 1, 1]})
    for ax, (rule, title) in zip(axs, panels):
        ks = kern[rule]
        off = np.linspace(-0.3, 0.3, len(ks)) if len(ks) > 2 else np.array([-0.12, 0.12])
        for o, k in zip(off, ks):
            c, mk, lab = STY[k]
            r = [s[f"{k}|{rule}|{h}|2000"] for h in ["0.5", "0.6", "0.7", "0.8", "0.9", "0.95"]]
            v = np.array([x["v"] for x in r]) * 100
            lo = np.array([x["lo"] for x in r]) * 100
            hi = np.array([x["hi"] for x in r]) * 100
            xs = np.arange(len(Hf)) + o
            ax.errorbar(xs, v, yerr=[v - lo, hi - v], color=c, marker=mk, ms=4.5, lw=1.1, elinewidth=0.8,
                        capsize=0, label=lab, **RING)
        five_line(ax)
        ax.set_title(title, loc="left")
        ax.set_xticks(range(6))
        ax.set_xticklabels([".5", ".6", ".7", ".8", ".9", ".95"] if len(ks) > 2 else [".5", "", ".7", "", ".9", ""])
        ax.grid(axis="x", visible=False)
        ax.set_xlabel("H (both series)")
        ax.set_ylim(0, None)
    axs[0].set_ylabel("false-alarm rate, % (nominal 5%)")
    axs[0].set_ylim(0, 22)
    axs[1].set_ylim(0, 22)
    axs[2].set_ylim(0, 70)
    axs[3].set_ylim(0, 22)
    for i in (0, 1, 3):
        axs[i].set_yticks([0, 5, 10, 15, 20])
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "fig3_size.png")


# ---------------------------------------------------------------- 4. power at matched size

def fig4():
    p = D["power"]
    Hs = ["0.5", "0.7", "0.8", "0.9"]
    Hf = np.array([float(h) for h in Hs])
    series = [("raw", "shift", "raw"), ("inc1", "shift", "inc1"), ("herm1", "shift", "herm1"),
              ("dcov_inc1", "shift", "dcov_inc1"), ("raw", "harris", "harris")]
    names = {"lin": "linear (lin)", "sq": "squared (sq)", "lag": "one-lag (lag)",
             "flat": "flat moving average", "factor": "slow common factor", "scale": "scale (volatility)"}
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6), sharex=True, sharey=True)
    for ax, a in zip(axs.flat, ["lin", "sq", "lag", "flat", "factor", "scale"]):
        off = np.linspace(-0.012, 0.012, len(series))
        for o, (k, r, sk) in zip(off, series):
            c, mk, lab = STY[sk]
            rr = [p[f"{a}|{k}|{r}|{h}"] for h in Hs]
            v = np.array([x["v"] for x in rr]) * 100
            lo = np.array([x["lo"] for x in rr]) * 100
            hi = np.array([x["hi"] for x in rr]) * 100
            ax.errorbar(Hf + o, v, yerr=[v - lo, hi - v], color=c, marker=mk, ms=4.5, lw=1.3, elinewidth=0.7,
                        capsize=0, label=lab, **RING)
        ax.axhline(5, color=INK2, lw=0.9, zorder=1)
        ax.set_title(names[a], loc="left")
        ax.set_ylim(0, 102)
        ax.set_xticks(Hf)
    for ax in axs[1]:
        ax.set_xlabel("H")
    for ax in axs[:, 0]:
        ax.set_ylabel("power at matched 5% size, %")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "fig4_power.png")


# ---------------------------------------------------------------- 5. seasonal and drift cells

def dot_panel(ax, cells, series, getter, title):
    """Horizontal dot plot on a log x axis. cells: list of (label, key-part). series: list of style keys."""
    ys = np.arange(len(cells))[::-1]
    off = np.linspace(0.16, -0.16, len(series)) if len(series) > 1 else [0]
    for o, sk in zip(off, series):
        c, mk, lab = STY[sk]
        for y, (_, cell) in zip(ys, cells):
            r = getter(sk, cell)
            if r is None:
                continue
            v, lo, hi = (max(r[q] * 100, 0.5) for q in ("v", "lo", "hi"))
            ax.plot([lo, hi], [y + o, y + o], color=c, lw=0.9)
            ax.plot([v], [y + o], color=c, marker=mk, ls="none", ms=5.5, label=lab, **RING)
    ax.axvline(5, color=INK2, lw=0.9, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.8, 130)
    ax.set_xticks([1, 5, 20, 100])
    ax.set_xticklabels(["1", "5", "20", "100"])
    ax.minorticks_off()
    ax.set_yticks(ys)
    ax.set_yticklabels([c[0] for c in cells])
    ax.set_title(title, loc="left")
    ax.grid(axis="y", visible=False)


def fig5():
    A = D["arms"]
    cells = [("shared season, H 0.5", ("SEAS-M-shared", "0.5")), ("shared season, H 0.8", ("SEAS-M-shared", "0.8")),
             ("slow drift, H 0.5", ("SEAS-V-slow", "0.5")), ("slow drift, H 0.8", ("SEAS-V-slow", "0.8"))]
    fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.5), sharey=True)
    for ax, (rule, title) in zip(axs, [("shift", "circular shift"), ("wild", "wild bootstrap"),
                                       ("sub", "subsampling"), ("harris", "Harris shift")]):
        dot_panel(ax, cells, ["raw", "inc1"],
                  lambda sk, cell, rule=rule: A.get(f"{cell[0]}|{sk}|{rule}|{cell[1]}"), title)
    fig.supxlabel("false-alarm rate, % (log scale, 5% line marked)", fontsize=8.5, color=INK2, y=0.13)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h[::4][:2], l[::4][:2], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    save(fig, "fig5_seasonal.png")


# ---------------------------------------------------------------- 6. Hermite split

def fig6():
    P, M, A = D["power"], D["mmd_power"], D["arms"]
    series = ["raw", "herm1", "hsplit1", "hsplit2"]
    pcells = [("lin, H 0.5", ("lin", "0.5")), ("lin, H 0.9", ("lin", "0.9")),
              ("sq, H 0.5", ("sq", "0.5")), ("factor, H 0.5", ("factor", "0.5")),
              ("factor, H 0.9", ("factor", "0.9")), ("split MMD, mean shift\n0.5, H 0.5 (wild)", ("mmd", "0.5")),
              ("split MMD, mean shift\n0.5, H 0.8 (wild)", ("mmd", "0.8"))]

    def pget(sk, cell):
        a, H = cell
        if a == "mmd":
            return M[f"mean0.5|{sk}|{H}"]
        return P[f"{a}|{sk}|shift|{H}"]

    scells = [("Gaussian, H 0.8", ("G0", "0.8")), ("Gaussian, H 0.95", ("G0", "0.95")),
              ("SV, H 0.8", ("NL-SV", "0.8")), ("season, H 0.5", ("SEAS-M", "0.5")),
              ("seasonal scale, H 0.5", ("SEAS-V", "0.5")), ("shared season, H 0.5", ("SEAS-M-shared", "0.5")),
              ("slow drift, H 0.5", ("SEAS-V-slow", "0.5"))]

    def sget(sk, cell):
        g, H = cell
        if g == "G0":
            return D["size"][f"{sk}|shift|{H}|2000"]
        return A.get(f"{g}|{sk}|shift|{H}")

    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
    # Power panel on a linear axis.
    ax = axs[0]
    ys = np.arange(len(pcells))[::-1]
    off = np.linspace(0.24, -0.24, len(series))
    for o, sk in zip(off, series):
        c, mk, lab = STY[sk]
        for y, (_, cell) in zip(ys, pcells):
            r = pget(sk, cell)
            ax.plot([r["lo"] * 100, r["hi"] * 100], [y + o, y + o], color=c, lw=0.9)
            ax.plot([r["v"] * 100], [y + o], color=c, marker=mk, ls="none", ms=5.5, label=lab, **RING)
    ax.axvline(5, color=INK2, lw=0.9)
    ax.set_xlim(0, 102)
    ax.set_yticks(ys)
    ax.set_yticklabels([c[0] for c in pcells])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("power at matched 5% size, %")
    ax.set_title("(a) power", loc="left")
    # Size panel.
    dot_panel(axs[1], scells, series, sget, "(b) false alarms, shift rule")
    axs[1].set_xlabel("false-alarm rate, % (log)")
    h, l = ax.get_legend_handles_labels()
    k = len(pcells)
    fig.legend(h[::k], l[::k], loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "fig6_hsplit.png")


# ---------------------------------------------------------------- 7. split-MMD false alarms

GREEN, RED = "#008300", "#e34948"


def fig7():
    """Split MMD, wild rule. Validated order (adjacent pairs): raw, herm1, datadiff, hsplit1, sym, hsplit2."""
    S = D["mmd_size"]
    Hs = ["0.5", "0.6", "0.7", "0.75", "0.8", "0.85", "0.9"]
    sty = [("raw", BLUE, "o"), ("herm1", AQUA, "^"), ("datadiff", YELLOW, "D"),
           ("hsplit1", MAGENTA, "v"), ("sym", GREEN, "s"), ("hsplit2", RED, "X")]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 4.6), sharex=True, sharey=True)
    for i, (g, gl) in enumerate([("G0", "Gaussian (G0)"), ("HT-S", "heavy tails (HT-S)")]):
        for j, n in enumerate(["500", "2000"]):
            ax = axs[i, j]
            off = np.linspace(-0.25, 0.25, len(sty))
            for o, (k, c, mk) in zip(off, sty):
                r = [S[f"{g}|{k}|{h}|{n}"] for h in Hs]
                v = np.array([x["v"] for x in r]) * 100
                lo = np.maximum(np.array([x["lo"] for x in r]) * 100, 0.6)
                hi = np.array([x["hi"] for x in r]) * 100
                ax.errorbar(np.arange(len(Hs)) + o, v, yerr=[v - lo, hi - v], color=c, marker=mk, ms=4.5, lw=1.2,
                            elinewidth=0.7, capsize=0, label=k, **RING)
            ax.axhline(5, color=INK2, lw=0.9, zorder=1)
            ax.set_yscale("log")
            ax.set_ylim(1.5, 50)
            ax.set_yticks([2, 5, 10, 20, 40])
            ax.set_yticklabels(["2", "5", "10", "20", "40"])
            ax.minorticks_off()
            ax.set_xticks(range(len(Hs)))
            ax.set_xticklabels([".5", ".6", ".7", ".75", ".8", ".85", ".9"])
            ax.grid(axis="x", visible=False)
            ax.set_title(f"{gl}, n = {n}", loc="left")
    for ax in axs[1]:
        ax.set_xlabel("H")
    for ax in axs[:, 0]:
        ax.set_ylabel("false-alarm rate, % (log)")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "fig7_mmd_size.png")


# ---------------------------------------------------------------- 8. Harris shift test

def fig8():
    """Harris on the raw kernel. Palette blue, orange, aqua: validated all pairs on white."""
    HC, P = D["harris_cells"], D["power"]
    groups = ["Gaussian memory", "process", "stress null", "non-stationary"]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    ax = axs[0]
    rng = np.random.default_rng(0)
    for gi, g in enumerate(groups):
        cells = [c for c in HC if c["group"] == g]
        y0 = len(groups) - 1 - gi
        for k, c, mk, o in [("raw", AQUA, "P", 0.12), ("inc1", ORANGE, "s", -0.12)]:
            v = np.array([x[k] for x in cells]) * 100
            ax.plot(v, y0 + o + rng.uniform(-0.06, 0.06, v.size), ls="none", marker=mk, color=c,
                    ms=6 if mk == "P" else 4.5, markeredgewidth=0 if mk == "P" else 0.9, markeredgecolor="white",
                    label=("raw kernel" if k == "raw" else "inc1") if gi == 0 else None)
    ax.axvline(5, color=INK2, lw=0.9)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([f"{g}\n({sum(c['group'] == g for c in HC)} cells)" for g in groups][::-1])
    ax.set_xscale("log")
    ax.set_xlim(0.8, 40)
    ax.set_xticks([1, 2, 5, 10, 20])
    ax.set_xticklabels(["1", "2", "5", "10", "20"])
    ax.minorticks_off()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("false-alarm rate, % (log)")
    ax.set_title("(a) false alarms, Harris rule", loc="left")
    ax.legend(loc="upper right", handletextpad=0.2, borderaxespad=0.2)
    alts = ["lin", "sq", "lag", "flat", "factor", "scale"]
    series = [("raw", "shift", BLUE, "o", "raw kernel, shift"), ("inc1", "shift", ORANGE, "s", "inc1, shift"),
              ("raw", "harris", AQUA, "P", "raw kernel, Harris")]
    for ax, H, t in [(axs[1], "0.5", "(b) power, H = 0.5"), (axs[2], "0.9", "(c) power, H = 0.9")]:
        ys = np.arange(len(alts))[::-1]
        for o, (k, r, c, mk, lab) in zip([0.2, 0, -0.2], series):
            for y, a in zip(ys, alts):
                q = P[f"{a}|{k}|{r}|{H}"]
                ax.plot([q["lo"] * 100, q["hi"] * 100], [y + o] * 2, color=c, lw=0.9)
                ax.plot([q["v"] * 100], [y + o], color=c, marker=mk, ls="none", ms=6 if mk == "P" else 5, label=lab,
                        markeredgecolor="white", markeredgewidth=0 if mk == "P" else 0.9)
        ax.set_xlim(0, 102)
        ax.set_yticks(ys)
        ax.set_yticklabels(alts)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("power at matched 5% size, %")
        ax.set_title(t, loc="left")
    axs[2].set_yticklabels([])
    h, l = axs[1].get_legend_handles_labels()
    fig.legend(h[::6], l[::6], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, "fig8_harris.png")


# ---------------------------------------------------------------- 9. talk slide: raw against inc1

def fig9():
    """Slide figure for the Gareth talk: false alarms against H, raw against inc1, three rules, n = 2000."""
    s = D["size"]
    Hs = ["0.5", "0.6", "0.7", "0.8", "0.9", "0.95"]
    with plt.rc_context({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 11}):
        fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.2))
        for ax, (rule, title, top) in zip(axs, [("wild", "wild bootstrap", 25), ("shift", "circular shift", 25),
                                                 ("sub", "subsampling", 70)]):
            for k, o in [("raw", -0.1), ("inc1", 0.1)]:
                c, mk, lab = STY[k]
                r = [s[f"{k}|{rule}|{h}|2000"] for h in Hs]
                v = np.array([x["v"] for x in r]) * 100
                lo = np.array([x["lo"] for x in r]) * 100
                hi = np.array([x["hi"] for x in r]) * 100
                ax.errorbar(np.arange(len(Hs)) + o, v, yerr=[v - lo, hi - v], color=c, marker=mk, ms=6, lw=2,
                            elinewidth=1, capsize=0, label=lab, **RING)
            ax.axhline(5, color=INK2, lw=1, zorder=1)
            ax.text(5.45, 5, " 5%", color=INK2, va="center", fontsize=10)
            ax.set_xticks(range(len(Hs)))
            ax.set_xticklabels(Hs)
            ax.grid(axis="x", visible=False)
            ax.set_ylim(0, top)
            ax.set_xlabel("H (both series)")
            ax.set_title(title, loc="left")
        axs[0].set_ylabel("false alarms, %")
        axs[0].legend(loc="upper left")
        fig.tight_layout()
        save(fig, "fig9_talk_size.png")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9):
        f()
