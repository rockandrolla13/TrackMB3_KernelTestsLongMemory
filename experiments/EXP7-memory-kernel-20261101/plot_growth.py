"""EXP7 headline picture: 95th percentile of the bare statistic against n, raw and inc1.

Usage: python plot_growth.py <tag>  -> growth_q95_<tag>.png, from results_<tag>.json. Plain python: matplotlib does not import under -I here.
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def write_figure(growth: dict, tag: str) -> None:
    hs = sorted({float(k.split('|')[1]) for k in growth["q95"]})
    ng = sorted({int(k.split('|')[2]) for k in growth["q95"]})
    # H is ordered, so one hue from light to dark. Text stays in ink colours. The table view is tables_<tag>.md.
    # Five steps is the most this ramp separates, so H = 0.6 stays in the table and not in the picture.
    ramp = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
    hs = [h for h in hs if h != 0.6]
    ink, muted = "#0b0b0b", "#52514e"
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=True)
    titles = {"raw": "Raw kernel", "inc1": "Increment kernel (inc1)"}
    handles = []
    for ax, kern in zip(axes, ("raw", "inc1")):
        ends = []
        for h, col in zip(hs, ramp):
            base = growth["q95"][f"{kern}|{h}|{ng[0]}"]["q95"]
            y = [growth["q95"][f"{kern}|{h}|{n}"]["q95"] / base for n in ng]
            lo = [growth["q95"][f"{kern}|{h}|{n}"]["lo"] / base for n in ng]
            hi = [growth["q95"][f"{kern}|{h}|{n}"]["hi"] / base for n in ng]
            ax.fill_between(ng, lo, hi, color=col, alpha=0.16, lw=0)
            line, = ax.plot(ng, y, marker="o", ms=5, lw=2, color=col, markeredgecolor="white", markeredgewidth=0.8, label=f"H = {h}")
            if kern == "raw":
                handles.append(line)
            ends.append((y[-1], h))
        # direct labels only where they do not collide: top down, at least 15% apart on the log scale
        last = None
        for v, h in sorted(ends, reverse=True):
            if last is None or last / v >= 1.15:
                ax.annotate(f"H = {h}", (ng[-1], v), xytext=(7, 0), textcoords="offset points", fontsize=8.5, color=muted, va="center")
                last = v
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(ng)
        ax.set_xticklabels([str(n) for n in ng], color=muted)
        ax.tick_params(axis="y", colors=muted)
        ax.tick_params(length=0)
        ax.minorticks_off()
        ticks = [0.5, 0.75, 1, 1.5, 2, 3, 4]
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t:g}" for t in ticks])
        ax.set_xlim(ng[0] * 0.88, ng[-1] * 1.85)
        ax.set_title(titles[kern], fontsize=11, loc="left", color=ink)
        ax.set_xlabel("sample size n", color=muted)
        ax.grid(True, axis="y", color="#ececec", lw=0.7)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color("#c9c8c3")
    axes[0].set_ylabel(f"95th percentile of T, relative to n = {ng[0]}", color=muted)
    fig.legend(handles=handles, loc="upper center", ncol=len(hs), frameon=False, fontsize=9, labelcolor=muted)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(HERE, f"growth_q95_{tag}.png"), dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    tag = sys.argv[1]
    with open(os.path.join(HERE, f"results_{tag}.json")) as f:
        write_figure(json.load(f)["growth"], tag)
