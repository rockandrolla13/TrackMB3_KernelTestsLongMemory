"""EXP8 chart for the talk: false-alarm rates on the real-data surrogate null, with the real-real placebo for
comparison. Numbers come from experiments/EXP8-stock-surrogate-20261108/results.json (pooled, n = 2000).
Run: python make_exp8_figure.py
"""
from __future__ import annotations

import json

import matplotlib.pyplot as plt

from make_talk_figures import INK, MUTED, OUT, REPO, SHIFT, SUB, WILD, save

ROWS = [("raw", "sub", "ordinary kernel, subsampling", SUB), ("raw", "wild", "ordinary kernel, wild", WILD),
        ("raw", "shift", "ordinary kernel, circular shift", SHIFT), ("inc1", "sub", "increment kernel, subsampling", SUB),
        ("inc1", "shift", "increment kernel, circular shift", SHIFT), ("raw", "harris", "ordinary kernel, Harris shift", MUTED)]


def main():
    res = json.load(open(REPO / "experiments/EXP8-stock-surrogate-20261108/results.json"))
    pooled = {(r["arm"], r["n"], r["kernel"], r["rule"]): r for r in res["pooled"]}
    fig, ax = plt.subplots(figsize=(10, 5.0))
    for i, (k, rule, label, col) in enumerate(ROWS):
        y = len(ROWS) - 1 - i
        s = pooled[("S1", 2000, k, rule)]
        a = pooled[("A", 2000, k, rule)]
        ax.barh(y, s["rate"], height=0.55, color=col, alpha=0.85)
        ax.plot([s["lo"], s["hi"]], [y, y], color=INK, lw=1.4)
        ax.scatter(a["rate"], y, marker="D", s=60, facecolor="white", edgecolor=INK, lw=1.4, zorder=3)
        ax.text(max(s["hi"], a["rate"]) + 0.25, y, f"{s['rate']:.1f}%", va="center", color=INK, fontsize=13)
    ax.axvline(5, color=INK, lw=1, ls="--")
    ax.set_ylim(-0.6, len(ROWS) - 0.2)
    ax.text(5.08, len(ROWS) - 0.45, "promised 5%", color=INK, fontsize=12)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([r[2] for r in ROWS][::-1])
    ax.set_xlim(0, 11)
    ax.set_xlabel("false alarms (%), n = 2000, four stock pairings pooled")
    ax.scatter([], [], marker="s", color=MUTED, label="bar: real window vs shuffled copy (known null)")
    ax.scatter([], [], marker="D", facecolor="white", edgecolor=INK, label="diamond: real vs real windows (placebo)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.35, 1.04), ncol=2, fontsize=11)
    save(fig, "EXP8-surrogate-null")


if __name__ == "__main__":
    main()
