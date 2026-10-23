"""Regenerate manuscript figures from the saved EXP1 cell summaries.

Run from any directory: python paper/long-memory-kernel-tests/make_figures.py
This script does not rerun simulations or change their summaries.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXP1 = REPO / "experiments" / "EXP1-20261016"
FIGURES = HERE / "figures"


def read_rows(name: str) -> list[dict[str, str]]:
    with (EXP1 / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def mmd_figure(rows: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    for process, label, color in (
        ("G1", "Gaussian fGn (rank one)", "#26659b"),
        ("G5", "Squared fGn (rank two)", "#b04c3f"),
        ("G6", "Absolute fGn (rank two)", "#729448"),
    ):
        subset = sorted(
            (row for row in rows if row["process"] == process
             and row["calibration"] == "wild_unc_05" and row["n"] == "2000"
             and row["power"] == "False"),
            key=lambda row: float(row["H"]),
        )
        if len(subset) != 7:
            raise ValueError(f"Expected seven MMD H cells for {process}, found {len(subset)}")
        h = [float(row["H"]) for row in subset]
        rate = [100 * float(row["size"]) for row in subset]
        low = [100 * (float(row["size"]) - float(row["ci_lo"])) for row in subset]
        high = [100 * (float(row["ci_hi"]) - float(row["size"])) for row in subset]
        ax.errorbar(h, rate, yerr=[low, high], label=label, color=color,
                    marker="o", markersize=4, capsize=2, linewidth=1.5)
    ax.axhline(5, color="black", linestyle="--", linewidth=1, label="Nominal 5%")
    ax.set(xlabel="Hurst parameter H", ylabel="Null rejection (%)", ylim=(0, 43))
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURES / "mmd_size.pdf")
    plt.close(fig)


def independence_figure(rows: list[dict[str, str]]) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 3.55), sharey=True)
    styles = (
        ("hsic_wb", "Wild multiplier", "#26659b"),
        ("hsic_shift", "Circular shift", "#729448"),
        ("hsic_sub", "Subsampling", "#b04c3f"),
    )
    for ax, prefix, title in ((left, "hsic", "HSIC"), (right, "dcov", "Distance covariance")):
        choices = styles if prefix == "hsic" else (
            ("dcov_wb", "Wild multiplier", "#26659b"),
            ("dcov_sub", "Subsampling", "#b04c3f"),
        )
        for test, label, color in choices:
            subset = [row for row in rows if row["test"] == test and row["n"] == "2000"]
            equal = sorted((row for row in subset if row["HX"] == row["HY"]),
                           key=lambda row: float(row["sumH"]))
            mixed = [row for row in subset if row["HX"] != row["HY"]]
            if len(equal) != 7 or len(mixed) != 3:
                raise ValueError(f"Unexpected E3 grid for {test}: {len(equal)} equal, {len(mixed)} mixed")
            ax.plot([float(r["sumH"]) for r in equal], [100 * float(r["size"]) for r in equal],
                    marker="o", markersize=4, linewidth=1.5, color=color, label=label)
            ax.scatter([float(r["sumH"]) for r in mixed], [100 * float(r["size"]) for r in mixed],
                       marker="x", s=30, linewidth=1.3, color=color)
        ax.axhline(5, color="black", linestyle="--", linewidth=1)
        ax.axvline(1.5, color="#666666", linestyle=":", linewidth=1)
        ax.set(title=title, xlabel=r"Memory sum $H_X+H_Y$", xlim=(1.27, 1.73), ylim=(0, 41))
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="upper left", fontsize=7, frameon=False)
    left.set_ylabel("Null rejection (%)")
    fig.tight_layout()
    fig.savefig(FIGURES / "independence_size.pdf")
    plt.close(fig)


def adjacent_slopes(x: list[float], y: list[float]) -> tuple[list[float], list[float]]:
    """Adjacent change in rejection percentage points per 0.1 memory unit."""
    if len(x) != len(y) or len(x) < 2 or any(b <= a for a, b in zip(x, x[1:])):
        raise ValueError("Expected aligned, strictly increasing memory grid")
    midpoint = [(a + b) / 2 for a, b in zip(x, x[1:])]
    slope = [0.1 * (b - a) / (xb - xa)
             for xa, xb, a, b in zip(x, x[1:], y, y[1:])]
    return midpoint, slope


def slope_figure(mmd_rows: list[dict[str, str]], ind_rows: list[dict[str, str]]) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 3.5))
    for process, label, color in (
        ("G1", "fGn, rank one", "#26659b"),
        ("G5", "Squared fGn, rank two", "#b04c3f"),
        ("G6", "Absolute fGn, rank two", "#729448"),
    ):
        cells = sorted(
            (row for row in mmd_rows if row["process"] == process
             and row["calibration"] == "wild_unc_05" and row["n"] == "2000"
             and row["power"] == "False"),
            key=lambda row: float(row["H"]),
        )
        if len(cells) != 7:
            raise ValueError(f"Unexpected MMD slope grid for {process}")
        x, y = adjacent_slopes(
            [float(row["H"]) for row in cells],
            [100 * float(row["size"]) for row in cells],
        )
        left.plot(x, y, marker="o", markersize=4, color=color, label=label)

    for test, label, color in (
        ("hsic_wb", "Wild multiplier", "#26659b"),
        ("hsic_shift", "Circular shift", "#729448"),
        ("hsic_sub", "Subsampling", "#b04c3f"),
    ):
        cells = sorted(
            (row for row in ind_rows if row["test"] == test and row["n"] == "2000"
             and row["HX"] == row["HY"]),
            key=lambda row: float(row["sumH"]),
        )
        if len(cells) != 7:
            raise ValueError(f"Unexpected HSIC slope grid for {test}")
        x, y = adjacent_slopes(
            [float(row["sumH"]) for row in cells],
            [100 * float(row["size"]) for row in cells],
        )
        right.plot(x, y, marker="o", markersize=4, color=color, label=label)

    for ax in (left, right):
        ax.axhline(0, color="black", linestyle="--", linewidth=1)
        ax.grid(axis="y", alpha=0.2)
        ax.legend(fontsize=7, frameon=False)
    left.set(title="Split MMD", xlabel="Midpoint H",
             ylabel="Change in rejection (pp per 0.1 H)")
    right.set(title="HSIC, equal-memory pairs", xlabel=r"Midpoint $H_X+H_Y$",
              ylabel="Change in rejection (pp per 0.1 sum)")
    fig.tight_layout()
    fig.savefig(FIGURES / "size_slopes.pdf")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    mmd_rows = read_rows("e2_summary.csv")
    ind_rows = read_rows("e3_summary.csv")
    mmd_figure(mmd_rows)
    independence_figure(ind_rows)
    slope_figure(mmd_rows, ind_rows)
    manifest = {
        "description": "Figures reproduce saved EXP1 summaries; no simulations rerun.",
        "inputs": {
            name: hashlib.sha256((EXP1 / name).read_bytes()).hexdigest()
            for name in ("e2_summary.csv", "e3_summary.csv")
        },
        "outputs": {
            name: hashlib.sha256((FIGURES / name).read_bytes()).hexdigest()
            for name in ("mmd_size.pdf", "independence_size.pdf", "size_slopes.pdf")
        },
    }
    (FIGURES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
