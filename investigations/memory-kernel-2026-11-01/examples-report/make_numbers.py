"""Write numbers.tex from data/report_numbers.json.

Each value becomes a key for the LaTeX command \\N{key}. An unknown key stops the build.
Rates print as percent with one decimal. Keys follow the JSON paths, joined by '|':
  size|<kernel>|<rule>|<H>|<n>, arms|<gen>|<kernel>|<rule>|<H>, power|<alt>|<kernel>|<rule>|<H>,
  mmd|<alt><param>|<kernel>|<H>. Add '|lo' or '|hi' for the 95% interval ends.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = json.loads((HERE / "data" / "report_numbers.json").read_text())
OUT = HERE / "numbers.tex"


def pct(x: float) -> str:
    return f"{100 * x:.1f}"


def main() -> None:
    m: dict[str, str] = {}

    def rate_block(prefix: str, table: dict) -> None:
        for k, r in table.items():
            m[f"{prefix}|{k}"] = pct(r["v"])
            m[f"{prefix}|{k}|lo"] = pct(r["lo"])
            m[f"{prefix}|{k}|hi"] = pct(r["hi"])

    rate_block("size", D["size"])
    rate_block("arms", D["arms"])
    rate_block("power", D["power"])
    rate_block("mmd", D["mmd_power"])

    for k, g in D["growth"].items():
        m[f"growth|{k}|ratio"] = f"{g['ratio_8000_500']:.2f}"
        m[f"growth|{k}|slope"] = f"{g['slope']['v']:.3f}"
        m[f"growth|{k}|slope|lo"] = f"{g['slope']['lo']:.3f}"
        m[f"growth|{k}|slope|hi"] = f"{g['slope']['hi']:.3f}"

    e = D["ex1"]
    for L, a in e["acf_at"].items():
        for kk, v in a.items():
            m[f"ex1|acf|{L}|{kk}"] = f"{v:.2f}"
    m["ex1|seed"] = str(e["seed"])
    m["ex1|pair|theory|raw"] = f"{e['pairsum_theory_raw'][-1]:.0f}"
    m["ex1|pair|theory|inc"] = f"{e['pairsum_theory_inc'][-1]:.2f}"
    m["ex1|pair|sample|raw"] = f"{e['pairsum_sample_raw'][-1]:.1f}"
    m["ex1|pair|sample|inc"] = f"{e['pairsum_sample_inc'][-1]:.2f}"
    m["ex1|pair|theory|raw|L100"] = f"{e['pairsum_theory_raw'][99]:.0f}"
    m["ex1|maxlag"] = str(e["maxlag"])

    for k, v in D["verdicts"].items():
        m[f"verdict|{k}|valid"] = v["valid"].replace("/", " of ")
        m[f"verdict|{k}|worst"] = f"{100 * v['worst_rate']:.1f}"

    rm = [k.split("|") for k, v in D["verdicts"].items() if v["removes_memory"]]
    rm = sorted(rm, key=lambda kr: (kr[1] != "shift", kr[0]))
    m["verdict|removes|list"] = ", ".join(r"\texttt{" + k.replace("_", r"\_") + "} with " + r for k, r in rm)
    m["verdict|removes|count"] = str(len(rm))
    m["verdict|total"] = str(len(D["verdicts"]))

    ws = [k for k, v in D["verdicts"].items() if k.split("|")[1] in ("wild", "sub")]
    m["verdict|wildsub|tested"] = str(len(ws))
    m["verdict|wildsub|pass"] = str(sum(D["verdicts"][k]["removes_memory"] for k in ws))

    # Split-MMD size (review C1).
    rate_block("mmdsize", {k: r for k, r in D["mmd_size"].items() if "v" in r})

    # Pooled rates and logistic H trend (review M1, M2).
    def pfmt(p: float) -> str:
        return "< 0.001" if p < 0.001 else f"= {p:.3f}" if p < 0.01 else f"= {p:.2f}"

    for k, t in D["trend"].items():
        m[f"trend|{k}|pooled"] = f"{100 * t['pooled']:.2f}"
        m[f"trend|{k}|lo"], m[f"trend|{k}|hi"] = f"{100 * t['lo']:.1f}", f"{100 * t['hi']:.1f}"
        m[f"trend|{k}|slope"] = f"{t['slope']:+.2f}"
        m[f"trend|{k}|p"] = pfmt(t["p"])
    memk = ["inc1", "herm1", "herm2", "sym", "datadiff", "dcov_inc1", "hp10", "hsplit1", "hsplit2"]
    sh = [D["trend"][f"{k}|shift"]["pooled"] for k in memk]
    m["trend|shift|min"], m["trend|shift|max"] = f"{100 * min(sh):.2f}", f"{100 * max(sh):.2f}"
    wl = [D["trend"][f"{k}|wild"]["pooled"] for k in memk]
    m["trend|wild|min"], m["trend|wild|max"] = f"{100 * min(wl):.1f}", f"{100 * max(wl):.1f}"
    flat6 = ["inc1", "herm1", "herm2", "sym", "datadiff", "hsplit1"]
    m["trend|flat6|minp"] = f"{min(D['trend'][f'{k}|{r}']['p'] for k in flat6 for r in ('wild', 'shift')):.2f}"
    m["trend|mem|m"] = f"{D['trend']['inc1|shift']['m']:,}".replace(",", "{,}")

    b = D["binom"]
    m["binom|cell"] = f"{100 * b['fail_cell_1000']:.1f}"
    m["binom|all24"] = f"{b['all24']:.2f}"

    # Median bandwidth and the falling inc1 null (review M4).
    bw = D["bandwidth"]
    m["bw|limit"] = f"{bw['sigma_limit']:.3f}"
    for H in ["0.5", "0.8", "0.9", "0.95"]:
        m[f"bw|lambda|{H}"] = f"{bw['lambda|' + H]:.3f}"
        for n in ["500", "1000", "2000", "4000", "8000"]:
            c = bw[f"{H}|{n}"]
            m[f"bw|{H}|{n}|sigma"] = f"{c['mean_sigma']:.3f}"
            m[f"bw|{H}|{n}|T"] = f"{c['mean_T_inc1']:.3f}"
            m[f"bw|{H}|{n}|formula"] = f"{c['formula_at_mean_sigma']:.3f}"
    for k, c in bw.items():
        if k.startswith("check|"):
            m[f"bw|{k}|fixed"] = f"{c['inc1_ratio_fixed']:.2f}"
            m[f"bw|{k}|sample"] = f"{c['inc1_ratio_sample']:.2f}"
    bw = {k: c for k, c in bw.items() if not k.startswith("check|")}
    err = [abs(c["formula_at_mean_sigma"] / c["mean_T_inc1"] - 1) for k, c in bw.items() if isinstance(c, dict)]
    m["bw|maxerr"] = f"{100 * max(err):.0f}"

    hv = [r["v"] for k, r in D["size"].items() if k.startswith("raw|harris|") and k.endswith("|2000")]
    m["size|raw|harris|2000|min"], m["size|raw|harris|2000|max"] = pct(min(hv)), pct(max(hv))

    hs = D["harris_stationary"]
    m["harris|min"], m["harris|max"], m["harris|cells"] = pct(hs["min"]), pct(hs["max"]), str(hs["cells"])
    m["sv|layer1"] = f"{D['sv_layers']['layer1']:.3f}"
    m["sv|layer2"] = f"{D['sv_layers']['layer2']:.3f}"
    lt = [abs(float(r["rate_est"]) - float(r["rate_true"])) for r in D["layer_true_vs_est"]]
    m["layer|maxdiff"] = f"{100 * max(lt):.1f}"
    m["layer|count"] = str(len(lt))
    run = D["run"]
    m["run|cells"], m["run|rows"] = str(run["cells"]), f"{run['rows']:,}".replace(",", "{,}")
    m["run|reps"] = f"{run['replications']['n<=2000']:,}".replace(",", "{,}")
    m["run|reps4000"] = str(run["replications"]["n=4000"])
    m["run|repsgrowth"] = f"{run['replications']['growth']:,}".replace(",", "{,}")
    m["run|errors"] = str(run["error_rows"])

    lines = ["% Generated by make_numbers.py from data/report_numbers.json. Do not edit.",
             r"\makeatletter",
             r"\newcommand{\N}[1]{\ifcsname rn@#1\endcsname\csname rn@#1\endcsname\else"
             r"\errmessage{Unknown number key: #1}\fi}"]
    for k, v in sorted(m.items()):
        lines.append(r"\expandafter\def\csname rn@" + k + r"\endcsname{" + v + "}")
    lines.append(r"\makeatother")
    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT} ({len(m)} keys)")


if __name__ == "__main__":
    main()
