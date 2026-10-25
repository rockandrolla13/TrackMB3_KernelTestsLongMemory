"""EXP2 multiple-testing re-analysis (no new simulation).

Rebuilds the test family behind verdicts_main.json, attaches p-values, and adjusts them with
Holm (primary), Bonferroni, Benjamini-Hochberg and Benjamini-Yekutieli. Also checks the
quantitative claims of paper section 4.3 ("Process features").

Usage: python -I multiple_testing/mt_analysis.py   ->  mt_results.csv, mt_claims.csv, mt_summary.json
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.path.insert(0, EXP)
from scipy import stats  # noqa: E402

import analyze as A  # noqa: E402  (frozen EXP2 code: load, verdict, interaction, RULES)
import gens  # noqa: E402

ALPHA = 0.05
FLOOR = 0.02  # the 2 pp floor of the pre-declared reading rule


# ---------- p-values ----------
def two_prop(k1, m1, k0, m0):
    """Pooled two-proportion z for p1 - p0. Returns (diff, z, p_one_sided_greater, p_two_sided)."""
    p1, p0 = k1 / m1, k0 / m0
    pp = (k1 + k0) / (m1 + m0)
    se = math.sqrt(pp * (1 - pp) * (1 / m1 + 1 / m0))
    z = (p1 - p0) / se if se > 0 else 0.0
    return p1 - p0, z, stats.norm.sf(z), 2 * stats.norm.sf(abs(z))


def fisher(k1, m1, k0, m0, alt):
    return stats.fisher_exact([[k1, m1 - k1], [k0, m0 - k0]], alternative=alt).pvalue


def interaction_test(d, fam, g, rule, n, h):
    """Same contrast and unpooled SE as analyze.interaction. Returns (I, z, p_two_sided)."""
    i, z = A.interaction(d, fam, g, rule, n, h)
    return i, z, 2 * stats.norm.sf(abs(z))


# ---------- adjustments ----------
def holm(p):
    m = len(p)
    order = sorted(range(m), key=lambda j: p[j])
    out, run = [0.0] * m, 0.0
    for r, j in enumerate(order):
        run = max(run, min(1.0, (m - r) * p[j]))
        out[j] = run
    return out


def bonf(p):
    return [min(1.0, len(p) * x) for x in p]


def bh(p, c=1.0):
    m = len(p)
    order = sorted(range(m), key=lambda j: p[j], reverse=True)
    out, run = [0.0] * m, 1.0
    for r, j in enumerate(order):
        rank = m - r
        run = min(run, min(1.0, c * m * p[j] / rank))
        out[j] = run
    return out


def by(p):
    return bh(p, c=sum(1.0 / k for k in range(1, len(p) + 1)))


ADJ = {"holm": holm, "bonf": bonf, "bh": bh, "by": by}


def adjust(rows, key, scope_fn, tag):
    """Adjust rows[key] within groups given by scope_fn; write rows[f'{tag}_{method}']."""
    groups = {}
    for r in rows:
        s = scope_fn(r)
        if s is not None:
            groups.setdefault(s, []).append(r)
    for g in groups.values():
        p = [r[key] for r in g]
        for name, f in ADJ.items():
            for r, q in zip(g, f(p)):
                r[f"{tag}_{name}"] = q


# ---------- build the family ----------
def build(d):
    rows = []
    for fam in ("mmd", "ind"):
        for rule in A.VERDICT_RULES[fam]:
            for n in A.N_GRID:
                for g in gens.MAIN[1:]:
                    a, b = d[(fam, g, 0.5, n)][rule], d[(fam, "G0", 0.5, n)][rule]
                    diff, z, p1, p2 = two_prop(a["k"], a["m"], b["k"], b["m"])
                    rows.append({"kind": "alone", "family": fam, "rule": rule, "n": n, "gen": g, "h": 0.5,
                                 "effect_pp": 100 * diff, "z": z, "p_raw": p1, "p_two_sided": p2,
                                 "p_fisher": fisher(a["k"], a["m"], b["k"], b["m"], "greater"),
                                 "flag_uncorrected": bool(diff > FLOOR and a["lo"] > b["hi"])})
                    for h in (0.7, 0.8):
                        i, zi, pi = interaction_test(d, fam, g, rule, n, h)
                        rows.append({"kind": "interaction", "family": fam, "rule": rule, "n": n, "gen": g, "h": h,
                                     "effect_pp": 100 * i, "z": zi, "p_raw": pi, "p_two_sided": pi,
                                     "p_fisher": "", "flag_uncorrected": bool(abs(i) > FLOOR and abs(zi) > 1.96)})
    return rows


def level_tests(d):
    """Supplementary family: F(h=0.8) vs G0(h=0.8), two-sided. Several section-4.3 claims are
    level differences at H = 0.8, which the verdict family does not test."""
    rows = []
    for fam in ("mmd", "ind"):
        for rule in A.VERDICT_RULES[fam]:
            for n in A.N_GRID:
                for g in gens.MAIN[1:]:
                    a, b = d[(fam, g, 0.8, n)][rule], d[(fam, "G0", 0.8, n)][rule]
                    diff, z, _, p2 = two_prop(a["k"], a["m"], b["k"], b["m"])
                    rows.append({"kind": "level_h0.8", "family": fam, "rule": rule, "n": n, "gen": g, "h": 0.8,
                                 "effect_pp": 100 * diff, "z": z, "p_raw": p2, "p_two_sided": p2,
                                 "p_fisher": fisher(a["k"], a["m"], b["k"], b["m"], "two-sided"),
                                 "flag_uncorrected": bool(abs(diff) > FLOOR and p2 < ALPHA)})
    return rows


# ---------- corrected verdicts ----------
def corrected_verdict(rs, col):
    """Re-apply the DESIGN.md reading rule with 'significant' = adjusted p < 0.05 (2 pp floor kept)."""
    al = rs[("alone", 0.5)]
    A_ = al["effect_pp"] > 100 * FLOOR and al[col] < ALPHA
    ints = [rs[("interaction", h)] for h in (0.7, 0.8)]
    ip = any(r["effect_pp"] > 100 * FLOOR and r["z"] > 0 and r[col] < ALPHA for r in ints)
    im = any(r["effect_pp"] < -100 * FLOOR and r["z"] < 0 and r[col] < ALPHA for r in ints)
    if A_ and ip:
        v = "breaks alone and amplifies memory"
    elif A_:
        v = "breaks alone"
    elif ip:
        v = "only with memory"
    elif im:
        v = "dampens memory"
    else:
        v = "no effect"
    if im and not ip and A_:
        v = "breaks alone; dampens memory"
    return v


def components(v):
    """Set of findings a verdict asserts, for partial-survival counting."""
    s = set()
    if v.startswith("breaks alone"):
        s.add("A")
    if "amplifies" in v or v == "only with memory":
        s.add("I+")
    if "dampens" in v:
        s.add("I-")
    return s


# ---------- section 4.3 claims ----------
def claims(rows_by_key):
    def R(kind, fam, rule, n, g, h):
        return rows_by_key[(kind, fam, rule, n, g, h)]

    C = []

    def add(claim, kind, fam, rule, g, h, n=2000, note=""):
        r = R(kind, fam, rule, n, g, h)
        C.append({"claim": claim, "test": f"{kind} {fam}/{rule} {g} H={h} n={n}", "effect_pp": r["effect_pp"],
                  "p_raw": r["p_raw"], "p_two_sided": r["p_two_sided"],
                  **{k: r[k] for k in r if k.startswith(("all_", "arm_", "aug_", "two_"))},
                  "note": note})

    for rule in ("hsic_sub", "dcov_sub"):
        add("HT-L adds ~8 (HSIC) / ~3 (dCov) points at H=0.8", "level_h0.8", "ind", rule, "HT-L", 0.8,
            note="claim as worded: level gap vs G0 at H=0.8")
        add("HT-L adds ~8 / ~3 points at H=0.8", "interaction", "ind", rule, "HT-L", 0.8,
            note="same claim read as interaction")
    add("HT-L: no excess at H=0.5", "alone", "ind", "dcov_sub", "HT-L", 0.5,
        note="one-sided p; HT-L is BELOW G0 here, see two-sided")
    add("HT-L: no excess at H=0.5", "alone", "ind", "hsic_sub", "HT-L", 0.5)
    for rule in ("hsic_sub", "dcov_sub"):
        add("NL-TAR adds 11-12 points to subsampling at H=0.8", "level_h0.8", "ind", rule, "NL-TAR", 0.8,
            note="claim as worded: level gap")
        add("NL-TAR adds 11-12 points to subsampling at H=0.8", "interaction", "ind", rule, "NL-TAR", 0.8)
    for rule in ("hsic_sub", "dcov_sub"):
        add("NL-SV keeps subsampling <=6.2% at H=0.8 (vs G0 17.5/22.8)", "level_h0.8", "ind", rule, "NL-SV", 0.8)
        add("NL-SV dampens memory effect", "interaction", "ind", rule, "NL-SV", 0.8)
    add("NL-SV breaks MMD sqrt subsampling alone (13.3% at H=0.5)", "alone", "mmd", "mmd_sub_sqrt", "NL-SV", 0.5)
    for g in ("SEAS-M", "SEAS-V"):
        for rule in ("hsic_wild", "hsic_sub", "dcov_wild", "dcov_sub"):
            add(f"{g} breaks kernel tests alone", "alone", "ind", rule, g, 0.5)
        add(f"{g}: circular shift not broken", "alone", "ind", "hsic_shift", g, 0.5, note="null claim")
    for g in ("SEAS-M", "SEAS-V"):
        for rule in ("mmd_sub_sqrt", "mmd_sn_block"):
            add(f"{g} makes split MMD conservative at H=0.5", "alone", "mmd", rule, g, 0.5,
                note="direction is BELOW G0; one-sided p is ~1, see two-sided")
    for g in ("HT-S", "HT-L", "NL-SV", "NL-TAR"):
        add("MMD unc wild 26-31% at H=0.8, no feature effect", "interaction", "mmd", "mmd_unc_wild", g, 0.8,
            note="null claim: no interaction expected")
    return C


def main():
    d = A.load("main")
    rows = build(d)
    lv = level_tests(d)
    # scopes: overall verdict family; per arm; n=2000 only (sensitivity); augmented = verdict + level family
    adjust(rows, "p_raw", lambda r: "all", "all")
    adjust(rows, "p_raw", lambda r: r["family"], "arm")
    adjust(rows, "p_raw", lambda r: "n2000" if r["n"] == 2000 else None, "n2000")
    two = [dict(r) for r in rows]
    adjust(two, "p_two_sided", lambda r: "all", "two")  # alone tested two-sided
    for r, t in zip(rows, two):
        r["two_holm"], r["two_by"] = t["two_holm"], t["two_by"]
    aug = rows + lv
    for r in aug:
        r["_aug_p"] = r["p_raw"]
    adjust(aug, "_aug_p", lambda r: "aug", "aug")
    adjust(lv, "p_raw", lambda r: "all", "all")
    adjust(lv, "p_raw", lambda r: r["family"], "arm")

    out = rows + lv
    cols = ["kind", "family", "rule", "n", "gen", "h", "effect_pp", "z", "p_raw", "p_two_sided", "p_fisher",
            "flag_uncorrected", "all_holm", "all_bonf", "all_bh", "all_by", "arm_holm", "arm_bonf", "arm_bh",
            "arm_by", "n2000_holm", "n2000_by", "two_holm", "two_by", "aug_holm", "aug_bonf", "aug_bh", "aug_by"]
    with open(os.path.join(HERE, "mt_results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in out:
            w.writerow([(f"{r[c]:.4g}" if isinstance(r.get(c), float) else r.get(c, "")) for c in cols])

    # verdict survival
    by_key = {(r["kind"], r["family"], r["rule"], r["n"], r["gen"], r["h"]): r for r in out}
    old = json.load(open(os.path.join(EXP, "verdicts_main.json")))
    surv = {}
    for col in ("all_holm", "arm_holm", "all_by", "all_bh", "two_holm"):
        tab = []
        for key, v in old.items():
            rule, ns, g = key.split("|")
            n = int(ns[1:])
            fam = "mmd" if rule.startswith("mmd") else "ind"
            assert A.verdict(d, fam, g, rule, n)["verdict"] == v["verdict"]  # family rebuild matches frozen file
            rs = {("alone", 0.5): by_key[("alone", fam, rule, n, g, 0.5)],
                  **{("interaction", h): by_key[("interaction", fam, rule, n, g, h)] for h in (0.7, 0.8)}}
            nv = corrected_verdict(rs, col)
            tab.append({"key": key, "n": n, "old": v["verdict"], "new": nv,
                        "full": nv == v["verdict"], "kept": components(v["verdict"]) & components(nv),
                        "lost": components(v["verdict"]) - components(nv),
                        "gained": components(nv) - components(v["verdict"])})
        nonnull = [t for t in tab if t["old"] != "no effect"]
        surv[col] = {"verdicts": len(tab), "non_null": len(nonnull),
                     "non_null_n2000": sum(t["n"] == 2000 for t in nonnull),
                     "survive_exact": sum(t["full"] for t in nonnull),
                     "survive_exact_n2000": sum(t["full"] for t in nonnull if t["n"] == 2000),
                     "survive_any_component": sum(bool(t["kept"]) for t in nonnull),
                     "new_nonnull": sum(t["old"] == "no effect" and t["new"] != "no effect" for t in tab),
                     "changed": [{k: (sorted(t[k]) if isinstance(t[k], set) else t[k]) for k in t}
                                 for t in tab if not t["full"]]}
    # uncorrected flag counts
    fam_counts = {"verdict_family_total": len(rows),
                  "alone": sum(r["kind"] == "alone" for r in rows),
                  "interaction": sum(r["kind"] == "interaction" for r in rows),
                  "mmd": sum(r["family"] == "mmd" for r in rows), "ind": sum(r["family"] == "ind" for r in rows),
                  "level_family": len(lv), "augmented": len(aug),
                  "flags_uncorrected": sum(bool(r["flag_uncorrected"]) for r in rows),
                  "raw_p_lt_05": int(sum(r["p_raw"] < ALPHA for r in rows)),
                  "holm_lt_05": int(sum(r["all_holm"] < ALPHA for r in rows)),
                  "by_lt_05": int(sum(r["all_by"] < ALPHA for r in rows)),
                  "bh_lt_05": int(sum(r["all_bh"] < ALPHA for r in rows)),
                  "bonf_lt_05": int(sum(r["all_bonf"] < ALPHA for r in rows)),
                  "flag_and_holm": int(sum(bool(r["flag_uncorrected"]) and r["all_holm"] < ALPHA for r in rows)),
                  "flag_and_by": int(sum(bool(r["flag_uncorrected"]) and r["all_by"] < ALPHA for r in rows)),
                  "alone_fisher_disagree": [f'{r["rule"]}|n{r["n"]}|{r["gen"]} z_p={r["p_raw"]:.3g} fisher_p={r["p_fisher"]:.3g}'
                                            for r in rows if r["kind"] == "alone"
                                            and (r["p_fisher"] < ALPHA) != (r["p_raw"] < ALPHA)]}
    cl = claims(by_key)
    ccols = ["claim", "test", "effect_pp", "p_raw", "p_two_sided", "all_holm", "all_by", "arm_holm", "arm_by",
             "two_holm", "two_by", "aug_holm", "aug_by", "note"]
    with open(os.path.join(HERE, "mt_claims.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(ccols)
        for c in cl:
            w.writerow([(f"{c[k]:.4g}" if isinstance(c.get(k), float) else c.get(k, "")) for k in ccols])
    json.dump({"family": fam_counts, "survival": surv}, open(os.path.join(HERE, "mt_summary.json"), "w"),
              indent=1, default=str)
    print(json.dumps(fam_counts, indent=1, default=str))
    for col, s in surv.items():
        print(col, {k: v for k, v in s.items() if k != "changed"})


if __name__ == "__main__":
    main()
