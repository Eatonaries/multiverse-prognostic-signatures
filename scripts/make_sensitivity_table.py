# -*- coding: utf-8 -*-
"""Supplementary Table S7: strict vs lenient sensitivity comparison.

Reads the per-cancer sensitivity JSON written by sens_compare.py and emits a
flat CSV that a reader can check by hand. No numbers are hand-entered.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\deepseek\jieru\_bio4")
ANA = ROOT / "results" / "analysis"
SUPP = ROOT / "results" / "supplementary"
SUPP.mkdir(parents=True, exist_ok=True)

CANCERS = [("hnsc", "TCGA-HNSC"), ("brca", "TCGA-BRCA"), ("luad", "TCGA-LUAD")]
KNOBS = ["norm", "depool", "fs", "ngenes", "cutoff"]

rows = []
for stem, label in CANCERS:
    p = ANA / f"sensitivity_{stem}.json"
    if not p.exists():
        print(f"MISSING {p.name}")
        continue
    d = json.load(open(p, encoding="utf-8"))
    for ver in ("strict", "lenient"):
        v = d["versions"][ver]
        for k in KNOBS:
            rows.append({
                "cohort": label,
                "version": ver,
                "decision_point": k,
                "marginal_swing_cindex": round(v["marginal_swing"][k], 6),
                "rank_by_swing": v["knob_order"].index(k) + 1,
            })
        rows.append({
            "cohort": label, "version": ver, "decision_point": "__SUMMARY__",
            "marginal_swing_cindex": "",
            "rank_by_swing": "",
            "n_specs": v["n_specs"],
            "cindex_min": round(v["cindex_min"], 6),
            "cindex_max": round(v["cindex_max"], 6),
            "range_cindex": round(v["range"], 6),
            "frac_ge_060": round(v["frac_ge_060"], 6),
            "knob_order": ">".join(v["knob_order"]),
        })
    a = d["agreement"]
    rows.append({
        "cohort": label, "version": "AGREEMENT", "decision_point": "__AGREEMENT__",
        "marginal_swing_cindex": "", "rank_by_swing": "",
        "range_diff_strict_minus_lenient": round(a["range_diff"], 6),
        "swing_rank_spearman_rho": round(a["swing_rank_spearman_rho"], 4),
        "swing_rank_p": round(a["swing_rank_p"], 4),
        "same_top_knob": a["same_top_knob"],
        "same_conclusion_range_ge_005": a["same_conclusion_range_ge_005"],
    })

fields = ["cohort", "version", "decision_point", "marginal_swing_cindex", "rank_by_swing",
          "n_specs", "cindex_min", "cindex_max", "range_cindex", "frac_ge_060", "knob_order",
          "range_diff_strict_minus_lenient", "swing_rank_spearman_rho", "swing_rank_p",
          "same_top_knob", "same_conclusion_range_ge_005"]

out = SUPP / "Table_S7_sensitivity_strict_vs_lenient.csv"
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)
print("wrote", out, out.stat().st_size, "bytes,", len(rows), "rows")

# console summary
print()
for stem, label in CANCERS:
    p = ANA / f"sensitivity_{stem}.json"
    if not p.exists():
        continue
    d = json.load(open(p, encoding="utf-8"))
    s, l = d["versions"]["strict"], d["versions"]["lenient"]
    a = d["agreement"]
    ng_top2 = all("ngenes" in v["knob_order"][:2] for v in (s, l))
    dep_bottom2 = all("depool" in v["knob_order"][-2:] for v in (s, l))
    print(f"{label}: range {s['range']:.4f}/{l['range']:.4f} (diff {a['range_diff']:+.4f}) | "
          f"top knob strict={s['knob_order'][0]} lenient={l['knob_order'][0]} same={a['same_top_knob']} | "
          f"ngenes in top-2 both={ng_top2} | depool in bottom-2 both={dep_bottom2} | "
          f"rho={a['swing_rank_spearman_rho']:.2f} p={a['swing_rank_p']:.3f}")
