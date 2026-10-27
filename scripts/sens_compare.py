# -*- coding: utf-8 -*-
"""Compare the strict and lenient sensitivity runs: does either fixed setting
change the conclusion?

Usage: python sens_compare.py <stem>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KNOBS = ["norm", "depool", "fs", "ngenes", "cutoff"]


def load(stem):
    d = ROOT / "results" / f"stageB_sens_{stem}"
    files = sorted(d.glob("*.csv")) if d.exists() else []
    if not files:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df["version"] = df["split"].str.split("__").str[0]
    # stageB derives 'split' from the stage-A split_kind; rebuild the version
    # label from the file name recorded in the payload instead
    return df


def load_v2(stem):
    d = ROOT / "results" / f"stageB_sens_{stem}"
    out = []
    for f in sorted(d.glob("*.csv")) if d.exists() else []:
        v = f.name.split("__")[0]
        t = pd.read_csv(f)
        t["version"] = v
        out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def spec_table(df):
    g = (df[df["extval"] == "internal"]
         .groupby(KNOBS + ["version"])["cindex"]
         .agg(["mean", "std", "count"]).reset_index())
    return g


def main():
    stem = sys.argv[1]
    df = load_v2(stem)
    if df.empty:
        print("no sensitivity data")
        return 1
    sp = spec_table(df)
    rep = {"stem": stem, "versions": {}}
    for v in sorted(sp["version"].unique()):
        s = sp[sp["version"] == v]
        swings = {k: float(s.groupby(k)["mean"].mean().max()
                           - s.groupby(k)["mean"].mean().min()) for k in KNOBS}
        # rank knobs by marginal swing
        rank = [k for k, _ in sorted(swings.items(), key=lambda kv: -kv[1])]
        rep["versions"][v] = {
            "n_specs": int(len(s)),
            "cindex_min": float(s["mean"].min()),
            "cindex_max": float(s["mean"].max()),
            "range": float(s["mean"].max() - s["mean"].min()),
            "frac_ge_060": float((s["mean"] >= 0.60).mean()),
            "marginal_swing": swings,
            "knob_order": rank,
        }
    # agreement
    if len(rep["versions"]) == 2:
        a, b = [rep["versions"][v] for v in ("strict", "lenient")
                if v in rep["versions"]]
        from scipy import stats
        ia = pd.Series(a["marginal_swing"])
        ib = pd.Series(b["marginal_swing"])
        rho, pv = stats.spearmanr(ia, ib)
        rep["agreement"] = {
            "knob_order_strict": a["knob_order"],
            "knob_order_lenient": b["knob_order"],
            "swing_rank_spearman_rho": float(rho),
            "swing_rank_p": float(pv),
            "range_diff": float(b["range"] - a["range"]),
            "same_conclusion_range_ge_005": bool(min(a["range"], b["range"]) >= 0.05),
            "same_top_knob": a["knob_order"][0] == b["knob_order"][0],
        }
    out = ROOT / "results" / "analysis" / f"sensitivity_{stem}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1), encoding="utf-8")
    print(json.dumps(rep, indent=1))
    sp.to_csv(ROOT / "results" / "analysis" / f"sensitivity_specs_{stem}.csv", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
