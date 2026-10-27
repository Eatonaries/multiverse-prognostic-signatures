# -*- coding: utf-8 -*-
"""Aggregate the per-fold stage-B output into a cut-off-effect table.

Purpose
-------
Figure S3 (manuscript/图注.md) promises: "Cut-off effects: hazard ratio and
log-rank p for the three cut-off rules, holding the signature fixed."

The main analysis table (results/analysis/specs_internal_<stem>.csv) already
carries the two cut-off-sensitive readouts at specification level --
`hr_median` and `sig_rate` (the fraction of folds with p < 0.05) -- but NOT
the log-rank p itself. This script rebuilds the same specification-level
table directly from the per-fold stage-B chunks and adds the p-value
aggregates, so the log-rank p can be plotted per cut-off rule.

The aggregation reproduces `agg()` in scripts/aggregate_analyze.py
(lines 50-69) EXACTLY for cindex_mean / cindex_sd / cindex_min /
cindex_max / n_folds / sig_rate / hr_median, and adds:
    p_median         median log-rank p over the folds of that specification
    p_min            smallest (most significant) log-rank p
    p_mean           arithmetic mean log-rank p
    neg_log10_p_med  median of -log10(p) (a scale on which "more
                     significant" is larger; this is the usual way to plot a
                     p-value distribution)

Writes: results/analysis/cutoff_effects_<stem>.csv
Usage:  python aggregate_cutoff_effects.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results"
OUT = SRC / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

CANCERS = {"TCGA-HNSC": "hnsc", "TCGA-BRCA": "brca", "TCGA-LUAD": "luad"}
KNOBS = ["norm", "depool", "fs", "ngenes", "cutoff", "split"]


def load_chunks(stem, prefix="stageB"):
    d = SRC / f"{prefix}_{stem}"
    files = sorted(d.glob("*.csv")) if d.exists() else []
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def agg_cutoff(df, group, extval="internal"):
    """Same aggregation as aggregate_analyze.agg(), plus p-value aggregates."""
    sub = df[df["extval"] == extval]
    if sub.empty:
        return pd.DataFrame()
    rows = []
    for keys, g in sub.groupby(group, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        rec = dict(zip(group, keys))
        ci = g["cindex"].dropna()
        rec["cindex_mean"] = float(ci.mean()) if len(ci) else np.nan
        rec["cindex_sd"] = float(ci.std(ddof=1)) if len(ci) > 1 else 0.0
        rec["cindex_min"] = float(ci.min()) if len(ci) else np.nan
        rec["cindex_max"] = float(ci.max()) if len(ci) else np.nan
        rec["n_folds"] = int(len(ci))
        pl = g["p_logrank"].dropna()
        rec["sig_rate"] = float((pl < 0.05).mean()) if len(pl) else np.nan
        hr = g["hr"].dropna()
        rec["hr_median"] = float(hr.median()) if len(hr) else np.nan
        # ---- added here -------------------------------------------------
        rec["p_median"] = float(pl.median()) if len(pl) else np.nan
        rec["p_min"] = float(pl.min()) if len(pl) else np.nan
        rec["p_mean"] = float(pl.mean()) if len(pl) else np.nan
        # -log10(p); clip at 1e-300 so log10 is finite
        with np.errstate(divide="ignore"):
            nl = -np.log10(np.clip(pl.to_numpy(), 1e-300, 1.0))
        rec["neg_log10_p_med"] = float(np.median(nl)) if len(nl) else np.nan
        rec["frac_p_lt_0.05"] = float((pl < 0.05).mean()) if len(pl) else np.nan
        rec["frac_p_lt_0.01"] = float((pl < 0.01).mean()) if len(pl) else np.nan
        rec["frac_p_lt_0.001"] = float((pl < 0.001).mean()) if len(pl) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    for proj, stem in CANCERS.items():
        df = load_chunks(stem)
        if df.empty:
            print(f"{proj}: no stageB chunks -- skipped")
            continue
        # identical scope to the main table: internal rows, and the splits
        # that actually hold out data ("full" trains on everything, so it
        # has no held-out set and its cindex/hr/p are NA)
        sub = df[df["split"] != "full"]
        spec = agg_cutoff(sub, KNOBS)
        spec.insert(0, "cancer", proj)
        path = OUT / f"cutoff_effects_{stem}.csv"
        spec.to_csv(path, index=False)

        # ---- cross-check against the main analysis table ----------------
        main_p = OUT / f"specs_internal_{stem}.csv"
        msg = []
        if main_p.exists():
            m = pd.read_csv(main_p)
            j = spec.merge(m, on=["cancer"] + KNOBS, suffixes=("_new", "_old"))
            for col in ("cindex_mean", "cindex_sd", "cindex_min", "cindex_max",
                        "n_folds", "sig_rate", "hr_median"):
                a = j[f"{col}_new"].to_numpy()
                b = j[f"{col}_old"].to_numpy()
                same = np.allclose(a, b, rtol=0, atol=1e-12, equal_nan=True)
                msg.append(f"{col}={'OK' if same else '*** DIFF'}")
            msg.append(f"rows={len(spec)} (main {len(m)})")
            msg.append(f"merged={len(j)}")
        print(f"{proj}: wrote {path.name}  " + "  ".join(msg))
    return 0


if __name__ == "__main__":
    sys.exit(main())
