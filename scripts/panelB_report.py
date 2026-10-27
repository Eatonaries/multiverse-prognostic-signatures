# -*- coding: utf-8 -*-
"""Panel B: variance decomposition restricted to the CV split (pure construction
freedom, with the validation-design lever held fixed). Also the real-vs-random
gap and the internal->external relationship, for all three cancers."""
import json

import numpy as np
import pandas as pd

R = "results/analysis"
KNOBS = ["norm", "depool", "fs", "ngenes", "cutoff", "split"]
KNOBS_B = ["norm", "depool", "fs", "ngenes", "cutoff"]


def eta2(df, col, knobs):
    grand = df[col].mean()
    ss_tot = float(((df[col] - grand) ** 2).sum())
    out = {}
    for k in knobs:
        ss_b = sum(len(g) * (g[col].mean() - grand) ** 2 for _, g in df.groupby(k))
        out[k] = float(ss_b / ss_tot) if ss_tot > 0 else 0.0
    return out


def swing(df, col, knobs):
    return {k: float(df.groupby(k)[col].mean().max() - df.groupby(k)[col].mean().min())
            for k in knobs}


out = {}
for stem in ("hnsc", "brca", "luad"):
    s = pd.read_csv(f"{R}/specs_internal_{stem}.csv").dropna(subset=["cindex_mean"])
    cv = s[s.split == "cv"]
    row = {
        "n_all": int(len(s)),
        "n_cv": int(len(cv)),
        "panelA_eta2": eta2(s, "cindex_mean", KNOBS),
        "panelA_swing": swing(s, "cindex_mean", KNOBS),
        "panelB_eta2": eta2(cv, "cindex_mean", KNOBS_B),
        "panelB_swing": swing(cv, "cindex_mean", KNOBS_B),
        "panelB_range": float(cv.cindex_mean.max() - cv.cindex_mean.min()),
        "split_marginals": {k: float(v) for k, v in s.groupby("split").cindex_mean.mean().items()},
    }
    # pair agreement: for CV specs, does the internal rank predict the external rank?
    a = pd.read_csv(f"{R}/attenuation_all.csv")
    a = a[a.cancer == f"TCGA-{stem.upper()}"]
    if not a.empty:
        row["att"] = {
            "n": int(len(a)),
            "ext_min": float(a.external_c.min()), "ext_med": float(a.external_c.median()),
            "ext_max": float(a.external_c.max()),
            "int_med": float(a.internal_cv_c.median()),
            # The manuscript reports the internal->external coupling as a RANK
            # correlation, so compute Spearman explicitly; keep Pearson too so the
            # two are never confused at write-up time.
            "corr_int_ext_spearman": float(
                a["internal_cv_c"].corr(a["external_c"], method="spearman")),
            "corr_int_ext_pearson": float(
                a[["internal_cv_c", "external_c"]].corr().iloc[0, 1]),
            # kept for backward compatibility with make_p4_report.py; equals the
            # Spearman value, which is the statistic the manuscript quotes.
            "corr_int_ext": float(
                a["internal_cv_c"].corr(a["external_c"], method="spearman")),
            "decayed_pct": float((a.delta < 0).mean() * 100),
            "delta_med": float(a.delta.median()),
            "delta_q05": float(a.delta.quantile(0.05)),
            "delta_q95": float(a.delta.quantile(0.95)),
            "delta_min": float(a.delta.min()),
            "delta_max": float(a.delta.max()),
            "n_cutoffs": int(a.cutoff.nunique()),
            "cutoffs": sorted(a.cutoff.unique().tolist()),
        }
    out[stem] = row

# real vs random gap
rc = pd.read_csv(f"{R}/random_control_summary.csv")
for stem, cname in (("hnsc", "TCGA-HNSC"), ("brca", "TCGA-BRCA"), ("luad", "TCGA-LUAD")):
    s = pd.read_csv(f"{R}/specs_internal_{stem}.csv")
    cv = s[s.split == "cv"]
    real = cv.groupby("ngenes")["cindex_mean"].mean()
    rnd = rc[rc.cancer == cname].groupby("ngenes")["mean"].mean()
    out[stem]["real_minus_random"] = {int(k): float(real[k] - rnd[k]) for k in real.index}
    out[stem]["real_cv_by_ngenes"] = {int(k): float(v) for k, v in real.items()}
    out[stem]["random_by_ngenes"] = {int(k): float(v) for k, v in rnd.items()}

json.dump(out, open(f"{R}/panelB_and_gaps.json", "w"), indent=1)
for stem, r in out.items():
    print("=" * 60)
    print(stem.upper(), f"(n_all={r['n_all']}, n_cv={r['n_cv']})")
    print(" Panel A eta2:", {k: round(v, 3) for k, v in r["panelA_eta2"].items()})
    print(" Panel B eta2:", {k: round(v, 3) for k, v in r["panelB_eta2"].items()})
    print(" Panel B range:", round(r["panelB_range"], 4))
    print(" split marginals:", {k: round(v, 4) for k, v in r["split_marginals"].items()})
    print(" real-random:", {k: f"{v:+.4f}" for k, v in r["real_minus_random"].items()})
    if "att" in r:
        a = r["att"]
        print(f" att: n={a['n']} ext {a['ext_min']:.3f}-{a['ext_max']:.3f} med {a['ext_med']:.3f} "
              f"| int med {a['int_med']:.3f} | spearman {a['corr_int_ext_spearman']:+.3f} "
              f"| pearson {a['corr_int_ext_pearson']:+.3f} | decayed {a['decayed_pct']:.0f}%"
              f" | cutoffs {a['cutoffs']}")
