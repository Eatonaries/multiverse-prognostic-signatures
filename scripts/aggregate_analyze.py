# -*- coding: utf-8 -*-
"""P4: aggregate the multiverse grid -> specification curve + variance
decomposition + external attenuation + random-signature comparison.

Spec space (per cancer, internal validation):
    norm(4) x depool(3) x fs(5) x ngenes(4) x cutoff(3) x split(4) = 2,880
  (cv folds are averaged within a spec; s50/s60/s70 average their 3 repeats;
   cutoff does not affect C-index, only HR/log-rank p - reported as such)
External validation (split=full, model trained on the whole cohort):
    same 3x5x4x3 signatures x the cancer's own external queue(s):
    HNSC -> geo1 (GSE41613) + geo2 (GSE65858); BRCA -> geo3 (GSE20685);
    LUAD -> geo4 (GSE31210). Each is compared against the SAME spec's
    internal cross-validated C-index (never an in-sample full-train value,
    which does not exist for the full split) -> attenuation delta.

Outputs (results/analysis/):
  specs_internal_<stem>.csv      one row per spec (all splits)
  spec_curve_<stem>.csv          sorted curve
  attenuation_<stem>.csv         internal(CV) vs external, per queue
  variance_decomp.json           3 methods + cross-cancer consistency
  random_control_summary.csv     real vs random C-index by ngenes
"""
from __future__ import annotations

import json
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
KNOBS_NOSPLIT = ["norm", "depool", "fs", "ngenes", "cutoff"]


def load_chunks(stem, prefix="stageB"):
    d = SRC / f"{prefix}_{stem}"
    files = sorted(d.glob("*.csv")) if d.exists() else []
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def agg(df, group, extval="internal"):
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
        rows.append(rec)
    return pd.DataFrame(rows)


# ------------------------------------------------------------ variance decomposition
def decomp_anova(spec, knobs=KNOBS):
    spec = spec.dropna(subset=["cindex_mean"])
    grand = spec["cindex_mean"].mean()
    ss_tot = float(((spec["cindex_mean"] - grand) ** 2).sum())
    out = {}
    for k in knobs:
        ss_b = 0.0
        for _, grp in spec.groupby(k):
            ss_b += len(grp) * (grp["cindex_mean"].mean() - grand) ** 2
        out[k] = float(ss_b / ss_tot) if ss_tot > 0 else 0.0
    return out


def decomp_swing(spec, knobs=KNOBS):
    spec = spec.dropna(subset=["cindex_mean"])
    return {k: float(spec.groupby(k)["cindex_mean"].mean().max()
                     - spec.groupby(k)["cindex_mean"].mean().min()) for k in knobs}


def decomp_rf(spec, knobs=KNOBS):
    """drop-one-knob loss in R^2 of a random forest fitted on one-hot knobs."""
    try:
        from sklearn.ensemble import RandomForestRegressor
    except Exception:  # noqa: BLE001
        return {}
    spec = spec.dropna(subset=["cindex_mean"])
    if spec.empty:
        return {}
    X = pd.get_dummies(spec[knobs].astype(str), columns=knobs)
    yv = spec["cindex_mean"].to_numpy()
    f = RandomForestRegressor(n_estimators=400, random_state=7, min_samples_leaf=2)
    f.fit(X, yv)
    r2 = f.score(X, yv)
    out = {}
    for k in KNOBS:
        cols = [c for c in X.columns if not c.startswith(k + "_")]
        m = RandomForestRegressor(n_estimators=400, random_state=7, min_samples_leaf=2)
        m.fit(X[cols], yv)
        out[k] = float(max(0.0, r2 - m.score(X[cols], yv)))
    s = sum(out.values()) or 1.0
    return {k: v / s for k, v in out.items()}


def main():
    report = {"generated": pd.Timestamp.now().isoformat()}
    internal, attenuations = {}, []
    for proj, stem in CANCERS.items():
        df = load_chunks(stem)
        if df.empty:
            report[proj] = {"status": "no data"}
            continue
        # internal spec space over the split levels that actually hold out data
        # ("full" trains on everything and has no internal validation set)
        spec = agg(df[df["split"] != "full"], KNOBS)
        spec.insert(0, "cancer", proj)
        spec.to_csv(OUT / f"specs_internal_{stem}.csv", index=False)
        internal[proj] = spec
        sc = spec.sort_values("cindex_mean").reset_index(drop=True)
        sc.to_csv(OUT / f"spec_curve_{stem}.csv", index=False)

        # ---- external attenuation: external queue vs the SAME spec's internal
        # cross-validated performance (the full-train units have no held-out
        # set, so an in-sample baseline does not exist and must not be used)
        internal_cv = agg(df[df["split"] == "cv"], KNOBS_NOSPLIT, "internal")
        frames = []
        for q in sorted(x for x in df["extval"].unique() if x != "internal"):
            g = agg(df[df["split"] == "full"], KNOBS_NOSPLIT, q)
            if g.empty or internal_cv.empty:
                continue
            g = g[g["cutoff"] == "median"]
            base = internal_cv[internal_cv["cutoff"] == "median"]
            att = g.merge(base[KNOBS_NOSPLIT + ["cindex_mean"]], on=KNOBS_NOSPLIT,
                          how="inner", suffixes=("", "_int"))
            if att.empty:
                continue
            att = att.rename(columns={"cindex_mean": "external_c",
                                      "cindex_mean_int": "internal_cv_c"})
            att["delta"] = att["external_c"] - att["internal_cv_c"]
            att.insert(0, "queue", q)
            att.insert(0, "cancer", proj)
            frames.append(att)
        if frames:
            a = pd.concat(frames)
            a.to_csv(OUT / f"attenuation_{stem}.csv", index=False)
            attenuations.append(a)

        vd = {"n_specs": int(len(spec)),
              "anova_eta2": decomp_anova(spec),
              "marginal_swing": decomp_swing(spec),
              "rf_drop_r2share": decomp_rf(spec),
              "overall": {
                  "cindex_min": float(spec["cindex_mean"].min()),
                  "cindex_max": float(spec["cindex_mean"].max()),
                  "range": float(spec["cindex_mean"].max() - spec["cindex_mean"].min()),
                  "frac_ge_060": float((spec["cindex_mean"] >= 0.60).mean()),
                  "frac_sig_gt50": float((spec["sig_rate"] > 0.5).mean()),
              },
              "n_events": None}
        # ---- two-panel decomposition -------------------------------------
        # Panel A (above) mixes two different kinds of freedom: choices that
        # BUILD the signature (norm/depool/fs/ngenes) and choices that define
        # how it is VALIDATED (split/cutoff - the latter changes the training
        # set size and nothing else). Panel B holds the validation scheme at
        # cross-validation so the construction hierarchy is not masked by it.
        cv_only = spec[spec["split"] == "cv"]
        if not cv_only.empty:
            vd["panelB_within_cv"] = {
                "n_specs": int(len(cv_only)),
                "anova_eta2": decomp_anova(cv_only, KNOBS_NOSPLIT),
                "marginal_swing": decomp_swing(cv_only, KNOBS_NOSPLIT),
                "rf_drop_r2share": decomp_rf(cv_only, KNOBS_NOSPLIT),
                "overall": {
                    "cindex_min": float(cv_only["cindex_mean"].min()),
                    "cindex_max": float(cv_only["cindex_mean"].max()),
                    "range": float(cv_only["cindex_mean"].max()
                                   - cv_only["cindex_mean"].min()),
                },
            }
        # validation-design effect: marginal means per split level
        vd["split_marginals"] = {
            k: float(v) for k, v in
            spec.groupby("split")["cindex_mean"].mean().items()}
        report[proj] = vd

    # random-signature control
    rc_all = []
    for proj, stem in CANCERS.items():
        rc = load_chunks(stem, prefix="stageC")
        if rc.empty:
            continue
        rc = rc[rc["split"].isin(["cv", "s50", "s60", "s70"])]
        g = (rc.groupby(["depool", "ngenes"])["cindex"]
               .agg(["mean", "std", "count"]).reset_index())
        g.insert(0, "cancer", proj)
        rc_all.append(g)
    if rc_all:
        rcs = pd.concat(rc_all)
        rcs.to_csv(OUT / "random_control_summary.csv", index=False)
        report["random_control"] = rcs.to_dict(orient="records")

    # cross-cancer framework: which knob matters most, consistently?
    if len(internal) >= 2:
        sw = pd.DataFrame({p: decomp_swing(s) for p, s in internal.items()})
        sw.to_csv(OUT / "knob_swing_by_cancer.csv")
        from scipy import stats
        corr = {}
        cols = list(sw.columns)
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                try:
                    rho, pv = stats.spearmanr(sw[cols[i]], sw[cols[j]])
                    corr[f"{cols[i]}~{cols[j]}"] = {"rho": float(rho), "p": float(pv)}
                except Exception:  # noqa: BLE001
                    pass
        report["cross_cancer_swing_rank_corr"] = corr
        report["knob_swing_by_cancer"] = sw.to_dict()
    if attenuations:
        att = pd.concat(attenuations)
        att.to_csv(OUT / "attenuation_all.csv", index=False)
        report["attenuation_overall"] = {
            "n_rows": int(len(att)),
            "median_delta": float(att["delta"].median()),
            "frac_decayed": float((att["delta"] < 0).mean()),
            "external_median": float(att["external_c"].median()),
            "frac_external_ge_060": float((att["external_c"] >= 0.60).mean()),
        }
    (OUT / "variance_decomp.json").write_text(
        json.dumps(report, indent=1), encoding="utf-8")
    print(json.dumps(report, indent=1)[:5000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
