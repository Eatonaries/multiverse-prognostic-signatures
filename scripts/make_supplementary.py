# -*- coding: utf-8 -*-
"""Build the supplementary tables from the finished grid.

Run after run_analysis_pipeline.py:
    python make_supplementary.py

Writes to results/supplementary/:
  Table_S1_clinical_provenance.csv   per-cancer censoring-time source counts (R2)
  Table_S2_specification_space.csv   the full design matrix (levels + counts)
  Table_S3_all_specs_<stem>.csv      every internal spec with C-index/HR/p
  Table_S4_external_validation.csv   internal(full-train) vs geo1/geo2 per spec
  Table_S5_random_control.csv        random-signature reference distribution
  Table_S6_independent_check.csv     Python-vs-R agreement (P5)
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results"
ANA = SRC / "analysis"
SUP = SRC / "supplementary"
SUP.mkdir(parents=True, exist_ok=True)
CANCERS = {"TCGA-HNSC": "hnsc", "TCGA-BRCA": "brca", "TCGA-LUAD": "luad"}
KNOBS = ["norm", "depool", "fs", "ngenes", "cutoff", "split"]


def t1_clinical():
    rows = []
    for proj, stem in CANCERS.items():
        meta = ROOT / "data" / "processed" / f"tcga-{stem}_full" / "analysis_meta.json"
        stats = ROOT / "data" / "processed" / f"tcga-{stem}_full" / "staging_stats.json"
        rec = {"cancer": proj}
        if meta.exists():
            m = json.loads(meta.read_text(encoding="utf-8"))
            cs = m.get("clinical_stats", {})
            rec.update({
                "dataset_n": m.get("dataset_n", m.get("n_samples")),
                "model_valid_n": m.get("model_valid_n", m.get("n_samples")),
                "events": m.get("n_events"),
                "clinical_n_raw": cs.get("n_raw"),
                "dropped_no_time": cs.get("n_dropped_no_time"),
                "n_used_clinical": cs.get("n_used"),
            })
            for k, v in (cs.get("time_source_counts") or {}).items():
                rec[f"time_source::{k}"] = v
        if stats.exists():
            s = json.loads(stats.read_text(encoding="utf-8"))
            rec["n_files_used"] = s.get("n_files_used")
            rec["n_cases_expression"] = s.get("n_cases_expression")
            rec["n_overlap_expression_clinical"] = s.get("n_overlap_expression_clinical")
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(SUP / "Table_S1_clinical_provenance.csv", index=False)
    return len(df)


def t2_space():
    levels = {
        "norm": ["log2tpm1", "voom", "deseq2vst", "quantile"],
        "depool": ["limma", "edger", "deseq2"],
        "fs": ["lasso", "ridge", "enet", "rsf", "unicox"],
        "ngenes": [5, 10, 20, 50],
        "cutoff": ["median", "bestlogrank", "tertile"],
        "split": ["cv", "s50", "s60", "s70"],
    }
    rows = []
    n = 1
    for k, v in levels.items():
        n *= len(v)
        rows.append({"decision_point": k, "n_levels": len(v),
                     "levels": " | ".join(map(str, v))})
    rows.append({"decision_point": "TOTAL internal specs per cancer",
                 "n_levels": n, "levels": f"{n} = product of the above"})
    rows.append({"decision_point": "cancers", "n_levels": len(CANCERS),
                 "levels": " | ".join(CANCERS)})
    rows.append({"decision_point": "TOTAL internal specs (all cancers)",
                 "n_levels": n * len(CANCERS), "levels": ""})
    df = pd.DataFrame(rows)
    df.to_csv(SUP / "Table_S2_specification_space.csv", index=False)
    return len(df)


def t3_4_5():
    made = []
    for proj, stem in CANCERS.items():
        f = ANA / f"specs_internal_{stem}.csv"
        if f.exists():
            pd.read_csv(f).to_csv(SUP / f"Table_S3_all_specs_{stem}.csv", index=False)
            made.append(f"S3_{stem}")
        a = ANA / f"attenuation_{stem}.csv"
        if a.exists():
            pd.read_csv(a).to_csv(SUP / f"Table_S4_external_{stem}.csv", index=False)
            made.append(f"S4_{stem}")
    rc = ANA / "random_control_summary.csv"
    if rc.exists():
        pd.read_csv(rc).to_csv(SUP / "Table_S5_random_control.csv", index=False)
        made.append("S5")
    made += t6_independent_check()
    return made


# The independent re-implementation (P5) was run in two of the three cohorts
# (HNSC and BRCA), each comparing 3 gene pools x 5 selectors x 4 sizes = 60
# specifications. Each per-cohort file is the authoritative source; the
# generic p5_spec_comparison.csv is only a copy of the HNSC one.
P5_SOURCES = [
    ("TCGA-HNSC", "p5_spec_comparison_hnsc.csv"),
    ("TCGA-BRCA", "p5_spec_comparison_brca.csv"),
]
P5_EXPECTED_FS = {"lasso", "ridge", "enet", "rsf", "unicox"}
P5_EXPECTED_ROWS = 60


def t6_independent_check():
    """Table S6: Python (NumPy/SciPy) vs R C-index for the same specs.

    Hard-fails on an incomplete or stale source instead of silently
    publishing a table that contradicts the manuscript.
    """
    frames = []
    for cohort, fname in P5_SOURCES:
        f = SRC / fname
        if not f.exists():
            print(f"  S6: MISSING SOURCE {fname} -- skipped")
            continue
        d = pd.read_csv(f)
        fs = set(d["fs"].unique())
        if fs != P5_EXPECTED_FS:
            raise SystemExit(
                f"S6 source {fname} has selectors {sorted(fs)}, "
                f"expected {sorted(P5_EXPECTED_FS)} -- refusing to publish a "
                f"table that contradicts the manuscript's '5 selectors'"
            )
        if len(d) != P5_EXPECTED_ROWS:
            raise SystemExit(
                f"S6 source {fname} has {len(d)} rows, expected "
                f"{P5_EXPECTED_ROWS} (3 pools x 5 selectors x 4 sizes)"
            )
        d = d.copy()
        d.insert(0, "cohort", cohort)
        frames.append(d)
    if not frames:
        return []
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(SUP / "Table_S6_independent_check.csv", index=False)
    print(f"  S6 rows={len(out)} cohorts={sorted(out.cohort.unique().tolist())} "
          f"selectors={sorted(out.fs.unique().tolist())}")
    return ["S6"]


def main():
    n1 = t1_clinical()
    n2 = t2_space()
    made = t3_4_5()
    print(f"S1 rows={n1}  S2 rows={n2}  others={made}")
    print("written to", SUP)


if __name__ == "__main__":
    main()
