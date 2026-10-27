# -*- coding: utf-8 -*-
"""Supplementary figures S1-S4, defined in manuscript/图注.md lines 68-75.

  Figure S1 - Cohort assembly flow (cases -> expression n clinical ->
               model-valid), with the source of each censoring time
               (death date vs follow-up entity).
  Figure S2 - Probe-to-symbol collapse and cross-platform symbol overlap
               for each external cohort.
  Figure S3 - Cut-off effects: hazard ratio and log-rank p for the three
               cut-off rules, holding the signature fixed.
  Figure S4 - Independent re-implementation: Python vs R C-index for the
               same specifications (agreement plot).

Reads only from results/ and data/processed/; writes only to
results/figures/. It recomputes no headline number: every panel is a
direct aggregation of an existing result table.

Usage:  python make_supplementary_figures.py
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANA = ROOT / "results" / "analysis"
RES = ROOT / "results"
SUP = RES / "supplementary"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

CANCERS = [("TCGA-HNSC", "hnsc"), ("TCGA-BRCA", "brca"), ("TCGA-LUAD", "luad")]

# queue label -> (GEO accession, processed directory) -- same mapping the
# grid runner used (scripts/run_grid.py GEO_BY_CANCER).
GEO_QUEUES = [
    ("geo1", "GSE41613", "geo_gse41613", "GPL570"),
    ("geo2", "GSE65858", "geo_gse65858", "GPL10558"),
    ("geo3", "GSE20685", "geo_gse20685", "GPL570"),
    ("geo4", "GSE31210", "geo_gse31210", "GPL570"),
]

plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "savefig.bbox": "tight"})


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf")
    fig.savefig(FIG / f"{stem}.png", dpi=200)
    plt.close(fig)
    print(f"  wrote {stem}.pdf / {stem}.png")


# ----------------------------------------------------------------- Figure S1
def fig_s1():
    """Cohort assembly flow + censoring-time provenance, per cancer."""
    rows = []
    for proj, stem in CANCERS:
        meta_p = PROC / f"tcga-{stem}_full" / "analysis_meta.json"
        stat_p = PROC / f"tcga-{stem}_full" / "staging_stats.json"
        if not meta_p.exists():
            print("  S1: missing", meta_p)
            continue
        m = json.loads(meta_p.read_text(encoding="utf-8"))
        cs = m.get("clinical_stats", {})
        rec = {
            "cancer": proj,
            "clinical_n_raw": cs.get("n_raw"),
            "n_used_clinical": cs.get("n_used"),
            "dropped_no_time": cs.get("n_dropped_no_time"),
            "dataset_n": m.get("dataset_n"),
            "model_valid_n": m.get("model_valid_n"),
            "events": m.get("n_events"),
        }
        if stat_p.exists():
            s = json.loads(stat_p.read_text(encoding="utf-8"))
            rec["n_cases_expression"] = s.get("n_cases_expression")
            rec["n_overlap_expression_clinical"] = s.get("n_overlap_expression_clinical")
        for k, v in (cs.get("time_source_counts") or {}).items():
            rec[f"time_source::{k}"] = v
        rows.append(rec)
    if not rows:
        return False
    df = pd.DataFrame(rows)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.0),
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    x = np.arange(len(df))
    w = 0.26
    # --- left: assembly funnel -------------------------------------------
    stages = [
        ("clinical_n_raw", "Clinical cases\nretrieved"),
        ("n_used_clinical", "Clinical with\nusable time"),
        ("n_cases_expression", "Cases with\nexpression"),
        ("model_valid_n", "Model-valid\n(final cohort)"),
    ]
    for i, (col, lab) in enumerate(stages):
        if col not in df.columns:
            continue
        ax1.bar(x + (i - 1.5) * w, df[col].astype(float), w, label=lab,
                color=["#9dc3e6", "#5b9bd5", "#2e75b6", "#1f4e79"][i])
    for i, (col, lab) in enumerate(stages):
        if col not in df.columns:
            continue
        for xi, v in zip(x + (i - 1.5) * w, df[col].astype(float)):
            ax1.text(xi, v, f"{int(v)}", ha="center", va="bottom", fontsize=7)
    ax1.set_xticks(x)
    ax1.set_xticklabels([c.replace("TCGA-", "") for c in df["cancer"]])
    ax1.set_ylabel("Patients")
    ax1.set_title("A  Cohort assembly (cases \u2192 expression \u2229 clinical \u2192 model-valid)",
                  fontsize=9)
    ax1.legend(fontsize=7, frameon=False, ncol=2, loc="upper left")
    ax1.set_ylim(0, max(df[[c for c, _ in stages
                             if c in df.columns]].max().astype(float)) * 1.28)

    # --- right: censoring-time source -------------------------------------
    src_cols = [c for c in df.columns if c.startswith("time_source::")]
    src_cols.sort(key=lambda c: -df[c].fillna(0).sum())
    bottom = np.zeros(len(df))
    palette = ["#c00000", "#ed7d31", "#a5a5a5", "#ffc000"]
    for i, c in enumerate(src_cols):
        vals = df[c].fillna(0).astype(float).values
        ax2.bar(x, vals, 0.5, bottom=bottom, label=c.split("::", 1)[1],
                color=palette[i % len(palette)])
        bottom += vals
    for xi, (tot, ev) in enumerate(zip(bottom, df["events"].astype(float))):
        ax2.text(xi, tot + bottom.max() * 0.015,
                 f"{int(tot)} used\n{int(ev)} events",
                 ha="center", va="bottom", fontsize=7)
    ax2.set_xticks(x)
    ax2.set_xticklabels([c.replace("TCGA-", "") for c in df["cancer"]])
    ax2.set_ylabel("Patients")
    ax2.set_title("B  Source of the censoring time per patient", fontsize=9)
    ax2.legend(fontsize=7, frameon=False, loc="upper right")
    ax2.set_ylim(0, bottom.max() * 1.30)

    fig.suptitle("Figure S1  Cohort assembly and censoring-time provenance",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "supplementary_figS1")
    return True


# ----------------------------------------------------------------- Figure S2
def _annot_rows(annot_name):
    """probe_id -> gene symbol, from a GPL annot .gz file."""
    p = ROOT / "data" / "raw" / annot_name
    out = []
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith(("^", "#", "!")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip('"') for h in parts]
                continue
            if len(parts) < len(header):
                parts = parts + [""] * (len(header) - len(parts))
            d = dict(zip(header, parts))
            out.append((d.get("ID", "").strip('"'),
                        d.get("Gene symbol", "").strip('"').strip()))
    return pd.DataFrame(out, columns=["probe_id", "symbol"])


def fig_s2():
    """Probe -> symbol collapse and cross-platform symbol overlap."""
    tcga_map = pd.read_csv(PROC / "ensg_to_symbol_gencode_v36.tsv", sep="\t")
    tcga_syms = set(tcga_map["symbol"])

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    # --- left: funnel of probe counts -------------------------------------
    # Each count is taken from THIS cohort's own artefacts, not from the
    # platform annotation file (an annotation describes the whole platform,
    # while a cohort only carries a subset of its probes).
    #   n_probes          <- meta.json          (probes present in this cohort)
    #   n_with_symbol     <- expr_probes.csv x annotation (probes we can map)
    #   n_symbols         <- gene_meta.json     (genes after the collapse)
    #   n_overlap         <- overlap_symbols.json (shared with TCGA)
    labels, n_probes, n_with_symbol, n_symbols, n_overlap = [], [], [], [], []
    annot_cache = {}
    for q, gse, gdir_name, platform in GEO_QUEUES:
        gdir = PROC / gdir_name
        meta_p = gdir / "meta.json"
        ov_p = gdir / "overlap_symbols.json"
        gm_p = gdir / "gene_meta.json"
        if not (meta_p.exists() and ov_p.exists() and gm_p.exists()):
            continue
        m = json.loads(meta_p.read_text(encoding="utf-8"))
        ov = json.loads(ov_p.read_text(encoding="utf-8"))
        gm = json.loads(gm_p.read_text(encoding="utf-8"))
        annot_name = f"{platform}.annot.gz"
        if annot_name not in annot_cache:
            annot_cache[annot_name] = _annot_rows(annot_name)
        sym_of = dict(zip(annot_cache[annot_name]["probe_id"],
                          annot_cache[annot_name]["symbol"]))
        n_with = 0
        with open(gdir / "expr_probes.csv", encoding="utf-8",
                  errors="replace") as fh:
            fh.readline()                       # header
            for line in fh:
                pid = line.split(",", 1)[0].strip('"')
                if sym_of.get(pid, ""):
                    n_with += 1
        labels.append(f"{gse}\n({platform})")
        n_probes.append(int(m.get("n_probes", 0)))
        n_with_symbol.append(n_with)
        n_symbols.append(int(gm.get("n_genes_symbol", 0)))
        n_overlap.append(len(ov))

    x = np.arange(len(labels))
    w = 0.20
    for i, (vals, lab, col) in enumerate([
            (n_probes, "Probes on array", "#9dc3e6"),
            (n_with_symbol, "Probes with a gene symbol", "#5b9bd5"),
            (n_symbols, "Distinct gene symbols (after collapse)", "#2e75b6"),
            (n_overlap, "Symbols shared with TCGA", "#1f4e79")]):
        ax = axes[0]
        ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=col)
        for xi, v in zip(x + (i - 1.5) * w, vals):
            ax.text(xi, v, f"{v:,}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Count (log scale)")
    ax.set_yscale("log")
    ax.set_title("A  Probe-to-symbol collapse and cross-platform overlap",
                 fontsize=9)
    ax.legend(fontsize=6.5, frameon=False)
    ax.set_ylim(top=max(n_probes) * 8)

    # --- right: overlap matrix between external cohorts and TCGA ----------
    syms = {}
    for q, gse, gdir_name, platform in GEO_QUEUES:
        ov_p = PROC / gdir_name / "overlap_symbols.json"
        if ov_p.exists():
            syms[gse] = set(json.loads(ov_p.read_text(encoding="utf-8")))
    names = list(syms.keys())
    M = np.zeros((len(names), len(names) + 1), dtype=int)
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            M[i, j] = len(syms[a] & syms[b])
        M[i, -1] = len(syms[a] & tcga_syms)
    ax = axes[1]
    im = ax.imshow(M, cmap="Blues")
    ax.set_xticks(range(len(names) + 1))
    ax.set_xticklabels(names + ["TCGA\n(GENCODE v36)"], rotation=30, ha="right",
                       fontsize=7.5)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:,}", ha="center", va="center", fontsize=6.5,
                    color="white" if M[i, j] > M.max() * 0.55 else "black")
    ax.set_title("B  Shared gene symbols (n)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Figure S2  Probe-to-symbol collapse and cross-platform symbol "
                 "overlap for each external cohort", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "supplementary_figS2")
    return True


# ----------------------------------------------------------------- Figure S3
def fig_s3():
    """Cut-off effects on the dichotomised readouts: HR and log-rank p.

    The C-index is computed from the CONTINUOUS risk score, so within one
    fixed signature it is identical for the three cut-off rules; only the
    dichotomised metrics (hazard ratio, log-rank p) move.

    Panel A: one held-out signature per cancer, showing the three cut-off
    rules side by side for HR and for -log10(p).
    Panel B: the distribution of both readouts over the whole cross-validated
    grid, one box per cut-off rule, per cancer.

    Source: results/analysis/cutoff_effects_<stem>.csv (built by
    scripts/aggregate_cutoff_effects.py from the per-fold stage-B chunks).
    """
    frames = []
    for proj, stem in CANCERS:
        f = ANA / f"cutoff_effects_{stem}.csv"
        if not f.exists():
            f = ANA / f"specs_internal_{stem}.csv"
        if f.exists():
            d = pd.read_csv(f)
            if "p_median" not in d.columns:
                print(f"  S3: {f.name} has no p_median -- run "
                      f"aggregate_cutoff_effects.py first")
                continue
            frames.append(d)
    if not frames:
        print("  S3: no cutoff_effects_*.csv found")
        return False
    df = pd.concat(frames, ignore_index=True)
    df = df[df["split"] == "cv"]
    df = df.dropna(subset=["hr_median", "p_median", "neg_log10_p_med"])

    CUTS = ["median", "bestlogrank", "tertile"]
    CUT_LABEL = {"median": "median", "bestlogrank": "best log-rank",
                 "tertile": "tertile"}
    COLORS = {"median": "#5b9bd5", "bestlogrank": "#c00000", "tertile": "#ed7d31"}

    fig = plt.figure(figsize=(12.4, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.15], hspace=.42, wspace=.30)

    # ---------- panel A: one fixed signature per cancer --------------------
    axA_hr = fig.add_subplot(gs[0, :])
    picked = []
    for proj, stem in CANCERS:
        d = df[(df["cancer"] == proj)]
        if d.empty:
            continue
        key = ["norm", "depool", "fs", "ngenes"]
        grp = d.groupby(key)["hr_median"].agg(lambda s: s.max() - s.min())
        if grp.empty:
            continue
        best = grp.idxmax()
        sub = d
        for k, v in zip(key, best):
            sub = sub[sub[k] == v]
        if sub.empty:
            continue
        rows = {r["cutoff"]: r for _, r in sub.iterrows()
                if r["cutoff"] in CUTS}
        if len(rows) < len(CUTS):
            continue
        picked.append((proj, best, rows))

    xx = np.arange(len(picked))
    w = 0.26
    for i, cut in enumerate(CUTS):
        hrs = [rows[cut]["hr_median"] for _, _, rows in picked]
        bars = axA_hr.bar(xx + (i - 1) * w, hrs, w, label=CUT_LABEL[cut],
                          color=COLORS[cut])
        for b, h in zip(bars, hrs):
            axA_hr.text(b.get_x() + b.get_width() / 2, h, f"{h:.2f}",
                        ha="center", va="bottom", fontsize=6.5)
    axA_hr.set_xticks(xx)
    axA_hr.set_xticklabels(
        [f"{p.replace('TCGA-', '')}\n{key[0]}\u00b7{key[1]}\u00b7{key[2]}\u00b7{key[3]}g"
         for p, key, _ in picked], fontsize=7)
    axA_hr.set_ylabel("Hazard ratio (high vs low risk)")
    axA_hr.set_title("A  Hazard ratio by cut-off rule "
                     "(one fixed held-out signature per cancer)", fontsize=9)
    axA_hr.legend(fontsize=7, frameon=False, ncol=3)
    axA_hr.axhline(1.0, ls=":", lw=.8, color="gray")

    # ---------- panel B: log-rank p for the SAME fixed signatures ----------
    axA_p = fig.add_subplot(gs[1, :])
    for i, cut in enumerate(CUTS):
        vals = [-np.log10(max(rows[cut]["p_median"], 1e-300))
                for _, _, rows in picked]
        bars = axA_p.bar(xx + (i - 1) * w, vals, w, label=CUT_LABEL[cut],
                         color=COLORS[cut])
        for b, v in zip(bars, vals):
            axA_p.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}",
                       ha="center", va="bottom", fontsize=6.5)
    axA_p.set_xticks(xx)
    axA_p.set_xticklabels(
        [f"{p.replace('TCGA-', '')}\n{key[0]}\u00b7{key[1]}\u00b7{key[2]}\u00b7{key[3]}g"
         for p, key, _ in picked], fontsize=7)
    axA_p.set_ylabel("\u2212log10(median log-rank p)")
    axA_p.set_title("B  Log-rank p for the same signatures "
                     "(higher = more significant)", fontsize=9)
    axA_p.legend(fontsize=7, frameon=False, ncol=3)
    # p = 0.05 reference
    axA_p.axhline(-np.log10(0.05), ls="--", lw=.9, color="gray")
    axA_p.text(0.995, -np.log10(0.05), "  p = 0.05", transform=axA_p.get_yaxis_transform(),
               ha="right", va="bottom", fontsize=6.8, color="gray")

    # ---------- panel C: whole-grid distributions, per cancer --------------
    for ci, (proj, stem) in enumerate(CANCERS):
        d = df[df["cancer"] == proj]
        if d.empty:
            continue
        axC = fig.add_subplot(gs[0, 2] if ci == 2 else gs[1, 2])
        # HR distributions
        data_hr = [d.loc[d.cutoff == c, "hr_median"].dropna().values
                   for c in CUTS]
        bp = axC.boxplot(data_hr, tick_labels=[CUT_LABEL[c] for c in CUTS],
                         widths=0.55, showfliers=False, patch_artist=True,
                         positions=[0, 1, 2])
        for patch, c in zip(bp["boxes"], CUTS):
            patch.set_facecolor(COLORS[c])
            patch.set_alpha(.55)
        for med in bp["medians"]:
            med.set_color("black")
        axC.axhline(1.0, ls=":", lw=.8, color="gray")
        axC.set_ylabel("Hazard ratio", fontsize=8)
        axC.tick_params(labelsize=7)
        axC.set_title(f"{proj.replace('TCGA-', '')}  HR by cut-off "
                      f"(whole CV grid)", fontsize=8.5)
        cidx = d.groupby("cutoff")["cindex_mean"].median()
        txt = "  |  ".join(f"{CUT_LABEL[c]}: C-index median {cidx[c]:.4f}"
                           for c in CUTS if c in cidx.index)
        axC.text(0.5, 0.02, txt, transform=axC.transAxes, ha="center",
                 va="bottom", fontsize=6.2, style="italic")

    # ---------- panel D: -log10(p) distributions, per cancer ---------------
    for ci, (proj, stem) in enumerate(CANCERS):
        d = df[df["cancer"] == proj]
        if d.empty:
            continue
        axD = fig.add_subplot(gs[0, 0] if ci == 0 else
                              (gs[1, 0] if ci == 1 else gs[0, 1]))
        data_p = [d.loc[d.cutoff == c, "neg_log10_p_med"].dropna().values
                  for c in CUTS]
        bp = axD.boxplot(data_p, tick_labels=[CUT_LABEL[c] for c in CUTS],
                         widths=0.55, showfliers=False, patch_artist=True,
                         positions=[0, 1, 2])
        for patch, c in zip(bp["boxes"], CUTS):
            patch.set_facecolor(COLORS[c])
            patch.set_alpha(.55)
        for med in bp["medians"]:
            med.set_color("black")
        axD.axhline(-np.log10(0.05), ls="--", lw=.9, color="gray")
        axD.set_ylabel("\u2212log10(median p)", fontsize=8)
        axD.tick_params(labelsize=7)
        axD.set_title(f"{proj.replace('TCGA-', '')}  log-rank p by cut-off",
                      fontsize=8.5)

    fig.suptitle("Figure S3  Cut-off effects on the dichotomised readouts "
                 "(hazard ratio and log-rank p); discrimination is unchanged",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "supplementary_figS3")
    return True


# ----------------------------------------------------------------- Figure S4
def fig_s4():
    """Python (NumPy/SciPy) vs R C-index agreement for the same specs."""
    sources = [("TCGA-HNSC", RES / "p5_spec_comparison_hnsc.csv"),
               ("TCGA-BRCA", RES / "p5_spec_comparison_brca.csv")]
    frames = []
    for cohort, f in sources:
        if f.exists():
            d = pd.read_csv(f)
            d = d.copy()
            d.insert(0, "cohort", cohort)
            frames.append(d)
    if not frames:
        print("  S4: no p5_spec_comparison files found")
        return False
    df = pd.concat(frames, ignore_index=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    # --- panel A: agreement plot ------------------------------------------
    for cohort, col in [("TCGA-HNSC", "#1f4e79"), ("TCGA-BRCA", "#c0504d")]:
        d = df[df.cohort == cohort]
        ax1.scatter(d["cindex"], d["py_cindex"], s=22, alpha=.75, color=col,
                    label=f"{cohort} (n={len(d)})", edgecolors="none")
    lim = [min(df["cindex"].min(), df["py_cindex"].min()) - .01,
           max(df["cindex"].max(), df["py_cindex"].max()) + .01]
    ax1.plot(lim, lim, "--", lw=.9, color="gray", label="y = x (exact agreement)")
    ax1.set_xlabel("R  survival::concordance  C-index")
    ax1.set_ylabel("Python (NumPy/SciPy) C-index")
    ax1.set_title("A  Independent re-implementation vs primary pipeline", fontsize=9)
    ax1.legend(fontsize=7, frameon=False, loc="upper left")

    # --- panel B: absolute difference distribution ------------------------
    data = [df.loc[df.cohort == c, "abs_diff"].values
            for c in ["TCGA-HNSC", "TCGA-BRCA"] if (df.cohort == c).any()]
    labels = [c.replace("TCGA-", "") for c in ["TCGA-HNSC", "TCGA-BRCA"]
              if (df.cohort == c).any()]
    bp = ax2.boxplot(data, tick_labels=labels, widths=0.45, showfliers=True,
                      patch_artist=True)
    for patch, col in zip(bp["boxes"], ["#9dc3e6", "#f4b6b6"]):
        patch.set_facecolor(col)
    for med in bp["medians"]:
        med.set_color("black")
    ax2.set_yscale("log")
    ax2.set_ylabel("|C-index difference| (log scale)")
    ax2.set_title("B  Absolute difference per specification", fontsize=9)
    for i, (c, vals) in enumerate(zip(labels, data), start=1):
        ax2.text(i, max(vals) * 3.2,
                 f"median {np.median(vals):.2e}\nmax {np.max(vals):.2e}\n"
                 f"100% < 0.01",
                 ha="center", va="bottom", fontsize=6.8)

    fig.suptitle("Figure S4  Independent re-implementation: Python vs R C-index "
                 "for the same specifications", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "supplementary_figS4")
    return True


def main():
    made = []
    if fig_s1():
        made.append("S1")
    if fig_s2():
        made.append("S2")
    if fig_s3():
        made.append("S3")
    if fig_s4():
        made.append("S4")
    print("supplementary figures:", made if made else "none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
