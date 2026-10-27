# -*- coding: utf-8 -*-
"""P4 figures: specification curve, variance decomposition, external attenuation.

Outputs (results/figures/):
  fig1_specification_curve_<cancer>.pdf/png   spec curve + spec matrix
  fig2_variance_decomposition.pdf/png           knob contribution (3 methods)
  fig3_external_attenuation.pdf/png             internal vs GEO scatter by spec
  fig4_cross_cancer_swing.pdf/png               swing magnitude per cancer

All figures are generated from the aggregated analysis tables only.
"""
from __future__ import annotations

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
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
CANCERS = [("TCGA-HNSC", "hnsc"), ("TCGA-BRCA", "brca"), ("TCGA-LUAD", "luad")]
KNOBS = ["norm", "depool", "fs", "ngenes", "cutoff"]
KNOB_LABEL = {"norm": "Normalization", "depool": "Gene pool", "fs": "Feature selection",
              "ngenes": "# genes", "cutoff": "Risk cutoff",
              "split": "Validation design"}
C_NS = {"log2tpm1": "log2(TPM+1)", "voom": "voom", "deseq2vst": "DESeq2-VST",
        "quantile": "quantile", "limma": "limma", "edger": "edgeR", "deseq2": "DESeq2",
        "lasso": "LASSO", "ridge": "Ridge", "enet": "ElasticNet", "rsf": "RSF",
        "unicox": "univCox", "median": "median", "bestlogrank": "best-logrank",
        "tertile": "tertile"}
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "savefig.bbox": "tight"})


def fig_spec_curve(stem, proj):
    f = ANA / f"spec_curve_{stem}.csv"
    if not f.exists():
        return False
    sc = pd.read_csv(f).sort_values("cindex_mean").reset_index(drop=True)
    if sc.empty:
        return False
    fig, (ax, axm) = plt.subplots(2, 1, figsize=(7.2, 5.2),
                                   gridspec_kw={"height_ratios": [3, 1.5]})
    x = np.arange(len(sc))
    ax.plot(x, sc["cindex_mean"], ".", ms=3, color="#1f4e79", alpha=.75)
    # rolling median
    if len(sc) > 20:
        w = max(5, len(sc) // 20)
        roll = sc["cindex_mean"].rolling(w, center=True, min_periods=1).median()
        ax.plot(x, roll, "-", lw=1.4, color="#c00000")
    ax.axhline(0.5, ls="--", lw=.9, color="gray")
    ax.axhline(0.6, ls=":", lw=.9, color="green")
    ax.set_ylabel("Validation C-index")
    ax.set_title(f"Specification curve — {proj}  "
                 f"(n specs = {len(sc)}, range = "
                 f"{sc['cindex_mean'].max()-sc['cindex_mean'].min():.3f})")
    ax.set_xlabel("Specification (ranked by performance)")
    # spec matrix
    order = {k: i for i, k in enumerate(KNOBS)}
    mat = []
    for k in KNOBS:
        u = sc[k].astype(str).unique().tolist()
        codes = {v: j for j, v in enumerate(u)}
        mat.append([codes[v] for v in sc[k].astype(str)])
    mat = np.array(mat)
    im = axm.imshow(mat, aspect="auto", cmap="tab10", interpolation="none",
                    vmin=0, vmax=9)
    # Level labels are drawn as a colour-keyed legend BELOW the band, not on top
    # of it: the levels are interleaved in the ranked order, so per-level text
    # at their mean position collapses into an unreadable pile.
    for j, k in enumerate(KNOBS):
        vals = sc[k].astype(str).unique().tolist()
        x0 = 0.055 + j * 0.19
        fig.text(x0, 0.175, KNOB_LABEL[k], fontsize=6.5, weight="bold")
        for m, v in enumerate(vals):
            col = plt.get_cmap("tab10")(m)
            fig.text(x0 + 0.012, 0.145 - m * 0.028, C_NS.get(v, v),
                     fontsize=5.8, color=col,
                     fontweight="bold" if m == 0 else "normal")
            fig.add_artist(plt.Rectangle((x0, 0.1405 - m * 0.028), 0.008, 0.012,
                                         transform=fig.transFigure,
                                         facecolor=col, edgecolor="none"))
    axm.set_yticks(range(len(KNOBS)))
    axm.set_yticklabels([KNOB_LABEL[k] for k in KNOBS], fontsize=7)
    axm.set_xticks([])
    for s in ("top", "right", "bottom"):
        axm.spines[s].set_visible(False)
    fig.subplots_adjust(bottom=0.30)
    fig.savefig(FIG / f"fig1_specification_curve_{stem}.pdf")
    fig.savefig(FIG / f"fig1_specification_curve_{stem}.png", dpi=200)
    plt.close(fig)
    return True


def fig_spec_curve_merged():
    """Merged Figure 1: three specification-curve panels side by side.

    A = TCGA-HNSC, B = TCGA-BRCA, C = TCGA-LUAD. Each panel shows the same
    content as the per-cancer figures (ranked scatter of all 2,880
    specifications, red centred rolling median, 0.5 and 0.60 reference
    lines, and a lower strip with the level of each decision point), so a
    reader can compare the three cancers directly.

    This is a re-layout only: it reads the same spec_curve_<stem>.csv that
    fig_spec_curve() uses and recomputes no number. The per-cancer files
    are left in place and are not modified.
    """
    panels = []
    for proj, stem in CANCERS:
        f = ANA / f"spec_curve_{stem}.csv"
        if not f.exists():
            continue
        sc = pd.read_csv(f).sort_values("cindex_mean").reset_index(drop=True)
        if sc.empty:
            continue
        panels.append((proj, stem, sc))
    if not panels:
        return False

    # 2 rows x N columns: row 0 = the curve, row 1 = the decision-level strip
    fig, axes = plt.subplots(2, len(panels), figsize=(7.0 * len(panels), 5.2),
                             squeeze=False,
                             gridspec_kw={"height_ratios": [3, 1.5]})
    for pi, (proj, stem, sc) in enumerate(panels):
        ax = axes[0][pi]
        axm = axes[1][pi]
        letter = chr(ord("A") + pi)

        x = np.arange(len(sc))
        ax.plot(x, sc["cindex_mean"], ".", ms=2.2, color="#1f4e79", alpha=.75)
        # centred rolling median
        if len(sc) > 20:
            w = max(5, len(sc) // 20)
            roll = sc["cindex_mean"].rolling(w, center=True, min_periods=1).median()
            ax.plot(x, roll, "-", lw=1.2, color="#c00000")
        ax.axhline(0.5, ls="--", lw=.9, color="gray")
        ax.axhline(0.6, ls=":", lw=.9, color="green")
        ax.set_ylabel("Validation C-index")
        ax.set_title(f"{letter}  {proj}  (n specs = {len(sc)}, range = "
                     f"{sc['cindex_mean'].max() - sc['cindex_mean'].min():.3f})",
                     fontsize=9.5)
        ax.set_xlabel("Specification (ranked by performance)")

        # lower strip: level of each decision point for the spec above it
        mat = []
        for k in KNOBS:
            u = sc[k].astype(str).unique().tolist()
            codes = {v: j for j, v in enumerate(u)}
            mat.append([codes[v] for v in sc[k].astype(str)])
        mat = np.array(mat)
        axm.imshow(mat, aspect="auto", cmap="tab10", interpolation="none",
                   vmin=0, vmax=9)
        axm.set_yticks(range(len(KNOBS)))
        axm.set_yticklabels([KNOB_LABEL[k] for k in KNOBS], fontsize=7)
        axm.set_xticks([])
        for s in ("top", "right", "bottom"):
            axm.spines[s].set_visible(False)

    # colour key for the decision-level strip, drawn once at figure level
    # (the per-cancer figures draw it per figure; here the three panels
    # share the same five decision points, so one key serves all three)
    key_levels = []
    for k in KNOBS:
        vals = panels[0][2][k].astype(str).unique().tolist()
        key_levels.append((k, vals))
    n_cols = len(KNOBS)
    cell_w = 1.0 / n_cols
    for j, (k, vals) in enumerate(key_levels):
        x0 = 0.02 + j * cell_w
        fig.text(x0, 0.055, KNOB_LABEL[k], fontsize=6.5, weight="bold")
        for m, v in enumerate(vals):
            col = plt.get_cmap("tab10")(m)
            fig.text(x0 + 0.008, 0.040 - m * 0.021, C_NS.get(v, v),
                     fontsize=5.6, color=col,
                     fontweight="bold" if m == 0 else "normal")
            fig.add_artist(plt.Rectangle(
                (x0 + 0.001, 0.0355 - m * 0.021), 0.004, 0.009,
                transform=fig.transFigure, facecolor=col, edgecolor="none"))

    fig.suptitle("Figure 1  Specification curve of prognostic signature "
                 "performance in three cohorts", fontsize=11)
    fig.tight_layout(rect=(0, 0.075, 1, 0.96))
    fig.savefig(FIG / "fig1_specification_curve.pdf")
    fig.savefig(FIG / "fig1_specification_curve.png", dpi=200)
    plt.close(fig)
    return True


def fig_variance():
    p = ANA / "variance_decomp.json"
    if not p.exists():
        return False
    rep = json.loads(p.read_text(encoding="utf-8"))
    methods = [("anova_eta2", "ANOVA η²"), ("marginal_swing", "Marginal swing"),
               ("rf_drop_r2share", "RF drop R² share")]
    cands = [c for c, _ in CANCERS if isinstance(rep.get(c), dict)
             and "anova_eta2" in rep.get(c, {})]
    if not cands:
        return False
    fig, axes = plt.subplots(1, len(methods), figsize=(9.6, 3.4), sharey=True)
    w = 0.8 / len(cands)
    for ax, (m, title) in zip(axes, methods):
        for i, c in enumerate(cands):
            v = [rep[c].get(m, {}).get(k, 0) for k in KNOBS]
            ax.bar(np.arange(len(KNOBS)) + i * w, v, width=w, label=c.replace("TCGA-", ""))
        ax.set_xticks(np.arange(len(KNOBS)) + 0.4 - w / 2)
        ax.set_xticklabels([KNOB_LABEL[k] for k in KNOBS], rotation=30, ha="right", fontsize=7)
        ax.set_title(title, fontsize=9)
    axes[0].set_ylabel("Variance share")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Which analytic choice moves the result? (variance decomposition)",
                 fontsize=10)
    fig.savefig(FIG / "fig2_variance_decomposition.pdf")
    fig.savefig(FIG / "fig2_variance_decomposition.png", dpi=200)
    plt.close(fig)
    return True


def fig_variance_panels():
    """Panel A (all decisions) vs Panel B (construction decisions only, within CV).

    Panel A mixes two kinds of freedom: how the signature is BUILT
    (norm/depool/fs/ngenes) and how it is VALIDATED (split, which changes the
    training-set size). Panel B holds validation fixed at CV so the
    construction hierarchy is visible.
    """
    p = ANA / "variance_decomp.json"
    if not p.exists():
        return False
    rep = json.loads(p.read_text(encoding="utf-8"))
    cands = [c for c, _ in CANCERS if isinstance(rep.get(c), dict)
             and "panelB_within_cv" in rep.get(c, {})]
    if not cands:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4))
    knobsB = ["norm", "depool", "fs", "ngenes"]
    for ax, (key, knobs, title) in zip(
            axes,
            [("anova_eta2", KNOBS, "Panel A — all 7 decisions"),
             (None, knobsB, "Panel B — construction only (CV fixed)")]):
        w = 0.8 / len(cands)
        for i, c in enumerate(cands):
            if key is None:
                vals = [rep[c]["panelB_within_cv"]["anova_eta2"].get(k, 0) for k in knobs]
            else:
                vals = [rep[c].get(key, {}).get(k, 0) for k in knobs]
            ax.bar(np.arange(len(knobs)) + i * w, vals, width=w,
                   label=c.replace("TCGA-", ""))
        ax.set_xticks(np.arange(len(knobs)) + 0.4 - w / 2)
        ax.set_xticklabels([KNOB_LABEL[k] for k in knobs], rotation=30,
                           ha="right", fontsize=7)
        ax.set_title(title, fontsize=9)
    axes[0].set_ylabel("Variance share (ANOVA η²)")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Separating construction choices from validation design",
                 fontsize=10)
    fig.savefig(FIG / "fig5_variance_two_panels.pdf")
    fig.savefig(FIG / "fig5_variance_two_panels.png", dpi=200)
    plt.close(fig)
    return True


def fig_random_control():
    """Random-signature null (pool-matched) against the real specification space.

    The null diamond for each (cancer, ngenes) is the mean over the three
    candidate gene pools (limma / edgeR / DESeq2). This is the same
    口径 as the "real minus random" gaps reported in the text, which are
    computed against the pool-averaged null
    (scripts/panelB_report.py: rnd = rc.groupby("ngenes")["mean"].mean()).
    Previously the diamond was taken from r["mean"].iloc[0], i.e. the FIRST
    pool only (DESeq2), which understated the null and therefore overstated
    the gap.
    """
    rc_path = ANA / "random_control_summary.csv"
    if not rc_path.exists():
        return False
    rc = pd.read_csv(rc_path)
    specs = []
    for proj, stem in CANCERS:
        f = ANA / f"spec_curve_{stem}.csv"
        if f.exists():
            s = pd.read_csv(f)
            s = s[s["split"] == "cv"] if "split" in s.columns else s
            s = s.assign(cancer=proj)
            specs.append(s)
    if not specs:
        return False
    specs = pd.concat(specs)
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    cancers = list(specs["cancer"].unique())
    pos = np.arange(len(cancers))
    for i, c in enumerate(cancers):
        d = specs[specs["cancer"] == c]["cindex_mean"]
        ax.boxplot([d], positions=[i], widths=0.5, showfliers=False,
                   boxprops=dict(color="#1f4e79"), medianprops=dict(color="#c00000"))
        for j, ng in enumerate(sorted(rc[rc["cancer"] == c]["ngenes"].unique())):
            r = rc[(rc["cancer"] == c) & (rc["ngenes"] == ng)]
            if r.empty:
                continue
            # pool-averaged null: mean of `mean` across the three depool rows
            ax.plot([i - 0.18 + 0.12 * j], [float(r["mean"].mean())], "D",
                    ms=4, color="#e07b00")
    ax.axhline(0.5, ls="--", lw=.8, color="gray")
    ax.set_xticks(pos)
    ax.set_xticklabels([c.replace("TCGA-", "") for c in cancers])
    ax.set_ylabel("Validation C-index (CV)")
    ax.set_title("Real specification space (box) vs random signatures drawn from\n"
                 "the same candidate gene pools (orange diamonds, sizes 5/10/20/50;\n"
                 "pool-averaged over limma, edgeR and DESeq2)", fontsize=9)
    fig.savefig(FIG / "fig6_random_control.pdf")
    fig.savefig(FIG / "fig6_random_control.png", dpi=200)
    plt.close(fig)
    return True


def fig_attenuation():
    """Internal (cross-validated) vs external C-index, one panel per
    cancer x external-cohort combination.

    There are FOUR such combinations, not two:
        HNSC -> geo1 (GSE41613) + geo2 (GSE65858)
        BRCA -> geo3 (GSE20685)
        LUAD -> geo4 (GSE31210)
    The previous version hard-coded a 1x2 grid and zipped it with
    sorted(queues), so the BRCA and LUAD panels were silently dropped even
    though the text cites LUAD external results throughout. The panel layout
    is now derived from the data (one axes per observed queue), so every
    cancer x cohort combination present in the results is shown.
    """
    frames = []
    for proj, stem in CANCERS:
        f = ANA / f"attenuation_{stem}.csv"
        if f.exists():
            frames.append(pd.read_csv(f))
    if not frames:
        return False
    att = pd.concat(frames, ignore_index=True)
    xcol = "internal_cv_c" if "internal_cv_c" in att.columns else "internal_full_c"

    # queue -> (discovery cancer, GEO accession) follows the mapping the grid
    # runner used (run_grid.py GEO_BY_CANCER).
    QUEUE_INFO = {
        "geo1": "GSE41613",
        "geo2": "GSE65858",
        "geo3": "GSE20685",
        "geo4": "GSE31210",
    }
    queues = sorted(att["queue"].unique().tolist())
    n_panels = len(queues)
    ncols = 2 if n_panels > 1 else 1
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.9 * nrows),
                             squeeze=False)

    for idx, q in enumerate(queues):
        ax = axes[idx // ncols][idx % ncols]
        d = att[att["queue"] == q]
        canc = d["cancer"].unique().tolist()
        gse = QUEUE_INFO.get(q, q)
        for ci, proj in enumerate(canc):
            dd = d[d["cancer"] == proj]
            ax.scatter(dd[xcol], dd["external_c"], s=11, alpha=.65,
                       color=["#1f4e79", "#c0504d", "#e07b00"][ci % 3],
                       label=proj.replace("TCGA-", ""), edgecolors="none")
        lim = [min(att[xcol].min(), att["external_c"].min()),
               max(att[xcol].max(), att["external_c"].max())]
        ax.plot(lim, lim, "--", lw=.8, color="gray")
        ax.set_xlabel("Internal (cross-validated) C-index", fontsize=8)
        ax.set_ylabel(f"External C-index ({gse})", fontsize=8)
        sub = ", ".join(f"{c.replace('TCGA-', '')} (n={len(d[d['cancer'] == c])})"
                        for c in canc)
        ax.set_title(f"{q} \u2014 {gse}\n{sub}", fontsize=8.5)
        ax.tick_params(labelsize=7)
        if len(canc) == 1:
            ax.legend(fontsize=7, frameon=False, loc="upper left")

    # blank out any unused axes in the grid
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Internal discrimination does not transfer "
                 "(all four cancer x external-cohort pairs)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG / "fig3_external_attenuation.pdf")
    fig.savefig(FIG / "fig3_external_attenuation.png", dpi=200)
    plt.close(fig)
    return True


def fig_cross_cancer():
    f = ANA / "knob_swing_by_cancer.csv"
    if not f.exists():
        return False
    sw = pd.read_csv(f, index_col=0)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    w = 0.8 / sw.shape[1]
    for i, c in enumerate(sw.columns):
        ax.bar(np.arange(len(sw)) + i * w, sw[c].to_numpy(), width=w,
               label=c.replace("TCGA-", ""))
    ax.set_xticks(np.arange(len(sw)) + 0.4 - w / 2)
    ax.set_xticklabels([KNOB_LABEL.get(k, k) for k in sw.index], rotation=25,
                       ha="right")
    ax.set_ylabel("Marginal swing in C-index")
    ax.set_title("Magnitude of each analytic choice, by cancer", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    fig.savefig(FIG / "fig4_cross_cancer_swing.pdf")
    fig.savefig(FIG / "fig4_cross_cancer_swing.png", dpi=200)
    plt.close(fig)
    return True


def main():
    made = []
    for proj, stem in CANCERS:
        if fig_spec_curve(stem, proj):
            made.append(f"fig1_{stem}")
    if fig_spec_curve_merged():
        made.append("fig1_merged")
    if fig_variance():
        made.append("fig2")
    if fig_attenuation():
        made.append("fig3")
    if fig_cross_cancer():
        made.append("fig4")
    if fig_variance_panels():
        made.append("fig5")
    if fig_random_control():
        made.append("fig6")
    print("figures:", made if made else "none (run aggregate_analyze.py first)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
