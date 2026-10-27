# Multiverse / specification-curve analysis of transcriptome-based prognostic signatures

This repository contains the full, fixed analysis pipeline behind the manuscript
*"One dataset, hundreds of prognostic signatures: a specification-curve analysis
of researcher degrees of freedom in transcriptome-based survival modelling."*

The pipeline enumerates a complete, pre-specified specification space of
construction choices for transcriptome-based overall-survival (OS) prognostic
signatures across three TCGA cancers (HNSC, BRCA, LUAD), and scores every
combination on internal (cross-validated) and external (GEO) cohorts.

## What is in here

| Script | Role |
|---|---|
| `stageA_grid2.R` | Stage A (expensive, per cancer×norm×split-fold): normalized matrices, gene pools (limma / edgeR / DESeq2), feature rankings (lasso / ridge / enet / rsf / unicox). Saves an RDS unit. |
| `stageB_grid.R` | Stage B (cheap): scores every pool×selector×signature-size×cutoff spec, internally (CV) and externally (GEO), into a CSV chunk. |
| `run_grid.py` | Single-cancer two-stage orchestrator (resumable; skips finished units). |
| `run_all_grids.py` | Watchdog over all three cancers (serial by cancer, adaptive workers, resumable). |
| `run_analysis_pipeline.py` | Stage-C aggregation → summary tables → figures → manuscript numbers. |
| `aggregate_analyze.py` | Specification curve, two-panel variance decomposition, external attenuation, cross-cancer rank correlations. |
| `aggregate_cutoff_effects.py` | Per-specification cut-off effects: hazard ratio, log-rank p and their aggregates, with the C-index shown to be cut-off-invariant. |
| `make_figures.py` | Figures 1–6 (Figure 1 merged 3-panel; Figure 3 all four cancer × external-cohort pairs). |
| `make_graphical_abstract.py` | Graphical abstract; every number is read from `results/analysis/*`. |
| `make_supplementary.py` | Supplementary Tables S1–S6. |
| `make_supplementary_figures.py` | Supplementary Figures S1–S4. |
| `panelB_report.py` | Within-CV panel + real-vs-random gaps + internal-vs-external (Spearman and Pearson both reported). |
| `make_p4_report.py` | Human-readable results report (filled only from `results/analysis/*`). |
| `sens_dual.py` | Dual-version (strict / lenient) sensitivity on three fixed settings. |
| `stageA_sens.R` | Sensitivity variant of Stage A. |
| `sens_compare.py` | Compares strict vs lenient spec distributions. |
| `build_jbi_package.py` | Assembles the journal submission package (self-contained manuscript + figures) from the markdown sources. |

## Environment

- **R** 4.6.1 with: `survival`, `limma`, `edgeR`, `DESeq2`, `glmnet`,
  `randomForestSRC`, `data.table`, `foreach`/`doParallel`.
- **Python** 3.x with: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`.
- R is invoked via `Rscript`. The path to `Rscript.exe` is hard-coded in
  `run_grid.py` / `run_all_grids.py` (`RSCRIPT`); change it to match your install.

## Paths

The R scripts resolve the project root from the environment variable
`BIO4_ROOT`, falling back to the original author's path:

```r
ROOT <- Sys.getenv("BIO4_ROOT", "D:/deepseek/jieru/_bio4")
```

Set `BIO4_ROOT` to the directory that contains `data/`, `scripts/`, `results/`
before running, e.g.:

```bash
export BIO4_ROOT=/path/to/multiverse-prog-sig
```

The Python scripts derive their root automatically from the script location
(`Path(__file__).resolve().parents[1]`).

## Input data (public, not bundled)

- TCGA HNSC / BRCA / LUAD: transcriptomic (TPM + counts) and clinical (OS time,
  event) retrieved via the GDC API.
- External validation: GSE41613, GSE65858, GSE20685, GSE31210 via GEO.

Place processed inputs under `data/processed/<project>_<tag>/slim.rds` as
consumed by `stageA_grid2.R` / `stageA_sens.R`.

## Reproducibility

- All random seeds are fixed (`set.seed(42)` in the R scripts; `seed=42` passed
  through the Python payloads). Re-running produces identical specification
  curves.
- The specification space is fully enumerated and pre-specified; no
  results-dependent choice was made after seeing the numbers.

## How to run (all three cancers)

```bash
python run_all_grids.py          # Stage A + Stage B for HNSC, BRCA, LUAD
python run_analysis_pipeline.py  # aggregate -> tables -> figures -> report
```

For the strict/lenient sensitivity:

```bash
python sens_dual.py TCGA-HNSC 5 4
python sens_dual.py TCGA-BRCA 5 4
python sens_dual.py TCGA-LUAD  5 4
```

## License

Research code, released for reproducibility. See the manuscript for the data
sources and the full method.
