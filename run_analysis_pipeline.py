# -*- coding: utf-8 -*-
"""After-grid pipeline: stageC (random control) -> aggregate -> figures -> summary.

Run once the cancer grids are complete:
    python run_analysis_pipeline.py [--skip-stageC]

Writes results/analysis/* and results/figures/*, plus
reports/P4_结果摘要.md (auto-generated numbers, no hand-typed numbers).
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANA = ROOT / "results" / "analysis"
CANCERS = {"TCGA-HNSC": "hnsc", "TCGA-BRCA": "brca", "TCGA-LUAD": "luad"}


def run(cmd):
    print(">>", " ".join(str(c) for c in cmd), flush=True)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print((p.stdout or "")[-1500:])
        print((p.stderr or "")[-1500:])
    return p.returncode


def stage_counts(stem, prefix):
    """Count stage-A .rds units and stage-B .csv units respectively.

    The two stages use different file types, so the glob must depend on the
    prefix; globbing "*.csv" for stageA made every cancer look empty and
    silently skipped its random control.
    """
    d = ROOT / "results" / f"{prefix}_{stem}"
    if not d.exists():
        return 0
    pat = "*.rds" if prefix == "stageA" else "*.csv"
    return len(list(d.glob(pat)))


def main():
    skip_c = "--skip-stageC" in sys.argv
    ready = {}
    for proj, stem in CANCERS.items():
        a = stage_counts(stem, "stageA")
        b = stage_counts(stem, "stageB")
        ready[proj] = (a, b)
        print(f"{proj}: stageA={a} stageB={b}")
    if not skip_c:
        for proj, stem in CANCERS.items():
            if ready[proj][0] == 0:
                continue
            run([sys.executable, str(ROOT / "scripts" / "run_stageC.py"), proj, "3", "20"])
    run([sys.executable, str(ROOT / "scripts" / "aggregate_analyze.py")])
    run([sys.executable, str(ROOT / "scripts" / "make_figures.py")])
    write_summary()


def write_summary():
    p = ANA / "variance_decomp.json"
    if not p.exists():
        print("no variance_decomp.json")
        return
    rep = json.loads(p.read_text(encoding="utf-8"))
    lines = [f"# P4 结果摘要（自动生成 {datetime.now():%Y-%m-%d %H:%M}）", ""]
    lines.append("> 本文件所有数字均由 `aggregate_analyze.py` 自动写入，禁止手工誊抄。\n")
    for proj, stem in CANCERS.items():
        v = rep.get(proj)
        if not isinstance(v, dict) or "anova_eta2" not in v:
            lines.append(f"## {proj}\n\n（无数据）\n")
            continue
        o = v["overall"]
        lines.append(f"## {proj}\n")
        lines.append(f"- 规格数：{v['n_specs']}")
        lines.append(f"- C-index 范围：{o['cindex_min']:.4f} – {o['cindex_max']:.4f}"
                     f"（**极差 {o['range']:.4f}**）")
        lines.append(f"- 达到 0.60 的比例：{o['frac_ge_060']:.1%}")
        lines.append(f"- 半数以上折 log-rank p<0.05 的规格比例：{o['frac_sig_gt50']:.1%}")
        lines.append("")
        lines.append("| 旋钮 | ANOVA η² | 边际摆动 | RF drop-R² 份额 |")
        lines.append("|---|---|---|---|")
        for k in ["norm", "depool", "fs", "ngenes", "cutoff", "split"]:
            lines.append(f"| {k} | {v['anova_eta2'].get(k, 0):.3f} | "
                         f"{v['marginal_swing'].get(k, 0):.4f} | "
                         f"{v['rf_drop_r2share'].get(k, 0):.3f} |")
        lines.append("")
    att = rep.get("attenuation_overall")
    if att:
        lines.append("## 外部衰减（全网格）\n")
        lines.append(f"- 规格×队列行数：{att['n_rows']}")
        lines.append(f"- ΔC-index 中位数：**{att['median_delta']:+.4f}**")
        lines.append(f"- 衰减（Δ<0）比例：{att['frac_decayed']:.1%}")
        lines.append(f"- 外部 C-index 中位数：{att['external_median']:.4f}")
        lines.append(f"- 外部达到 0.60 的比例：{att['frac_external_ge_060']:.2%}")
        lines.append("")
    cc = rep.get("cross_cancer_swing_rank_corr")
    if cc:
        lines.append("## 跨癌种一致性（旋钮摆动幅度的秩相关）\n")
        for k, v in cc.items():
            lines.append(f"- {k}: ρ = {v['rho']:.3f} (p = {v['p']:.3f})")
        lines.append("")
    rc = rep.get("random_control")
    if rc:
        lines.append("## 随机签名对照\n")
        lines.append("| 癌种 | 基因池 | 基因数 | 随机 C-index 均值 | sd | n |")
        lines.append("|---|---|---|---|---|---|")
        for r in rc[:24]:
            lines.append(f"| {r['cancer']} | {r['depool']} | {r['ngenes']} | "
                         f"{r['mean']:.4f} | {r['std']:.4f} | {r['count']} |")
        lines.append("")
    out = ROOT / "reports" / "P4_结果摘要.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print("written:", out)
    write_fill_in(rep)


def write_fill_in(rep):
    """Map every manuscript placeholder to its auto-generated value so the
    Results section can be filled in mechanically and audited."""
    rows = ["# 稿件占位符 → 数字 映射表（自动生成，禁止手抄）", ""]
    rows.append("| 稿件位置 | 占位符 | 数值 |")
    rows.append("|---|---|---|")
    att = rep.get("attenuation_overall") or {}
    for proj in CANCERS:
        v = rep.get(proj)
        if not isinstance(v, dict) or "anova_eta2" not in v:
            continue
        o = v["overall"]
        tag = proj.replace("TCGA-", "")
        rows.append(f"| Results 3.1 | `[[{tag}: 规格数]]` | {v['n_specs']} |")
        rows.append(f"| Results 3.1 | `[[{tag}: C-index range]]` | "
                    f"{o['cindex_min']:.4f} – {o['cindex_max']:.4f}（极差 {o['range']:.4f}） |")
        rows.append(f"| Results 3.1 | `[[{tag}: fraction >= 0.60]]` | {o['frac_ge_060']:.1%} |")
        rows.append(f"| Results 3.1 | `[[{tag}: fraction significant]]` | {o['frac_sig_gt50']:.1%} |")
        for k, val in (v.get("anova_eta2") or {}).items():
            rows.append(f"| Results 3.2 | `[[{tag}: anova {k}]]` | {val:.3f} |")
        for k, val in (v.get("rf_drop_r2share") or {}).items():
            rows.append(f"| Results 3.2 | `[[{tag}: rf {k}]]` | {val:.3f} |")
        for k, val in (v.get("marginal_swing") or {}).items():
            rows.append(f"| Results 3.2 | `[[{tag}: swing {k}]]` | {val:.4f} |")
    if att:
        rows.append(f"| Results 3.4 | `[[external rows]]` | {att['n_rows']} |")
        rows.append(f"| Results 3.4 | `[[median delta C-index]]` | {att['median_delta']:+.4f} |")
        rows.append(f"| Results 3.4 | `[[fraction decayed]]` | {att['frac_decayed']:.1%} |")
        rows.append(f"| Results 3.4 | `[[external median]]` | {att['external_median']:.4f} |")
        rows.append(f"| Results 3.4 | `[[fraction external >= 0.60]]` | {att['frac_external_ge_060']:.2%} |")
    cc = rep.get("cross_cancer_swing_rank_corr") or {}
    for k, v in cc.items():
        rows.append(f"| Results 3.3 | `[[rank corr {k}]]` | ρ = {v['rho']:.3f} (p = {v['p']:.3f}) |")
    rc = rep.get("random_control") or []
    for r in rc[:12]:
        rows.append(f"| Results 3.5 | `[[random {r['cancer']} {r['depool']} ng={r['ngenes']}]]` | "
                    f"{r['mean']:.4f} ± {r['std']:.4f} (n={r['count']}) |")
    p5 = ROOT / "results" / "p5_independent_check.json"
    if p5.exists():
        d = json.loads(p5.read_text(encoding="utf-8"))
        l3 = d.get("layer3_spec_comparison", {})
        if "median_abs_diff" in l3:
            rows.append(f"| Methods 2.5 | `[[P5 median abs diff]]` | {l3['median_abs_diff']:.2e} |")
            rows.append(f"| Methods 2.5 | `[[P5 max abs diff]]` | {l3['max_abs_diff']:.4f} |")
            rows.append(f"| Methods 2.5 | `[[P5 n specs]]` | {l3['n_specs']} |")
    out = ROOT / "reports" / "P4_占位符映射表.md"
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print("written:", out)


if __name__ == "__main__":
    sys.exit(main())
