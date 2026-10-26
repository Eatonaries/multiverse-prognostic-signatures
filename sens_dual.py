# -*- coding: utf-8 -*-
"""P2-3: strict vs lenient dual-version sensitivity analysis (R3).

The multiverse varies the choices INSIDE the specification space; this checks
whether three choices held FIXED in that space change the conclusions:

    setting      strict (= production)      lenient
    expression filter   TPM>1 in >=20%       none
    gene-pool size      2000                 5000
    glmnet lambda       lambda.1se           lambda.min

Runs the cross-validation units under both versions and compares the resulting
specification distributions (range, knob marginals) so the paper can state
whether the conclusion is an artefact of those fixed settings.

Usage: python sens_dual.py <TCGA-XXX> [folds] [workers]
"""
from __future__ import annotations

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
VERSIONS = {
    "strict": {"min_frac": 0.20, "pool_size": 2000, "lambda_rule": "1se"},
    "lenient": {"min_frac": 0.0, "pool_size": 5000, "lambda_rule": "min"},
}
NORMS = ["log2tpm1", "voom"]


def log(m):
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {m}", flush=True)
    with open(ROOT / "logs" / "sens_dual.log", "a", encoding="utf-8") as fh:
        fh.write(m + "\n")


def run_unit(project, stem, version, norm, fold):
    # stage-A writes its RDS named from (split, norm, fold) only, so the two
    # versions must live in separate directories or they would overwrite
    # each other.
    adir = ROOT / "results" / f"stageA_sens_{stem}" / version
    bdir = ROOT / "results" / f"stageB_sens_{stem}"
    adir.mkdir(parents=True, exist_ok=True)
    bdir.mkdir(parents=True, exist_ok=True)
    tag = f"{version}__cv__{norm}__f{fold:02d}"
    rds = adir / f"cv__{norm}__f{fold:02d}.rds"
    csv = bdir / f"{tag}.csv"
    if not rds.exists():
        p = {"project": project, "tag": "full", "norm": norm, "kind": "cv",
             "fold": fold, "seed": 42, "outdir": str(adir).replace("\\", "/"),
             **VERSIONS[version]}
        pin = ROOT / "checkpoints" / f"_sens_{stem}_{tag}.json"
        pin.write_text(json.dumps(p), encoding="utf-8")
        r = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "scripts" / "stageA_sens.R"),
                            str(pin)], capture_output=True, text=True, timeout=21600)
        if r.returncode != 0 or not rds.exists():
            log(f"FAIL stageA {tag}: {(r.stderr or r.stdout)[-200:]}")
            return False
        log(f"stageA {tag} ok")
    if not csv.exists():
        # stageB names its output from the stageA file stem
        pin = ROOT / "checkpoints" / f"_sensb_{stem}_{tag}.json"
        p = {"stageA_rds": str(rds).replace("\\", "/"), "project": project,
             "tag": "full", "out_csv": str(csv).replace("\\", "/"), "geo_dirs": {}}
        pin.write_text(json.dumps(p), encoding="utf-8")
        r = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "scripts" / "stageB_grid.R"),
                            str(pin)], capture_output=True, text=True, timeout=7200)
        if r.returncode != 0 or not csv.exists():
            log(f"FAIL stageB {tag}: {(r.stderr or r.stdout)[-200:]}")
            return False
    return True


def main():
    project = sys.argv[1] if len(sys.argv) > 1 else "TCGA-HNSC"
    folds = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    stem = project.lower().replace("tcga-", "")
    units = [(v, n, f) for v in VERSIONS for n in NORMS for f in range(1, folds + 1)]
    log(f"{project}: {len(units)} sensitivity units (folds={folds})")
    ok = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_unit, project, stem, v, n, f) for v, n, f in units]
        for i, fu in enumerate(as_completed(futs), 1):
            ok += int(bool(fu.result()))
            if i % 5 == 0:
                log(f"{project}: {i}/{len(units)} sensitivity units done")
    log(f"{project}: DONE ok={ok}/{len(units)}")
    subprocess.run([sys.executable, str(ROOT / "scripts" / "sens_compare.py"), stem])


if __name__ == "__main__":
    sys.exit(main())
