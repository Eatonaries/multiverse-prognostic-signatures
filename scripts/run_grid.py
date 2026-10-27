# -*- coding: utf-8 -*-
"""run_grid.py: orchestrate the two-stage multiverse grid.

Stage A (expensive, shared): one (cancer x norm x split-fold) unit ->
  RDS with train/test matrices, gene pools (depool levels), fs rankings.
Stage B (cheap): per unit -> score all depool x fs x ngenes x cutoff specs
  (+ external queue from the "full" units) into a CSV chunk.

Resumable: units with existing outputs are skipped.

Usage: python run_grid.py <TCGA-HNSC|TCGA-BRCA|TCGA-LUAD> [workers]
"""
from __future__ import annotations

import json
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
NORMS = ["log2tpm1", "voom", "deseq2vst", "quantile"]
SPLITS = {"cv": 10, "s50": 3, "s60": 3, "s70": 3, "full": 1}

# External cohorts are matched to the DISCOVERY cancer type: evaluating a
# BRCA-derived signature in a head-and-neck cohort measures cross-cancer
# transfer, not transportability within the disease. Each entry maps the
# queue label used in the results to its processed GEO directory.
PROC = ROOT / "data" / "processed"
GEO_BY_CANCER = {
    "hnsc": {"geo1": PROC / "geo_gse41613",   # OSCC, GPL570, OS + cancer-specific
             "geo2": PROC / "geo_gse65858"},  # HNSCC, GPL10558, OS + PFS
    "brca": {"geo3": PROC / "geo_gse20685"},  # BRCA, GPL570, OS (83 events)
    "luad": {"geo4": PROC / "geo_gse31210"},  # LUAD, GPL570, OS + RFS
}

_lock = threading.Lock()


def log(msg):
    with _lock:
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
        print(line, flush=True)
        with open(ROOT / "logs" / "run_grid.log", "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def stem_of(project):
    return project.lower().replace("tcga-", "")


def ensure_staged(project):
    """stage_expression + build_analysis_set + make_slim if missing."""
    cdir = ROOT / "data" / "processed" / f"{project.lower()}_full"
    if (cdir / "slim.rds").exists():
        return
    log(f"staging {project} ...")
    subprocess.run([sys.executable, str(ROOT / "scripts" / "pipeline_core.py"),
                    "--project", project, "--tag", "full"],
                   check=True, capture_output=True, text=True, timeout=7200)
    subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "scripts" / "make_slim.R"),
                    project, "full"], check=True, capture_output=True, text=True,
                   timeout=1800)
    log(f"staged {project}")


def unit_done(project, norm, split, fold):
    stem = stem_of(project)
    rds = ROOT / "results" / f"stageA_{stem}" / f"{split}__{norm}__f{fold:02d}.rds"
    csv = ROOT / "results" / f"stageB_{stem}" / f"{split}__{norm}__f{fold:02d}.csv"
    return rds.exists() and csv.exists()


def free_gb():
    """Available physical memory in GB (Windows)."""
    import ctypes

    class GS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    gs = GS(); gs.dwLength = ctypes.sizeof(gs)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(gs))
    return gs.ullAvailPhys / 1024 ** 3


def wait_for_memory(min_gb=3.0, limit_s=3600):
    """Never start a unit when RAM is tight: R units peak at 2-3 GB each."""
    t0 = time.time()
    while free_gb() < min_gb and time.time() - t0 < limit_s:
        log(f"  memory guard: {free_gb():.1f} GB free (<{min_gb}) - waiting")
        time.sleep(60)


def run_unit(project, norm, split, fold):
    stem = stem_of(project)
    adir = ROOT / "results" / f"stageA_{stem}"
    bdir = ROOT / "results" / f"stageB_{stem}"
    adir.mkdir(parents=True, exist_ok=True)
    bdir.mkdir(parents=True, exist_ok=True)
    rds = adir / f"{split}__{norm}__f{fold:02d}.rds"
    csv = bdir / f"{split}__{norm}__f{fold:02d}.csv"
    seed = 42
    if not rds.exists():
        wait_for_memory()
        payload = {"project": project, "tag": "full", "norm": norm,
                   "kind": split, "fold": fold, "seed": seed,
                   "outdir": str(adir).replace("\\", "/")}
        pin = ROOT / "checkpoints" / f"_u_{stem}_{split}_{norm}_{fold}.json"
        pin.write_text(json.dumps(payload), encoding="utf-8")
        t0 = time.time()
        p = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "scripts" / "stageA_grid2.R"),
                            str(pin)], capture_output=True, text=True, timeout=21600)
        if p.returncode != 0 or not rds.exists():
            log(f"FAIL stageA {stem}|{norm}|{split}f{fold} rc={p.returncode} "
                f"{(p.stderr or p.stdout)[-300:]}")
            return None
        log(f"stageA {stem}|{norm}|{split}f{fold} done ({time.time()-t0:.0f}s)")
    payload = {"stageA_rds": str(rds).replace("\\", "/"), "project": project,
               "tag": "full", "out_csv": str(csv).replace("\\", "/"),
               "geo_dirs": {k: str(p).replace("\\", "/")
                            for k, p in GEO_BY_CANCER.get(stem, {}).items()
                            if p.exists()}}
    pin = ROOT / "checkpoints" / f"_ub_{stem}_{split}_{norm}_{fold}.json"
    pin.write_text(json.dumps(payload), encoding="utf-8")
    p = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "scripts" / "stageB_grid.R"),
                        str(pin)], capture_output=True, text=True, timeout=7200)
    if p.returncode != 0 or not csv.exists():
        log(f"FAIL stageB {stem}|{norm}|{split}f{fold} rc={p.returncode} "
            f"{(p.stderr or p.stdout)[-300:]}")
        return None
    return csv


def main():
    project = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    ensure_staged(project)
    units = [(n, s, f) for n in NORMS for s, reps in SPLITS.items() for f in range(1, reps + 1)]
    todo = [u for u in units if not unit_done(project, *u)]
    log(f"{project}: units total={len(units)} remaining={len(todo)} workers={workers}")
    t0 = time.time()
    n_done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_unit, project, *u): u for u in todo}
        for fu in as_completed(futs):
            r = fu.result()
            n_done += 1
            if n_done % 5 == 0 or r is None:
                el = time.time() - t0
                log(f"{project}: {n_done}/{len(todo)} units ({el/60:.0f} min elapsed)")
    log(f"{project}: ALL DONE {n_done} units in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    sys.exit(main())
