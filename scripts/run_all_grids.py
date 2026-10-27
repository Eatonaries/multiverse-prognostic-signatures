# -*- coding: utf-8 -*-
"""run_all_grids.py - sequential multiverse grid across cancers with watchdog.

For each cancer: repeatedly launch run_grid.py (resumable; skips finished
units) until all 80 units have both stageA .rds and stageB .csv, or the
restart budget is exhausted. Guards memory (>=4 GB free before launch).

Usage: python run_all_grids.py [TCGA-HNSC TCGA-BRCA TCGA-LUAD ...]
"""
import ctypes
import subprocess
import sys
import time
from datetime import datetime
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs" / "run_all_grids.log"
NORMS = ["log2tpm1", "voom", "deseq2vst", "quantile"]
SPLITS = {"cv": 10, "s50": 3, "s60": 3, "s70": 3, "full": 1}
TOTAL_UNITS = len(list(product(NORMS, [(s, i) for s, r in SPLITS.items()
                                       for i in range(1, r + 1)])))
MAX_ROUNDS = 30  # restarts per cancer


def log(m):
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {m}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def free_gb():
    class GS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    gs = GS()
    gs.dwLength = ctypes.sizeof(gs)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(gs))
    return gs.ullAvailPhys / 1024 ** 3


def counts(project):
    stem = project.lower().replace("tcga-", "")
    a = ROOT / "results" / f"stageA_{stem}"
    b = ROOT / "results" / f"stageB_{stem}"
    return (len(list(a.glob("*.rds"))) if a.exists() else 0,
            len(list(b.glob("*.csv"))) if b.exists() else 0)


def commit_free_gb():
    """Free commit charge (commit limit - committed) in GB."""
    class GS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    gs = GS(); gs.dwLength = ctypes.sizeof(gs)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(gs))
    return gs.ullAvailPageFile / 1024 ** 3


def pick_workers():
    """Adaptive parallelism: R units peak at 2-3 GB and the commit ceiling on
    this box is ~42 GB. 4 workers drove commit to 96% and starved each process
    (CPU per worker fell from ~65% to ~17%), so only use 4 when there is real
    headroom."""
    cf = commit_free_gb()
    n = 4 if cf > 10 else 3
    log(f"workers={n} (free commit {cf:.1f} GB)")
    return n


def run_cancer(project):
    stall = 0
    for rnd in range(1, MAX_ROUNDS + 1):
        a, b = counts(project)
        if a >= TOTAL_UNITS and b >= TOTAL_UNITS:
            log(f"{project}: COMPLETE ({a}/{TOTAL_UNITS})")
            return True
        while free_gb() < 4:
            log(f"{project}: memory {free_gb():.1f} GB < 4, waiting")
            time.sleep(120)
        log(f"{project} round {rnd}: stageA={a}/{TOTAL_UNITS} "
            f"stageB={b}/{TOTAL_UNITS} launching run_grid (free {free_gb():.1f} GB)")
        t0 = time.time()
        subprocess.run([sys.executable, str(ROOT / "scripts" / "run_grid.py"),
                        project, str(pick_workers())])
        a2, b2 = counts(project)
        log(f"{project} round {rnd}: ended wall={(time.time()-t0)/60:.0f}min "
            f"stageA={a2}/{TOTAL_UNITS} stageB={b2}/{TOTAL_UNITS}")
        if a2 == a and b2 == b:
            stall += 1
            log(f"{project}: no progress in round {rnd} (stall {stall}/3)")
        else:
            stall = 0
        if stall >= 3:
            log(f"{project}: 3 consecutive stalled rounds -> abort cancer")
            return False
    a, b = counts(project)
    return a >= TOTAL_UNITS and b >= TOTAL_UNITS


def main():
    projects = sys.argv[1:] or ["TCGA-HNSC", "TCGA-BRCA", "TCGA-LUAD"]
    results = {}
    for p in projects:
        results[p] = run_cancer(p)
    log("SUMMARY: " + ", ".join(f"{k}={'OK' if v else 'INCOMPLETE'}"
                                for k, v in results.items()))


if __name__ == "__main__":
    sys.exit(main())
