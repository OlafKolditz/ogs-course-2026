#!/usr/bin/env python3
"""
run_ogs_all.py  -  run all (or selected) OGS projects of the STIMTEC DFN folder

Each project <name>.prj writes its output into its own subfolder
    results/<name>/        (*.pvd, *.vtu and <name>.log)
via the OGS command line option -o. The prj files stay unchanged; the .pvd files
keep plain file names, so ParaView and compare_dfn_injection.py read them directly.

Examples
    python run_ogs_all.py                          # all *.prj in this folder
    python run_ogs_all.py -p stimtec_dfn_40_kS_t   # one project (name or part of it)
    python run_ogs_all.py -p 40 -j 4               # all 40.6 projects, 4 in parallel
    python run_ogs_all.py --clean                  # delete old results/<name>/ first
    python run_ogs_all.py --ogs C:/ogs/bin/ogs.exe # other OGS executable (or env OGS)
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))


def run(prj, results, ogs, clean):
    name = os.path.splitext(os.path.basename(prj))[0]
    out = os.path.join(results, name)
    if clean and os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)
    log = os.path.join(out, name + ".log")
    t0 = time.time()
    with open(log, "w") as fh:
        r = subprocess.run([ogs, prj, "-o", out], cwd=os.path.dirname(prj),
                           stdout=fh, stderr=subprocess.STDOUT)
    n = len(glob.glob(os.path.join(out, "*.vtu")))
    return name, r.returncode, time.time() - t0, n, log


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-d", "--dir", default=HERE, help="folder with the *.prj files")
    ap.add_argument("-r", "--results", default=None, help="results folder (default: <dir>/results)")
    ap.add_argument("-p", "--projects", nargs="*", help="project names or substrings (default: all)")
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="number of OGS runs in parallel")
    ap.add_argument("--ogs", default=os.environ.get("OGS", "ogs"), help="OGS executable")
    ap.add_argument("--clean", action="store_true", help="delete results/<name>/ before the run")
    a = ap.parse_args()

    prjs = sorted(glob.glob(os.path.join(os.path.abspath(a.dir), "*.prj")))
    if a.projects:
        prjs = [p for p in prjs
                if any(w == os.path.basename(p)[:-4] or w.lower() in os.path.basename(p).lower()
                       for w in a.projects)]
    if not prjs:
        sys.exit("no matching *.prj found")
    if shutil.which(a.ogs) is None and not os.path.isfile(a.ogs):
        sys.exit(f"OGS executable not found: {a.ogs}  (use --ogs or set OGS)")
    results = os.path.abspath(a.results or os.path.join(a.dir, "results"))
    print(f"{len(prjs)} project(s) -> {results}  ({a.jobs} parallel)")

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run, p, results, a.ogs, a.clean) for p in prjs]
        failed = 0
        for f in futs:
            name, rc, dt, n, log = f.result()
            ok = "ok    " if rc == 0 else "FAILED"
            failed += rc != 0
            print(f"  [{ok}] {name:28s} {dt:6.1f} s  {n:4d} vtu   log: {os.path.relpath(log, a.dir)}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
