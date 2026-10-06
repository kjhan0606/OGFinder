#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.lsbg`` with the shared core on PYTHONPATH (Desktop shell).

Flags of the LSBG tool that the plugin manifest passes through unchanged (listed here so that
tools/validate_manifests.py can check the manifest against this driver; the real definitions are
the argparse options of ``ds10core/tools/lsbg.py``):  --core-dir --out --segmap --mask-out --summary
--detect-thresh --detect-minarea --smooth-fwhm --mu-eff-min --mu-eff-max --r-eff-min --r-eff-max
--ellipticity-max --mag-zeropoint --r --i --classify --classify-keep-all --classify-threshold
(--classify*: optional RBF-SVM purity stage of the core, model builtin/ds10_lsbg/lsbg_classifier_v1.json)
"""
from __future__ import annotations
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ds10 import find_core, alias_env  # reuse core discovery

alias_env()
args = sys.argv[1:]
core_dir = None
# optional --core-dir DIR at the front
if len(args) >= 2 and args[0] == "--core-dir":
    core_dir = args[1]
    args = args[2:]
d, why = find_core(core_dir)
if not d:
    sys.stderr.write(why + "\n")
    sys.exit(2)
env = os.environ.copy()
env["PYTHONPATH"] = d + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
py = env.get("DS10_PYTHON") or env.get("ASTRAFEX_PYTHON") or sys.executable
# the Desktop step declares output mode "set": the catalogue table is filled from this process's stdout, so run the
# tool (it writes --out and prints nothing on stdout) and then echo the catalogue it wrote
rc = subprocess.call([py, "-m", "ds10core.tools.lsbg", *args], env=env, stdout=sys.stderr)
out = args[args.index("--out") + 1] if "--out" in args[:-1] else None
if rc == 0 and out and os.path.exists(out):
    with open(out, encoding="utf-8") as fh:
        sys.stdout.write(fh.read())
sys.exit(rc)
