#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.lsbg`` with the shared core on PYTHONPATH (Desktop shell).

Flags of the LSBG tool that the plugin manifest passes through unchanged (listed here so that
tools/validate_manifests.py can check the manifest against this driver; the real definitions are
the argparse options of ``ds10core/tools/lsbg.py``):  --core-dir --out --segmap --mask-out --summary
--detect-thresh --detect-minarea --smooth-fwhm --mu-eff-min --mu-eff-max --r-eff-min --r-eff-max
--ellipticity-max --mag-zeropoint --r --i
"""
from __future__ import annotations
import os, sys
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
os.execvpe(py, [py, "-m", "ds10core.tools.lsbg", *args], env)
