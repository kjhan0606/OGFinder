#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.userplugin run`` with the shared core on PYTHONPATH (Desktop shell).

Flags: --core-dir DIR (plugin directory) --work --catalog --image --param --params-json --timeout --mem-mb --summary
"""
from __future__ import annotations
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ds10 import find_core, alias_env
alias_env()
args = sys.argv[1:]
core_dir = None
if len(args) >= 2 and args[0] == "--core-dir":
    core_dir, args = args[1], args[2:]
d, why = find_core(core_dir)
if not d:
    sys.stderr.write(why + "\n"); sys.exit(2)
env = os.environ.copy()
env["PYTHONPATH"] = d + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
py = env.get("DS10_PYTHON") or env.get("ASTRAFEX_PYTHON") or sys.executable
# args: <plugin_dir> --work ... --catalog ...
rc = subprocess.call([py, "-m", "ds10core.tools.userplugin", "run", *args], env=env, stdout=sys.stderr)
# echo catalogue if written
work = args[args.index("--work") + 1] if "--work" in args[:-1] else None
cat = os.path.join(work, "catalog.tsv") if work else None
if rc == 0 and cat and os.path.exists(cat):
    sys.stdout.write(open(cat, encoding="utf-8").read())
sys.exit(rc)
