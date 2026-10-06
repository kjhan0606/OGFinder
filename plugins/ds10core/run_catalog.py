#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.catalog`` with the shared core on PYTHONPATH (Desktop shell).

Flags: --core-dir --survey --image --out --reg --summary --max-rows --pad-arcsec --epoch --preset --mag-limit
--columns --index --no-cache --refresh --crop --catalog --ref --radius --prefix --column-style --xmatch-out --out
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
# first remaining arg is the catalog subcommand (query / xmatch / ...)
rc = subprocess.call([py, "-m", "ds10core.tools.catalog", *args], env=env, stdout=sys.stderr)
out = args[args.index("--out") + 1] if "--out" in args[:-1] else (args[args.index("--xmatch-out") + 1] if "--xmatch-out" in args[:-1] else None)
if rc == 0 and out and os.path.exists(out):
    # the catalogue table of the Desktop takes the first line as the header: drop the "# provenance" comment lines of the TSV
    sys.stdout.write("".join(ln for ln in open(out, encoding="utf-8") if not ln.startswith("#")))
sys.exit(rc)
