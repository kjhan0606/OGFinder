#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.pixmask`` with the shared core on PYTHONPATH (Desktop shell).

OGFinder does not import the core. Flags the plugin manifest passes through unchanged (listed here so that
tools/validate_manifests.py can check the manifest against this driver; the real definitions are the argparse
options of ``ds10core/tools/pixmask.py``):  --core-dir --out --summary --bleed --persist --bleed-axis
--bleed-nsig --bleed-max --persist-nsig --persist-grow --dq-bits --saturate --previous --dq
"""
from __future__ import annotations
import os, subprocess, sys


def main(argv: list[str]) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from ds10 import find_core, alias_env

    alias_env()
    args = list(argv)
    core_dir = None
    if len(args) >= 2 and args[0] == "--core-dir":
        core_dir = args[1]
        args = args[2:]
    d, why = find_core(core_dir)
    if not d:
        sys.stderr.write(why + "\n")
        return 2
    env = os.environ.copy()
    env["PYTHONPATH"] = d + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    py = env.get("DS10_PYTHON") or env.get("ASTRAFEX_PYTHON") or sys.executable
    return subprocess.call([py, "-m", "ds10core.tools.pixmask", *args], env=env)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
