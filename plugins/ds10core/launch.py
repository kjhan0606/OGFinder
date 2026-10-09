#!/usr/bin/env python3
"""Start one shared-core tool in a separate process.

OGFinder does not import the core. ``call_core`` finds the Astrafex Web checkout and runs
``python -m <module>`` with that directory on PYTHONPATH. ``echo`` is for steps whose output mode is
``set``: the tool's status line goes to stderr and the ``--out`` file is printed, which is what the
catalogue loader reads.
"""
from __future__ import annotations

import os
import subprocess
import sys


def call_core(module: str, argv: list[str], *, echo: bool = False) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from ds10 import alias_env, find_core

    alias_env()
    args = list(argv)
    core_dir = None
    if len(args) >= 2 and args[0] == "--core-dir":
        core_dir = args[1] or None
        args = args[2:]
    found, why = find_core(core_dir)
    if not found:
        sys.stderr.write(why + "\n")
        return 2
    env = os.environ.copy()
    env["PYTHONPATH"] = found + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    py = env.get("DS10_PYTHON") or env.get("ASTRAFEX_PYTHON") or sys.executable
    rc = subprocess.call([py, "-m", module, *args], env=env, stdout=(sys.stderr if echo else None))
    if echo:
        out = args[args.index("--out") + 1] if "--out" in args[:-1] else None
        if rc == 0 and out and os.path.isfile(out):
            with open(out, encoding="utf-8") as fh:
                sys.stdout.write(fh.read())
    return rc


def drop_empty(argv: list[str], flag: str, *also: str) -> list[str]:
    """Drop ``flag`` and its path when that path is missing or empty, and drop each paired flag too."""
    out = list(argv)
    if flag not in out[:-1]:
        return out
    i = out.index(flag)
    path = out[i + 1]
    if path and not path.startswith("--") and os.path.isfile(path) and os.path.getsize(path) > 0:
        return out
    del out[i:i + 2]
    for name in also:
        if name in out[:-1]:
            j = out.index(name)
            del out[j:j + 2]
    return out
