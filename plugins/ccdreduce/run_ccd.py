#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.ccd_reduce`` for the desktop shell.

``--frames`` is a list of FITS paths joined by ``|`` or newlines. Each path becomes one ``--frame``
JSON spec so the classifier sees the original name and the assignment can name the full path.
An empty ``--core-dir`` asks the shared core discovery for the copy next to OGFinder.

Flags the manifest and the dialog pass through, listed here so ``tools/validate_manifests.py`` can see them:
--core-dir --frames --assignment --combine --out --summary --op --frame
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "ds10core"))
from ds10 import alias_env, find_core  # noqa: E402


def frame_args(text: str) -> list[str]:
    """One JSON ``--frame`` value per path. A value that already starts with ``{`` is kept."""
    out = []
    for part in (text or "").replace("\n", "|").split("|"):
        part = part.strip()
        if not part:
            continue
        if part.startswith("{"):
            out.append(part)
        else:
            out.append(json.dumps({"file": part, "name": os.path.basename(part), "id": part}, separators=(",", ":")))
    return out


def main(argv: list[str] | None = None) -> int:
    alias_env()
    args = list(sys.argv[1:] if argv is None else argv)
    core = None
    frames = ""
    rest: list[str] = []
    i = 0
    while i < len(args):
        if args[i] == "--core-dir" and i + 1 < len(args):
            core = args[i + 1] or None
            i += 2
        elif args[i] == "--frames" and i + 1 < len(args):
            frames = args[i + 1]
            i += 2
        else:
            rest.append(args[i])
            i += 1
    if not rest:
        sys.stderr.write("usage: run_ccd.py --core-dir DIR inspect|reduce|op [--assignment A] [--combine median|mean] [--op OP] --out FILE --summary FILE --frames a.fits|b.fits\n")
        return 2
    d, why = find_core(core)
    if not d:
        sys.stderr.write(why + "\n")
        return 2
    cmd, tail = rest[0], rest[1:]
    expanded = [cmd]
    for spec in frame_args(frames):
        expanded += ["--frame", spec]
    expanded += tail
    env = os.environ.copy()
    env["PYTHONPATH"] = d + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    py = env.get("DS10_PYTHON") or env.get("ASTRAFEX_PYTHON") or sys.executable
    return subprocess.call([py, "-m", "ds10core.tools.ccd_reduce", *expanded], env=env)


if __name__ == "__main__":
    sys.exit(main())
