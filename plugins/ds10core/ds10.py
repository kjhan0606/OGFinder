#!/usr/bin/env python3
"""Stand-alone shell of the shared ds10 compute core (regions, pixel table, map association + calibration, forced photometry, ds10-script / replay runner).

The compute code is NOT part of OGFinder (and not GPL): it is the Python package ``ds10core`` of the ds10-web repository (original code, numpy + scipy only).  This
launcher finds that package and runs ``python -m ds10core ARGS...`` in a SEPARATE PROCESS; OGFinder code never imports it.  Nothing of ds9/AST/SEP is used by the core.

Where the core is looked for (first hit wins; the directory must contain ``ds10core/__init__.py``):
  1. ``--core-dir DIR``                      (the plugin parameter "Core directory")
  2. environment ``DS10_CORE``
  3. ``~/.ds9/ogf_params/ds10core.json``    {"core-dir": "..."}
  4. ``<OGFinder root>/../ds10-web/server``, ``/workspace/ds10-web/server``

The interpreter that runs the core: ``--python PATH``, else env ``DS10_PYTHON``, else the interpreter running this launcher (needs numpy and scipy).

  ds10.py [--core-dir D] [--python P] [--where] COMMAND ARGS...      (COMMAND = info, calib, maps, regions, pixtab, bands, validate, run, forced, compare, export-script, list)

Words of the form ``split:A B C`` (used for the list of band images typed into a dialog) are expanded to the separate arguments A B C.

Flags of the core command line that the plugin manifest passes through unchanged (listed here so that tools/validate_manifests.py can check the manifest
against this driver; the real definitions are the argparse options of ``ds10core/cli.py``):  --aperture-radii --background --bkg-stat --cog-radius --core-dir --dir --files --frame --image --minarea --mode --out --size --smooth-fwhm --text --thresh --to --work --xy --zp
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def candidates(explicit=None):
    if explicit:
        yield explicit, "--core-dir"
    if os.environ.get("DS10_CORE"):
        yield os.environ["DS10_CORE"], "DS10_CORE"
    cfg = os.path.expanduser("~/.ds9/ogf_params/ds10core.json")
    try:
        with open(cfg, encoding="utf-8") as f:
            d = json.load(f).get("core-dir")
        if d:
            yield d, cfg
    except (OSError, ValueError, AttributeError):
        pass
    yield os.path.join(os.path.dirname(ROOT), "ds10-web", "server"), "next to OGFinder"
    yield "/workspace/ds10-web/server", "default location"


def find_core(explicit=None):
    """(directory containing ds10core, where it came from) or (None, message)."""
    tried = []
    for d, why in candidates(explicit):
        d = os.path.abspath(os.path.expanduser(d))
        if os.path.isfile(os.path.join(d, "ds10core", "__init__.py")):
            return d, why
        tried.append(f"{d} ({why})")
    return None, "ds10core not found; looked in: " + "; ".join(tried)


def main(argv):
    explicit = python = None
    where = False
    args = list(argv)
    while args and args[0].startswith("--") and args[0] in ("--core-dir", "--python", "--where"):
        o = args.pop(0)
        if o == "--where":
            where = True
        else:
            if not args:
                print(f"{o} needs a value", file=sys.stderr)
                return 2
            v = args.pop(0)
            if o == "--core-dir":
                explicit = v or None
            else:
                python = v or None
    d, why = find_core(explicit)
    if d is None:
        print(why + "\nSet DS10_CORE to the 'server' directory of ds10-web (or the plugin parameter 'Core directory').", file=sys.stderr)
        return 3
    if where:
        print(json.dumps({"core_dir": d, "found_by": why, "python": python or os.environ.get("DS10_PYTHON") or sys.executable}))
        return 0
    py = python or os.environ.get("DS10_PYTHON") or sys.executable
    env = dict(os.environ)
    env["PYTHONPATH"] = d + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.setdefault("DS10_OGF_ROOT", ROOT)
    words = []
    for w in args:
        words += w[6:].split() if w.startswith("split:") else [w]
    return subprocess.call([py, "-m", "ds10core", *words], env=env)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
