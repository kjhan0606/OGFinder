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

  ds10.py [--core-dir D] [--python P] [--where] COMMAND ARGS...      (COMMAND = info, calib, maps, regions, pixtab, bands, validate, run, forced, compare, export-script, list, open, status)

Words of the form ``split:LIST`` (the list of band images typed into a dialog) are expanded to separate arguments.  LIST is separated by ``|`` or by new
lines when it contains one of them (so file names may contain spaces: ``/my data/a.fits|/my data/b.fits``); otherwise by spaces, where a name with
spaces can be quoted (``"/my data/a.fits" /b.fits``); an unbalanced quote falls back to plain spaces.

``open BUNDLE --dest DIR`` restores an offline bundle of the web app into the folder DIR (see docs/offline_workflow.md); with an empty or missing ``--dest`` the
folder is ``~/ds10-offline/<bundle file name without .zip>``.

Flags of the core command line that the plugin manifest passes through unchanged (listed here so that tools/validate_manifests.py can check the manifest
against this driver; the real definitions are the argparse options of ``ds10core/cli.py``):  --aperture-radii --background --bkg-stat --cog-radius --core-dir --dir --dest --files --force --frame --image --minarea --mode --out --size --smooth-fwhm --text --thresh --to --work --xy --zp
"""
import json
import os
import shlex
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


def split_words(v):
    """The arguments of a ``split:`` word (see the module docstring)."""
    if "|" in v or "\n" in v:
        return [x.strip() for x in v.replace("\n", "|").split("|") if x.strip()]
    try:
        return shlex.split(v)
    except ValueError:                                              # unbalanced quote (an apostrophe in a name): plain spaces
        return v.split()


def default_open_dest(words):
    """``open BUNDLE [--dest D]``: D empty or missing -> ~/ds10-offline/<bundle name>."""
    if not words or words[0] != "open" or len(words) < 2:
        return words
    w = list(words)
    if "--dest" in w:
        i = w.index("--dest")
        if i + 1 < len(w) and w[i + 1].strip() and not w[i + 1].startswith("--"):
            return w
        del w[i:i + (2 if i + 1 < len(w) and not w[i + 1].strip() else 1)]
    stem = os.path.splitext(os.path.basename(w[1]))[0] or "bundle"
    return w + ["--dest", os.path.join(os.path.expanduser("~"), "ds10-offline", stem)]


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
        words += split_words(w[6:]) if w.startswith("split:") else [w]
    words = default_open_dest(words)
    return subprocess.call([py, "-m", "ds10core", *words], env=env)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
