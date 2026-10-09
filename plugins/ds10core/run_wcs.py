#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.wcs_solve`` with the shared core on PYTHONPATH (Desktop shell).

OGFinder does not import the core. The status line stays on stdout. Flags the plugin manifest passes through
unchanged (listed here so that tools/validate_manifests.py can check the manifest against this driver):
--core-dir --out --summary --cat --survey --ra --dec --scale --pa --mag-limit --radius
"""
from __future__ import annotations

import sys

from launch import call_core


def main(argv: list[str]) -> int:
    return call_core("ds10core.tools.wcs_solve", argv, echo=False)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
