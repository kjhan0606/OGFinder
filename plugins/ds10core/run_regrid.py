#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.regrid`` with the shared core on PYTHONPATH (Desktop shell).

OGFinder does not import the core. The status line stays on stdout. Flags the plugin manifest passes through
unchanged (listed here so that tools/validate_manifests.py can check the manifest against this driver):
--core-dir --ref --out --summary
"""
from __future__ import annotations

import sys

from launch import call_core


def main(argv: list[str]) -> int:
    return call_core("ds10core.tools.regrid", argv, echo=False)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
