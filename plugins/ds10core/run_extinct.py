#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.galactic_ext`` with the shared core on PYTHONPATH (Desktop shell).

OGFinder does not import the core. Output mode is ``set``: the tool writes ``--out`` and this driver prints that
catalogue on stdout so the table is replaced. Flags the plugin manifest passes through unchanged (listed here so
that tools/validate_manifests.py can check the manifest against this driver):  --core-dir --ebv --rv --filter
--wavelength-nm --catalog --mag-column --out --summary
"""
from __future__ import annotations

import sys

from launch import call_core


def main(argv: list[str]) -> int:
    return call_core("ds10core.tools.galactic_ext", argv, echo=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
