#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.detect`` with the shared core on PYTHONPATH (Desktop shell).

OGFinder does not import the core. Output mode is ``set``: the tool writes ``--out`` and this driver prints that
catalogue on stdout so the table is replaced. Flags the plugin manifest passes through unchanged (listed here so
that tools/validate_manifests.py can check the manifest against this driver):  --core-dir --out --segmap --summary
--thresh --minarea --smooth-fwhm --back-size --mag-zeropoint --max-sources --no-smooth --deblend --deblend-nlevels
--deblend-mincont --no-noise-correction --noise-map --noise-kind --bad-mask --bad-mask-sense --mask
"""
from __future__ import annotations

import sys

from launch import call_core


def main(argv: list[str]) -> int:
    return call_core("ds10core.tools.detect", argv, echo=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
