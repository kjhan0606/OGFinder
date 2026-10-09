#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.counts`` with the shared core on PYTHONPATH (Desktop shell).

OGFinder does not import the core. Output mode is ``set``: the tool writes ``--out`` and this driver prints that
catalogue on stdout so the table is replaced. An empty ``--catalog`` file (the catalogue table is absent) is
dropped together with ``--catalog-out``, so the core bins the image alone. Flags the plugin manifest passes
through unchanged (listed here so that tools/validate_manifests.py can check the manifest against this driver):
--core-dir --out --summary --kind --n-bins --per-bin --per-image --mag-zeropoint --psf-fwhm --match-radius
--thresh --minarea --smooth-fwhm --back-size --seed --mag-min --mag-max --pixel-scale --deblend --deblend-nlevels
--deblend-mincont --re-bins --catalog --catalog-out
"""
from __future__ import annotations

import sys

from launch import call_core, drop_empty


def main(argv: list[str]) -> int:
    argv = drop_empty(list(argv), "--catalog", "--catalog-out")
    return call_core("ds10core.tools.counts", argv, echo=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
