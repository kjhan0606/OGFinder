#!/usr/bin/env python3
"""OGFinder AI bridge CLI (thin driver; the code lives in the `ai_bridge` package).

  python ds9_ai_bridge.py --mode {list-services,check-service,run,dry-run,validate-profile} --service NAME \\
      --task TASK --catalog CAT.tsv --image A.fits[,B.fits] --output OUT.tsv

See docs/ai_services.md.  The package is found next to the OGFinder tree (<root>/ai_bridge), via $OGF_AI_BRIDGE_DIR,
or from sys.path.
"""
import os
import sys


def _find_package_parent():
    cands = []
    if os.environ.get('OGF_AI_BRIDGE_DIR'):
        cands.append(os.path.dirname(os.path.abspath(os.environ['OGF_AI_BRIDGE_DIR'])))
    here = os.path.dirname(os.path.abspath(__file__))
    cands += [os.path.dirname(os.path.dirname(here)),                    # <root>/ds9/library -> <root>
              os.path.dirname(here), here]
    exe = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv else ''
    cands.append(os.path.dirname(exe))
    for c in cands:
        if c and os.path.isfile(os.path.join(c, 'ai_bridge', '__init__.py')):
            return c
    return None


_p = _find_package_parent()
if _p and _p not in sys.path:
    sys.path.insert(0, _p)
try:
    from ai_bridge.cli import main
except ImportError as e:
    sys.stderr.write('ERROR: cannot import the ai_bridge package (%s). Set OGF_AI_BRIDGE_DIR=<OGFinder>/ai_bridge\n' % e)
    sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())
