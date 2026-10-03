import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
for p in (ROOT, os.path.join(HERE, '..')):
    if p not in sys.path:
        sys.path.insert(0, p)
