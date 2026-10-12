import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
for p in (ROOT, os.path.dirname(HERE), os.path.join(ROOT, 'ds9', 'library')):
    if p not in sys.path:
        sys.path.insert(0, p)
