"""pytest plugin: force the spawn start method (what Windows/macOS use).  PYTHONPATH=tools python -m pytest -p spawnpatch <tests>  (round 2, item 5)"""
import multiprocessing as mp
_orig = mp.get_context
def get_context(method=None):
    return _orig('spawn' if method in (None, 'fork', 'forkserver') else method)
mp.get_context = get_context
mp.set_start_method('spawn', force=True)
