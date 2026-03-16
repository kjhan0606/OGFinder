"""Parallel map using multiprocessing with fork context."""

import multiprocessing as mp
import os
import sys


def resolve_n_workers(n_workers):
    """Resolve n_workers value: 0 → cpu_count()-1 (max 64), 1 → sequential."""
    if n_workers == 0:
        return min(max(1, os.cpu_count() - 1), 64)
    return max(1, n_workers)


def parallel_map(worker_fn, tasks, n_workers=0, label="Processing"):
    """Apply worker_fn to each task, optionally in parallel.

    Parameters
    ----------
    worker_fn : callable
        Function applied to each task. Must accept a single argument (tuple).
    tasks : list
        List of task arguments to map over.
    n_workers : int
        0 = auto (cpu_count-1, max 64), 1 = sequential, N = N workers.
    label : str
        Label for progress messages on stderr.

    Returns
    -------
    list : Results in the same order as tasks.
    """
    from .progress import ProgressTracker

    n = len(tasks)
    if n == 0:
        return []

    actual_workers = resolve_n_workers(n_workers)
    tracker = ProgressTracker(total=n, label=label)

    if actual_workers == 1 or n == 1:
        # Sequential — no pool overhead, easier to debug
        results = []
        for task in tasks:
            results.append(worker_fn(task))
            tracker.update()
        tracker.finish()
        return results

    # Parallel via fork context (safe for numpy/sep on Linux)
    ctx = mp.get_context('fork')
    chunksize = max(1, n // (actual_workers * 4))

    print(f"{label}: {n} tasks, {actual_workers} workers", file=sys.stderr)

    results = []
    with ctx.Pool(processes=actual_workers) as pool:
        for result in pool.imap(worker_fn, tasks, chunksize=chunksize):
            results.append(result)
            tracker.update()

    tracker.finish()
    return results
