"""Progress tracking for parallel and sequential pipelines."""

import sys
import time


class ProgressTracker:
    """Stderr progress reporter with rate and ETA.

    Prints at most every 2 seconds to avoid flooding output.
    Supports an optional MPI rank prefix.
    """

    def __init__(self, total, label="Processing", rank=None, interval=2.0):
        self.total = total
        self.label = label
        self.rank = rank
        self.interval = interval
        self.count = 0
        self.start_time = time.monotonic()
        self._last_print = 0.0

    def update(self, n=1):
        self.count += n
        now = time.monotonic()
        if now - self._last_print >= self.interval or self.count >= self.total:
            self._print(now)
            self._last_print = now

    def _print(self, now):
        elapsed = now - self.start_time
        rate = self.count / elapsed if elapsed > 0 else 0
        pct = 100.0 * self.count / self.total if self.total > 0 else 0
        remaining = (self.total - self.count) / rate if rate > 0 else 0

        prefix = f"[rank {self.rank}] " if self.rank is not None else ""
        msg = (f"{prefix}{self.label}: {self.count}/{self.total} "
               f"({pct:.0f}%) [{rate:.1f}/s, {remaining:.0f}s]")
        print(f"\r{msg}", end='', file=sys.stderr, flush=True)

    def finish(self):
        elapsed = time.monotonic() - self.start_time
        prefix = f"[rank {self.rank}] " if self.rank is not None else ""
        print(f"\r{prefix}{self.label}: {self.count}/{self.total} done "
              f"in {elapsed:.1f}s", file=sys.stderr, flush=True)
        print(file=sys.stderr)  # newline
