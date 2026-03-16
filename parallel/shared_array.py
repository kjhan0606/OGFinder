"""Shared memory array for zero-copy image sharing across worker processes."""

import atexit
import uuid
from dataclasses import dataclass
from multiprocessing import shared_memory

import numpy as np

# Track all shared memory names created by this process for cleanup
_active_shm_names = set()


def _cleanup_shm():
    """atexit handler: unlink any leaked shared memory blocks."""
    for name in list(_active_shm_names):
        try:
            shm = shared_memory.SharedMemory(name=name, create=False)
            shm.close()
            shm.unlink()
        except FileNotFoundError:
            pass
        except Exception:
            pass
    _active_shm_names.clear()


atexit.register(_cleanup_shm)


@dataclass(frozen=True)
class SharedArraySpec:
    """Lightweight picklable descriptor for a shared memory array.

    Workers receive this instead of the full array (~100 bytes vs megabytes).
    """
    shm_name: str
    shape: tuple
    dtype: str

    def attach(self):
        """Attach to the shared memory block and return (ndarray view, shm handle).

        The caller MUST call shm.close() when done (but NOT unlink — only the
        parent SharedArray context manager does that).
        """
        shm = shared_memory.SharedMemory(name=self.shm_name, create=False)
        arr = np.ndarray(self.shape, dtype=self.dtype, buffer=shm.buf)
        return arr, shm


class SharedArray:
    """Context manager that publishes a numpy array to shared memory.

    Usage::

        with SharedArray(image_data) as spec:
            # spec is a SharedArraySpec, safe to pickle to workers
            results = parallel_map(worker_fn, [(spec, src) for src in catalog])
        # shared memory is unlinked on exit
    """

    def __init__(self, arr):
        self._arr = np.ascontiguousarray(arr)
        name = f"ds9_shm_{uuid.uuid4().hex[:16]}"
        self._shm = shared_memory.SharedMemory(name=name, create=True,
                                                size=self._arr.nbytes)
        # Copy data into shared memory
        view = np.ndarray(self._arr.shape, dtype=self._arr.dtype,
                          buffer=self._shm.buf)
        view[:] = self._arr
        self.spec = SharedArraySpec(
            shm_name=name,
            shape=self._arr.shape,
            dtype=str(self._arr.dtype),
        )
        _active_shm_names.add(name)

    def __enter__(self):
        return self.spec

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self._shm.close()
            self._shm.unlink()
        except Exception:
            pass
        _active_shm_names.discard(self.spec.shm_name)
        return False
