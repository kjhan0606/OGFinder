"""Own an Xvfb the way scripts/run_all_checks.sh does.

The display is left alone when a server is already up (run_all_checks.sh starts
:77 itself and kills that PID with its own trap).  A server this process starts
is killed on a normal return, an exception, and on INT, TERM and HUP.  Only that
PID is signalled.  The child is started with fds above 2 closed, so it cannot
inherit a lock fd.
"""
import atexit
import os
import shutil
import signal
import subprocess
import time


def display_up(display):
    """True when something is already serving this display.  Never kills a process."""
    if shutil.which('xdpyinfo'):
        r = subprocess.run(['xdpyinfo'], env=dict(os.environ, DISPLAY=display),
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return r.returncode == 0
    r = subprocess.run(['pgrep', '-f', 'Xvfb %s' % display], stdout=subprocess.DEVNULL)
    return r.returncode == 0


class Guard:
    def __init__(self):
        self.proc = None
        self._installed = False

    def install(self):
        if self._installed:
            return
        self._installed = True
        atexit.register(self.stop)
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(sig, self._on_signal)

    def ensure(self, display, wait=1.5):
        """Start Xvfb when the display is down.  Return the PID we started, else None."""
        if self.proc is not None and self.proc.poll() is None:
            return self.proc.pid
        if display_up(display):
            return None
        # close_fds is the default; set it so a future change cannot inherit lock fds.
        self.proc = subprocess.Popen(
            ['Xvfb', display, '-screen', '0', '1400x1000x24'],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            close_fds=True)
        if wait:
            time.sleep(wait)
        return self.proc.pid

    def stop(self):
        proc = self.proc
        if proc is None:
            return
        self.proc = None
        if proc.poll() is not None:
            return
        try:
            proc.terminate()  # SIGTERM to this PID only, not the process group
        except OSError:
            return
        try:
            proc.wait(timeout=2)
        except subprocess.TimeoutExpired:
            proc.kill()
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass

    def _on_signal(self, signum, frame):
        self.stop()
        signal.signal(signum, signal.SIG_DFL)
        os.kill(os.getpid(), signum)
