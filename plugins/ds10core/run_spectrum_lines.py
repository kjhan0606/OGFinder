#!/usr/bin/env python3
"""Original shared-core isolated spectral line adapter.
Forwarded flags: --core-dir --summary --wave-frame --flux-unit --exclude
--lines --redshift --window --polarity --min-snr
"""
import sys
from launch import call_core
if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.spectrum_lines', sys.argv[1:]))
