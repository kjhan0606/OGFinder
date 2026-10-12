#!/usr/bin/env python3
"""Original shared-core spectral peak proposal adapter.
Forwarded flags: --core-dir --summary --lines --wave-frame --flux-unit --exclude
--polarity --min-snr --z-min --z-max --tolerance --baseline
"""
import sys
from launch import call_core
if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.spectrum_search',sys.argv[1:]))
