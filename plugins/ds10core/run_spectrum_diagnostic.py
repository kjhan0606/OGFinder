#!/usr/bin/env python3
"""Original shared-core joint-report ratio adapter.
Forwarded flags: --core-dir --summary --ratios --min-snr
"""
import sys
from launch import call_core
if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.spectrum_diagnostic',sys.argv[1:]))
