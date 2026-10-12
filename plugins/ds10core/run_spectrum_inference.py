#!/usr/bin/env python3
"""Original shared-core bundled inference adapter.
Forwarded flags: --core-dir --summary --config
"""
import sys
from launch import call_core
if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.spectrum_inference',sys.argv[1:]))
