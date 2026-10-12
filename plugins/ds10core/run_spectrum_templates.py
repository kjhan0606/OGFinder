#!/usr/bin/env python3
"""Original shared-core fixed-template decomposition adapter.
Forwarded flags: --core-dir --summary --wave-frame --flux-unit --exclude
--fit-range --agn-index --pivot --template-note
"""
import sys
from launch import call_core
if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.spectrum_templates',sys.argv[1:]))
