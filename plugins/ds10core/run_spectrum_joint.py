#!/usr/bin/env python3
"""Original shared-core joint emission Gaussian adapter.
Forwarded flags: --core-dir --summary --components --fit-range --shift-limit
--continuum-model --agn-index --pivot --template-note
--redshift --wave-frame --flux-unit --exclude --min-snr
"""
import sys
from launch import call_core
if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.spectrum_joint',sys.argv[1:]))
