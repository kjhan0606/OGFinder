#!/usr/bin/env python3
"""Original shared-core spectrum preparation across a process boundary.

Forwarded flags: --core-dir --summary --out --wave-frame --flux-unit --exclude
--degree --clip --iterations
"""
import sys
from launch import call_core

if __name__ == '__main__':
    args = sys.argv[1:]
    task_args = args[2:] if args[:1] == ['--core-dir'] else args
    sys.exit(call_core('ds10core.tools.spectrum_prepare', args, echo=bool(task_args and task_args[0] == 'continuum')))
