#!/usr/bin/env python3
"""Run the original shared-core audit or profile stack across a process boundary.

Forwarded CLI flags: --core-dir --summary --out --bootstrap --seed
"""
import sys
from launch import call_core

if __name__ == '__main__':
    sys.exit(call_core('ds10core.tools.workflow_support', sys.argv[1:], echo='stack' in sys.argv[1:]))
