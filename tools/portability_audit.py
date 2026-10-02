#!/usr/bin/env python3
"""Static audit of POSIX-only constructs in the OGFinder Tcl / Python / shell sources (docs/windows_macos_build.md).

Nothing here proves a port works: it lists WHERE the code assumes a POSIX system, so a person with a Windows or macOS machine knows what
to look at first.  Usage:  tools/portability_audit.py [--list]   (counts; --list also prints every hit)
"""
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAN = [('Tcl', ['ds9/library', 'plugins'], ('.tcl',)),
        ('Python', ['ds9/library', 'plugins', 'moving', 'ai_bridge', 'tools', 'parallel', 'lsbg', 'icl', 'ai_merge', 'morphology', 'photoz'], ('.py',)),
        ('shell', ['scripts'], ('.sh',))]
SKIP_DIRS = {'__pycache__', 'tests', 'validation', 'golden'}
# (name, regex, applies to, why it matters)
RULES = [
    ('kill', r'\bexec\s+kill\b|OGFsess_exec\s+kill\b', 'Tcl', 'no kill on Windows (taskkill)'),
    ('tmp', r'["\' ]/tmp/', 'Tcl Python', 'hard-coded /tmp (use tempfile / ::ogf temp dir)'),
    ('shebang-exec', r'exec\s+(?:\./)?\S+\.py\b', 'Tcl', 'exec of a .py file directly (needs the interpreter on Windows)'),
    ('python3-name', r'["\']python3["\']', 'Tcl Python', 'literal "python3" (Windows uses python / py -3; OGFINDER_PYTHON overrides)'),
    ('fork', r'os\.fork\(|get_context\(["\']fork', 'Python', 'fork start method (not on Windows, not default on macOS)'),
    ('signal', r'signal\.(SIGKILL|SIGALRM|alarm)|os\.killpg|os\.setsid|preexec_fn', 'Python', 'POSIX signals / sessions'),
    ('posix-path', r'["\']/(?:usr|bin|etc|home|proc|dev)/', 'Tcl Python shell', 'absolute POSIX path'),
    ('shell-true', r'shell\s*=\s*True|os\.system\(|os\.popen\(', 'Python', 'goes through a shell'),
    ('bash', r'#!/bin/(?:ba)?sh|\bexec\s+(?:bash|sh)\b', 'Tcl shell', 'needs a POSIX shell (Cygwin / MSYS2 / WSL on Windows)'),
    ('xdotool', r'\bxdotool\b|\bxdpyinfo\b|\bXvfb\b', 'Tcl shell', 'X11-only test tooling'),
    ('pool', r'ProcessPoolExecutor|multiprocessing\.Pool|mp\.Pool|\bctx\.Pool|get_context\(', 'Python', 'process pool (spawn start method on Windows/macOS: entry points must be importable and guarded)'),
    ('ctypes-so', r'\.so["\']|CDLL\(|cdll\.', 'Python', 'shared library loading (.so / .dll / .dylib)'),
]

def files():
    for lang, dirs, exts in SCAN:
        for d in dirs:
            for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
                dn[:] = [x for x in dn if x not in SKIP_DIRS and not x.startswith('.')]
                for f in fn:
                    if f.endswith(exts) and f != 'portability_audit.py':
                        yield lang, os.path.join(dp, f)

def main():
    show = '--list' in sys.argv
    counts = {r[0]: [] for r in RULES}
    nfiles = 0
    for lang, path in files():
        nfiles += 1
        try:
            text = open(path, encoding='utf-8', errors='replace').read().splitlines()
        except OSError:
            continue
        for n, line in enumerate(text, 1):
            s = line.strip()
            if s.startswith('#') and not s.startswith('#!'):
                continue
            for name, rx, applies, why in RULES:
                if lang in applies.split() and re.search(rx, line):
                    counts[name].append('%s:%d: %s' % (os.path.relpath(path, ROOT), n, s[:110]))
    print('scanned %d files' % nfiles)
    print('%-14s %5s  %s' % ('construct', 'hits', 'why it matters'))
    for name, rx, applies, why in RULES:
        print('%-14s %5d  %s' % (name, len(counts[name]), why))
        if show:
            for h in counts[name]:
                print('      ' + h)
    return 0

if __name__ == '__main__':
    sys.exit(main())
