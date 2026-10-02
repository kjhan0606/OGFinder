#!/usr/bin/env python3
"""Run the agent_cli adapter against REAL agent-CLI binaries that happen to be installed - WITHOUT a login.

What this proves (and what not):
  * proves  : every flag the profile passes is accepted by the real binary's argument parser (all five reject unknown
              flags), that the binary starts, and that its *error* envelope for 'not logged in' is parsed into a clear
              per-object error with the login hint.
  * NOT     : it never obtains a model answer, so the success envelope of each CLI is only known from its documentation
              (and, for grok, not at all).  No data leaves the machine: with no credentials the CLIs fail before sending.

HOME is an empty temporary directory and no *_API_KEY variable is passed, so a user's real login is never used.
Exit 0 when every installed CLI behaved as described (not-installed CLIs are reported as SKIP; none installed -> 77).
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ai_bridge import agent_cli, cli  # noqa: E402

CAT = "NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\tMAG_AUTO\n1\t10\t10\t202.4\t47.1\t20.1\n2\t20\t20\t202.5\t47.2\t21.0\n"


def main():
    tmp = tempfile.mkdtemp(prefix='agent_real_')
    home = os.path.join(tmp, 'home')
    os.makedirs(home)
    cat = os.path.join(tmp, 'cat.tsv')
    open(cat, 'w').write(CAT)
    for k in list(os.environ):
        if k.endswith(('_API_KEY', '_TOKEN')) or k in ('GOOGLE_GENAI_USE_VERTEXAI', 'GOOGLE_GENAI_USE_GCA', 'CODEX_HOME', 'GROK_HOME'):
            del os.environ[k]
    os.environ['HOME'] = home
    rows = agent_cli.detect_all()
    exes = {}
    for name in ('codex', 'claude', 'agy', 'gemini', 'grok'):
        w = shutil.which(name)
        if w:
            exes[name] = w
    if not exes:
        print('no agent CLI installed: SKIP')
        return 77
    bad = 0
    for svc, exe in (('agent_codex', 'codex'), ('agent_claude', 'claude'), ('agent_agy', 'agy'), ('agent_agy', 'gemini'), ('agent_grok', 'grok')):
        if exe not in exes:
            print('SKIP  %-12s %s not installed' % (svc, exe))
            continue
        if svc == 'agent_agy' and exe == 'gemini' and 'agy' in exes:
            pass                                       # still test the gemini argv via an explicit executable
        prof = dict(agent_cli.builtin_agent(svc), name='real_' + exe, executable=exes[exe], retries=0, timeout_s=90)
        sf = os.path.join(tmp, 'sf_%s.json' % exe)
        json.dump({'services': [prof]}, open(sf, 'w'))
        out = os.path.join(tmp, 'out_%s.tsv' % exe)
        ver = agent_cli.version_of(exes[exe])
        rc = cli.main(['--mode', 'run', '--service', prof['name'], '--task', 'star_galaxy', '--catalog', cat,
                       '--services-file', sf, '--allow-agent-cli', '--no-cache', '--output', out, '-q'])
        txt = open(out).read() if os.path.exists(out) else ''
        err = txt.split('\n')[1].split('\t')[-1] if txt.count('\n') > 1 else ''
        arg_error = any(w in err.lower() for w in ('unknown option', 'unexpected argument', 'unknown arguments',
                                                   'flags provided but not defined', 'unrecognized'))
        login_hint = 'log in with the CLI itself' in err or 'not logged in' in err.lower() or 'not signed in' in err.lower()
        ok = rc == 2 and not arg_error and login_hint
        bad += not ok
        print('%s  %-12s %-8s %s\n        -> %s' % ('PASS' if ok else 'FAIL', svc, exe, ver, err[:260]))
    shutil.rmtree(tmp, ignore_errors=True)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
