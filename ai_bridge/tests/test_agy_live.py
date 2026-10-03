"""LIVE test of the agy profile (`agent_agy`) against the real `agy` CLI on a tiny star/galaxy task.

Skipped unless `agy` is installed and logged in (an OAuth token exists after the one-time interactive `agy` login)
and OGF_AI_OFFLINE is not 1.  Costs a few hundred tokens.
Run:  python -m pytest ai_bridge/tests/test_agy_live.py -v -rs
"""
import json
import os
import shutil
import subprocess
import tempfile
import unittest

from ai_bridge import agent_cli, cli

CAT = ("NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\tMAG_AUTO\tCLASS_STAR\tFWHM_IMAGE\n"
       "1\t10\t10\t202.4\t47.1\t18.0\t0.98\t2.1\n"
       "2\t20\t20\t202.5\t47.2\t17.5\t0.02\t12.5\n")


# a live service that is busy, rate-limited, slow or unreachable is not a failure of the code under test: SKIP with the reason
TRANSIENT = ('timed out', 'timeout', 'rate limit', 'rate-limit', '429', '502', '503', '504', 'overloaded', 'temporarily', 'try again', 'connection', 'unavailable', 'busy')


def logged_in():
    exe = shutil.which('agy')
    if not exe or os.environ.get('OGF_AI_OFFLINE') == '1':
        return None
    tok = os.path.expanduser('~/.gemini/antigravity-cli/antigravity-oauth-token')
    return exe if os.path.exists(tok) else None


EXE = logged_in()


@unittest.skipUnless(EXE, 'agy CLI not installed or not logged in (run `agy (interactive login once)`)')
class AgyLive(unittest.TestCase):
    def test_star_galaxy_roundtrip(self):
        tmp = tempfile.mkdtemp(prefix='grok_live_')
        try:
            cat = os.path.join(tmp, 'cat.tsv'); open(cat, 'w').write(CAT)
            prof = dict(agent_cli.builtin_agent('agent_agy'), name='live_agy', executable=EXE, retries=0, timeout_s=180)
            sf = os.path.join(tmp, 'sf.json'); json.dump({'services': [prof]}, open(sf, 'w'))
            out = os.path.join(tmp, 'out.tsv')
            rc = cli.main(['--mode', 'run', '--service', 'live_agy', '--task', 'star_galaxy', '--catalog', cat, '--services-file', sf,
                           '--allow-agent-cli', '--no-cache', '--output', out, '-q'])
            txt = open(out).read() if os.path.exists(out) else ''
            low = txt.lower()
            if 'not signed in' in low or '401' in low or 'unauthor' in low or 'invalid api key' in low or any(m in low for m in TRANSIENT):
                self.skipTest('agy is not really authenticated: %s' % txt.strip().splitlines()[-1][:200])
            self.assertEqual(rc, 0, txt[:500])
            rows = [l.split('\t') for l in txt.strip().splitlines()]
            self.assertEqual(len(rows), 3, txt)                       # header + 2 objects
            self.assertTrue(all('error' not in c.lower() for r in rows[1:] for c in r[-3:]), txt)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
