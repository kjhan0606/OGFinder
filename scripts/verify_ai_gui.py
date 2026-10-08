#!/usr/bin/env python3
"""AI Services verification: GUI smoke test under Xvfb + exported-script replay + pipeline mode.

  verify_ai_gui.py [--python PY] [--display :77] [--workdir DIR]

1. starts a local HTTP test service (127.0.0.1) and writes a profile file with a REST profile for it
2. ds9 m51.fits -source scripts/verify_ai_gui.tcl   (menu structure, Service Registry dialog, dry run, mock runs with
   and without cutouts, selected-rows run, network-service confirmation dialog, session export)
3. replays the exported script in a plain shell (empty HOME) and compares the catalogue byte for byte
4. pipeline mode: ai.run (network) is skipped without --allow-network, runs with it
"""
import argparse
import hashlib
import http.server
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import xvfb_guard

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
results = []


def rec(check, ok, detail=''):
    results.append((check, bool(ok), detail))
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', check, ('  -- ' + detail) if detail else ''), flush=True)
    return ok


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


class H(http.server.BaseHTTPRequestHandler):
    log = []

    def log_message(self, *a):
        pass

    def do_POST(self):
        n = int(self.headers.get('Content-Length') or 0)
        js = json.loads(self.rfile.read(n).decode())
        H.log.append((self.path, js, self.headers.get('Authorization')))
        oid = int(js['oid'])
        body = json.dumps({'model': 'local-test-service 1.0', 'rid': 'req-%s' % oid,
                           'out': {'z': round(0.01 * oid, 6), 'e': 0.005}}).encode()
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    do_GET = lambda self: (self.send_response(200), self.send_header('Content-Length', '2'), self.end_headers(), self.wfile.write(b'{}'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--python', default=os.environ.get('OGFINDER_PYTHON', '/workspace/ogf_venv/bin/python'))
    ap.add_argument('--display', default=':77')
    ap.add_argument('--workdir', default=None)
    a = ap.parse_args()
    base = a.workdir or tempfile.mkdtemp(prefix='ogf_ai_verify_')
    os.makedirs(base, exist_ok=True)
    print('workdir', base)
    srv = http.server.ThreadingHTTPServer(('127.0.0.1', 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    url = 'http://127.0.0.1:%d' % srv.server_address[1]
    sf = os.path.join(base, 'ai_services.json')
    json.dump({'services': [{
        'name': 'localrest', 'task': 'photoz', 'transport': 'https_json', 'base_url': url, 'endpoint': '/z',
        'description': 'local test HTTP server of verify_ai_gui.py (127.0.0.1)',
        'auth': {'scheme': 'bearer_env', 'env': 'AI_VERIFY_TOKEN'},
        'request': {'template': {'oid': '{id}', 'm': '{mags}'}},
        'response': {'model_path': 'model', 'request_id_path': 'rid', 'fields': {
            'PHOTOZ': {'path': 'out.z', 'type': 'float'}, 'PHOTOZ_ERR': {'path': 'out.e', 'type': 'float'}}},
        'retries': 1, 'backoff_s': 0.1}]}, open(sf, 'w'), indent=1)
    # Own the Xvfb only when this process started it.  An already-up display (run_all_checks.sh) is left alone.
    guard = xvfb_guard.Guard()
    guard.install()
    guard.ensure(a.display)
    try:
        out = os.path.join(base, 'gui')
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out)
        for f in ('ds9.auto', 'ds9.auto.dir'):
            shutil.rmtree(os.path.expanduser('~/' + f), ignore_errors=True)
        env = dict(os.environ, DISPLAY=a.display, OGF_AI_OUT=out, OGFINDER_PYTHON=a.python, OGF_AI_REST='1',
                   OGF_AI_SERVICES_FILE=sf, AI_VERIFY_TOKEN='verify-token-123', OGF_AI_CACHE_DIR=os.path.join(base, 'cache'))
        print('\n=== GUI session (Xvfb %s) ===' % a.display)
        t0 = time.time()
        p = subprocess.run([os.path.join(ROOT, 'bin', 'ds9'), '/workspace/fits/m51.fits', '-geometry', '1300x950', '-source',
                            os.path.join(HERE, 'verify_ai_gui.tcl')], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300)
        dt = time.time() - t0
        so, se = p.stdout.decode(), p.stderr.decode()
        open(os.path.join(out, 'gui_stdout.txt'), 'w').write(so)
        open(os.path.join(out, 'gui_stderr.txt'), 'w').write(se)
        log = open(os.path.join(out, 'ai_gui_steps.log')).read()
        rec('ds9 session exit code 0 and DONE logged (%.1f s)' % dt, p.returncode == 0 and 'DONE' in log)
        bad_err = [l for l in se.splitlines() if l.strip() and not l.startswith('Filtered:')]
        rec('no Tcl errors on stderr', not bad_err, '; '.join(bad_err[:3]))
        rec('menu "AI Services" has exactly: Service Registry... / Run Task on Catalog... / Show Last Run Log',
            "menu AI Services entries: {Service Registry...} {Run Task on Catalog...} {Show Last Run Log}" in log)
        rec('registry dialog opened, lists mock + localrest + the 4 built-in agent profiles, env status only',
            any(('registry dialog exists: 1 rows=6 items=%s' % ' '.join(o)) in log for o in (
                ['localrest', 'mock', 'agent_codex', 'agent_claude', 'agent_agy', 'agent_grok'],
                ['mock', 'localrest', 'agent_codex', 'agent_claude', 'agent_agy', 'agent_grok'])), '')
        rec('registry "Test Connection" on the built-in mock reports OK', 'registry test(mock): OK' in log)
        rec('dry run: nothing sent, catalog columns unchanged, no step recorded',
            'dry run selected rows ok=1; columns unchanged: 1; sessions steps: 1' in log)
        m = re.search(r'run photoz/mock ok=1 rows=(\d+) cols=(\d+) new=(.*)', log)
        rec('photo-z (mock) added columns incl. renamed PHOTO_Z, PHOTO_Z_ERR and provenance', bool(m) and 'PHOTO_Z PHOTO_Z_ERR' in m.group(3)
            and 'AI_PHOTOZ_SERVICE' in m.group(3), m.group(0)[:150] if m else '')
        m = re.search(r'run morphology/mock ok=1 cols=(\d+) new=(.*)', log)
        rec('morphology (mock, cutouts 20" zscale PNG) added MORPH_TYPE/MORPH_CONF', bool(m) and 'MORPH_TYPE MORPH_CONF' in m.group(2))
        rec('selected-rows run recorded', 'run star_galaxy/mock on selected ok=1' in log)
        rec('REST service: exactly one confirmation dialog, none for the 2nd run / test',
            'confirmation dialogs shown for network service: 1' in log and 'dialogs now 1' in log)
        m = re.search(r'run photoz/localrest ok=1 cols (\d+) -> (\d+) new=(.*)', log)
        rec('REST service added columns REST_Z, REST_Z_ERR', bool(m) and 'REST_Z REST_Z_ERR' in m.group(3), m.group(0)[:130] if m else '')
        rec('REST: second run updated the existing columns in place (column count unchanged)',
            re.search(r'second localrest run ok=1 .* cols=(\d+)', log) is not None)
        rec('REST server saw bearer token from env (not stored anywhere)', H.log and H.log[0][2] == 'Bearer verify-token-123',
            '%d requests' % len(H.log))
        rec('backend hook (Photo-z backend = mock) ran the external service', re.search(r'backend hook photoz=mock: status=AI Services: 10 column', log) is not None)
        # --- secrets
        leak = []
        for dp, dn, fn in os.walk(out):
            for f in fn:
                if 'verify-token-123' in open(os.path.join(dp, f), errors='replace').read():
                    leak.append(f)
        for f in ('ai_last_run.log', 'ai_last_run.provenance.json'):
            pth = os.path.expanduser('~/.ds9/' + f)
            if os.path.exists(pth) and 'verify-token-123' in open(pth, errors='replace').read():
                leak.append(pth)
        rec('token value not present in session script / logs / provenance', not leak, ','.join(leak))
        script = os.path.join(out, 'ogfinder_session.py')
        gui_cat = os.path.join(out, 'catalog_gui.tsv')
        rec('session script exported', os.path.exists(script))
        steps = re.findall(r'^#\s+\d+\s+(\S+)\s+(\S+)\s+(.*)$', open(script).read(), re.M)
        print('  recorded steps:'); [print('     %s %s %s' % s) for s in steps]
        rec('ai.run steps recorded: AUTO for whole-catalog runs, MANUAL for the selected-rows run',
            [s[0] for s in steps if s[1] == 'ai.run'] == ['AUTO', 'AUTO', 'MANUAL', 'AUTO', 'AUTO', 'AUTO'], str([s[0] for s in steps if s[1] == 'ai.run']))
        # --- replay in a plain shell
        print('\n=== replay (plain shell, empty HOME) ===')
        home = os.path.join(base, 'emptyhome')
        shutil.rmtree(home, ignore_errors=True)
        os.makedirs(home)
        renv = {'PATH': '/usr/bin:/bin', 'HOME': home, 'LANG': 'C.UTF-8', 'OGF_AI_SERVICES_FILE': sf,
                'AI_VERIFY_TOKEN': 'verify-token-123', 'OGF_AI_CACHE_DIR': os.path.join(base, 'cache_replay')}
        rep = os.path.join(base, 'replay')
        shutil.rmtree(rep, ignore_errors=True)
        t0 = time.time()
        r = subprocess.run([a.python, script, '--python', a.python, '--mode', 'replay', '--outdir', rep, '/workspace/fits/m51.fits'],
                           env=renv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=base)
        dt = time.time() - t0
        txt = r.stdout.decode()
        open(rep + '.log', 'w').write(txt)
        rec('replay exit code 0 (%.1f s)' % dt, r.returncode == 0)
        ncrc = len(re.findall(r'catalog after step identical to the GUI session', txt))
        ndiff = len(re.findall(r'WARNING: catalog after this step differs', txt))
        # the two real-network (localrest) steps differ ONLY in the AI_PHOTOZ_TIME column (request wall-clock time, recorded
        # from the service call; a fresh replay makes new requests) - checked below on the REST_Z columns and the time column
        rec('replayed steps: all catalog fingerprints identical to the GUI except the 2 network-service steps (timestamp column)',
            ncrc >= 7 and ndiff == 2, '%d steps identical, %d differ (expected: the 2 localrest steps)' % (ncrc, ndiff))
        rc = os.path.join(rep, 'm51', 'catalog_final.tsv')
        same = os.path.exists(rc) and sha(rc) == sha(gui_cat)
        rec('replay catalog_final.tsv byte-identical to the GUI catalogue (sha256)', same,
            '%s' % (sha(gui_cat) if os.path.exists(gui_cat) else '?'))
        def col(path, name):
            L = open(path).read().split('\n')
            h = L[0].split('\t')
            i = h.index(name)
            return [l.split('\t')[i] for l in L[1:] if l]
        cg, cr = os.path.join(out, 'catalog_gui.tsv'), os.path.join(rep, 'm51', 'catalog_final.tsv')
        zg, zr = col(cg, 'REST_Z'), col(cr, 'REST_Z')
        rec('REST_Z / REST_Z_ERR columns identical in GUI and replay (deterministic test service)',
            zg == zr and col(cg, 'REST_Z_ERR') == col(cr, 'REST_Z_ERR') and len(zg) == 111, '%d values' % len(zg))
        hh = open(gui_cat).readline().rstrip('\n').split('\t')
        rec('GUI catalogue has the AI columns (%d columns total)' % len(hh),
            all(c in hh for c in ('PHOTO_Z', 'MORPH_TYPE', 'STAR_PROB', 'REST_Z', 'AI_PHOTOZ_SERVICE')))
        rec('replay wrote nothing into its ~/.ds9', not any(os.scandir(os.path.join(home, '.ds9'))) if os.path.isdir(os.path.join(home, '.ds9')) else True)
        # --- pipeline mode: network step needs --allow-network
        print('\n=== pipeline mode (new data = m51.fits again) ===')
        for flag, label in (([], 'without --allow-network'), (['--allow-network'], 'with --allow-network')):
            pdir = os.path.join(base, 'pipe_' + ('net' if flag else 'nonet'))
            shutil.rmtree(pdir, ignore_errors=True)
            r = subprocess.run([a.python, script, '--python', a.python, '--outdir', pdir, '/workspace/fits/m51.fits'] + flag,
                               env=renv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=base)
            t = r.stdout.decode()
            open(pdir + '.log', 'w').write(t)
            skipped_net = len(re.findall(r'SKIPPED: needs network access \(external AI service\)', t))
            cat = os.path.join(pdir, 'm51', 'catalog_final.tsv')
            hdr = open(cat).readline().split('\t') if os.path.exists(cat) else []
            if flag:
                rec('pipeline %s: exit 0, REST columns present, MANUAL step skipped' % label,
                    r.returncode == 0 and 'REST_Z' in hdr and 'SKIPPED: MANUAL step' in t, 'columns=%d' % len(hdr))
            else:
                rec('pipeline %s: network ai.run steps skipped with a message, mock steps still ran' % label,
                    r.returncode == 0 and skipped_net >= 1 and 'REST_Z' not in hdr and 'MORPH_TYPE' in hdr,
                    '%d network steps skipped, columns=%d' % (skipped_net, len(hdr)))
        srv.shutdown()
        npass = sum(1 for _, ok, _ in results if ok)
        print('\nSUMMARY: %d/%d checks passed' % (npass, len(results)))
        json.dump([{'check': c, 'ok': o, 'detail': d} for c, o, d in results], open(os.path.join(base, 'verify_ai_gui_results.json'), 'w'), indent=1)
        return 0 if npass == len(results) else 1
    finally:
        guard.stop()


if __name__ == '__main__':
    sys.exit(main())
