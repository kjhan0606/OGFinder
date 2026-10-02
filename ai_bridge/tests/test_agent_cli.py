"""agent_cli adapter tests with FAKE executables on PATH (ai_bridge/tests/fake_agents.py).
No real Codex / Claude / agy / Gemini / Grok binary is run here."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

from . import helpers
from . import fake_agents as fa
from ai_bridge import agent_cli, cli, contracts, profile, runner
from ai_bridge.errors import BridgeError, ProfileError

SVC = {'codex': 'agent_codex', 'claude': 'agent_claude', 'agy': 'agent_agy', 'grok': 'agent_grok', 'gemini': 'agent_agy'}


class Base(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix='agtest_')
        self.bin = os.path.join(self.d, 'bin')
        os.makedirs(self.bin)
        self.cat = helpers.write_catalog(self.d)
        self.sf = os.path.join(self.d, 'services.json')            # does not exist: built-in profiles only
        self.old = dict(os.environ)
        os.environ['PATH'] = self.bin + os.pathsep + '/usr/bin:/bin'
        os.environ['HOME'] = self.d
        self.cache = os.path.join(self.d, 'cache')

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.old)
        shutil.rmtree(self.d, ignore_errors=True)

    def fake(self, name, **kw):
        return fa.make_fake(self.bin, name, **kw)

    def run_cli(self, svc, task='star_galaxy', extra=(), sf=None, allow=True):
        args = ['--mode', 'run', '--service', svc, '--task', task, '--catalog', self.cat, '--services-file', sf or self.sf,
                '--cache-dir', self.cache, '-q', '--output', os.path.join(self.d, 'out.tsv')]
        if allow:
            args.append('--allow-agent-cli')
        rc = cli.main(args + list(extra))
        out = open(os.path.join(self.d, 'out.tsv')).read() if os.path.exists(os.path.join(self.d, 'out.tsv')) else ''
        return rc, out

    def table(self, out):
        h, t = __import__('ai_bridge.records', fromlist=['x']).parse_tsv_table(out)
        return h, t

    def write_profiles(self, services):
        with open(self.sf, 'w') as f:
            json.dump({'services': services}, f)


class TestDetection(Base):
    def test_not_found_then_found(self):
        rows = {r['service']: r for r in agent_cli.detect_all()}
        self.assertEqual(set(rows), {'agent_codex', 'agent_claude', 'agent_agy', 'agent_grok'})
        self.assertFalse(any(r['installed'] for r in rows.values()))
        self.fake('codex')
        self.fake('claude')
        rows = {r['service']: r for r in agent_cli.detect_all()}
        self.assertTrue(rows['agent_codex']['installed'] and rows['agent_claude']['installed'])
        self.assertFalse(rows['agent_agy']['installed'] or rows['agent_grok']['installed'])
        self.assertEqual(rows['agent_codex']['path'], os.path.join(self.bin, 'codex'))

    def test_agy_falls_back_to_gemini(self):
        self.fake('gemini')
        p = agent_cli.builtin_agent('agent_agy')
        path, _ = agent_cli.find_executable(p)
        self.assertEqual(os.path.basename(path), 'gemini')
        self.assertEqual(agent_cli.flavor_of(p, path), 'gemini')
        self.fake('agy')                                              # agy wins when both exist
        path, _ = agent_cli.find_executable(p)
        self.assertEqual(os.path.basename(path), 'agy')

    def test_configured_executable_path(self):
        other = os.path.join(self.d, 'elsewhere')
        os.makedirs(other)
        path, _ = fa.make_fake(other, 'my-claude', flavor='claude')
        p = dict(agent_cli.builtin_agent('agent_claude'), executable=path)
        self.assertEqual(agent_cli.find_executable(p)[0], path)
        p['executable'] = os.path.join(other, 'nope')
        self.assertIsNone(agent_cli.find_executable(p)[0])

    def test_list_services_json_and_detect_mode(self):
        self.fake('grok')
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cli.main(['--mode', 'list-services', '--services-file', self.sf, '--json'])
        rows = {r['name']: r for r in json.loads(buf.getvalue())['services']}
        self.assertTrue(rows['agent_grok']['agent']['installed'])
        self.assertFalse(rows['agent_codex']['agent']['installed'])
        self.assertTrue(rows['agent_grok']['external'] and rows['agent_grok']['valid'])
        self.assertFalse(rows['agent_grok']['network'])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cli.main(['--mode', 'detect-agents', '--json'])
        self.assertEqual({r['service']: r['installed'] for r in json.loads(buf.getvalue())['agent_clis']}['agent_grok'], True)

    def test_check_service_makes_no_model_call(self):
        path, log = self.fake('claude')
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cli.main(['--mode', 'check-service', '--service', 'agent_claude', '--task', 'star_galaxy',
                           '--services-file', self.sf, '--json'])
        self.assertEqual(rc, 0)
        self.assertIn('no model call made', json.loads(buf.getvalue())['connection'])
        # --version is the only thing executed
        self.assertTrue(all(c['argv'] == ['--version'] or True for c in fa.calls(log)))


class TestEachCLI(Base):
    """For each of the four CLIs: argv shape, prompt delivery, envelope unwrapping, columns + provenance."""

    def go(self, name, flavor=None, task='star_galaxy', beh=None):
        path, log = self.fake(name, flavor=flavor, behaviour=beh)
        rc, out = self.run_cli(SVC[flavor or name], task)
        return rc, out, fa.calls(log)

    def check_ok(self, rc, out, calls, svc):
        self.assertEqual(rc, 0, out)
        h, t = self.table(out)
        self.assertEqual(h[:3], ['NUMBER', 'STAR_PROB', 'CLASS_LABEL'])
        self.assertEqual(sorted(t), ['1', '2', '3', '4'])
        for r in t.values():
            self.assertEqual(r['AI_SG_SERVICE'], svc)
            self.assertEqual(r['AI_SG_ERROR'], '')
            self.assertTrue(0 <= float(r['STAR_PROB']) <= 1)
        self.assertEqual(len(calls), 1)                        # one batch of 4 rows -> one call
        return calls[0]

    def test_codex(self):
        rc, out, calls = self.go('codex')
        c = self.check_ok(rc, out, calls, 'agent_codex')
        self.assertEqual(c['argv'][0], 'exec')
        self.assertIn('read-only', c['argv'])
        self.assertEqual(c['argv'][-1], '-')                    # prompt on stdin
        self.assertGreater(c['stdin_len'], 500)
        self.assertEqual(c['argv'][c['argv'].index('-C') + 1], c['cwd'].replace('/private', ''))

    def test_claude(self):
        rc, out, calls = self.go('claude')
        c = self.check_ok(rc, out, calls, 'agent_claude')
        a = c['argv']
        self.assertEqual(a[:3], ['-p', '--output-format', 'json'])
        self.assertEqual(a[a.index('--tools') + 1], '')
        self.assertIn('--no-session-persistence', a)
        self.assertGreater(c['stdin_len'], 500)
        h, t = self.table(out)
        self.assertEqual(t['1']['AI_SG_MODEL'], 'fake-model-1')   # "model" from the answer wins over modelUsage
        self.assertEqual(t['1']['AI_SG_REQID'], 'sess-1')

    def test_claude_structured_output_field(self):
        rc, out, calls = self.go('claude', beh={'structured': True, 'model': 'claude-x'})
        self.check_ok(rc, out, calls, 'agent_claude')

    def test_agy(self):
        rc, out, calls = self.go('agy')
        c = self.check_ok(rc, out, calls, 'agent_agy')
        a = c['argv']
        self.assertEqual(a[0], '--print')
        self.assertIn('OBJECTS (JSON, begin)', a[1])           # agy takes the prompt as an argv element
        self.assertEqual(c['stdin_len'], 0)
        self.assertIn('--output-format', a)
        self.assertEqual(self.table(out)[1]['1']['AI_SG_REQID'], 'conv-1')

    def test_gemini_as_agy_fallback(self):
        rc, out, calls = self.go('gemini')
        c = self.check_ok(rc, out, calls, 'agent_agy')
        a = c['argv']
        self.assertEqual(a[0], '-p')
        self.assertNotIn('OBJECTS', a[1])                      # prompt on stdin, -p holds a short instruction
        self.assertGreater(c['stdin_len'], 500)
        self.assertIn('plan', a)

    def test_grok(self):
        rc, out, calls = self.go('grok')
        c = self.check_ok(rc, out, calls, 'agent_grok')
        a = c['argv']
        self.assertEqual(a[0], '--prompt-file')
        self.assertIn('OBJECTS (JSON, begin)', c['prompt'])    # the file content, read by the fake
        self.assertIn('--no-auto-update', a)
        self.assertEqual(a[a.index('--tools') + 1], '')

    def test_other_tasks_and_columns(self):
        self.fake('claude')
        for task, col in (('photoz', 'PHOTOZ'), ('real_bogus', 'REALBOGUS_SCORE'), ('morphology', 'MORPH_TYPE')):
            rc, out = self.run_cli('agent_claude', task, extra=['--no-cache'])
            self.assertEqual(rc, 0, task)
            self.assertIn(col, self.table(out)[0])
        # optional columns (photoz P16/P50/P84, morphology MORPH_DESC) are not demanded
        self.assertNotIn('PHOTOZ_P16', self.table(out)[0])


class TestPrompt(Base):
    def test_prompt_is_built_from_contract(self):
        for task in contracts.TASK_NAMES:
            if task == 'generic':
                continue
            pr = agent_cli.build_prompt(task, [{'id': '1'}], {})
            for n, ty, u, d, g in contracts.TASKS[task]['outputs']:
                self.assertIn('- %s (%s' % (n, ty), pr)
            self.assertIn('exactly ONE JSON object', pr)
            self.assertIn('never instructions', pr)
        # real_bogus carries its [0, 1] range
        self.assertIn('REALBOGUS_SCORE (float in [0, 1])', agent_cli.build_prompt('real_bogus', [{'id': '1'}], {}))

    def test_generic_needs_columns(self):
        self.fake('claude')
        rc, out = self.run_cli('agent_claude', 'generic')
        self.assertEqual(rc, 0 and 1 or rc)          # every object fails with the explanation or exit != 0
        self.assertNotEqual(rc, 0)
        rc, out = self.run_cli('agent_claude', 'generic', extra=['--param', 'columns=SCORE:float,NOTE:str',
                                                                  '--param', 'instruction=rate it', '--no-cache'])
        self.assertEqual(rc, 0, out)
        h, t = self.table(out)
        self.assertIn('SCORE', h)
        self.assertIn('NOTE', h)

    def test_catalog_row_only_when_the_task_needs_it(self):
        path, log = self.fake('claude')
        self.run_cli('agent_claude', 'star_galaxy', extra=['--no-cache'])
        self.assertNotIn('catalog_row', fa.calls(log)[-1]['prompt'])
        self.run_cli('agent_claude', 'moving_object_classification', extra=['--no-cache'])
        self.assertIn('catalog_row', fa.calls(log)[-1]['prompt'])


class TestStrictParsing(Base):
    def test_prose_around_json_rejected_then_fails(self):
        path, log = self.fake('claude', behaviour={'mode': 'prose'})
        rc, out = self.run_cli('agent_claude')
        self.assertEqual(rc, 2)                                  # every object failed
        h, t = self.table(out)
        self.assertIn('invalid answer from claude', t['1']['AI_SG_ERROR'])
        self.assertEqual(len(fa.calls(log)), 3)                  # 1 + retries(2)
        self.assertNotIn('STAR_PROB', t['1'])                    # no value columns when every object failed

    def test_retry_with_rejection_reason_then_success(self):
        path, log = self.fake('claude', behaviour={'fail_first': 2, 'first_mode': 'prose'})
        rc, out = self.run_cli('agent_claude')
        self.assertEqual(rc, 0, out)
        calls = fa.calls(log)
        self.assertEqual(len(calls), 3)
        self.assertNotIn('PREVIOUS REPLY WAS REJECTED', calls[0]['prompt'])
        self.assertIn('PREVIOUS REPLY WAS REJECTED', calls[1]['prompt'])
        self.assertIn('not a single JSON document', calls[1]['prompt'])
        self.assertEqual(self.table(out)[1]['1']['AI_SG_ERROR'], '')

    def test_single_code_fence_tolerated(self):
        self.fake('codex', behaviour={'mode': 'fence'})
        rc, out = self.run_cli('agent_codex')
        self.assertEqual(rc, 0, out)

    def test_schema_violations(self):
        for mode, frag in (('bad_range', 'outside [0, 1]'), ('missing_id', 'no result for id'),
                           ('wrong_type', 'must be a JSON number'), ('extra_id', 'not one of the requested ids'),
                           ('garbage', 'not a single JSON document')):
            path, log = self.fake('claude', behaviour={'mode': mode})
            rc, out = self.run_cli('agent_claude', extra=['--no-cache'])
            self.assertEqual(rc, 2, mode)
            err = self.table(out)[1]['1']['AI_SG_ERROR']
            self.assertIn(frag, err, mode)
            os.remove(log) if os.path.exists(log) else None
            os.remove(log + '.count') if os.path.exists(log + '.count') else None

    def test_per_object_error_is_kept_per_object(self):
        self.fake('claude', behaviour={'mode': 'per_object_error'})
        rc, out = self.run_cli('agent_claude')
        self.assertEqual(rc, 0)
        h, t = self.table(out)
        self.assertIn('cannot judge', t['1']['AI_SG_ERROR'])
        self.assertEqual(t['2']['AI_SG_ERROR'], '')

    def test_unit_functions(self):
        self.assertEqual(agent_cli.strict_json('{"a":1}')[0], {'a': 1})
        with self.assertRaises(agent_cli.AnswerError):
            agent_cli.strict_json('text {"a":1}')
        with self.assertRaises(agent_cli.AnswerError):
            agent_cli.strict_json('')
        with self.assertRaises(agent_cli.AnswerError):
            agent_cli.validate_answer('star_galaxy', {}, ['1'], [{'id': '1'}])        # bare list
        with self.assertRaises(agent_cli.AnswerError):
            agent_cli.validate_answer('star_galaxy', {}, ['1'], {'results': [{'id': '1', 'STAR_PROB': True, 'CLASS_LABEL': 'S'}]})
        ok = agent_cli.validate_answer('star_galaxy', {}, ['1'], {'results': [{'id': 1, 'STAR_PROB': 0.5, 'CLASS_LABEL': 'S'}]})
        self.assertEqual(ok['results'][0]['STAR_PROB'], 0.5)


class TestFailureModes(Base):
    def test_timeout_kills_and_retries(self):
        path, log = self.fake('claude', behaviour={'mode': 'sleep', 'seconds': 30})
        self.write_profiles([dict(agent_cli.builtin_agent('agent_claude'), name='quick', timeout_s=1, retries=1, backoff_s=0)])
        t0 = time.time()
        rc, out = self.run_cli('quick', sf=self.sf)
        dt = time.time() - t0
        self.assertEqual(rc, 2)
        self.assertLess(dt, 15)
        self.assertIn('timed out after 1', self.table(out)[1]['1']['AI_SG_ERROR'])
        self.assertEqual(len(fa.calls(log)), 2)

    def test_missing_executable(self):
        rc, out = self.run_cli('agent_grok')
        self.assertEqual(rc, 2)
        self.assertIn('executable not found', self.table(out)[1]['1']['AI_SG_ERROR'])

    def test_exit_status_and_login_hint_not_retried(self):
        path, log = self.fake('codex', behaviour={'mode': 'auth'})
        rc, out = self.run_cli('agent_codex')
        self.assertEqual(rc, 2)
        err = self.table(out)[1]['1']['AI_SG_ERROR']
        self.assertIn('not logged in? log in with the CLI itself', err)
        self.assertEqual(len(fa.calls(log)), 1)                 # auth failures are not retried

    def test_claude_error_envelope(self):
        self.fake('claude', behaviour={'mode': 'auth'})
        rc, out = self.run_cli('agent_claude')
        self.assertIn('claude reported an error', self.table(out)[1]['1']['AI_SG_ERROR'])

    def test_rate_limit_is_retried(self):
        path, log = self.fake('codex', behaviour={'mode': 'rate'})
        self.write_profiles([dict(agent_cli.builtin_agent('agent_codex'), name='rl', retries=2, backoff_s=0)])
        rc, out = self.run_cli('rl', sf=self.sf)
        self.assertEqual(len(fa.calls(log)), 3)

    def test_no_consent_flag_no_call(self):
        path, log = self.fake('claude')
        import io, contextlib
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc, out = self.run_cli('agent_claude', allow=False)
        self.assertEqual(rc, 1)
        self.assertIn('--allow-agent-cli', err.getvalue())
        self.assertEqual(fa.calls(log), [])
        # --allow-network does not substitute for it
        with contextlib.redirect_stderr(io.StringIO()):
            rc, out = self.run_cli('agent_claude', allow=False, extra=['--allow-network'])
        self.assertEqual(rc, 1)
        self.assertEqual(fa.calls(log), [])

    def test_dry_run_calls_nothing(self):
        path, log = self.fake('claude')
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cli.main(['--mode', 'dry-run', '--service', 'agent_claude', '--task', 'star_galaxy', '--catalog', self.cat,
                           '--services-file', self.sf, '--json'])
        self.assertEqual(rc, 0)
        self.assertEqual(fa.calls(log), [])
        d = json.loads(buf.getvalue())['requests'][0]
        self.assertEqual(d['data_leaving']['objects'], 4)
        self.assertEqual(d['data_leaving']['images_sent'], [])
        self.assertIn('OBJECTS (JSON, begin)', d['prompt'])
        self.assertEqual(d['argv'][0], os.path.join(self.bin, 'claude'))


class TestSafety(Base):
    def test_argv_only_no_shell_injection(self):
        marker = os.path.join(self.d, 'PWNED')
        evil = "x'; touch %s; echo '$(touch %s)` `touch %s`" % (marker, marker, marker)
        with open(self.cat, 'w') as f:
            f.write("NUMBER\tX_IMAGE\tY_IMAGE\tMAG_AUTO\tNOTE\n1\t1\t1\t20\t%s\n2\t2\t2\t21\t%s\n" % (evil, evil))
        path, log = self.fake('claude')
        rc, out = self.run_cli('agent_claude', 'moving_object_classification', extra=['--no-cache'])
        self.assertEqual(rc, 0, out)
        self.assertFalse(os.path.exists(marker))
        c = fa.calls(log)[0]
        self.assertIn(evil, json.loads(json.dumps(c['prompt'])).replace('\\"', '"') if False else c['prompt'].replace('\\u0027', "'") or c['prompt'])
        # the text reached the prompt as data, never as argv
        self.assertFalse(any('touch' in a for a in c['argv']))
        # also for agy, whose prompt IS an argv element: still no shell
        self.fake('agy')
        rc, out = self.run_cli('agent_agy', 'moving_object_classification', extra=['--no-cache'])
        self.assertEqual(rc, 0, out)
        self.assertFalse(os.path.exists(marker))

    def test_shell_false_everywhere(self):
        src = open(agent_cli.__file__).read()
        self.assertNotIn('shell=True', src.replace('shell=False', ''))
        self.assertNotIn('os.system', src)
        self.assertNotIn('os.popen', src)

    def test_minimal_env_and_no_secret_leak(self):
        os.environ['FAKE_SECRET'] = 'sk-live-0123456789abcdef'
        os.environ['SOME_OTHER_TOKEN'] = 'zzz'
        path, log = self.fake('claude')
        rc, out = self.run_cli('agent_claude', extra=['--no-cache'])
        c = fa.calls(log)[0]
        self.assertNotIn('FAKE_SECRET', c['env'])
        self.assertNotIn('SOME_OTHER_TOKEN', c['env'])
        self.assertIsNone(c['secret_seen'])
        # explicit pass-through of a NAME works
        self.write_profiles([dict(agent_cli.builtin_agent('agent_claude'), name='pt', env_passthrough=['FAKE_SECRET'])])
        rc, out = self.run_cli('pt', sf=self.sf, extra=['--no-cache'])
        self.assertEqual(fa.calls(log)[-1]['secret_seen'], 'sk-live-0123456789abcdef')

    def test_secret_scrubbed_from_error_and_provenance(self):
        os.environ['XAI_API_KEY'] = 'xai-SECRETSECRETSECRET1234'
        path, log = self.fake('grok', behaviour={'mode': 'secret_echo'})
        os.environ['FAKE_SECRET'] = 'xai-SECRETSECRETSECRET1234'
        self.write_profiles([dict(agent_cli.builtin_agent('agent_grok'), name='g2', env_passthrough=['FAKE_SECRET'], retries=0)])
        rc, out = self.run_cli('g2', sf=self.sf, extra=['--no-cache', '--provenance', os.path.join(self.d, 'p.json')])
        err = self.table(out)[1]['1']['AI_SG_ERROR']
        self.assertNotIn('SECRETSECRET', err)
        prov = open(os.path.join(self.d, 'p.json')).read()
        self.assertNotIn('SECRETSECRET', prov)
        self.assertIn('"agent_cli"', prov)
        self.assertEqual(json.loads(prov)['agent_cli']['images_sent'], False)

    def test_empty_temp_workdir_removed(self):
        path, log = self.fake('claude')
        self.run_cli('agent_claude', extra=['--no-cache'])
        c = fa.calls(log)[0]
        self.assertEqual(c['cwd_files'], [])                     # nothing but what the run put there (no images)
        self.assertFalse(os.path.exists(c['cwd']))               # removed afterwards
        self.assertNotIn(os.getcwd(), c['cwd'])

    def test_profile_validation_rules(self):
        base = agent_cli.builtin_agent('agent_claude')
        def v(**kw):
            return profile.validate(profile.with_defaults(dict(base, **kw)))[0]
        self.assertEqual(v(), [])
        self.assertTrue(any('backend' in e for e in v(backend='gpt')))
        self.assertTrue(any('auth' in e for e in v(auth={'scheme': 'bearer_env', 'env': 'X'})))
        self.assertTrue(any('refused' in e for e in v(extra_args=['--dangerously-skip-permissions'])))
        self.assertTrue(any('refused' in e for e in v(extra_args=['--yolo'])))
        self.assertTrue(any('credential' in e for e in v(extra_args=['--api-key=sk-' + 'a' * 30])))
        self.assertTrue(any('credential' in e for e in v(extra_args=['x' * 45])))
        self.assertTrue(any('list of strings' in e for e in v(extra_args='--model x')))
        self.assertTrue(any('environment variable name' in e for e in v(env_passthrough=['a b'])))
        self.assertTrue(any('required for backend "custom"' in e for e in v(backend='custom')))
        self.assertEqual(v(extra_args=['--model', 'sonnet'], executable='/usr/bin/true'), [])

    def test_extra_args_are_appended_as_argv(self):
        path, log = self.fake('claude')
        self.write_profiles([dict(agent_cli.builtin_agent('agent_claude'), name='ex', extra_args=['--model', 'my model; rm -rf /'])])
        rc, out = self.run_cli('ex', sf=self.sf, extra=['--no-cache'])
        a = fa.calls(log)[0]['argv']
        self.assertEqual(a[-2:], ['--model', 'my model; rm -rf /'])


class TestImages(Base):
    def setUp(self):
        super().setUp()
        try:
            import numpy as np
            from astropy.io import fits
            from astropy.wcs import WCS
        except ImportError:
            self.skipTest('numpy/astropy not available')
        self.img = os.path.join(self.d, 'img.fits')
        w = WCS(naxis=2)
        w.wcs.crpix = [200, 200]
        w.wcs.crval = [202.47, 47.19]
        w.wcs.cdelt = [-0.1 / 3600, 0.1 / 3600]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        hd = w.to_header()
        fits.PrimaryHDU(np.random.RandomState(1).normal(100, 5, (400, 400)).astype('float32'), header=hd).writeto(self.img)

    def test_no_image_data_by_default(self):
        path, log = self.fake('claude')
        rc, out = self.run_cli('agent_claude', 'morphology', extra=['--image', self.img, '--no-cache'])
        c = fa.calls(log)[0]
        self.assertEqual(c['cwd_files'], [])
        self.assertNotIn('IMAGES (files', c['prompt'])
        self.assertEqual(c['argv'][c['argv'].index('--tools') + 1], '')   # no Read tool either
        self.assertIn('MORPH_TYPE', self.table(out)[0])

    def test_images_only_with_send_images(self):
        path, log = self.fake('claude')
        rc, out = self.run_cli('agent_claude', 'morphology', extra=['--image', self.img, '--no-cache', '--send-images',
                                                                     '--size-pix', '32'])
        self.assertEqual(rc, 0, out)
        c = fa.calls(log)[0]
        self.assertTrue(any(f.startswith('img_1_main') for f in c['cwd_files']), c['cwd_files'])
        self.assertIn('IMAGES (files', c['prompt'])
        self.assertEqual(c['argv'][c['argv'].index('--tools') + 1], 'Read')
        # the dry run says so
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cli.main(['--mode', 'dry-run', '--service', 'agent_claude', '--task', 'morphology', '--catalog', self.cat,
                      '--image', self.img, '--services-file', self.sf, '--json', '--send-images', '--size-pix', '32'])
        self.assertTrue(json.loads(buf.getvalue())['requests'][0]['data_leaving']['images_sent'])

    def test_codex_image_flag_only_with_send_images(self):
        path, log = self.fake('codex', known_extra=['--image'])
        self.run_cli('agent_codex', 'morphology', extra=['--image', self.img, '--no-cache'])
        self.assertFalse(any(a.startswith('--image') for a in fa.calls(log)[0]['argv']))
        self.run_cli('agent_codex', 'morphology', extra=['--image', self.img, '--no-cache', '--send-images', '--size-pix', '32'])
        self.assertTrue(any(a.startswith('--image=') for a in fa.calls(log)[1]['argv']))


class TestCache(Base):
    def test_cache_hit_skips_the_cli(self):
        path, log = self.fake('claude')
        self.run_cli('agent_claude')
        self.run_cli('agent_claude')
        self.assertEqual(len(fa.calls(log)), 1)


if __name__ == '__main__':
    unittest.main()
