import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from . import helpers
from ai_bridge import adapters, cache as cachemod, cli, contracts, jsonpath, mapping, profile, records, runner, templating, units
from ai_bridge.errors import ProfileError, TemplateError
from .helpers import Server, rest_profile, write_catalog


class Tmp(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix='aibr_test_')
        self.cat = write_catalog(self.d)
        self.cache_dir = os.path.join(self.d, 'cache')
        self.sleeps = []
        self._p = mock.patch.object(runner, '_sleep', lambda s: self.sleeps.append(s))
        self._p.start()
        self.logs = []

    def tearDown(self):
        self._p.stop()
        shutil.rmtree(self.d, ignore_errors=True)

    def run_it(self, prof, task='photoz', **opts):
        o = dict(allow_network=True, cache_dir=self.cache_dir)
        o.update(opts)
        return runner.run(prof, task, self.cat, [], None, o, log=self.logs.append)


class TestTemplating(unittest.TestCase):
    def lk(self, n):
        return {'id': '7', 'ra': 1.5, 'mags': {'AUTO': 20.0}, 'none': None, 's': 'a"b'}[n]

    def test_typed_and_text(self):
        t = {'a': '{ra}', 'b': 'id={id}!', 'c': ['{mags}', '{none}'], 'd': '{{literal}}', 'e': 'm={mags}'}
        r = templating.render(t, self.lk)
        self.assertEqual(r['a'], 1.5)
        self.assertEqual(r['b'], 'id=7!')
        self.assertEqual(r['c'], [{'AUTO': 20.0}, None])
        self.assertEqual(r['d'], '{literal}')
        self.assertEqual(r['e'], 'm={"AUTO":20.0}')

    def test_unknown_placeholder(self):
        with self.assertRaises(TemplateError):
            templating.render({'a': '{nope}'}, self.lk)
        with self.assertRaises(TemplateError):
            templating.render({'a': 'x {nope} y'}, self.lk)

    def test_placeholders_listing(self):
        self.assertEqual(templating.placeholders({'a': '{ra}', 'b': ['{dec} {{x}}']}), ['ra', 'dec'])


class TestMapping(unittest.TestCase):
    def test_paths(self):
        o = {'a': {'b': [{'c': 1}, {'c': 2}]}, 'x': 3}
        self.assertEqual(jsonpath.get_path(o, 'a.b[1].c'), 2)
        self.assertEqual(jsonpath.get_path(o, '$.a.b[-1].c'), 2)
        self.assertEqual(jsonpath.get_path(o, 'a.b[*].c'), [1, 2])
        self.assertIs(jsonpath.get_path(o, 'a.z'), jsonpath.MISSING)
        self.assertIs(jsonpath.get_path(o, 'a.b[5]'), jsonpath.MISSING)

    def test_types_units_missing(self):
        item = {'z': '0.5', 'zerr': -99, 'n': 3.4, 'f': 'true', 'angle': 2.0, 'lab': 7, 'nanv': 'nan'}
        f = {'Z': {'path': 'z', 'type': 'float'}, 'ZE': {'path': 'zerr', 'type': 'float', 'missing_values': [-99]},
             'N': {'path': 'n', 'type': 'int'}, 'F': {'path': 'f', 'type': 'bool'},
             'A': {'path': 'angle', 'type': 'float', 'unit': 'arcsec', 'to_unit': 'arcmin'},
             'L': {'path': 'lab', 'type': 'str'}, 'NV': {'path': 'nanv', 'type': 'float'}, 'MISS': {'path': 'nothere', 'type': 'float'}}
        r = mapping.apply_fields(item, f)
        self.assertEqual(r['Z'], 0.5)
        self.assertIsNone(r['ZE'])
        self.assertEqual(r['N'], 3)
        self.assertIs(r['F'], True)
        self.assertAlmostEqual(r['A'], 2.0 / 60.0)
        self.assertEqual(r['L'], '7')
        self.assertIsNone(r['NV'])
        self.assertIsNone(r['MISS'])

    def test_scale_offset_and_errors(self):
        r = mapping.apply_fields({'v': 10}, {'X': {'path': 'v', 'type': 'float', 'scale': 0.1, 'offset': 1}})
        self.assertAlmostEqual(r['X'], 2.0)
        with self.assertRaises(mapping.MapError):
            mapping.apply_fields({'v': 'abc'}, {'X': {'path': 'v', 'type': 'float'}})
        with self.assertRaises(mapping.MapError):          # probability outside [0,1]
            mapping.apply_fields({'p': 1.7}, {'STAR_PROB': {'path': 'p', 'type': 'float'}})
        with self.assertRaises(mapping.MapError):          # nothing found at all
            mapping.apply_fields({'q': 1}, {'X': {'path': 'v', 'type': 'float'}})
        with self.assertRaises(ValueError):
            units.factor('arcsec', 'Jy')


class TestProfile(unittest.TestCase):
    def test_valid_and_invalid(self):
        good = rest_profile('https://example.invalid')
        e, w = profile.validate(profile.with_defaults(good))
        self.assertEqual(e, [])
        self.assertTrue(any('example.invalid' in x for x in w))
        bad = rest_profile('https://x.org', auth={'scheme': 'bearer_env', 'env': 'a b'}, transport='nope')
        e, _ = profile.validate(profile.with_defaults(bad))
        self.assertTrue(any('auth.env' in x for x in e) and any('transport' in x for x in e))

    def test_no_secrets_in_profile(self):
        p = rest_profile('https://x.org', auth={'scheme': 'bearer_env', 'env': 'TOK', 'value': 'abc'})
        e, _ = profile.validate(profile.with_defaults(p))
        self.assertTrue(any('auth.value' in x for x in e))
        p = rest_profile('https://x.org', headers={'Authorization': 'Bearer ' + 'a' * 30})
        e, w = profile.validate(profile.with_defaults(p))
        self.assertTrue(any('literal secret' in x for x in e))

    def test_http_refused_except_localhost(self):
        e, _ = profile.validate(profile.with_defaults(rest_profile('http://example.org')))
        self.assertTrue(any('http://' in x for x in e))
        e, _ = profile.validate(profile.with_defaults(rest_profile('http://127.0.0.1:1234')))
        self.assertEqual(e, [])

    def test_unknown_placeholder_in_template(self):
        p = rest_profile('https://x.org', request={'template': {'a': '{bogus}'}})
        e, _ = profile.validate(profile.with_defaults(p))
        self.assertTrue(any('{bogus}' in x for x in e))

    def test_mapping_must_hit_contract(self):
        p = rest_profile('https://x.org')
        p['response']['fields'] = {'FOO': {'path': 'a', 'type': 'float'}}
        e, _ = profile.validate(profile.with_defaults(p))
        self.assertTrue(any('standard columns' in x for x in e))

    def test_hash_ignores_notes(self):
        a = rest_profile('https://x.org')
        b = dict(a, _note='hello')
        self.assertEqual(profile.profile_hash(a), profile.profile_hash(b))
        c = dict(a, timeout_s=6)
        self.assertNotEqual(profile.profile_hash(a), profile.profile_hash(c))

    def test_example_file_validates(self):
        path = os.path.join(helpers.ROOT, 'ai_services.example.json')
        lst, _, err = profile.load_file(path)
        self.assertIsNone(err)
        self.assertGreaterEqual(len(lst), 6)
        for p in lst:
            e, _ = profile.validate(profile.with_defaults(p))
            self.assertEqual(e, [], p['name'])
            if p.get('template'):
                self.assertFalse(p.get('enabled', True))
                self.assertTrue(p['transport'] in ('local_command', 'python_callable', 'agent_cli') or 'example.invalid' in p['base_url'])


class TestRecords(unittest.TestCase):
    def test_records(self):
        d = tempfile.mkdtemp()
        try:
            h, rows = records.read_catalog(write_catalog(d))
            recs = records.build_records(h, rows)
            self.assertEqual([r['id'] for r in recs], ['1', '2', '3', '4'])
            self.assertEqual(recs[0]['mags'], {'AUTO': 20.1, 'F160W': 20.0})
            self.assertIsNone(recs[1]['mags']['F160W'])            # 99 -> none
            self.assertEqual(recs[1]['mag_errs']['F160W'], None)
            sel = records.build_records(h, rows, numbers=['3', '1'])
            self.assertEqual([r['id'] for r in sel], ['1', '3'])   # catalogue order, not request order
        finally:
            shutil.rmtree(d)


class TestMockAndDeterminism(Tmp):
    def test_mock_tags_and_determinism(self):
        r1 = self.run_it(profile.builtin_mock(), no_cache=True)
        r2 = self.run_it(profile.builtin_mock(), no_cache=True)
        self.assertEqual(r1['table'], r2['table'])
        for oid, row in r1['table'].items():
            self.assertEqual(row['AI_PHOTOZ_SERVICE'], 'mock')
            self.assertIn('FAKE', row['AI_PHOTOZ_MODEL'])
            self.assertLessEqual(row['PHOTOZ_P16'], row['PHOTOZ_P84'])
        self.assertEqual(r1['summary']['objects_ok'], 4)
        self.assertEqual(r1['ids'], ['1', '2', '3', '4'])

    def test_every_task_has_mock_output(self):
        for t in contracts.TASK_NAMES:
            r = self.run_it(profile.builtin_mock(), task=t, no_cache=True) if t != 'generic' else None
            if r:
                cols = [c for c, _ in r['columns']]
                for n in contracts.output_names(t):
                    self.assertIn(n, cols, t)


class TestRest(Tmp):
    def test_request_template_and_auth_from_env(self):
        def beh(h, rec):
            return 200, {'model': 'm-1.2', 'rid': 'r-%d' % len(srv.requests), 'out': {'z': 0.5, 'e': 0.1}}, {}
        srv = Server(beh)
        try:
            prof = rest_profile(srv.url, auth={'scheme': 'bearer_env', 'env': 'AIB_TEST_TOKEN'},
                                headers={'X-Static': 'yes'})
            with mock.patch.dict(os.environ, {'AIB_TEST_TOKEN': 'sekret-123'}):
                r = self.run_it(prof, max_objects=2)
            self.assertEqual(len(srv.requests), 2)
            q = srv.requests[0]
            self.assertEqual(q['path'], '/z')
            self.assertEqual(q['headers']['Authorization'], 'Bearer sekret-123')
            self.assertEqual(q['headers']['X-Static'], 'yes')
            self.assertEqual(q['json'], {'oid': '1', 'pos': {'ra': 202.4696, 'dec': 47.1952}, 'm': {'AUTO': 20.1, 'F160W': 20.0}})
            row = r['table']['1']
            self.assertEqual((row['PHOTOZ'], row['PHOTOZ_ERR']), (0.5, 0.1))
            self.assertEqual(row['AI_PHOTOZ_MODEL'], 'm-1.2')          # model string comes from the service
            self.assertTrue(row['AI_PHOTOZ_REQID'].startswith('r-'))
            self.assertTrue(row['AI_PHOTOZ_TIME'])
            # the secret appears neither in the profile hash material nor in provenance
            with mock.patch.dict(os.environ, {'AIB_TEST_TOKEN': 'sekret-123'}):
                prov, text = runner.write_outputs(r, os.path.join(self.d, 'o.tsv'), {}, self.cat, [])
            self.assertNotIn('sekret-123', json.dumps(prov) + text)
            self.assertEqual(prov['auth_env_var'], 'AIB_TEST_TOKEN')
            self.assertEqual(prov['auth_env_status_at_run'], 'set')
        finally:
            srv.close()

    def test_auth_header_and_query_schemes(self):
        srv = Server(lambda h, rec: (200, {'out': {'z': 1.0, 'e': 0.1}}, {}))
        try:
            with mock.patch.dict(os.environ, {'AIB_K': 'k-val'}):
                self.run_it(rest_profile(srv.url, auth={'scheme': 'header_env', 'header': 'X-Api-Key', 'env': 'AIB_K'}), max_objects=1)
                self.assertEqual(srv.requests[-1]['headers']['X-Api-Key'], 'k-val')
                self.run_it(rest_profile(srv.url, auth={'scheme': 'query_env', 'param': 'key', 'env': 'AIB_K'}), max_objects=1, no_cache=True)
                self.assertIn('key=k-val', srv.requests[-1]['path'])
        finally:
            srv.close()

    def test_unset_env_var_fails_cleanly_without_sending(self):
        srv = Server(lambda h, rec: (200, {'out': {'z': 1.0, 'e': 0.1}}, {}))
        try:
            env = {k: v for k, v in os.environ.items() if k != 'AIB_MISSING'}
            with mock.patch.dict(os.environ, env, clear=True):
                r = self.run_it(rest_profile(srv.url, auth={'scheme': 'bearer_env', 'env': 'AIB_MISSING'}))
            self.assertEqual(len(srv.requests), 0)
            self.assertEqual(r['summary']['objects_ok'], 0)
            self.assertIn('AIB_MISSING', r['table']['1']['AI_PHOTOZ_ERROR'])
            self.assertEqual(r['exit_code'], 2)
        finally:
            srv.close()

    def test_retry_backoff_then_success(self):
        state = {'n': 0}

        def beh(h, rec):
            state['n'] += 1
            if state['n'] <= 3:
                return 503, {'error': 'busy'}, {}
            return 200, {'out': {'z': 0.9, 'e': 0.1}}, {}
        srv = Server(beh)
        try:
            r = self.run_it(rest_profile(srv.url), max_objects=1)
            self.assertEqual(len(srv.requests), 4)
            self.assertEqual(self.sleeps, [0.5, 1.0, 2.0])           # exponential backoff
            self.assertEqual(r['stats'].retries, 3)
            self.assertEqual(r['table']['1']['PHOTOZ'], 0.9)
        finally:
            srv.close()

    def test_retry_after_and_backoff_cap(self):
        state = {'n': 0}

        def beh(h, rec):
            state['n'] += 1
            return (429, {}, {'Retry-After': '7'}) if state['n'] == 1 else (200, {'out': {'z': 0.9, 'e': 0.1}}, {})
        srv = Server(beh)
        try:
            self.run_it(rest_profile(srv.url, backoff_max_s=5), max_objects=1)
            self.assertEqual(self.sleeps, [5])                      # Retry-After 7 capped by backoff_max_s
        finally:
            srv.close()

    def test_no_retry_on_4xx_and_gives_up_after_retries(self):
        srv = Server(lambda h, rec: (400, {'error': 'bad'}, {}))
        try:
            r = self.run_it(rest_profile(srv.url), max_objects=1)
            self.assertEqual(len(srv.requests), 1)
            self.assertIn('HTTP 400', r['table']['1']['AI_PHOTOZ_ERROR'])
        finally:
            srv.close()
        srv = Server(lambda h, rec: (500, {'error': 'x'}, {}))
        try:
            r = self.run_it(rest_profile(srv.url, retries=2), max_objects=1)
            self.assertEqual(len(srv.requests), 3)
            self.assertEqual(r['stats'].retries, 2)
            self.assertEqual(r['exit_code'], 2)
        finally:
            srv.close()

    def test_timeout_is_retryable(self):
        import time as _t
        srv = Server(lambda h, rec: (_t.sleep(1.0), (200, {'out': {'z': 1, 'e': 1}}, {}))[1])
        try:
            r = self.run_it(rest_profile(srv.url, timeout_s=0.2, retries=1), max_objects=1)
            self.assertIn('timeout', r['table']['1']['AI_PHOTOZ_ERROR'])
            self.assertEqual(r['stats'].retries, 1)
        finally:
            srv.close()

    def test_partial_failure_and_resume(self):
        def beh(h, rec):
            if rec['json']['oid'] == '3':
                return 200, {'out': {'oops': 1}}, {}                 # maps nothing -> per-object error
            return 200, {'out': {'z': 0.2, 'e': 0.1}}, {}
        srv = Server(beh)
        out = os.path.join(self.d, 'out.tsv')
        try:
            r = self.run_it(rest_profile(srv.url))
            self.assertEqual(r['exit_code'], 0)
            self.assertEqual((r['summary']['objects_ok'], r['summary']['objects_failed']), (3, 1))
            self.assertIn('mapping', r['table']['3']['AI_PHOTOZ_ERROR'])
            prov, text = runner.write_outputs(r, out, {}, self.cat, [])
            lines = text.strip().split('\n')
            self.assertEqual(len(lines), 5)
            self.assertEqual(prov['objects_failed'], 1)
            self.assertEqual(prov['failures'][0]['id'], '3')
            # resume: only object 3 is requested again (cache disabled to prove it is the resume logic)
            n0 = len(srv.requests)
            r2 = self.run_it(rest_profile(srv.url), resume=True, no_cache=True) if False else \
                runner.run(rest_profile(srv.url), 'photoz', self.cat, [], out, dict(allow_network=True, no_cache=True, resume=True), log=self.logs.append)
            self.assertEqual(len(srv.requests) - n0, 1)
            self.assertEqual(r2['summary']['objects_failed'], 1)
        finally:
            srv.close()

    def test_cache_hit_and_key_sensitivity(self):
        srv = Server(lambda h, rec: (200, {'out': {'z': 0.3, 'e': 0.1}}, {}))
        try:
            p = rest_profile(srv.url)
            r1 = self.run_it(p)
            self.assertEqual((len(srv.requests), r1['stats'].cache_hits), (4, 0))
            r2 = self.run_it(p)
            self.assertEqual((len(srv.requests), r2['stats'].cache_hits), (4, 4))
            self.assertEqual(r1['table'], r2['table'])
            r3 = self.run_it(p, no_cache=True)
            self.assertEqual(len(srv.requests), 8)                    # --no-cache ignores the cache
            # a changed profile (different template) must not hit the old entries
            p2 = rest_profile(srv.url, request={'template': {'oid': '{id}', 'extra': 1}})
            self.run_it(p2)
            self.assertEqual(len(srv.requests), 12)
            # a changed request (different parameters) neither
            p3 = rest_profile(srv.url, request={'template': {'oid': '{id}', 'p': '{param_k}'}})
            self.run_it(p3, params={'k': 1})
            n = len(srv.requests)
            self.run_it(p3, params={'k': 2})
            self.assertEqual(len(srv.requests), n + 4)
            # a changed timeout is part of the profile hash too
            self.run_it(rest_profile(srv.url, timeout_s=6))
            self.assertEqual(len(srv.requests), n + 8)
            # cache key never contains credentials: file contents do not mention the token
            with mock.patch.dict(os.environ, {'AIB_T': 'tok-xyz'}):
                self.run_it(rest_profile(srv.url, auth={'scheme': 'bearer_env', 'env': 'AIB_T'}))
            blob = ''
            for dp, dn, fn in os.walk(self.cache_dir):
                for f in fn:
                    blob += open(os.path.join(dp, f)).read()
            self.assertNotIn('tok-xyz', blob)
        finally:
            srv.close()

    def test_garbage_response_not_cached(self):
        srv = Server(lambda h, rec: (200, {'unexpected': True}, {}))
        try:
            r = self.run_it(rest_profile(srv.url), max_objects=1)
            self.assertEqual(r['summary']['objects_ok'], 0)
            self.assertEqual(r['cache'].stores, 0)
        finally:
            srv.close()

    def test_batching_and_id_matching(self):
        def beh(h, rec):
            items = rec['json']['sources']
            return 200, {'model': 'bm', 'results': [{'id': i['id'], 'z': float(i['id']) / 10} for i in reversed(items) if i['id'] != '2']}, {}
        srv = Server(beh)
        try:
            prof = rest_profile(srv.url, batch_size=3, request={'template': {'sources': '{items}'},
                                                                'item_template': {'id': '{id}', 'm': '{mags}'}},
                                response={'items_path': 'results', 'id_path': 'id', 'model_path': 'model',
                                          'fields': {'PHOTOZ': {'path': 'z', 'type': 'float'}}})
            r = self.run_it(prof)
            self.assertEqual([len(q['json']['sources']) for q in srv.requests], [3, 1])
            self.assertEqual(r['table']['3']['PHOTOZ'], 0.3)             # matched by id although the order differs
            self.assertIn('not in response', r['table']['2']['AI_PHOTOZ_ERROR'])
            self.assertEqual(r['summary']['objects_failed'], 1)
        finally:
            srv.close()

    def test_max_payload_and_rate_limit(self):
        srv = Server(lambda h, rec: (200, {'out': {'z': 0.3, 'e': 0.1}}, {}))
        try:
            r = self.run_it(rest_profile(srv.url, max_payload_bytes=20), max_objects=1)
            self.assertEqual(len(srv.requests), 0)
            self.assertIn('max_payload_bytes', r['table']['1']['AI_PHOTOZ_ERROR'])
            self.sleeps.clear()
            with mock.patch.object(runner, '_clock', lambda: 100.0):
                self.run_it(rest_profile(srv.url, rate_limit_per_s=2), no_cache=True)
            self.assertEqual(self.sleeps, [0.5, 0.5, 0.5])           # 4 requests at 2/s with a frozen clock
        finally:
            srv.close()

    def test_redirect_to_other_host_refused(self):
        srv = Server(lambda h, rec: (302, b'', {'Location': 'http://localhost:9/other'}))
        try:
            with mock.patch.dict(os.environ, {'AIB_T': 'tok'}):
                r = self.run_it(rest_profile(srv.url, auth={'scheme': 'bearer_env', 'env': 'AIB_T'}, retries=0), max_objects=1)
            self.assertIn('redirect to a different host refused', r['table']['1']['AI_PHOTOZ_ERROR'])
        finally:
            srv.close()

    def test_network_needs_allow_network(self):
        from ai_bridge.errors import BridgeError
        with self.assertRaises(BridgeError):
            runner.run(rest_profile('https://example.invalid'), 'photoz', self.cat, [], None, {}, log=self.logs.append)

    def test_dry_run_sends_nothing_and_redacts_secret(self):
        srv = Server(lambda h, rec: (200, {}, {}))
        try:
            with mock.patch.dict(os.environ, {'AIB_T': 'sekret-xyz'}):
                r = runner.run(rest_profile(srv.url, auth={'scheme': 'bearer_env', 'env': 'AIB_T'}), 'photoz', self.cat, [], None,
                               {'dry_run': True}, log=self.logs.append)
            self.assertEqual(len(srv.requests), 0)
            self.assertEqual(len(r['previews']), 4)
            txt = json.dumps(r['previews'])
            self.assertNotIn('sekret-xyz', txt)
            self.assertIn('AIB_T', txt)
            self.assertEqual(r['previews'][0]['body']['oid'], '1')
            self.assertEqual(r['previews'][0]['method'], 'POST')
        finally:
            srv.close()


class TestMultipart(Tmp):
    def test_multipart_upload(self):
        try:
            import numpy  # noqa
            from astropy.io import fits  # noqa
        except ImportError:
            self.skipTest('numpy/astropy needed for cutouts')
        import numpy as np
        from astropy.io import fits
        img = os.path.join(self.d, 'im.fits')
        fits.PrimaryHDU(np.random.RandomState(1).normal(100, 5, (400, 400)).astype('float32')).writeto(img)
        srv = Server(lambda h, rec: (200, {'label': 'E', 'score': 0.8}, {}))
        try:
            prof = {'name': 'up', 'task': 'morphology', 'transport': 'https_multipart', 'base_url': srv.url, 'endpoint': '/c',
                    'cutouts': {'size_pix': 32, 'format': 'png'}, 'retries': 0,
                    'request': {'multipart': {'fields': {'oid': '{id}'}, 'files': {'image': '{cutout_png_path}'}}},
                    'response': {'fields': {'MORPH_TYPE': {'path': 'label', 'type': 'str'}, 'MORPH_CONF': {'path': 'score', 'type': 'float'}}}}
            r = runner.run(prof, 'morphology', self.cat, [img], None, dict(allow_network=True, no_cache=True, cutout_dir=os.path.join(self.d, 'cut'),
                                                                           max_objects=2), log=self.logs.append)
            q = srv.requests[0]
            self.assertEqual(len(srv.requests), 2)
            self.assertIn('multipart/form-data', q['headers']['Content-Type'])
            self.assertIn(b'name="image"; filename="obj1_main.png"', q['raw'])
            self.assertTrue(q['raw'].count(b'\x89PNG') == 1)
            self.assertEqual(r['table']['1']['MORPH_TYPE'], 'E')
            self.assertTrue(os.path.exists(os.path.join(self.d, 'cut', 'obj1_main.png')))
        finally:
            srv.close()


class TestLocalAndCallable(Tmp):
    def test_local_command(self):
        ex = os.path.join(helpers.ROOT, 'ai_bridge', 'examples', 'echo_local_command.py')
        prof = {'name': 'lc', 'task': 'star_galaxy', 'transport': 'local_command', 'command': ['{python}', ex], 'batch_size': 3,
                'timeout_s': 30}
        r = self.run_it(prof, task='star_galaxy')
        self.assertEqual(r['summary']['objects_ok'], 4)
        self.assertEqual(r['table']['3']['CLASS_LABEL'], 'GALAXY')        # min mag 18.2 -> galaxy in the toy rule
        self.assertEqual(r['table']['3']['STAR_PROB'], 0.1)
        self.assertIn('echo-example', r['table']['1']['AI_SG_MODEL'])
        self.assertEqual(r['stats'].requests_sent, 2)

    def test_local_command_failure_modes(self):
        prof = {'name': 'lc', 'task': 'star_galaxy', 'transport': 'local_command', 'command': ['{python}', '-c', 'import sys; sys.exit(3)'],
                'retries': 0}
        r = self.run_it(prof, task='star_galaxy')
        self.assertEqual(r['exit_code'], 2)
        self.assertIn('status 3', r['table']['1']['AI_SG_ERROR'])
        prof['command'] = ['/no/such/exe']
        r = self.run_it(prof, task='star_galaxy', no_cache=True)
        self.assertIn('not found', r['table']['1']['AI_SG_ERROR'])
        prof['command'] = ['{python}', '-c', 'print("not json")']
        r = self.run_it(prof, task='star_galaxy', no_cache=True)
        self.assertIn('not JSON', r['table']['1']['AI_SG_ERROR'])

    def test_python_callable(self):
        prof = {'name': 'pc', 'task': 'real_bogus', 'transport': 'python_callable', 'callable': 'example_callable:score',
                'sys_path': [os.path.join(helpers.ROOT, 'ai_bridge', 'examples')], 'batch_size': 10, 'params': {'scale': 2.0}}
        r = self.run_it(prof, task='real_bogus')
        self.assertEqual(r['summary']['objects_ok'], 4)
        self.assertAlmostEqual(r['table']['1']['REALBOGUS_SCORE'], min(1.0, 2.0 * 100.5 / 1000.0))
        prof['callable'] = 'example_callable:nope'
        r = self.run_it(prof, task='real_bogus', no_cache=True)
        self.assertIn('cannot import', r['table']['1']['AI_RB_ERROR'])


class TestCutouts(Tmp):
    def setUp(self):
        super().setUp()
        try:
            import numpy as np
            from astropy.io import fits
        except ImportError:
            self.skipTest('numpy/astropy needed')
        self.np, self.fits = np, fits
        rs = np.random.RandomState(0)
        data = rs.normal(100, 5, (300, 300)).astype('float32')
        data[148:153, 148:153] += 500
        self.img = os.path.join(self.d, 'a.fits')
        h = fits.Header()
        for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CRVAL1=10.0, CRVAL2=20.0, CRPIX1=150.5, CRPIX2=150.5,
                         CDELT1=-0.2 / 3600, CDELT2=0.2 / 3600).items():
            h[k] = v
        fits.PrimaryHDU(data, h).writeto(self.img)

    def test_cutout_formats_sizes_and_edges(self):
        from ai_bridge import cutouts
        im = cutouts.ImageSet({'main': self.img})
        self.assertAlmostEqual(im.pixel_scale_arcsec(), 0.2, places=6)
        rec = {'id': '1', 'x': 150.0, 'y': 150.0, 'ra': 10.0, 'dec': 20.0}
        for fmt in ('png', 'fits', 'npy'):
            mk = cutouts.CutoutMaker(im, os.path.join(self.d, 'c_' + fmt), size_arcsec=6.4, norm='asinh', fmt=fmt)
            self.assertEqual(mk.size_for('main'), 32)             # 6.4" / 0.2"/pix
            info = mk.make(rec, 'main')
            self.assertTrue(os.path.exists(info['path']))
            self.assertEqual(info['size'], 32)
        a = self.np.load(os.path.join(self.d, 'c_npy', 'obj1_main.npy'))
        self.assertEqual(a.shape, (32, 32))
        self.assertTrue(0.0 <= a.min() and a.max() <= 1.0)
        self.assertEqual(a[16, 16], a.max())                       # the bright source sits at the centre
        with self.fits.open(os.path.join(self.d, 'c_fits', 'obj1_main.fits')) as hl:
            self.assertEqual(hl[0].data.shape, (32, 32))
            self.assertAlmostEqual(float(hl[0].data[16, 16]), float(im.data['main'][149, 149]))
        with open(os.path.join(self.d, 'c_png', 'obj1_main.png'), 'rb') as f:
            self.assertEqual(f.read(8), b'\x89PNG\r\n\x1a\n')
        # edge object: padded with NaN (npy raw) and still produced; outside object -> None
        mk = cutouts.CutoutMaker(im, os.path.join(self.d, 'edge'), size_pix=20, norm='none', fmt='npy')
        info = mk.make({'id': '9', 'x': 2.0, 'y': 2.0}, 'main')
        arr = self.np.load(info['path'])
        self.assertTrue(self.np.isnan(arr[0, 0]) and self.np.isfinite(arr[-1, -1]))
        self.assertIsNone(mk.make({'id': '10', 'x': 5000.0, 'y': 5000.0}, 'main'))

    def test_deterministic_bytes_and_norm_modes(self):
        from ai_bridge import cutouts
        im = cutouts.ImageSet({'main': self.img})
        rec = {'id': '1', 'x': 150.0, 'y': 150.0}
        h = []
        for _ in range(2):
            mk = cutouts.CutoutMaker(im, os.path.join(self.d, 'det'), size_pix=40, norm='zscale', fmt='png')
            h.append(mk.make(rec, 'main')['sha256'])
        self.assertEqual(h[0], h[1])
        for how in cutouts.NORMALIZATIONS:
            n = cutouts.normalize(im.data['main'][100:140, 100:140], how)
            self.assertEqual(n.shape, (40, 40))
        with self.assertRaises(ValueError):
            cutouts.normalize(im.data['main'], 'bogus')

    def test_mock_with_images_uses_cutouts_and_cli_end_to_end(self):
        out = os.path.join(self.d, 'o.tsv')
        rc = cli.main(['--mode', 'run', '--service', 'mock', '--task', 'morphology', '--catalog', self.cat, '--image', self.img,
                       '--size-pix', '24', '--output', out, '--services-file', os.path.join(self.d, 'none.json'), '-q'])
        self.assertEqual(rc, 0)
        prov = json.load(open(out + '.provenance.json'))
        self.assertEqual(prov['service'], 'mock')
        self.assertEqual(prov['images'][0]['sha256'], runner.sha256_file(self.img))
        self.assertEqual(prov['objects_ok'], 4)


class TestCli(Tmp):
    def sf(self, profiles):
        p = os.path.join(self.d, 'svc.json')
        json.dump({'services': profiles}, open(p, 'w'))
        return p

    def call(self, *args):
        return cli.main(['--services-file', self.sf_path] + list(args))

    def test_modes_and_exit_codes(self):
        srv = Server(lambda h, rec: (200, {'out': {'z': 0.3, 'e': 0.1}}, {}))
        try:
            self.sf_path = self.sf([rest_profile(srv.url), rest_profile('https://example.invalid', name='tmpl', template=True, enabled=False)])
            out = os.path.join(self.d, 'o.tsv')
            base = ['--service', 'testrest', '--task', 'photoz', '--catalog', self.cat, '--cache-dir', self.cache_dir]
            # network service without --allow-network -> usage error, nothing sent
            self.assertEqual(self.call('--mode', 'run', *base, '--output', out, '-q'), 1)
            self.assertEqual(len(srv.requests), 0)
            self.assertEqual(self.call('--mode', 'dry-run', *base, '-q'), 0)
            self.assertEqual(len(srv.requests), 0)
            self.assertEqual(self.call('--mode', 'run', *base, '--allow-network', '--output', out, '-q'), 0)
            self.assertEqual(len(srv.requests), 4)
            lines = open(out).read().strip().split('\n')
            self.assertEqual(lines[0].split('\t')[:3], ['NUMBER', 'PHOTOZ', 'PHOTOZ_ERR'])
            self.assertEqual([l.split('\t')[0] for l in lines[1:]], ['1', '2', '3', '4'])
            prov = json.load(open(out + '.provenance.json'))
            self.assertEqual(prov['requests']['sent'], 4)
            self.assertEqual(len(prov['profile_sha256']), 64)
            self.assertEqual(self.call('--mode', 'validate-profile', '--service', 'testrest'), 0)
            self.assertEqual(self.call('--mode', 'validate-profile', '--service', 'tmpl'), 0)
            self.assertEqual(self.call('--mode', 'check-service', '--service', 'testrest', '--allow-network'), 0)  # GET base_url answers 200
            self.assertEqual(self.call('--mode', 'check-service', '--service', 'testrest'), 1)   # network needs --allow-network
        finally:
            srv.close()

    def test_total_failure_exit_code_and_disabled(self):
        srv = Server(lambda h, rec: (400, {}, {}))
        try:
            self.sf_path = self.sf([rest_profile(srv.url, retries=0), rest_profile(srv.url, name='off', enabled=False)])
            base = ['--task', 'photoz', '--catalog', self.cat, '--allow-network', '--no-cache', '-q']
            self.assertEqual(self.call('--mode', 'run', '--service', 'testrest', *base, '--output', os.path.join(self.d, 'o.tsv')), 2)
            self.assertTrue(os.path.exists(os.path.join(self.d, 'o.tsv')))                # still written, with errors
            self.assertEqual(self.call('--mode', 'run', '--service', 'off', *base), 1)
            self.assertEqual(self.call('--mode', 'run', '--service', 'nope', *base), 1)
        finally:
            srv.close()

    def test_strict_partial(self):
        self.sf_path = self.sf([])
        base = ['--service', 'mock', '--task', 'photoz', '--catalog', self.cat, '-q', '--mock-fail-ids', '2']
        self.assertEqual(self.call('--mode', 'run', *base), 0)
        self.assertEqual(self.call('--mode', 'run', *base, '--strict'), 3)

    def test_list_services_env_status_never_value(self):
        self.sf_path = self.sf([rest_profile('https://example.invalid', auth={'scheme': 'bearer_env', 'env': 'AIB_LIST'})])
        import io, contextlib
        buf = io.StringIO()
        with mock.patch.dict(os.environ, {'AIB_LIST': 'topsecret'}), contextlib.redirect_stdout(buf):
            self.call('--mode', 'list-services')
        self.assertIn('AIB_LIST\tset', buf.getvalue())
        self.assertNotIn('topsecret', buf.getvalue())


class TestSetEnabled(Tmp):
    def test_set_enabled_roundtrip(self):
        p = os.path.join(self.d, 'svc.json')
        json.dump({'_README': 'keep me', 'services': [rest_profile('https://example.invalid'), rest_profile('https://example.invalid', name='b')]}, open(p, 'w'))
        self.assertEqual(cli.main(['--services-file', p, '--mode', 'set-enabled', '--service', 'b', '--enabled', 'no']), 0)
        d = json.load(open(p))
        self.assertEqual(d['_README'], 'keep me')
        self.assertEqual([s.get('enabled', True) for s in d['services']], [True, False])
        self.assertTrue(os.path.exists(p + '.bak'))
        self.assertEqual(cli.main(['--services-file', p, '--mode', 'set-enabled', '--service', 'zzz', '--enabled', 'no']), 1)


class TestCheckService(Tmp):
    def test_check_service_http(self):
        srv = Server(lambda h, rec: (200, {'ok': True}, {}))
        try:
            p = rest_profile(srv.url, health={'path': '/health'})
            ad = adapters.make_adapter(runner.effective_profile(p, 'photoz'), {'params': {}, 'bands': []})
            ok, msg = ad.check()
            self.assertTrue(ok, msg)
            self.assertEqual(srv.requests[0]['path'], '/health')
        finally:
            srv.close()


class TestTap(Tmp):
    def test_tap_adql_and_csv(self):
        def beh(h, rec):
            return 200, b'main_id,sep\r\nSomething,1.5\r\n', {'Content-Type': 'text/csv'}
        srv = Server(beh)
        try:
            prof = {'name': 't', 'task': 'generic', 'transport': 'tap_query', 'base_url': srv.url + '/tap',
                    'params': {'radius_arcsec': 3.6}, 'provenance_prefix': 'T',
                    'request': {'adql': "SELECT main_id, 1.5 AS sep FROM basic WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{radius_deg}))"},
                    'response': {'fields': {'T_ID': {'path': 'rows[0].main_id', 'type': 'str'}, 'T_SEP': {'path': 'rows[0].sep', 'type': 'float'}}}}
            r = self.run_it(prof, task='generic', max_objects=1)
            q = srv.requests[0]
            self.assertEqual(q['path'], '/tap/sync')
            from urllib.parse import parse_qs
            f = parse_qs(q['raw'].decode())
            self.assertEqual(f['LANG'], ['ADQL'])
            self.assertEqual(f['REQUEST'], ['doQuery'])
            self.assertIn('CIRCLE(\'ICRS\',202.4696,47.1952,0.001)', f['QUERY'][0])
            self.assertEqual(r['table']['1']['T_ID'], 'Something')
            self.assertEqual(r['table']['1']['T_SEP'], 1.5)
            self.assertEqual(r['table']['1']['AI_T_SERVICE'], 't')
        finally:
            srv.close()


class TestContracts(unittest.TestCase):
    def test_contract_tables(self):
        for t, d in contracts.TASKS.items():
            self.assertIn('prefix', d)
            for o in d['outputs']:
                self.assertIn(o[1], ('float', 'int', 'str', 'bool', 'json'))
        md = contracts.contract_table_markdown()
        for t in contracts.TASK_NAMES:
            self.assertIn('`%s`' % t, md)
        self.assertEqual(len(contracts.TASK_NAMES), 12)


if __name__ == '__main__':
    unittest.main()
