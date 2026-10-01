"""LIVE test of the tap_query adapter against real public IVOA TAP services (non-ML).  Needs internet.
Skipped when OGF_AI_OFFLINE=1.  If the network is unreachable the tests are SKIPPED (not passed) with the reason.
Run:  python -m pytest ai_bridge/tests/test_tap_live.py -v -rs
"""
import json
import os
import shutil
import tempfile
import unittest
import urllib.request

from . import helpers
from ai_bridge import profile, runner
from .helpers import write_catalog

# Messier 51 core region objects; positions are the catalogue's own (ALPHA_J2000/DELTA_J2000 in degrees)
CAT = ("NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\tMAG_AUTO\tMAGERR_AUTO\n"
       "1\t1\t1\t202.469575\t47.195258\t9.0\t0.01\n"      # M51 nucleus (NGC 5194)
       "2\t2\t2\t10.684708\t41.268750\t3.4\t0.01\n"       # M31 nucleus
       "3\t3\t3\t0.0\t-89.0\t20.0\t0.1\n")                # empty-ish field near the south pole


def reachable(url):
    try:
        urllib.request.urlopen(url, timeout=15).read(10)
        return True, ''
    except Exception as e:                                         # noqa
        return False, '%s: %s' % (type(e).__name__, e)


@unittest.skipIf(os.environ.get('OGF_AI_OFFLINE') == '1', 'OGF_AI_OFFLINE=1')
class TestLiveTap(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix='aibr_tap_')
        self.cat = write_catalog(self.d, CAT)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _run(self, prof, task='generic'):
        return runner.run(prof, task, self.cat, [], None, dict(allow_network=True, cache_dir=os.path.join(self.d, 'c')),
                          log=lambda m: None)

    def test_simbad_cone_search(self):
        ok, why = reachable('https://simbad.cds.unistra.fr/simbad/sim-tap/availability')
        if not ok:
            self.skipTest('SIMBAD TAP unreachable: ' + why)
        lst, _, err = profile.load_file(os.path.join(helpers.ROOT, 'ai_services.example.json'))
        prof = next(p for p in lst if p['name'] == 'simbad_cone')
        prof['params'] = {'radius_arcsec': 20.0}
        r = self._run(prof)
        t = r['table']
        self.assertEqual(r['summary']['objects'], 3)
        # M51 nucleus: SIMBAD must return an object within 20" (name text is whatever SIMBAD answers; we only check non-empty + sep)
        self.assertTrue(t['1']['SIMBAD_ID'], t['1'])
        self.assertLess(t['1']['SIMBAD_SEP'], 20.0)
        self.assertTrue(t['2']['SIMBAD_ID'], t['2'])
        self.assertLess(t['2']['SIMBAD_SEP'], 20.0)
        print('\nSIMBAD live: #1 ->', t['1']['SIMBAD_ID'], t['1']['SIMBAD_TYPE'], '%.2f"' % t['1']['SIMBAD_SEP'],
              '| #2 ->', t['2']['SIMBAD_ID'], t['2']['SIMBAD_TYPE'], '%.2f"' % t['2']['SIMBAD_SEP'], '| #3 ->', t['3'])
        # the third position has no SIMBAD object within 20": a valid empty answer -> empty cells, no error
        self.assertEqual(t['3']['AI_SIMBAD_ERROR'], '')
        self.assertIsNone(t['3']['SIMBAD_ID'])

    def test_vizier_tap_gaia_dr3_cone(self):
        ok, why = reachable('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/availability')
        if not ok:
            self.skipTest('VizieR TAP unreachable: ' + why)
        prof = {'name': 'vizier_gaia_dr3', 'task': 'generic', 'transport': 'tap_query', 'provenance_prefix': 'GAIA',
                'base_url': 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap', 'params': {'radius_arcsec': 30.0},
                'request': {'adql': 'SELECT TOP 1 "Source", "Gmag", "RA_ICRS", "DE_ICRS" FROM "I/355/gaiadr3" WHERE 1=CONTAINS('
                                    'POINT(\'ICRS\', "RA_ICRS", "DE_ICRS"), CIRCLE(\'ICRS\', {ra}, {dec}, {radius_deg})) ORDER BY "Gmag"',
                            'format': 'csv'},
                'response': {'fields': {'GAIA_SOURCE': {'path': 'rows[0].Source', 'type': 'str'},
                                        'GAIA_G': {'path': 'rows[0].Gmag', 'type': 'float'}}},
                'timeout_s': 90, 'retries': 1, 'backoff_s': 2.0}
        r = self._run(prof)
        t = r['table']
        if r['summary']['objects_ok'] == 0:
            self.fail('VizieR TAP returned nothing usable: %s' % t['1'])
        print('\nVizieR TAP live: #2 (M31 nucleus, r=30") brightest Gaia DR3 source:', t['2'].get('GAIA_SOURCE'), 'G=', t['2'].get('GAIA_G'))
        self.assertTrue(t['2']['GAIA_SOURCE'])
        self.assertTrue(isinstance(t['2']['GAIA_G'], float))

    def test_tap_error_is_reported_per_object(self):
        ok, why = reachable('https://simbad.cds.unistra.fr/simbad/sim-tap/availability')
        if not ok:
            self.skipTest('SIMBAD TAP unreachable: ' + why)
        prof = {'name': 'bad_adql', 'task': 'generic', 'transport': 'tap_query', 'base_url': 'https://simbad.cds.unistra.fr/simbad/sim-tap',
                'request': {'adql': 'SELECT nonexistent_column FROM basic WHERE ra = {ra}'},
                'response': {'fields': {'X': {'path': 'rows[0].nonexistent_column', 'type': 'str'}}}, 'retries': 0}
        r = self._run(prof)
        self.assertEqual(r['summary']['objects_ok'], 0)
        self.assertIn('TAP error', r['table']['1']['AI_GEN_ERROR'])
        print('\nTAP error text from the real service:', r['table']['1']['AI_GEN_ERROR'][:160])
