"""Standalone static review entry points; no submitted code execution."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'tools/user_tool_import.py'
EXAMPLE = ROOT / 'tools/user_tool_kit/examples/user_tool_submission'


class UserToolImportTests(unittest.TestCase):
    def test_print_submission_registration_audit(self):
        example = ROOT / 'examples/user_tool_print'
        baseline = subprocess.run([sys.executable, str(example / 'original.py')], capture_output=True, text=True, check=True)
        self.assertEqual(baseline.stdout.strip(), 'ASTRAFEX_USER_PRINT_OK')
        with tempfile.TemporaryDirectory(prefix='print audit ') as folder:
            env = dict(os.environ, UTI_ROOT=str(ROOT), UTI_PYTHON=sys.executable, UTI_EXAMPLE=str(example), UTI_TMP=folder)
            script = r'''
source [file join $::env(UTI_ROOT) ds9 library ogf_core.tcl]
source [file join $::env(UTI_ROOT) ds9 library ogf_ui.tcl]
proc OGFPython {} {return $::env(UTI_PYTHON)}
proc OGFSessRoot {} {return $::env(UTI_ROOT)}
proc tk_chooseDirectory {args} {return $::env(UTI_EXAMPLE)}
proc tk_getSaveFile {args} {return [file join $::env(UTI_TMP) report.json]}
proc OGFTextWindow {title text} {puts $text}
OGFUIUserToolReview
puts "REGISTERED=[dict exists $::ogf::plugins user_print_smoke]"
'''
            result = subprocess.run(['tclsh'], input=script, text=True, env=env, capture_output=True, check=True)
            self.assertEqual(result.stderr, '')
            report = json.loads((Path(folder) / 'report.json').read_text())
            self.assertEqual(report['status'], 'static_pass_needs_review')
            self.assertFalse(report['accepted_for_execution'])
            self.assertIn('REGISTERED=0', result.stdout)
            self.assertNotIn('ASTRAFEX_USER_PRINT_OK', result.stdout)
        print('PRINT AUDIT STANDALONE: original print works; menu static check passed; automatic registration absent; script not executed by Import')

    def test_review(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--review', str(EXAMPLE)], capture_output=True, text=True, check=True)
        report = json.loads(result.stdout)
        self.assertEqual(report['status'], 'static_pass_needs_review')
        self.assertFalse(report['accepted_for_execution'])

    def test_kit(self):
        with tempfile.TemporaryDirectory() as folder:
            dest = Path(folder) / 'kit.zip'
            subprocess.run([sys.executable, str(SCRIPT), '--kit', str(dest)], check=True, capture_output=True)
            with zipfile.ZipFile(dest) as archive:
                self.assertIsNone(archive.testzip())
                self.assertEqual(len(archive.namelist()), 9)
                self.assertEqual(json.loads(archive.read('examples/user_tool_submission/submission.json'))['target'], 'standalone')
            result = subprocess.run([sys.executable, str(SCRIPT), '--kit', str(dest)], capture_output=True)
            self.assertNotEqual(result.returncode, 0)  # never overwrite existing archive

    def test_validator_copy_matches_web(self):
        web = ROOT.parent / 'ds10-web/tools/validate_user_tool.py'
        if not web.exists():
            self.skipTest('Sibling Web checkout unavailable')
        self.assertEqual(web.read_bytes(), (ROOT / 'tools/user_tool_validation.py').read_bytes())

    def test_tcl_menu_commands(self):
        with tempfile.TemporaryDirectory(prefix='user tool ') as folder:
            env = dict(os.environ, UTI_ROOT=str(ROOT), UTI_PYTHON=sys.executable, UTI_EXAMPLE=str(EXAMPLE), UTI_TMP=folder)
            script = r'''
source [file join $::env(UTI_ROOT) ds9 library ogf_ui.tcl]
proc OGFPython {} {return $::env(UTI_PYTHON)}
proc OGFSessRoot {} {return $::env(UTI_ROOT)}
proc tk_chooseDirectory {args} {return $::env(UTI_EXAMPLE)}
proc tk_getSaveFile {args} {
    if {[lsearch -exact $args .zip] >= 0} {return [file join $::env(UTI_TMP) kit.zip]}
    return [file join $::env(UTI_TMP) report.json]
}
proc OGFTextWindow {title text} {puts $title}
OGFUIUserToolKit
OGFUIUserToolReview
'''
            result = subprocess.run(['tclsh'], input=script, text=True, env=env, capture_output=True, check=True)
            self.assertEqual(result.stderr, '')
            self.assertTrue((Path(folder) / 'kit.zip').exists())
            report = json.loads((Path(folder) / 'report.json').read_text())
            self.assertFalse(report['accepted_for_execution'])


if __name__ == '__main__':
    unittest.main()
