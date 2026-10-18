import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT))
from user_tool_registry import register, run
from user_tool_checkpoint import capture, recover, definition
from ogfkit.science_workflow import recipes, checkpoint, restore

class PersonalPipelineTests(unittest.TestCase):
    def test_register_run_checkpoint_restore_and_integrity(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp);registry = root / 'registry'
            old = os.environ.get('ASTRAFEX_USER_TOOLS');os.environ['ASTRAFEX_USER_TOOLS'] = str(registry)
            try:
                info = register(ROOT / 'examples/user_tool_print', registry)
                plugin = json.loads(Path(info['plugin']).read_text())
                self.assertEqual(plugin['steps'][0]['cli'][2:4], ['--run','user_print_smoke'])
                work = root / 'work';result = run(info['id'], registry, work, {}, {})
                self.assertEqual(result['status'], 'succeeded')
                self.assertIn('ASTRAFEX_USER_PRINT_OK', (work / 'personal-log.txt').read_text())
                with self.assertRaises(ValueError): run(info['id'], registry, work, {}, {})
                nodes = [dict(id='personal_123', after='detect', tool=info['id'], digest=info['digest'])]
                state = dict(recipe='lsbg',personal_steps=nodes)
                state['recipe_definition'] = definition(next(r for r in recipes() if r['id']=='lsbg'), nodes)
                output = capture(state, [str(work/'personal-result.json'),str(work/'personal-log.txt')])
                manifest = checkpoint(root/'checkpoint', state, output)
                data = restore(manifest, root/'resumed');recover(data['state'])
                self.assertEqual(data['state']['recipe_definition']['steps'][5]['id'], 'personal_123')
                package = registry / info['id'] / info['digest']
                (package/'adapter.py').write_text('raise RuntimeError("tampered")')
                with self.assertRaises(ValueError): run(info['id'], registry, root/'tamper', {}, {})
            finally:
                if old is None: os.environ.pop('ASTRAFEX_USER_TOOLS',None)
                else: os.environ['ASTRAFEX_USER_TOOLS'] = old

    def test_native_input_and_artifact(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'source';source.mkdir()
            template=ROOT/'tools/user_tool_kit/examples/user_tool_submission'
            m=json.loads((template/'submission.json').read_text());m['target']='standalone'
            for name in m['files']: (source/name).write_bytes((template/name).read_bytes())
            (source/'submission.json').write_text(json.dumps(m))
            register(source,root/'registry')
            incoming=root/'values.json';incoming.write_text('[1,2,3]')
            result=run(m['id'],root/'registry',root/'work',{'values':str(incoming)}, {})
            self.assertEqual(result['status'],'succeeded')
            self.assertEqual(json.loads((root/'work/result.json').read_text())['mean'],2.0)
            old=os.environ.get('ASTRAFEX_USER_TOOLS');os.environ['ASTRAFEX_USER_TOOLS']=str(root/'registry')
            try:
                plugin=json.loads((root/'registry'/m['id']/'plugin.json').read_text())
                state={'personal_steps':[{'id':'personal_123','after':'detect','tool':m['id'],'digest':plugin['personal_digest']}], 'parameters':{m['id']:{'inputs-json':json.dumps({'values':str(incoming)})}}}
                outputs=capture(state,[str(root/'work/personal-result.json')])
                checkpoint_file=checkpoint(root/'checkpoint',state,outputs)
                incoming.unlink()
                data=restore(checkpoint_file,root/'resumed')
                from user_tool_checkpoint import remap_inputs
                remap_inputs(data)
                params=json.loads(data['state']['parameters'][m['id']]['inputs-json'])
                self.assertEqual(Path(params['values']).read_text(),'[1,2,3]')
                self.assertEqual(run(m['id'],root/'registry',root/'again',params,{})['status'],'succeeded')
            finally:
                if old is None: os.environ.pop('ASTRAFEX_USER_TOOLS',None)
                else: os.environ['ASTRAFEX_USER_TOOLS']=old

    def test_timeout_and_native_exception(self):
        for code in ('while True: pass', 'raise ValueError("EXPECTED_FAILURE")'):
            with tempfile.TemporaryDirectory() as temp:
                root=Path(temp);source=root/'source';source.mkdir()
                template=ROOT/'examples/user_tool_print'
                m=json.loads((template/'submission.json').read_text());m['limits']['seconds']=1
                for name in m['files']: (source/name).write_bytes((template/name).read_bytes())
                (source/'adapter.py').write_text('def run(inputs, params, output_dir):\n    '+code+'\n')
                m['files']['adapter.py']=hashlib.sha256((source/'adapter.py').read_bytes()).hexdigest()
                (source/'submission.json').write_text(json.dumps(m))
                register(source,root/'registry')
                result=run(m['id'],root/'registry',root/'work',{}, {})
                self.assertEqual(result['status'],'failed')
                self.assertFalse(result['artifacts'])

    def test_tcl_registration_and_inline_attachment(self):
        with tempfile.TemporaryDirectory() as temp:
            env=dict(os.environ,UTI_ROOT=str(ROOT),UTI_PYTHON=sys.executable,ASTRAFEX_USER_TOOLS=temp)
            script=r'''
source [file join $::env(UTI_ROOT) ds9 library ogf_json.tcl]
source [file join $::env(UTI_ROOT) ds9 library ogf_core.tcl]
source [file join $::env(UTI_ROOT) ds9 library ogf_ui.tcl]
source [file join $::env(UTI_ROOT) plugins workflow_designer workflow_designer.tcl]
source [file join $::env(UTI_ROOT) ds9 library ogf_session.tcl]
proc OGFPython {} {return $::env(UTI_PYTHON)}
proc OGFSessRoot {} {return $::env(UTI_ROOT)}
proc tk_chooseDirectory {args} {return [file join $::env(UTI_ROOT) examples user_tool_print]}
proc OGFTextWindow {args} {puts $args}
proc ::ogf::workflow::draw {} {}
proc ::ogf::job::busy {} {return 0}
OGFUIUserToolRegister
set ::ogf::workflow::recipe [dict create id lsbg steps [list [dict create id detect]]]
set ::ogf::workflow::selected detect
::ogf::workflow::attach_personal user_print_smoke
puts [::ogf::workflow::personal_json]
puts "REGISTERED=[dict exists $::ogf::plugins user_print_smoke]"
'''
            result=subprocess.run(['tclsh'],input=script,text=True,env=env,capture_output=True,check=True)
            self.assertEqual(result.stderr,'', result.stdout)
            self.assertIn('REGISTERED=1',result.stdout)
            self.assertIn('"after": "detect"',result.stdout)

if __name__=='__main__': unittest.main()
