import json
from pathlib import Path
import sys
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from ogfkit.science_workflow import checkpoint, recipes, restore, validate_catalog


def test_recipe_actions_exist_and_both_products_use_same_data():
    validate_catalog(recipes(), ROOT)
    desktop = json.loads((ROOT / 'plugins/workflow_designer/recipes.json').read_text())
    assert desktop['recipes'] == recipes()
    web = ROOT.parent / 'ds10-web/server/ds10web/science_recipes.json'
    if web.exists():
        assert json.loads(web.read_text()) == desktop
    assert len(recipes()) == 11


def test_checkpoint_preserves_outputs_and_resume_uses_working_copies(tmp_path):
    source = tmp_path / 'input.fits'
    source.write_bytes(b'input')
    output = tmp_path / 'result.tsv'
    output.write_text('flux\n5\n')
    state = {'selected': 'fit', 'progress': {'fit': 'reviewed'}}
    manifest = checkpoint(tmp_path / 'sessions', state, [output], [source])
    output.unlink()
    resumed = restore(manifest, tmp_path / 'resume')
    assert resumed['state'] == state
    assert not resumed['input_changes']
    Path(resumed['artifacts'][0]['resolved_path']).write_text('working edit')
    assert restore(manifest)['artifacts'][0]['sha256'] == resumed['artifacts'][0]['sha256']
    source.unlink()
    assert restore(manifest)['input_changes'][0]['reason'] == 'missing'


def test_tampered_output_is_refused(tmp_path):
    output = tmp_path / 'r'
    output.write_bytes(b'ab')
    manifest = checkpoint(tmp_path / 'saved', {}, [output])
    entry = json.loads(manifest.read_text())['artifacts'][0]
    (manifest.parent / entry['path']).write_bytes(b'cd')
    with pytest.raises(ValueError, match='damaged'):
        restore(manifest)


def test_failure_keeps_previous_snapshot_and_no_partial_success(tmp_path):
    dest = tmp_path / 'saved'
    manifest = checkpoint(dest, {'selected': 'load'})
    with pytest.raises(FileNotFoundError):
        checkpoint(dest, {}, [tmp_path / 'missing'])
    assert manifest.exists()
    assert len(list(dest.iterdir())) == 1
    assert restore(manifest)['state']['selected'] == 'load'


def test_snapshot_escape_is_rejected(tmp_path):
    outside = tmp_path / 'outside'
    outside.write_bytes(b'x')
    manifest = checkpoint(tmp_path / 'saved', {}, [outside])
    data = json.loads(manifest.read_text())
    data['artifacts'][0]['path'] = str(outside)
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='path'):
        restore(manifest)


def test_writer_lock_is_not_removed_by_second_writer(tmp_path):
    dest = tmp_path / 'saved'
    dest.mkdir()
    (dest / '.saving').write_text('owner')
    with pytest.raises(FileExistsError):
        checkpoint(dest, {})
    assert (dest / '.saving').read_text() == 'owner'


def test_desktop_tracks_actual_records_and_resume_does_not_restart_processes():
    script = r'''
source ds9/library/ogf_json.tcl
source plugins/workflow_designer/workflow_designer.tcl
proc check {actual expected} {if {$actual ne $expected} {puts stderr "Expected $expected, got $actual"; exit 1}}
proc record {seq ms failed note} {dict create seq $seq ms $ms failed $failed note $note undone 0 argv {python tool.py}}
set ok [record 1 20 0 {}]
set pending [record 2 {} 0 {}]
check [::ogf::workflow::execution_status {} {}] opened
check [::ogf::workflow::execution_status [list $ok $pending] 2] running
check [::ogf::workflow::execution_status [list $ok $pending] {}] opened
check [::ogf::workflow::execution_status [list $ok] {}] succeeded
check [::ogf::workflow::execution_status [list [record 3 4 1 error]] {}] failed
check [::ogf::workflow::execution_status [list [record 4 4 1 {cancelled by the user}]] {}] cancelled
namespace eval ::ogf::session {}
namespace eval ::ogf::job {variable cur; array set cur {seq 2 busy 1}}
proc ::ogf::session::steps {} {return $::records}
proc ::ogf::job::busy {} {return $::ogf::job::cur(busy)}
set ::records [list $pending]
set ::ogf::workflow::bindings [dict create detect {2}]
::ogf::workflow::reconcile
check [dict get $::ogf::workflow::progress detect] running
set ::ogf::job::cur(busy) 0
set ::records [list [record 2 20 0 {}]]
::ogf::workflow::reconcile
check [dict get $::ogf::workflow::progress detect] succeeded
dict set ::ogf::workflow::progress detect reviewed
::ogf::workflow::reconcile
check [dict get $::ogf::workflow::progress detect] reviewed
rename ::ogf::workflow::draw ::ogf::workflow::original_draw
proc ::ogf::workflow::draw {} {}
set state [dict create recipe_definition [dict create id icl steps {}] progress [dict create stack running] selected stack notes {} records {} parameters {} catalog_source /no/catalog]
::ogf::workflow::apply_restore [dict create state $state artifacts {} inputs {} input_changes {}]
check [dict get $::ogf::workflow::progress stack] needs_review
check $::ogf::workflow::bindings {}
puts OK
'''
    result = subprocess.run(['tclsh'], input=script, text=True, cwd=ROOT, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert not result.stderr, result.stderr
    assert result.stdout.strip() == 'OK'


def test_declared_reports_enter_desktop_recorder_without_a_popup():
    script = r'''
source ds9/library/ogf_json.tcl
source ds9/library/ogf_core.tcl
source plugins/workflow_designer/workflow_designer.tcl
set fd [open plugins/ds10core/plugin.json]; set ::manifest [::ogf::json::parse [read $fd]]; close $fd
rename ::ogf::reg::get ::ogf::reg::original_get
proc ::ogf::reg::get {id} {return $::manifest}
rename ::ogf::reg::missing ::ogf::reg::original_missing
proc ::ogf::reg::missing {id} {return {}}
rename ::ogf::step::check_needs ::ogf::step::original_check_needs
proc ::ogf::step::check_needs {step} {return {}}
rename ::ogf::step::context ::ogf::step::original_context
proc ::ogf::step::context {id {catname {}}} {return [dict create python python3 plugin_dir /plugins work /work image /input.fits catalog /profiles.tsv core-dir /core workflow-bootstrap 200 workflow-seed 0]}
rename ::ogf::job::run ::ogf::job::original_run
proc ::ogf::job::run {argv args} {set ::options $args; return 1}
proc check {actual expected} {if {$actual ne $expected} {puts stderr "Expected $expected, got $actual"; exit 1}}
check [::ogf::step::run ds10core workflow-audit] 1
check [dict get $::options -outputs] /work/audit.json
check [::ogf::step::run ds10core profile-stack] 1
check [dict get $::options -outputs] {/work/profile_stack.json /work/profile_stack.tsv /work/profile_stack_input.tsv}
rename ::ogf::workflow::keep_note ::ogf::workflow::original_keep_note
proc ::ogf::workflow::keep_note {} {}
rename ::ogf::workflow::save ::ogf::workflow::original_save
proc ::ogf::workflow::save {} {set ::saved 1}
set ::ogf::workflow::recipe [dict create steps [list [dict create id results action checkpoint]]]
set ::ogf::workflow::selected results
::ogf::workflow::launch
check $::saved 1
puts OK
'''
    result = subprocess.run(['tclsh'], input=script, text=True, cwd=ROOT, capture_output=True)
    assert result.returncode == 0 and not result.stderr, result.stderr
    assert result.stdout.strip() == 'OK'


def test_desktop_shared_core_audit_runs_as_a_separate_process(tmp_path):
    import struct
    core = ROOT.parent / 'ds10-web/server'
    if not (core / 'ds10core/tools/workflow_support.py').exists():
        pytest.skip('shared core checkout is unavailable')
    path = tmp_path / 'image.fits'
    cards = ["SIMPLE  =                    T", "BITPIX  =                  -32", "NAXIS   =                    2", "NAXIS1  =                    2", "NAXIS2  =                    2", "END"]
    header = ''.join(card.ljust(80) for card in cards).ljust(2880).encode('ascii')
    path.write_bytes(header + struct.pack('>4f', 1, 2, 3, 4) + bytes(2864))
    report = tmp_path / 'audit.json'
    result = subprocess.run([sys.executable, str(ROOT/'plugins/ds10core/run_workflow.py'), '--core-dir', str(core), 'audit', str(path), '--summary', str(report)], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(report.read_text())['mean'] == 2.5


def test_desktop_spectral_table_load_uses_catalogue_without_fits_or_popup(tmp_path):
    import os
    path=tmp_path/'spectrum.tsv';path.write_text('WAVE_A\tFLUX\tFLUX_ERR\n5000\t100\t1\n5001\t101\t1\n')
    script=r'''
source ds9/library/ogf_json.tcl
source ds9/library/ogf_core.tcl
source plugins/workflow_designer/workflow_designer.tcl
proc tk_getOpenFile {args} {return [list $::env(SPECTRUM_TEST_FILE)]}
rename ::ogf::workflow::keep_note ::ogf::workflow::original_keep_note
proc ::ogf::workflow::keep_note {} {}
rename ::ogf::workflow::draw ::ogf::workflow::original_draw
proc ::ogf::workflow::draw {} {}
rename ::ogf::session::steps ::ogf::session::original_steps
proc ::ogf::session::steps {} {return {}}
rename ::ogf::cat::load_tsv ::ogf::cat::original_load_tsv
proc ::ogf::cat::load_tsv {data name} {::set ::loaded $data}
set ::ogf::workflow::window .mock
set ::ogf::workflow::recipe [dict create id stellar_spectrum steps [list [dict create id load action load]]]
set ::ogf::workflow::selected load
::ogf::workflow::launch
if {![info exists ::loaded] || [string first WAVE_A $::loaded] != 0} {puts stderr "TSV was not loaded: $::ogf::workflow::message";exit 1}
if {$::ogf::workflow::inputs ne [list $::env(SPECTRUM_TEST_FILE)]} {puts stderr "Input reference missing";exit 1}
puts OK
'''
    result=subprocess.run(['tclsh'],input=script,text=True,cwd=ROOT,capture_output=True,env={**os.environ,'SPECTRUM_TEST_FILE':str(path)})
    assert result.returncode==0 and not result.stderr,result.stderr
    assert result.stdout.strip()=='OK'


def test_desktop_spectral_continuum_products_survive_checkpoint(tmp_path):
    core=ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/spectrum_prepare.py').exists():pytest.skip('shared core checkout unavailable')
    source=tmp_path/'spectrum.tsv'
    source.write_text('WAVE_A\tFLUX\tFLUX_ERR\n'+''.join(f'{4900+i}\t{100+.01*i}\t1\n' for i in range(101)))
    out=tmp_path/'continuum.tsv';report=tmp_path/'continuum.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum.py'),'--core-dir',str(core),'continuum',str(source),'--out',str(out),'--summary',str(report),'--wave-frame','vacuum','--flux-unit','counts','--degree','1'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert result.stdout.startswith('WAVE_A\tFLUX\tFLUX_ERR\tCONTINUUM')
    meta=json.loads(report.read_text());assert not meta['science_validated'] and meta['n_fit']==101
    snapshot=tmp_path/'continuum_input.tsv';assert snapshot.read_bytes()==source.read_bytes()
    manifest=checkpoint(tmp_path/'saved',{'recipe':'stellar_spectrum','selected':'continuum'},[out,report,snapshot],[source])
    original=out.read_bytes();out.unlink()
    resumed=restore(manifest,tmp_path/'resume')
    assert resumed['state']['selected']=='continuum'
    assert any(Path(item['resolved_path']).read_bytes()==original for item in resumed['artifacts'])


def test_desktop_isolated_line_products_survive_checkpoint(tmp_path):
    import math
    core=ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/spectrum_lines.py').exists():pytest.skip('shared core checkout unavailable')
    source=tmp_path/'s.tsv'
    source.write_text('WAVE_A\tFLUX\tFLUX_ERR\n'+''.join(f'{4900+i*.5}\t{100+20*math.exp(-.5*((4900+i*.5-5000)/3)**2)}\t1\n' for i in range(401)))
    report=tmp_path/'lines.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_lines.py'),'--core-dir',str(core),'lines',str(source),'--summary',str(report),'--wave-frame','vacuum','--flux-unit','counts/Angstrom','--lines','synthetic:5000'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    meta=json.loads(result.stdout);assert meta['lines'][0]['status']=='accepted'
    assert meta['redshift_candidate'] is None
    files=[tmp_path/name for name in ('lines.json','lines_input.tsv','lines_measurements.tsv','lines_model.tsv')]
    manifest=checkpoint(tmp_path/'saved',{'recipe':'stellar_spectrum','selected':'science'},files,[source])
    resumed=restore(manifest,tmp_path/'resume')
    assert len(resumed['artifacts'])==4
    assert resumed['state']['selected']=='science'


def test_desktop_spectral_search_products_survive_checkpoint(tmp_path):
    import math
    core=ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/spectrum_search.py').exists():pytest.skip('shared core checkout unavailable')
    source=tmp_path/'s.tsv'
    source.write_text('WAVE_A\tFLUX\tFLUX_ERR\n'+''.join(f'{5400+i*.5}\t{20+15*math.exp(-.5*((5400+i*.5-5500)/3)**2)+15*math.exp(-.5*((5400+i*.5-5720)/3)**2)}\t.2\n' for i in range(820)))
    report=tmp_path/'search.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_search.py'),'--core-dir',str(core),str(source),'--summary',str(report),'--wave-frame','vacuum','--flux-unit','counts/Angstrom','--lines','A:5000;B:5200'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    meta=json.loads(result.stdout)
    assert meta['candidates'][0]['redshift']==pytest.approx(.1,abs=.001)
    assert meta['requires_review'] and not meta['science_validated']
    files=[tmp_path/name for name in ('search.json','search_input.tsv','search_peaks.tsv','search_candidates.tsv')]
    manifest=checkpoint(tmp_path/'saved',{'recipe':'galaxy_spectrum','selected':'search'},files,[source])
    resumed=restore(manifest,tmp_path/'resume')
    assert len(resumed['artifacts'])==4
    assert resumed['state']['selected']=='search'


def test_desktop_joint_component_products_survive_checkpoint(tmp_path):
    import math
    core=ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/spectrum_joint.py').exists():pytest.skip('shared core checkout unavailable')
    source=tmp_path/'s.tsv'
    source.write_text('WAVE_A\tFLUX\tFLUX_ERR\n'+''.join(f'{4870+i*.5}\t{50+20*math.exp(-.5*((4870+i*.5-5000)/3)**2)+8*math.exp(-.5*((4870+i*.5-5000)/15)**2)}\t.2\n' for i in range(540)))
    report=tmp_path/'joint.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_joint.py'),'--core-dir',str(core),str(source),'--summary',str(report),'--wave-frame','vacuum','--flux-unit','counts/Angstrom','--components','narrow:5000:1:6;broad:5000:8:25','--fit-range','4880:5130'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    meta=json.loads(result.stdout)
    assert meta['status']=='accepted_model' and len(meta['components'])==2
    files=[tmp_path/name for name in ('joint.json','joint_input.tsv','joint_components.tsv','joint_model.tsv')]
    manifest=checkpoint(tmp_path/'saved',{'recipe':'agn_spectrum','selected':'lines'},files,[source])
    resumed=restore(manifest,tmp_path/'resume')
    assert len(resumed['artifacts'])==4
    assert resumed['state']['selected']=='lines'


def test_desktop_template_decomposition_checkpoint(tmp_path):
    import math
    core=ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/spectrum_templates.py').exists():pytest.skip('shared core checkout unavailable')
    source=tmp_path/'s.tsv';lines=['WAVE_A\tFLUX\tFLUX_ERR\tTEMPLATE_HOST\tTEMPLATE_FEII']
    for i in range(1000):
        w=4500+i;host=1+.4*math.sin((w-4500)/90);iron=math.exp(-.5*((w-4750)/25)**2)
        lines.append(f'{w}\t{30*(w/5000)**-1.5+12*host+7*iron}\t.1\t{host}\t{iron}')
    source.write_text('\n'.join(lines)+'\n');report=tmp_path/'decompose.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_templates.py'),'--core-dir',str(core),str(source),'--summary',str(report),'--wave-frame','vacuum','--flux-unit','counts/Angstrom','--fit-range','4510:5480','--template-note','synthetic observed-grid templates'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    meta=json.loads(result.stdout);assert meta['coefficients']==pytest.approx([30,12,7],rel=.001)
    files=[tmp_path/name for name in ('decompose.json','decompose_input.tsv','decomposition.tsv')]
    saved=checkpoint(tmp_path/'saved',{'recipe':'agn_spectrum','selected':'decompose'},files,[source])
    assert len(restore(saved,tmp_path/'resume')['artifacts'])==3


def test_desktop_joint_diagnostic_report_and_checkpoint(tmp_path):
    core=ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/spectrum_diagnostic.py').exists():pytest.skip('shared core checkout unavailable')
    source=tmp_path/'joint.json'
    source.write_text(json.dumps(dict(tool='spectrum-joint',components=[dict(name=n,integrated_flux=f,integrated_flux_error=e,status='accepted',flags=[]) for n,f,e in [('a',100.,2.),('b',50.,3.)]],component_flux_covariance=[[4,3.6],[3.6,9]])))
    report=tmp_path/'diagnostic.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_diagnostic.py'),'--core-dir',str(core),str(source),'--summary',str(report),'--ratios','ab=a/b'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert json.loads(result.stdout)['ratios'][0]['ratio']==pytest.approx(2)
    files=[tmp_path/name for name in ('diagnostic.json','diagnostic_joint_input.json','diagnostic_ratios.tsv')]
    saved=checkpoint(tmp_path/'saved',{'recipe':'agn_spectrum','selected':'diagnostic'},files,[source])
    assert len(restore(saved,tmp_path/'resume')['artifacts'])==3


def test_desktop_diagnostic_source_resolves_restored_artifact_and_rejects_latest_failure(tmp_path):
    import os
    original=tmp_path/'joint.json';resumed=tmp_path/'0000_joint.json';resumed.write_text('{}')
    rec=dict(step='ds10core.spectrum-joint',seq=5,outputs=[str(original)],ms='100',failed='0',undone='0')
    script=r'''
source ds9/library/ogf_json.tcl
source plugins/workflow_designer/workflow_designer.tcl
namespace eval ::ogf::session {}
proc ::ogf::session::steps {} {return $::records}
set ::records {}
set ::ogf::workflow::recipe [dict create steps [list [dict create id lines action {ds10core spectrum-joint}]]]
set ::ogf::workflow::saved_records [list $::env(SAVED_RECORD)]
set ::ogf::workflow::restored_artifacts [dict create $::env(ORIGINAL) $::env(RESUMED)]
set path [::ogf::workflow::input_artifact lines joint.json]
if {$path ne $::env(RESUMED)} {puts stderr "Restored file was not connected";exit 1}
set ::records [list [dict create step ds10core.spectrum-joint seq 6 outputs {} ms 100 failed 1 undone 0]]
if {![catch {::ogf::workflow::input_artifact lines joint.json}]} {puts stderr "Failed latest execution reused old success";exit 1}
puts OK
'''
    run=subprocess.run(['tclsh'],input=script,text=True,cwd=ROOT,capture_output=True,env={**os.environ,'SAVED_RECORD':json.dumps(rec),'ORIGINAL':str(original),'RESUMED':str(resumed)})
    assert run.returncode==0 and not run.stderr,run.stderr
    assert run.stdout.strip()=='OK'


def test_desktop_bundled_inference_bridge_and_checkpoint(tmp_path):
    import math
    core=ROOT.parent/'ds10-web/server'
    source=tmp_path/'s.tsv'
    source.write_text('WAVE_A\tFLUX\tFLUX_ERR\n'+''.join(f'{6400+i}\t{100+30*math.exp(-.5*((6400+i-6563)/12)**2)}\t.1\n' for i in range(330)))
    fitted=tmp_path/'joint.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_joint.py'),'--core-dir',str(core),str(source),'--summary',str(fitted),'--wave-frame','air','--flux-unit','1e-17 erg/s/cm2/Angstrom','--components','broad:6563:8:20','--fit-range','6410:6720'],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    cfg=dict(extinction=dict(mode='external',ebv=.1,ebv_error=.02,k={'broad':2.5},law_note='Synthetic independent screen'),luminosity=dict(distance_mpc=100,distance_error_mpc=2,distance_note='Synthetic luminosity distance',flux_to_cgs=1e-17,flux_scale_note='Synthetic source numeric unit conversion',components=['broad'],dust_corrected=True),black_hole=dict(component='broad',broad_halpha_confirmed=True,calibration='greene_ho_2005_eq6',applicability_note='Synthetic broad Halpha test',scatter_dex=.4,instrument_fwhm_kms=100,instrument_error_kms=5,resolution_note='Synthetic Gaussian resolution'))
    report=tmp_path/'inference.json'
    result=subprocess.run([sys.executable,str(ROOT/'plugins/ds10core/run_spectrum_inference.py'),'--core-dir',str(core),str(fitted),'--summary',str(report),'--config',json.dumps(cfg)],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert json.loads(result.stdout)['results']['black_hole']['status']=='conditional_mass_estimate'
    files=[tmp_path/name for name in ('inference.json','inference_joint_input.json','inference_config.json','inference_summary.tsv')]
    saved=checkpoint(tmp_path/'saved',{'recipe':'agn_spectrum','selected':'inference','parameters':cfg},files,[fitted])
    resumed=restore(saved,tmp_path/'resume')
    assert len(resumed['artifacts'])==4 and resumed['state']['parameters']==cfg
