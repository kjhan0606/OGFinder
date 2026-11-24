import json
import struct
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from ogfkit.workflow_cli import inventory, run_step


def test_inventory_matches_recipe_actions():
    data = inventory()
    assert len(data) == 11
    icl = next(r for r in data if r['id'] == 'icl')
    steps = {s['id']: s for s in icl['steps']}
    assert steps['load']['status'] == 'manual_or_unavailable'
    assert steps['audit']['status'] == 'direct_cli'
    assert steps['audit']['action'] == ['ds10core', 'workflow-audit']
    assert steps['validate']['status'] == 'manual_or_unavailable'
    replay = steps['s23b_evidence']
    assert replay['optional'] and not replay['network']
    assert replay['status'] == 'direct_cli'
    assert replay['needs'] == []
    assert replay['output']['files'] == ['s23b_icl_evidence/summary.json', 's23b_icl_evidence/assessment.md']


def test_audit_runs_and_records_inputs_outputs(tmp_path):
    core = ROOT.parent / 'ds10-web/server'
    if not (core / 'ds10core/tools/workflow_support.py').exists():
        pytest.skip('shared core checkout is unavailable')
    image = tmp_path / 'image.fits'
    cards = ['SIMPLE  =                    T', 'BITPIX  =                  -32', 'NAXIS   =                    2', 'NAXIS1  =                    2', 'NAXIS2  =                    2', 'END']
    image.write_bytes(''.join(c.ljust(80) for c in cards).ljust(2880).encode('ascii') + struct.pack('>4f', 1, 2, 3, 4) + bytes(2864))
    work = tmp_path / 'audit'
    receipt = run_step('icl', 'audit', work, image=str(image), params={'core-dir': str(core)})
    assert receipt['status'] == 'completed', (work / '.astrafek.stderr.log').read_text()
    assert json.loads((work / 'audit.json').read_text())['mean'] == 2.5
    assert receipt['input_sha256']['image']
    assert receipt['outputs'][0]['sha256']
    assert receipt['science_validated'] is False
    with pytest.raises(FileExistsError):
        run_step('icl', 'audit', work, image=str(image), params={'core-dir': str(core)})


def test_missing_input_does_not_create_work(tmp_path):
    work = tmp_path / 'not-created'
    with pytest.raises(ValueError, match='existing file'):
        run_step('icl', 'audit', work, image=str(tmp_path / 'missing.fits'))
    assert not work.exists()


def test_plan_executes_two_audits_and_keeps_separate_receipts(tmp_path):
    from ogfkit.workflow_cli import run_plan
    core = ROOT.parent / 'ds10-web/server'
    if not (core / 'ds10core/tools/workflow_support.py').exists():
        pytest.skip('shared core checkout is unavailable')
    image = tmp_path / 'image.fits'
    cards = ['SIMPLE  =                    T', 'BITPIX  =                  -32', 'NAXIS   =                    2', 'NAXIS1  =                    2', 'NAXIS2  =                    2', 'END']
    image.write_bytes(''.join(c.ljust(80) for c in cards).ljust(2880).encode('ascii') + struct.pack('>4f', 1, 2, 3, 4) + bytes(2864))
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': name, 'recipe': 'icl', 'step': 'audit', 'image': str(image), 'params': {'core-dir': str(core)}}
        for name in ('first', 'second')]}
    work = tmp_path / 'batch'
    result = run_plan(plan, work)
    assert result['status'] == 'completed'
    assert [step['status'] for step in result['steps']] == ['completed', 'completed']
    assert json.loads((work / 'first/audit.json').read_text())['mean'] == 2.5
    assert json.loads((work / 'second/step_receipt.json').read_text())['input_sha256']['image']
    assert result['science_validated'] is False
    with pytest.raises(FileExistsError):
        run_plan(plan, work)


def test_plan_preflights_all_steps_before_work(tmp_path):
    from ogfkit.workflow_cli import run_plan
    image = tmp_path / 'input.fits'
    image.touch()
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': 'audit', 'recipe': 'icl', 'step': 'audit', 'image': str(image)},
        {'id': 'bad', 'recipe': 'icl', 'step': 'audit', 'image': {'step': 'audit', 'output': 'unlisted.fits'}}]}
    work = tmp_path / 'batch'
    with pytest.raises(ValueError, match='prior output'):
        run_plan(plan, work)
    assert not work.exists()


def test_plan_connects_prior_declared_product_and_stops_on_failure(tmp_path, monkeypatch):
    import ogfkit.workflow_cli as cli
    image = tmp_path / 'input.fits'
    image.touch()
    seen = []
    def fake_step(recipe, step, work, image=None, catalog=None, params=None, python=None, mask=None):
        seen.append(image)
        work.mkdir()
        product = work / 'audit.json'
        product.write_text('{}')
        result = {'status': 'completed' if len(seen) == 1 else 'failed',
                  'outputs': [{'path': str(product), 'sha256': cli.digest(product)}]}
        (work / 'step_receipt.json').write_text(json.dumps(result))
        return result
    monkeypatch.setattr(cli, 'run_step', fake_step)
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': 'first', 'recipe': 'icl', 'step': 'audit', 'image': str(image)},
        {'id': 'second', 'recipe': 'icl', 'step': 'audit', 'image': {'step': 'first', 'output': 'audit.json'}},
        {'id': 'third', 'recipe': 'icl', 'step': 'audit', 'image': str(image)}]}
    work = tmp_path / 'batch'
    result = cli.run_plan(plan, work)
    assert result['status'] == 'failed'
    assert seen == [str(image), str(work / 'first/audit.json')]
    assert not (work / 'third').exists()
    assert json.loads((work / 'plan_result.json').read_text())['status'] == 'failed'


def test_plan_real_spectrum_audit_to_continuum(tmp_path):
    from ogfkit.workflow_cli import run_plan
    core = ROOT.parent / 'ds10-web/server'
    if not (core / 'ds10core/tools/spectrum_prepare.py').exists():
        pytest.skip('shared core checkout is unavailable')
    source = tmp_path / 'spectrum.tsv'
    source.write_text('WAVE_A\tFLUX\tFLUX_ERR\n' + ''.join(f'{4900+i}\t{100+.01*i}\t1\n' for i in range(101)))
    params = {'core-dir': str(core), 'spectrum-wave-frame': 'vacuum',
              'spectrum-flux-unit': 'counts', 'spectrum-degree': 1}
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': 'audit', 'recipe': 'stellar_spectrum', 'step': 'audit', 'catalog': str(source), 'params': params},
        {'id': 'continuum', 'recipe': 'stellar_spectrum', 'step': 'continuum',
         'catalog': {'step': 'audit', 'output': 'spectrum_audit_input.tsv'}, 'params': params}]}
    work = tmp_path / 'spectral-plan'
    result = run_plan(plan, work)
    assert result['status'] == 'completed', result
    assert (work / 'continuum/continuum_input.tsv').read_bytes() == source.read_bytes()
    assert json.loads((work / 'continuum/continuum.json').read_text())['n_fit'] == 101
    assert json.loads((work / 'continuum/step_receipt.json').read_text())['science_validated'] is False


def test_icl_profile_batch_excludes_masked_pixels_and_requires_center(tmp_path):
    import csv
    import numpy as np
    from astropy.io import fits
    image = tmp_path / 'image.fits'
    mask = tmp_path / 'mask.fits'
    data = np.full((21, 21), 10.0, dtype=np.float32)
    data[10, 11] = 1000
    masked = np.zeros(data.shape, dtype=np.uint8)
    masked[10, 11] = 1
    fits.writeto(image, data)
    fits.writeto(mask, masked)
    params = {'center-x': 10, 'center-y': 10, 'rmin': 2,
              'rmax': 8, 'nsteps': 3, 'spacing': 'linear', 'pixel-scale': 1}
    work = tmp_path / 'profile'
    receipt = run_step('icl', 'profile', work, image=str(image), mask=str(mask), params=params)
    assert receipt['status'] == 'completed', (work / '.astrafek.stderr.log').read_text()
    with (work / 'icl_profile.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert int(rows[0]['NPIX']) == 8
    assert float(rows[0]['FLUX']) == pytest.approx(80)
    assert receipt['input_sha256']['mask']
    assert receipt['science_validated'] is False
    with pytest.raises(ValueError, match='center'):
        run_step('icl', 'profile', tmp_path / 'missing-center', image=str(image), mask=str(mask))


def test_icl_profile_rejects_wrong_shape_mask(tmp_path):
    import numpy as np
    from astropy.io import fits
    image, mask = tmp_path / 'image.fits', tmp_path / 'mask.fits'
    fits.writeto(image, np.ones((21, 21)))
    fits.writeto(mask, np.zeros((10, 10)))
    receipt = run_step('icl', 'profile', tmp_path / 'profile', image=str(image), mask=str(mask),
                       params={'center-x': 10, 'center-y': 10})
    assert receipt['status'] == 'failed'
    assert not receipt['outputs']


@pytest.mark.parametrize('selected_frame', ['source', 'bgsub'])
def test_icl_gui_passes_original_image_mask_when_using_background_subtracted_image(tmp_path, selected_frame):
    import os
    import subprocess
    image = tmp_path / 'source.fits'
    bgsub = tmp_path / 'background_subtracted.fits'
    mask = tmp_path / 'source_mask.fits'
    for path in (image, bgsub, mask):
        path.touch()
    script = r'''
namespace eval ::ogf::cat {
    variable v
    proc get {key} {variable v; return $v($key)}
    proc set {key value} {variable v; set v($key) $value}
}
source plugins/icl/icl.tcl
rename CatalogPanelICLUpdateFiles oldUpdateFiles
proc CatalogPanelICLUpdateFiles {fn} {
    if {$fn ne $::env(TEST_IMAGE)} {error "resolved paths against the wrong image"}
    ::ogf::cat::set icl,mask_file $::env(TEST_MASK)
    ::ogf::cat::set icl,bgsub_file $::env(TEST_BGSUB)
    ::ogf::cat::set icl,profile_file /tmp/profile-test.tsv
}
proc CatalogPanelGetFITS {} {return $::env(TEST_ACTIVE)}
set ::ogf::icl_bgsub_sources([file normalize $::env(TEST_BGSUB)]) [list $::env(TEST_IMAGE) $::env(TEST_MASK)]
proc CatalogPanelGetScript {x} {return ds9/library/ds9_icl.py}
proc OGFPython {} {return python3}
proc CatalogPanelCmdLog {args} {}
proc CatalogPanelLoadTSV {args} {}
proc CatalogPanelICLDrawAnnuli {} {}
rename exec original_exec
proc exec {args} {set ::cmd $args; return "header\nvalue"}
foreach {key value} {center_x 10 center_y 10 param,rmin 1 param,rmax 10 param,nsteps 3 param,spacing linear param,ellipticity 0 param,pa 0 param,mag-zeropoint 25 param,pixel-scale 1} {
    ::ogf::cat::set icl,$key $value
}
CatalogPanelICLProfile
if {[lsearch -exact $::cmd $::env(TEST_BGSUB)] < 0} {error "bgsub missing"}
set i [lsearch -exact $::cmd --mask]
if {$i < 0 || [lindex $::cmd [expr {$i + 1}]] ne $::env(TEST_MASK)} {error "original mask missing"}
puts OK
'''
    env = {**os.environ, 'TEST_IMAGE': str(image), 'TEST_BGSUB': str(bgsub), 'TEST_MASK': str(mask),
           'TEST_ACTIVE': str(image if selected_frame == 'source' else bgsub)}
    result = subprocess.run(['tclsh'], input=script, text=True, capture_output=True, cwd=ROOT, env=env)
    assert result.returncode == 0 and result.stdout.strip() == 'OK', result.stderr


def test_headless_background_records_refined_mask_for_result_frame(tmp_path):
    import os
    import subprocess
    image = tmp_path / 'source.fits'
    bgsub = tmp_path / 'icl_bgsub_source.fits'
    mask = tmp_path / 'mask_source_bool.fits'
    refined = tmp_path / 'mask_source_iclrefined.fits'
    for path in (image, bgsub, mask, refined):
        path.touch()
    script = r'''
namespace eval ::ogf::cat {
    variable v
    proc get {key} {variable v; return $v($key)}
    proc set {key value} {variable v; set v($key) $value}
}
namespace eval ::ogf::params {proc get {plugin key} {if {$key eq "bkg-iterative"} {return 1}; return polynomial}}
source plugins/icl/icl.tcl
rename CatalogPanelICLUpdateFiles oldUpdateFiles
proc CatalogPanelICLUpdateFiles {fn} {::ogf::cat::set icl,mask_file $::env(TEST_MASK)}
proc CatalogPanelGetFITS {} {return $::env(TEST_IMAGE)}
proc OGFMaskRefinedPath {id} {return $::env(TEST_REFINED)}
proc OGFSessWorkDir {} {return $::env(TEST_DIR)}
proc CatalogPanelFitsBaseName {fn} {return source}
proc CreateFrame {} {}
proc LoadFitsFile {args} {}
proc ChangeScaleMode {} {}
OGFIclBkgHeadlessAfter
set entry $::ogf::icl_bgsub_sources([file normalize $::env(TEST_BGSUB)])
if {[lindex $entry 0] ne $::env(TEST_IMAGE) || [lindex $entry 1] ne $::env(TEST_REFINED)} {error "source/refined mask missing"}
puts OK
'''
    env = {**os.environ, 'TEST_IMAGE': str(image), 'TEST_BGSUB': str(bgsub), 'TEST_MASK': str(mask),
           'TEST_REFINED': str(refined), 'TEST_DIR': str(tmp_path)}
    result = subprocess.run(['tclsh'], input=script, text=True, capture_output=True, cwd=ROOT, env=env)
    assert result.returncode == 0 and result.stdout.strip() == 'OK', result.stderr


def test_plan_icl_auto_mask_feeds_real_profile(tmp_path):
    import numpy as np
    from astropy.io import fits
    from ogfkit.workflow_cli import run_plan
    rng = np.random.default_rng(12)
    yy, xx = np.mgrid[:96, :96]
    data = 100 + rng.normal(0, 1, (96, 96)) + 900*np.exp(-((xx-18)**2+(yy-18)**2)/8)
    image = tmp_path / 'sky_checked.fits'
    fits.writeto(image, data.astype('f4'))
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': 'mask', 'recipe': 'icl', 'step': 'mask', 'image': str(image),
         'params': {'mag-zeropoint': 25}},
        {'id': 'profile', 'recipe': 'icl', 'step': 'profile', 'image': str(image),
         'mask': {'step': 'mask', 'output': 'mask_bool.fits'},
         'params': {'center-x': 55, 'center-y': 55, 'rmin': 4, 'rmax': 20,
                    'nsteps': 4, 'spacing': 'linear', 'pixel-scale': 1}}]}
    work = tmp_path / 'pipeline'
    result = run_plan(plan, work, stop_after='mask')
    assert result['status'] == 'paused' and not (work/'profile').exists()
    saved_receipt = work/'mask/step_receipt.json'
    saved_mtime = saved_receipt.stat().st_mtime_ns
    result = run_plan(plan, work, resume=True)
    assert saved_receipt.stat().st_mtime_ns == saved_mtime
    assert result['status'] == 'completed', result
    mask = fits.getdata(work/'mask/mask_bool.fits')
    assert mask[18, 18] != 0 and mask[55, 55] == 0
    assert (work/'profile/icl_profile.tsv').is_file()
    profile = json.loads((work/'profile/step_receipt.json').read_text())
    mask_receipt = json.loads((work/'mask/step_receipt.json').read_text())
    assert profile['input_sha256']['mask'] == next(p['sha256'] for p in mask_receipt['outputs'] if p['path'].endswith('mask_bool.fits'))
    assert 'BCG' in mask_receipt['warnings'][0]
    assert result['science_validated'] is False


@pytest.mark.parametrize('iterative', [False, True])
def test_plan_icl_background_chains_image_and_effective_mask(tmp_path, iterative):
    import numpy as np
    from astropy.io import fits
    from ogfkit.workflow_cli import run_plan
    rng = np.random.default_rng(41)
    yy, xx = np.mgrid[:64, :64]
    data = 100 + .01*xx + .02*yy + rng.normal(0, .1, (64, 64))
    data += 20*np.exp(-((xx-32)**2+(yy-32)**2)/100)
    image, mask = tmp_path/'image.fits', tmp_path/'mask.fits'
    fits.writeto(image, data.astype('f4'))
    m = np.zeros((64, 64), dtype='u1'); m[5:9, 5:9] = 1
    fits.writeto(mask, m)
    params = {'bkg-method': 'polynomial', 'bkg-order': 1,
              'bkg-iterative': int(iterative), 'bkg-n-iterations': 1,
              'interp-method': 'nearest', 'bkg-refine-thresh': 1000000}
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': 'sky', 'recipe': 'icl', 'step': 'sky', 'image': str(image), 'mask': str(mask), 'params': params},
        {'id': 'profile', 'recipe': 'icl', 'step': 'profile',
         'image': {'step': 'sky', 'output': 'icl_bgsub.fits'},
         'mask': {'step': 'sky', 'output': 'background_mask.fits'},
         'params': {'center-x': 32, 'center-y': 32, 'rmin': 2, 'rmax': 12,
                    'nsteps': 4, 'spacing': 'linear', 'pixel-scale': 1}}]}
    work = tmp_path / 'pipeline'
    result = run_plan(plan, work, stop_after='sky')
    assert result['status'] == 'paused' and not (work/'profile').exists()
    sky_receipt = work/'sky/step_receipt.json'
    sky_mtime = sky_receipt.stat().st_mtime_ns
    result = run_plan(plan, work, resume=True)
    assert sky_receipt.stat().st_mtime_ns == sky_mtime
    assert result['status'] == 'completed', result
    background = fits.getdata(work/'sky/icl_background.fits')
    bgsub = fits.getdata(work/'sky/icl_bgsub.fits')
    assert np.allclose(fits.getdata(image)-background, bgsub, atol=1e-5)
    assert np.all(fits.getdata(work/'sky/background_mask.fits')[5:9, 5:9] != 0)
    receipt = json.loads((work/'sky/step_receipt.json').read_text())
    assert 'remove ICL' in receipt['warnings'][0]
    assert (work/'profile/icl_profile.tsv').is_file()


def test_iterative_background_full_mask_stops_profile(tmp_path):
    import numpy as np
    from astropy.io import fits
    from ogfkit.workflow_cli import run_plan
    yy, xx = np.mgrid[:64, :64]
    rng = np.random.default_rng(41)
    data = 100 + rng.normal(0, .1, (64, 64)) + 20*np.exp(-((xx-32)**2+(yy-32)**2)/100)
    image, mask = tmp_path/'image.fits', tmp_path/'mask.fits'
    fits.writeto(image, data.astype('f4'))
    fits.writeto(mask, np.zeros((64, 64), dtype='u1'))
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': 'sky', 'recipe': 'icl', 'step': 'sky', 'image': str(image), 'mask': str(mask),
         'params': {'bkg-method': 'polynomial', 'bkg-order': 1, 'bkg-iterative': 1,
                    'bkg-n-iterations': 1, 'interp-method': 'nearest'}},
        {'id': 'profile', 'recipe': 'icl', 'step': 'profile',
         'image': {'step': 'sky', 'output': 'icl_bgsub.fits'},
         'mask': {'step': 'sky', 'output': 'background_mask.fits'},
         'params': {'center-x': 32, 'center-y': 32, 'rmin': 2, 'rmax': 12,
                    'nsteps': 4, 'spacing': 'linear', 'pixel-scale': 1}}]}
    work = tmp_path/'pipeline'
    result = run_plan(plan, work)
    assert result['status'] == 'failed' and result['error'] == 'profile: failed'
    assert fits.getdata(work/'sky/background_mask.fits')[32, 32] != 0
    assert 'No valid profile points' in (work/'profile/.astrafek.stderr.log').read_text()


def test_icl_background_requires_method_and_rejects_wrong_shape_mask(tmp_path):
    import numpy as np
    from astropy.io import fits
    image, mask = tmp_path/'image.fits', tmp_path/'mask.fits'
    fits.writeto(image, np.ones((32, 32)))
    fits.writeto(mask, np.zeros((8, 8)))
    with pytest.raises(ValueError, match='explicit bkg-method'):
        run_step('icl', 'sky', tmp_path/'no-method', image=str(image), mask=str(mask))
    result = run_step('icl', 'sky', tmp_path/'bad-mask', image=str(image), mask=str(mask),
                      params={'bkg-method': 'polynomial'})
    assert result['status'] == 'failed'


def test_plan_pause_resume_checks_saved_outputs_and_inputs(tmp_path):
    import fcntl
    from ogfkit.workflow_cli import run_plan
    core = ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/workflow_support.py').exists():
        pytest.skip('shared core checkout unavailable')
    image = tmp_path/'image.fits'
    cards = ['SIMPLE  =                    T', 'BITPIX  =                  -32', 'NAXIS   =                    2', 'NAXIS1  =                    2', 'NAXIS2  =                    2', 'END']
    image.write_bytes(''.join(c.ljust(80) for c in cards).ljust(2880).encode('ascii') + struct.pack('>4f', 1, 2, 3, 4) + bytes(2864))
    params = {'core-dir': str(core)}
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': name, 'recipe': 'icl', 'step': 'audit', 'image': str(image), 'params': params}
        for name in ('first', 'second')]}
    work = tmp_path/'pipeline'
    paused = run_plan(plan, work, stop_after='first')
    assert paused['status'] == 'paused' and not (work/'second').exists()
    with (work/'.plan.lock').open('a') as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(ValueError, match='another process'):
            run_plan(plan, work, resume=True)
    receipt = work/'first/step_receipt.json'
    original_time = receipt.stat().st_mtime_ns
    resumed = run_plan(plan, work, resume=True)
    assert resumed['status'] == 'completed' and receipt.stat().st_mtime_ns == original_time
    assert (work/'second/audit.json').is_file()
    (work/'first/audit.json').write_text('{}')
    with pytest.raises(ValueError, match='output missing or changed'):
        run_plan(plan, work, resume=True)
    assert json.loads((work/'plan_result.json').read_text())['status'] == 'completed'


def test_plan_resume_rejects_changed_input_or_code(tmp_path, monkeypatch):
    import ogfkit.workflow_cli as cli
    core = ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/workflow_support.py').exists():
        pytest.skip('shared core checkout unavailable')
    image = tmp_path/'image.fits'
    cards = ['SIMPLE  =                    T', 'BITPIX  =                  -32', 'NAXIS   =                    2', 'NAXIS1  =                    2', 'NAXIS2  =                    2', 'END']
    image.write_bytes(''.join(c.ljust(80) for c in cards).ljust(2880).encode('ascii') + struct.pack('>4f', 1, 2, 3, 4) + bytes(2864))
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': name, 'recipe': 'icl', 'step': 'audit', 'image': str(image), 'params': {'core-dir': str(core)}}
        for name in ('first', 'second')]}
    work = tmp_path/'pipeline'
    cli.run_plan(plan, work, stop_after='first')
    original = image.read_bytes()
    image.write_bytes(original + b'changed')
    with pytest.raises(ValueError, match='input missing or changed'):
        cli.run_plan(plan, work, resume=True)
    assert json.loads((work/'plan_result.json').read_text())['status'] == 'paused'
    image.write_bytes(original)
    monkeypatch.setattr(cli, 'execution_fingerprint', lambda *a: 'different-code')
    with pytest.raises(ValueError, match='code or execution environment differs'):
        cli.run_plan(plan, work, resume=True)
    assert not (work/'second').exists()


def test_export_python_runs_and_resumes_exact_plan_from_other_directory(tmp_path):
    import subprocess
    from ogfkit.workflow_cli import export_python
    core = ROOT.parent/'ds10-web/server'
    if not (core/'ds10core/tools/workflow_support.py').exists():
        pytest.skip('shared core checkout unavailable')
    image = tmp_path/'image.fits'
    cards = ['SIMPLE  =                    T', 'BITPIX  =                  -32', 'NAXIS   =                    2', 'NAXIS1  =                    2', 'NAXIS2  =                    2', 'END']
    image.write_bytes(''.join(c.ljust(80) for c in cards).ljust(2880).encode('ascii') + struct.pack('>4f', 1, 2, 3, 4) + bytes(2864))
    plan = {'schema': 'astrafek-workflow-plan-v1', 'steps': [
        {'id': name, 'recipe': 'icl', 'step': 'audit', 'image': str(image), 'params': {'core-dir': str(core)}}
        for name in ('first', 'second')]}
    script, work = tmp_path/'workflow.py', tmp_path/'output'
    receipt = export_python(plan, script)
    assert receipt['science_validated'] is False and receipt['sha256']
    with pytest.raises(FileExistsError):
        export_python(plan, script)
    cmd = [sys.executable, str(script), '--work-dir', str(work)]
    p1 = subprocess.run(cmd + ['--stop-after', 'first'], cwd=tmp_path, text=True, capture_output=True)
    assert p1.returncode == 0, p1.stderr
    assert json.loads(p1.stdout)['status'] == 'paused'
    p2 = subprocess.run(cmd + ['--resume'], cwd=tmp_path, text=True, capture_output=True)
    assert p2.returncode == 0, p2.stderr
    assert json.loads(p2.stdout)['status'] == 'completed'
    assert json.loads((work/'second/audit.json').read_text())['mean'] == 2.5
