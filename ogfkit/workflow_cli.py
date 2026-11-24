"""Run supported desktop science-workflow steps from a Unix shell.

The desktop recipe and plugin manifests are the source of truth. Steps that
need GUI hooks, catalog merging, or manual actions are reported as unsupported
instead of silently pretending to execute a complete pipeline.
"""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from ogfkit.cliexpand import build_argv, defaults, find_step, load_manifests
from ogfkit.science_workflow import recipes, validate_catalog

ROOT = Path(__file__).resolve().parents[1]


def inspect_step(step, manifests):
    action = step.get('action') or []
    record = {'id': step['id'], 'title': step['title'], 'action': action,
              'optional': bool(step.get('optional')), 'status': 'manual_or_unavailable'}
    if len(action) != 2 or action[0] not in manifests:
        return record
    try:
        definition = find_step(manifests[action[0]], action[1])
    except KeyError:
        record['status'] = 'missing_plugin_step'
        return record
    record.update(plugin=action[0], plugin_step=action[1],
                  needs=definition.get('needs', []), output=definition.get('output', {}),
                  session=definition.get('session'), network=bool(definition.get('network')))
    if action == ['icl', 'mask']:
        # Explicitly create a fresh mask with the same engine as the GUI's
        # automatic mask. Existing manual edits are deliberately not implied.
        record.update(status='headless_adapter', adapter='icl_auto_mask', needs=['image'],
                      output={'mode': 'capture', 'files': ['mask.fits', 'mask_bool.fits']})
        return record
    if action == ['icl', 'background']:
        record.update(status='headless_adapter', adapter='icl_background', needs=['image', 'mask'],
                      output={'mode': 'capture', 'files': ['icl_background.fits',
                                                            'icl_bgsub.fits', 'background_mask.fits']})
        return record
    if action == ['icl', 'profile']:
        # The GUI after-hook only updates display state. The batch adapter
        # requires an explicit mask and centre and stores a stable product name.
        record.update(status='headless_adapter', adapter='icl_profile', needs=['image', 'mask'],
                      output={'mode': 'capture', 'files': ['icl_profile.tsv']},
                      cli_template=definition['cli'])
        return record
    if not definition.get('cli'):
        record['status'] = 'no_cli_template'
    elif definition.get('before') or definition.get('after'):
        record['status'] = 'requires_gui_hooks'
    elif definition.get('output', {}).get('mode') in ('add_columns', 'append', 'merge'):
        record['status'] = 'requires_catalog_merge'
    elif not set(definition.get('needs', [])).issubset({'image', 'catalog'}):
        record['status'] = 'requires_gui_context'
    elif re.search(r'\{(?:cat:[^}]+|mask|psf)\}', json.dumps(definition['cli'])):
        # cliexpand silently substitutes empty GUI catalog/PSF/mask context.
        record['status'] = 'requires_gui_context'
    else:
        record['status'] = 'direct_cli'
        record['cli_template'] = definition['cli']
    return record


def inventory():
    catalog = recipes()
    validate_catalog(catalog, ROOT)
    manifests = load_manifests(str(ROOT))
    return [{'id': r['id'], 'title': r['title'], 'inputs': r['inputs'],
             'outputs': r['outputs'], 'steps': [inspect_step(s, manifests) for s in r['steps']]}
            for r in catalog]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def execution_fingerprint(tasks, python):
    """Fingerprint in-tree workflow code, known shared cores and environment.

    This is an integrity guard, not an assertion of complete reproducibility of
    external binaries, shared libraries or site-specific environment state.
    """
    executable = shutil.which(str(python or sys.executable))
    if not executable:
        raise ValueError('Python interpreter not found')
    executable = str(Path(executable).resolve())
    inspect = ("import importlib.metadata as m,json,sys; "
               "print(json.dumps({'executable':sys.executable,'version':sys.version,"
               "'packages':sorted((str(x.metadata.get('Name','')),x.version) "
               "for x in m.distributions())},sort_keys=True))")
    try:
        env_result = subprocess.run([executable, '-c', inspect], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise ValueError('cannot inspect Python execution environment') from exc
    locations = [ROOT / name for name in ('ogfkit', 'icl', 'lsbg', 'ogfmeas', 'plugins', 'ds9/library')]
    core_candidates = {ROOT.parent / 'ds10-web/server/ds10core', Path('/workspace/ds10-web/server/ds10core')}
    for var in ('DS10_CORE', 'ASTRAFEX_CORE'):
        if os.environ.get(var):
            core_candidates.add(Path(os.environ[var]).expanduser() / 'ds10core')
    config = Path.home() / '.ds9/ogf_params/ds10core.json'
    config_hash = digest(config) if config.is_file() else None
    if config_hash:
        try:
            config_core = json.loads(config.read_text()).get('core-dir')
            if config_core:
                core_candidates.add(Path(config_core).expanduser() / 'ds10core')
        except (ValueError, AttributeError):
            pass  # The desktop discovery code ignores malformed config too.
    for task in tasks:
        core = task.get('params', {}).get('core-dir')
        if core:
            core_candidates.add(Path(core).expanduser() / 'ds10core')
    locations.extend(sorted(core_candidates, key=str))
    source_files = {}
    for folder in locations:
        if not folder.is_dir():
            continue
        for path in sorted(folder.rglob('*')):
            if path.is_file() and path.suffix in ('.py', '.tcl', '.json', '.so') and '__pycache__' not in path.parts:
                source_files[str(path.resolve())] = digest(path)
    if any(task.get('recipe') == 'icl' and task.get('step') == 's23b_evidence' for task in tasks):
        report = ROOT / 'analysis/icl_wgl_membership/hsc_debug/build_s23b_icl_wl_evidence_report.py'
        source_files[str(report)] = digest(report)
    if any(task.get('recipe') == 'icl' and task.get('step') == 's16a_calibration' for task in tasks):
        report = ROOT / 'analysis/icl_wgl_membership/hsc_debug/build_s16a_calibration_evidence.py'
        source_files[str(report)] = digest(report)
    if not source_files:
        raise ValueError('no workflow code found to fingerprint')
    env = {k: os.environ.get(k) for k in ('DS10_CORE', 'ASTRAFEX_CORE', 'PYTHONPATH', 'PATH',
             'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}
    data = {'python_sha256': digest(executable), 'python_environment': json.loads(env_result.stdout),
            'environment': env, 'core_config_sha256': config_hash, 'sources': source_files}
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run_step(recipe_id, step_id, work, image=None, catalog=None, params=None, python=None, mask=None):
    entries = {r['id']: r for r in inventory()}
    if recipe_id not in entries:
        raise ValueError(f'unknown recipe: {recipe_id}')
    steps = {s['id']: s for s in entries[recipe_id]['steps']}
    if step_id not in steps:
        raise ValueError(f'unknown step: {recipe_id}/{step_id}')
    step = steps[step_id]
    if step['status'] not in ('direct_cli', 'headless_adapter'):
        raise ValueError(f'step is not directly executable: {step["status"]}')
    if step['network']:
        raise ValueError('network steps require a dedicated adapter')
    inputs = {'image': image, 'catalog': catalog, 'mask': mask}
    for need in step['needs']:
        if not inputs.get(need) or not Path(inputs[need]).is_file():
            raise ValueError(f'{need} must name an existing file')
    paths = {k: str(Path(v).resolve()) if v else '' for k, v in inputs.items()}
    manifests = load_manifests(str(ROOT))
    params = params or {}
    if not isinstance(params, dict):
        raise ValueError('parameters must be a JSON object')
    allowed = set(defaults(manifests[step['plugin']]))
    unknown = set(params) - allowed
    if unknown:
        raise ValueError('unknown plugin parameters: ' + ', '.join(sorted(unknown)))
    if step.get('adapter') == 'icl_profile':
        try:
            center = [float(params[k]) for k in ('center-x', 'center-y')]
        except (KeyError, TypeError, ValueError):
            raise ValueError('ICL profile requires explicit center-x and center-y')
        if not all(math.isfinite(v) for v in center):
            raise ValueError('ICL profile center must be finite')
    if step.get('adapter') == 'icl_background' and params.get('bkg-method') not in ('polynomial', 'chebyshev', 'sep_large'):
        raise ValueError('ICL background requires an explicit bkg-method')
    work = Path(work).resolve()
    if work.exists():
        raise FileExistsError(f'refusing to overwrite work directory: {work}')
    context = {'root': str(ROOT), 'work': str(work), 'python': python or sys.executable, **paths}
    if step.get('adapter') == 'icl_auto_mask':
        effective = defaults(manifests['icl'])
        effective.update(params)
        argv = [python or sys.executable, str(ROOT / 'ds9/library/ds9_mask.py'), paths['image'],
                '--mode', 'auto', '--mask', str(work / 'mask.fits'), '--fresh',
                '--minarea', '5', '--mag-threshold', '99', '--n-workers', '1']
        for name in ('detect-thresh', 'expand-factor', 'max-dilate-radius',
                     'bright-star-mag-limit', 'bright-star-radius-scale', 'mag-zeropoint'):
            argv.extend(['--' + name, str(effective[name])])
    elif step.get('adapter') == 'icl_background':
        effective = defaults(manifests['icl'])
        effective.update(params)
        argv = [python or sys.executable, str(ROOT / 'ds9/library/ds9_icl.py'), paths['image'],
                '--mode', 'background', '--mask', paths['mask'],
                '--bkg-output', str(work / 'icl_background.fits'),
                '--bgsub-output', str(work / 'icl_bgsub.fits')]
        for name in ('bkg-method', 'bkg-order', 'bkg-sigma-clip', 'bkg-sep-mesh'):
            argv.extend(['--' + name, str(effective[name])])
        iterative = str(effective['bkg-iterative']).lower() not in ('', '0', 'false', 'no', 'off')
        if iterative:
            argv.extend(['--iterative', '--mask-output', str(work / 'background_mask.fits')])
            for name in ('interp-method', 'bkg-n-iterations', 'bkg-convergence-tol', 'bkg-refine-thresh'):
                argv.extend(['--' + name, str(effective[name])])
    else:
        argv = build_argv(manifests, step['plugin'], step['plugin_step'], context, params)
    if step.get('adapter') == 'icl_profile':
        argv[argv.index('--profile-output') + 1] = str(work / 'icl_profile.tsv')
        argv.extend(['--mask', paths['mask']])
    input_hashes = {k: digest(v) for k, v in paths.items() if v and Path(v).is_file()}
    work.mkdir(parents=True, exist_ok=False)
    stdout_path, stderr_path = work / '.astrafek.stdout.log', work / '.astrafek.stderr.log'
    error = ''
    try:
        with stdout_path.open('xb') as stdout, stderr_path.open('xb') as stderr:
            completed = subprocess.run(argv, cwd=work, stdout=stdout, stderr=stderr, check=False)
        exit_code = completed.returncode
        if step.get('adapter') == 'icl_background' and exit_code == 0 and not iterative:
            shutil.copyfile(paths['mask'], work / 'background_mask.fits')
    except OSError as exc:
        exit_code = 127
        error = f'{type(exc).__name__}: {exc}'
    declared = list(step['output'].get('files', []))
    if step['output'].get('catalog'):
        declared.append(step['output']['catalog'])
    products, missing = [], []
    for name in dict.fromkeys(declared):
        path = (work / name).resolve()
        if not path.is_relative_to(work):
            raise ValueError('declared output escapes the work directory')
        if path.is_file():
            products.append({'path': str(path), 'sha256': digest(path), 'bytes': path.stat().st_size})
        else:
            missing.append(name)
    status = 'failed' if exit_code != 0 else ('incomplete_outputs' if missing else 'completed')
    receipt = {'schema': 'astrafek-workflow-step-receipt-v1', 'recipe': recipe_id,
               'step': step_id, 'action': step['action'], 'status': status,
               'exit_code': exit_code, 'error': error, 'argv': argv,
               'input_sha256': input_hashes,
               'params': params, 'outputs': products, 'missing_declared_outputs': missing,
               'warnings': (['Fresh automatic mask can exclude BCG and diffuse light; inspect before science use.']
                            if step.get('adapter') == 'icl_auto_mask' else
                            ['Local background fitting can remove ICL; validate sky transfer before science use.']
                            if step.get('adapter') == 'icl_background' else []),
               'stdout': str(stdout_path), 'stderr': str(stderr_path),
               'science_validated': False}
    (work / 'step_receipt.json').write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + '\n')
    return receipt


def validate_plan(plan, stop_after=None):
    """Preflight selected offline steps, input references and parameters."""
    if not isinstance(plan, dict) or plan.get('schema') != 'astrafek-workflow-plan-v1':
        raise ValueError('expected astrafek-workflow-plan-v1')
    tasks = plan.get('steps')
    if not isinstance(tasks, list) or not tasks:
        raise ValueError('plan steps must be a non-empty list')
    recipes_by_id = {r['id']: r for r in inventory()}
    manifests = load_manifests(str(ROOT))
    known = {}
    for task in tasks:
        if not isinstance(task, dict):
            raise ValueError('each plan step must be an object')
        name = task.get('id')
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', name) or name in known:
            raise ValueError('step IDs must be unique, simple names')
        recipe_id, step_id = task.get('recipe'), task.get('step')
        recipe = recipes_by_id.get(recipe_id)
        step = next((s for s in recipe['steps'] if s['id'] == step_id), None) if recipe else None
        if not step or step['status'] not in ('direct_cli', 'headless_adapter') or step['network']:
            raise ValueError(f'plan step {name} is not a supported offline CLI action')
        if not isinstance(task.get('params', {}), dict):
            raise ValueError(f'plan step {name}: params must be an object')
        if set(task.get('params', {})) - set(defaults(manifests[step['plugin']])):
            raise ValueError(f'plan step {name}: unknown plugin parameters')
        for key in ('image', 'catalog', 'mask'):
            value = task.get(key)
            if key in step['needs'] and value is None:
                raise ValueError(f'plan step {name}: missing {key}')
            if isinstance(value, dict):
                previous, filename = value.get('step'), value.get('output')
                if (set(value) != {'step', 'output'} or not isinstance(previous, str)
                        or not isinstance(filename, str) or previous not in known
                        or filename not in known[previous]):
                    raise ValueError(f'plan step {name}: invalid prior output reference')
            elif value is not None and (not isinstance(value, str) or not Path(value).is_file()):
                raise ValueError(f'plan step {name}: {key} must name an existing file')
        output = step['output']
        known[name] = set(output.get('files', [])) | ({output['catalog']} if output.get('catalog') else set())

    if stop_after is not None and stop_after not in known:
        raise ValueError('stop-after must name a step in the plan')
    return tasks


def run_plan(plan, work, python=None, resume=False, stop_after=None):
    """Execute an explicit sequence; omitted/manual stages remain unperformed."""
    tasks = validate_plan(plan, stop_after)

    work = Path(work).resolve()
    if resume:
        if not work.is_dir():
            raise ValueError('resume requires an existing work directory')
    else:
        work.mkdir(parents=True, exist_ok=False)
    lock = (work / '.plan.lock').open('a')
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise ValueError('another process is running this plan')
    try:
        return _run_plan_locked(plan, tasks, work, python, resume, stop_after)
    finally:
        lock.close()


def export_python(plan, destination):
    """Write a standalone launcher for this exact reviewed plan; never overwrite."""
    validate_plan(plan)
    payload = json.dumps(plan, ensure_ascii=False, sort_keys=True)
    script = '''#!/usr/bin/env python3
"""Astr afek selected workflow plan; omitted/manual stages are not performed."""
import argparse
import json
from pathlib import Path
import sys

PLAN = json.loads(%s)
DEFAULT_OGFINDER = %s

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ogfinder-root', type=Path, default=Path(DEFAULT_OGFINDER))
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--stop-after')
    args = parser.parse_args()
    sys.path.insert(0, str(args.ogfinder_root.resolve()))
    from ogfkit.workflow_cli import run_plan
    result = run_plan(PLAN, args.work_dir, args.python, args.resume, args.stop_after)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result['status'] in ('completed', 'paused') else 1

if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print('workflow launcher: ' + str(exc), file=sys.stderr)
        sys.exit(2)
''' % (repr(payload), repr(str(ROOT)))
    destination = Path(destination)
    with destination.open('x', encoding='utf-8') as stream:
        stream.write(script)
    return {'status': 'exported', 'path': str(destination.resolve()), 'sha256': digest(destination),
            'plan': plan, 'science_validated': False}


def _run_plan_locked(plan, tasks, work, python, resume, stop_after):
    plan_file = work / 'plan.json'
    summary_path = work / 'plan_result.json'
    plan_text = json.dumps(plan, indent=2, ensure_ascii=False) + '\n'
    fingerprint = execution_fingerprint(tasks, python)
    if resume:
        if not plan_file.is_file() or not summary_path.is_file() or plan_file.read_text() != plan_text:
            raise ValueError('saved plan is missing or differs')
        summary = json.loads(summary_path.read_text())
        if (summary.get('schema') != 'astrafek-workflow-plan-result-v2'
                or summary.get('plan_sha256') != digest(plan_file)
                or summary.get('execution_fingerprint') != fingerprint):
            raise ValueError('saved plan, code or execution environment differs')
        if summary.get('status') not in ('paused', 'running', 'completed'):
            raise ValueError('failed plan cannot be resumed; preserve it and use a new work directory')
        prior = summary.get('steps', [])
        if len(prior) > len(tasks) or any(x.get('id') != task['id'] or x.get('status') != 'completed'
                                          for x, task in zip(prior, tasks)):
            raise ValueError('completed-step history is inconsistent with the plan')
    else:
        plan_file.write_text(plan_text)
        summary = {'schema': 'astrafek-workflow-plan-result-v2', 'plan_sha256': digest(plan_file),
                   'execution_fingerprint': fingerprint, 'status': 'running', 'steps': [],
                   'science_validated': False}
    if resume and stop_after in {entry['id'] for entry in summary['steps']}:
        raise ValueError('stop-after step is already complete')

    def save():
        tmp = work / '.plan_result.json.tmp'
        tmp.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
        os.replace(tmp, summary_path)

    if not resume:
        save()
    completed = {}
    prior_count = len(summary['steps'])
    for index, task in enumerate(tasks):
        name = task['id']
        inputs = {}
        for key in ('image', 'catalog', 'mask'):
            value = task.get(key)
            if isinstance(value, dict):
                previous = completed[value['step']]
                path = (work / value['step'] / value['output']).resolve()
                product = next((p for p in previous['outputs'] if p['path'] == str(path)), None)
                if not product or not path.is_file() or digest(path) != product.get('sha256'):
                    summary.update(status='failed', error=f'{name}: prior output missing or changed')
                    save()
                    return summary
                inputs[key] = str(path)
            else:
                inputs[key] = value
        if index < prior_count:
            saved = summary['steps'][index]
            receipt_path = work / name / 'step_receipt.json'
            if (saved.get('receipt') != str(receipt_path) or not receipt_path.is_file()
                    or digest(receipt_path) != saved.get('receipt_sha256')):
                raise ValueError(f'{name}: saved receipt missing or changed')
            receipt = json.loads(receipt_path.read_text())
            if (receipt.get('status') != 'completed' or receipt.get('recipe') != task['recipe']
                    or receipt.get('step') != task['step'] or receipt.get('params') != task.get('params', {})):
                raise ValueError(f'{name}: saved receipt disagrees with plan')
            hashes = {k: digest(v) for k, v in inputs.items() if v and Path(v).is_file()}
            if hashes != receipt.get('input_sha256'):
                raise ValueError(f'{name}: input missing or changed')
            for product in receipt.get('outputs', []):
                path = Path(product['path'])
                if (not path.resolve().is_relative_to((work / name).resolve()) or not path.is_file()
                        or digest(path) != product.get('sha256') or path.stat().st_size != product.get('bytes')):
                    raise ValueError(f'{name}: saved output missing or changed')
            completed[name] = receipt
            continue
        summary['status'] = 'running'
        save()
        try:
            receipt = run_step(task['recipe'], task['step'], work / name,
                               params=task.get('params'), python=python, **inputs)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            summary.update(status='failed', error=f'{name}: {type(exc).__name__}: {exc}')
            save()
            return summary
        receipt_path = work / name / 'step_receipt.json'
        summary['steps'].append({'id': name, 'status': receipt['status'],
                                 'receipt': str(receipt_path), 'receipt_sha256': digest(receipt_path)})
        save()
        if receipt['status'] != 'completed':
            summary.update(status='failed', error=f'{name}: {receipt["status"]}')
            save()
            return summary
        completed[name] = receipt
        if stop_after == name and index < len(tasks) - 1:
            summary['status'] = 'paused'
            save()
            return summary
    summary['status'] = 'completed'
    save()
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('list', help='show all recipes and direct CLI availability')
    show = sub.add_parser('show', help='show one recipe and step mappings')
    show.add_argument('recipe')
    run = sub.add_parser('run-step', help='run one direct CLI stage in a new work directory')
    run.add_argument('recipe'); run.add_argument('step'); run.add_argument('--work-dir', type=Path, required=True)
    run.add_argument('--image'); run.add_argument('--catalog'); run.add_argument('--mask'); run.add_argument('--params-json', type=Path)
    run.add_argument('--python', default=sys.executable)
    batch = sub.add_parser('run-plan', help='run an explicit sequence of offline CLI steps')
    batch.add_argument('plan', type=Path)
    batch.add_argument('--work-dir', type=Path, required=True)
    batch.add_argument('--python', default=sys.executable)
    batch.add_argument('--resume', action='store_true', help='verify and continue an unchanged paused plan')
    batch.add_argument('--stop-after', help='save after the named step for later review')
    export = sub.add_parser('export-python', help='write a Python launcher for an explicit plan')
    export.add_argument('plan', type=Path)
    export.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'list':
        result = inventory()
    elif args.command == 'show':
        result = next((r for r in inventory() if r['id'] == args.recipe), None)
        if result is None:
            raise ValueError(f'unknown recipe: {args.recipe}')
    elif args.command == 'run-step':
        params = json.loads(args.params_json.read_text()) if args.params_json else {}
        result = run_step(args.recipe, args.step, args.work_dir, args.image, args.catalog, params, args.python, args.mask)
    elif args.command == 'run-plan':
        result = run_plan(json.loads(args.plan.read_text()), args.work_dir, args.python, args.resume, args.stop_after)
    else:
        result = export_python(json.loads(args.plan.read_text()), args.output)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 1 if args.command in ('run-step', 'run-plan') and result['status'] not in ('completed', 'paused') else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f'workflow-cli: {exc}', file=sys.stderr)
        sys.exit(2)
