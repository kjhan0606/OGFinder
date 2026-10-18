"""Personal native tools. Explicit local execution, never an upload or safety certificate."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from user_tool_validation import validate, read_member


def registry_path(value=None):
    return Path(value or os.environ.get('ASTRAFEX_USER_TOOLS', Path.home() / '.local/share/astrafex/user-tools')).expanduser().resolve()


def identifier(value):
    if not re.fullmatch(r'user_[a-z0-9_]+', value):
        raise ValueError('Invalid user tool id')
    return value


def register(source, registry):
    source = Path(source).resolve()
    report = validate(source)
    if report['errors']:
        raise ValueError(str(report['errors']))
    raw = read_member(source, 'submission.json')
    m = json.loads(raw)
    if m['target'] not in ('standalone', 'both'):
        raise ValueError('Use target standalone or both for native registration')
    pid = identifier(m['id'])
    digest = hashlib.sha256(raw).hexdigest()
    target = registry / pid
    if target.is_symlink():
        raise ValueError('Registry symlink refused')
    target.mkdir(parents=True, exist_ok=True)
    package = target / digest
    if not package.exists():
        stage = Path(tempfile.mkdtemp(prefix='staging-', dir=target))
        try:
            for name in ['submission.json', *m['files']]:
                data = read_member(source, name)
                dest = stage / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
            checked = validate(stage)
            if checked['errors']:
                raise ValueError('Source changed during registration')
            os.replace(stage, package)
        finally:
            if stage.exists(): shutil.rmtree(stage)
    plugin = dict(schema=1, id=pid, name=pid + ' (personal local)', tab='Measure', order=990,
                  description=m['purpose'] + ' — personal native Python; runs with your local user permissions.',
                  license='user-provided, not certified', primary='run', personal_digest=digest,
                  steps=[dict(id='run', label='Run personal tool locally', stage='measure', session='AUTO', needs=[],
                              output=dict(mode='capture', files=['personal-result.json', 'personal-log.txt']),
                              cli=['{python}', str(Path(__file__).with_name('user_tool_import.py')), '--run', pid,
                                   '--registry', str(registry), '--digest', digest, '--work', '{work}',
                                   '--inputs-json', '{inputs-json}', '--params-json', '{params-json}'])],
                  params=[dict(name='inputs-json', label='Input names to local paths (JSON)', type='string', default='{}'),
                          dict(name='params-json', label='Parameters (JSON)', type='string', default='{}')])
    temp = target / 'plugin.json.tmp'
    temp.write_text(json.dumps(plugin, ensure_ascii=False, indent=2))
    os.replace(temp, target / 'plugin.json')
    return dict(id=pid, digest=digest, plugin=str(target / 'plugin.json'), execution='standalone-native', accepted_for_execution=False)


def run(pid, registry, work, inputs, params, digest=None):
    pid = identifier(pid)
    if not digest:
        digest = json.loads((registry / pid / 'plugin.json').read_text())['personal_digest']
    if not re.fullmatch(r'[a-f0-9]{64}', digest):
        raise ValueError('Invalid registered version')
    package = registry / pid / digest
    report = validate(package)
    if report['errors'] or report.get('submission_sha256') != digest:
        raise ValueError('Registered package integrity failed')
    m = json.loads((package / 'submission.json').read_text())
    if not isinstance(inputs, dict) or not isinstance(params, dict) or any(k not in m['inputs'] for k in inputs) or any(k not in m['parameters'] for k in params):
        raise ValueError('Invalid input or parameter names')
    if any(spec.get('required', True) and name not in inputs for name, spec in m['inputs'].items()):
        raise ValueError('Required input missing')
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    if any((work / name).exists() for name in ('personal-result.json', 'personal-log.txt', *m['outputs'])):
        raise ValueError('Use a new work directory for every run')
    run_dir = Path(tempfile.mkdtemp(prefix=pid + '-', dir=work))
    incoming = run_dir / 'inputs';incoming.mkdir()
    output = run_dir / 'output';output.mkdir()
    mapped = {};total = 0
    for index, (name, filename) in enumerate(inputs.items()):
        path = Path(filename).resolve()
        if not path.is_file(): raise ValueError('Input is not a file')
        total += path.stat().st_size
        if total > 128 * 1024 * 1024: raise ValueError('Inputs exceed 128 MiB')
        dest = incoming / (str(index) + path.suffix)
        shutil.copyfile(path, dest);mapped[name] = str(dest)
    request = dict(package=str(package), output=str(output), inputs=mapped, original_inputs=inputs, params=params, manifest=m)
    request_path = run_dir / 'request.json';request_path.write_text(json.dumps(request))
    log_path = run_dir / 'log.txt'
    env = {'PATH': os.defpath, 'HOME': str(run_dir), 'PYTHONNOUSERSITE': '1', 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
    seconds = min(m['limits']['seconds'], 600)
    status = 'failed'
    class Cancelled(Exception): pass
    def cancel(_sig, _frame): raise Cancelled()
    old_handler = signal.signal(signal.SIGTERM, cancel)
    with log_path.open('wb') as log:
        process = subprocess.Popen([sys.executable, '-I', str(Path(__file__).with_name('user_tool_native_worker.py')), str(request_path)], cwd=run_dir, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=seconds)
            if code == 0: status = 'succeeded'
        except Cancelled:
            status = 'cancelled';code = -15
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL);process.wait();code = -9
        finally:
            # Clean up the child process group; this is not an OS security sandbox.
            try: os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            process.wait()
            signal.signal(signal.SIGTERM, old_handler)
    text = log_path.read_bytes()[:1024*1024].decode('utf-8', errors='replace')
    if code == -9: text += '\nTime limit exceeded'
    artifacts = []
    result_path = output / 'result-index.json'
    if status == 'succeeded':
        try:
            names = json.loads(result_path.read_text())['artifacts'];size = 0
            if not isinstance(names, list) or len(set(names)) != len(names): raise ValueError('Invalid artifact list')
            for name in names:
                if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', name) or name not in m['outputs'] or name in ('personal-result.json', 'personal-log.txt', 'result-index.json'): raise ValueError('Invalid artifact')
                path = output / name
                if path.is_symlink() or not path.is_file(): raise ValueError('Invalid artifact file')
                size += path.stat().st_size
                if size > min(m['limits']['output_mib'], 64)*1024*1024: raise ValueError('Output limit exceeded')
                dest = work / name
                if dest.exists(): raise ValueError('Refusing to overwrite output: ' + name)
            for name in names:
                path = output / name;dest = work / name
                shutil.copyfile(path, dest)
                artifacts.append(dict(name=name, path=str(dest), sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
        except Exception as exc:
            status='failed';text += '\n' + str(exc);artifacts=[]
    result = dict(id=pid, digest=digest, status=status, execution='standalone-native', science_validated=False,
                  artifacts=artifacts, run_directory=str(run_dir), time=datetime.now(timezone.utc).isoformat())
    for name, data in [('personal-result.json', json.dumps(result, indent=2)), ('personal-log.txt', text)]:
        path = work / name
        if path.exists(): raise ValueError('Use a new work directory for every run')
        path.write_text(data)
    return result
