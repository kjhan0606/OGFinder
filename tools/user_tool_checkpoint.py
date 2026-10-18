"""Trusted checkpoint bridge; reconstruct personal nodes without loading saved commands."""
import base64
import hashlib
import json
from pathlib import Path
import re
import tempfile
from user_tool_registry import registry_path, identifier, register
from user_tool_validation import validate, read_member


def definition(base, nodes):
    result = json.loads(json.dumps(base))
    steps = result['steps']
    for item in nodes:
        if not re.fullmatch(r'personal_[0-9]+', item['id']) or any(s['id'] == item['id'] for s in steps):
            raise ValueError('Invalid personal node id')
        pid = identifier(item['tool'])
        if not re.fullmatch(r'[a-f0-9]{64}', item['digest']): raise ValueError('Invalid digest')
        index = next(i for i, s in enumerate(steps) if s['id'] == item['after'])
        steps.insert(index + 1, dict(id=item['id'], title=pid + ' · 사용자 단계', detail='Personal native Python; local user privileges.',
                                    action=[pid, 'run'], optional=False, personal_digest=item['digest']))
    result['edges'] = [[a['id'], b['id']] for a,b in zip(steps, steps[1:])]
    return result


def capture(state, outputs):
    packages = {}
    for node in state.get('personal_steps', []):
        pid = identifier(node['tool']);digest = node['digest']
        if not re.fullmatch(r'[a-f0-9]{64}', digest): raise ValueError('Invalid digest')
        source = registry_path() / pid / digest
        report = validate(source)
        if report['errors'] or report['submission_sha256'] != digest: raise ValueError('Personal package changed')
        raw = read_member(source, 'submission.json');m = json.loads(raw)
        packages[digest] = {name: base64.b64encode(read_member(source, name)).decode() for name in ['submission.json', *m['files']]}
    state['personal_packages'] = packages
    state['personal_input_copies'] = {}
    # Native session records enumerate the report; expand its explicitly declared artifacts and copied inputs.
    for filename in list(outputs):
        path = Path(filename)
        if path.name != 'personal-result.json': continue
        report = json.loads(path.read_text())
        for artifact in report.get('artifacts', []):
            candidate = Path(artifact['path']).resolve()
            if candidate.parent != path.resolve().parent or hashlib.sha256(candidate.read_bytes()).hexdigest() != artifact['sha256']:
                raise ValueError('Personal output changed')
            outputs.append(str(candidate))
        request_path = Path(report['run_directory']) / 'request.json'
        if not request_path.resolve().is_relative_to(path.resolve().parent): raise ValueError('Invalid run directory')
        request = json.loads(request_path.read_text())
        for name, filename in request['inputs'].items():
            if name in request.get('original_inputs', {}):
                state['personal_input_copies'][request['original_inputs'][name]] = filename
            candidate = Path(filename).resolve()
            if not candidate.is_relative_to(request_path.parent): raise ValueError('Invalid copied input')
            outputs.append(str(candidate))
    return outputs


def recover(state):
    # Verify every package before any registration; restoring never executes adapters.
    with tempfile.TemporaryDirectory(prefix='personal-restore-') as tmp:
        sources = []
        for index, (digest, files) in enumerate(state.get('personal_packages', {}).items()):
            source = Path(tmp) / str(index);source.mkdir()
            if not re.fullmatch(r'[a-f0-9]{64}', digest): raise ValueError('Invalid package digest')
            total = 0
            for name, encoded in files.items():
                path = Path(name)
                if path.is_absolute() or '..' in path.parts or '\\' in name: raise ValueError('Invalid package member')
                raw = base64.b64decode(encoded, validate=True);total += len(raw)
                if total > 16*1024*1024: raise ValueError('Package too large')
                target = source / path;target.parent.mkdir(parents=True, exist_ok=True);target.write_bytes(raw)
            report = validate(source)
            if report['errors'] or report['submission_sha256'] != digest: raise ValueError('Saved package damaged')
            sources.append(source)
        for node in state.get('personal_steps', []):
            if node['digest'] not in state.get('personal_packages', {}): raise ValueError('Missing saved personal package')
            manifest = json.loads(base64.b64decode(state['personal_packages'][node['digest']]['submission.json']))
            if manifest['id'] != node['tool']: raise ValueError('Saved tool mismatch')
        for source in sources: register(source, registry_path())

def remap_inputs(data):
    state = data['state']
    mapping = {a['source']: a['resolved_path'] for a in data['artifacts']}
    copies = state.get('personal_input_copies', {})
    for pid, params in state.get('parameters', {}).items():
        if not pid.startswith('user_') or 'inputs-json' not in params: continue
        inputs = json.loads(params['inputs-json'])
        params['inputs-json'] = json.dumps({key: mapping.get(copies.get(path, path), path) for key,path in inputs.items()})
