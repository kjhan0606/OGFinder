"""Read-only submission preflight. Never imports, installs or executes submitted code."""
import argparse
import ast
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys

SCHEMA = 'astrafex-user-tool/1'
MAX_FILE = 2 * 1024 * 1024
MAX_TOTAL = 16 * 1024 * 1024
SHA = re.compile(r'[a-f0-9]{64}')


def read_member(root, name):
    if not isinstance(name, str) or not name or '\\' in name:
        raise ValueError('Invalid member path')
    parts = name.split('/')
    if PurePosixPath(name).is_absolute() or any(
        not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', p) or p in ('.', '..') for p in parts
    ):
        raise ValueError('Unsafe member path: ' + name)
    path = root
    for part in parts:
        path = path / part
        if path.is_symlink():
            raise ValueError('Symlink member: ' + name)
    if not path.is_file():
        raise ValueError('Missing regular file: ' + name)
    with path.open('rb') as stream:
        data = stream.read(MAX_FILE + 1)
    if len(data) > MAX_FILE:
        raise ValueError('Member exceeds 2 MiB: ' + name)
    return data


def unique_object(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise ValueError('Duplicate JSON key: ' + key)
        obj[key] = value
    return obj


def validate(directory):
    report = {'schema': 'astrafex-user-tool-review/1', 'status': 'needs_revision',
              'accepted_for_execution': False, 'errors': [], 'imports': [],
              'review_flags': [], 'files_checked': [],
              'scope': 'Declared files only; no code execution, dependency resolution or safety certification'}
    errors = report['errors']
    try:
        root = Path(directory).resolve(strict=True)
        raw = read_member(root, 'submission.json')
        report['submission_sha256'] = hashlib.sha256(raw).hexdigest()
        m = json.loads(raw, object_pairs_hook=unique_object)
        if not isinstance(m, dict):
            raise ValueError('Manifest must be an object')
        if m.get('schema') != SCHEMA:
            errors.append('Unsupported schema')
        if not re.fullmatch(r'user_[a-z0-9_]+', str(m.get('id', ''))):
            errors.append('id must start with user_ and contain lowercase letters/digits/underscores')
        for key in ('version', 'purpose', 'provenance'):
            if not isinstance(m.get(key), str) or not m[key].strip():
                errors.append('Missing nonempty text: ' + key)
        if m.get('target') not in ('browser', 'standalone', 'both', 'server'):
            errors.append('Invalid target')
        report['requested_target'] = m.get('target')
        if m.get('target') == 'server':
            report['review_flags'].append('User tools are not installed on the Web server; revise target to browser or standalone')
        if m.get('entrypoint') != 'adapter.py:run':
            errors.append('v1 entrypoint must be adapter.py:run')
        for key in ('inputs', 'parameters', 'outputs'):
            if not isinstance(m.get(key), dict):
                errors.append('Expected object: ' + key)
        permissions = m.get('permissions', {})
        if not isinstance(permissions, dict):
            raise ValueError('permissions must be an object')
        for key in ('network', 'subprocess', 'read_outside_inputs', 'write_outside_output'):
            if type(permissions.get(key)) is not bool:
                errors.append('Explicit boolean permission required: ' + key)
            elif permissions[key]:
                report['review_flags'].append('Requested permission: ' + key)
        limits = m.get('limits', {})
        if not isinstance(limits, dict):
            raise ValueError('limits must be an object')
        for key in ('seconds', 'memory_mib', 'output_mib'):
            if type(limits.get(key)) is not int or limits[key] <= 0:
                errors.append('Positive integer limit required: ' + key)
        dependencies = m.get('dependencies')
        if not isinstance(dependencies, list):
            raise ValueError('dependencies must be a list')
        declared = set()
        for dep in dependencies:
            if not isinstance(dep, dict):
                raise ValueError('Dependency must be an object')
            for key in ('name', 'version', 'source', 'license', 'license_source', 'sha256', 'kind'):
                if not isinstance(dep.get(key), str) or not dep[key].strip():
                    errors.append('Missing dependency field: ' + key)
            version = str(dep.get('version', ''))
            if not re.fullmatch(r'[0-9][A-Za-z0-9.!+_-]*', version):
                errors.append('Dependency must have an exact version')
            if not SHA.fullmatch(str(dep.get('sha256', ''))):
                errors.append('Dependency SHA-256 required')
            if dep.get('kind') not in ('python', 'wasm', 'executable'):
                errors.append('Invalid dependency kind')
            for key in ('source', 'license_source'):
                if not str(dep.get(key, '')).startswith('https://'):
                    errors.append('HTTPS dependency source required: ' + key)
            names = dep.get('imports')
            if not isinstance(names, list) or any(not isinstance(n, str) or not n.isidentifier() for n in names):
                errors.append('Dependency imports must be a list of top-level module names')
            else:
                declared.update(names)
        members = m.get('files')
        if not isinstance(members, dict) or not 1 <= len(members) <= 100:
            raise ValueError('files must contain 1..100 members')
        if 'submission.json' in members:
            errors.append('Do not self-hash submission.json')
        if 'adapter.py' not in members or m.get('test_plan') not in members:
            errors.append('adapter.py and test_plan must be listed')
        local_modules = {n[:-3] for n in members if '/' not in n and n.endswith('.py')}
        stdlib = set(getattr(sys, 'stdlib_module_names', ()))
        report['inspector_python'] = sys.version.split()[0]
        if not stdlib:
            errors.append('Inspector requires Python 3.10+ for standard-library inventory')
        total = 0
        for name, expected in members.items():
            data = read_member(root, name)
            total += len(data)
            if total > MAX_TOTAL:
                raise ValueError('Submission exceeds 16 MiB')
            actual = hashlib.sha256(data).hexdigest()
            report['files_checked'].append({'path': name, 'sha256': actual})
            if not isinstance(expected, str) or not SHA.fullmatch(expected) or actual != expected:
                errors.append('SHA-256 mismatch: ' + name)
            if not name.endswith('.py'):
                continue
            tree = ast.parse(data, filename=name)
            if name == 'adapter.py':
                runs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'run']
                if len(runs) != 1 or [a.arg for a in runs[0].args.args] != ['inputs', 'params', 'output_dir'] or runs[0].args.posonlyargs or runs[0].args.kwonlyargs or runs[0].args.vararg or runs[0].args.kwarg:
                    errors.append('Expected def run(inputs, params, output_dir)')
            for node in ast.walk(tree):
                roots = []
                if isinstance(node, ast.Import):
                    roots = [a.name.split('.')[0] for a in node.names]
                elif isinstance(node, ast.ImportFrom):
                    if node.level:
                        report['review_flags'].append(name + ': relative import requires manual resolution')
                    else:
                        roots = [(node.module or '').split('.')[0]]
                for module in roots:
                    report['imports'].append({'file': name, 'line': node.lineno, 'module': module})
                    if module not in stdlib | declared | local_modules:
                        errors.append(name + ': undeclared import ' + module)
                    if module in {'subprocess', 'socket', 'urllib', 'http', 'ctypes', 'js', 'pyodide', 'micropip', 'importlib', 'os', 'pickle'}:
                        report['review_flags'].append(name + ': inspect privileged/dynamic module ' + module)
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {'eval', 'exec', '__import__', 'compile'}:
                    report['review_flags'].append(name + ': dynamic code call ' + node.func.id)
        if local_modules & (stdlib | declared):
            errors.append('Local module shadows standard/dependency module')
        report['review_flags'].append('Manual review required: import side effects, dynamic access, dependency closure, licences, sandbox and numerical tests')
    except (OSError, ValueError, TypeError, SyntaxError, RecursionError) as exc:
        errors.append(str(exc))
    report['status'] = 'needs_revision' if errors else 'static_pass_needs_review'
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    args = parser.parse_args()
    report = validate(args.directory)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 1 if report['errors'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
