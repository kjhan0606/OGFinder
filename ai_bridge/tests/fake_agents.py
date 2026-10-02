"""Fake agent CLIs for the unit tests: small executables (python scripts with a shebang) named like the real CLIs
that mimic the *documented / observed* invocation and output format of each one:

  codex   `codex exec ... -`           prompt on stdin, final message printed to stdout (plain text)
  claude  `claude -p --output-format json ...`  prompt on stdin, JSON envelope {"type":"result","result":"<text>",...}
  agy     `agy --print PROMPT --output-format json ...`   envelope {"conversation_id","status","response",...}
  gemini  `gemini -p TEXT --output-format json ...`  prompt on stdin, envelope {"response","stats":{"models":{..}}}
  grok    `grok --prompt-file F --output-format json ...`   envelope shape UNVERIFIED -> {"result": "<text>"} here

They are NOT the real programs: they only prove that the bridge builds the argv / prompt, parses these shapes and
enforces its rules.  `behaviour` (a dict baked into the script) selects failure modes.
"""
import json
import os
import stat
import sys

TEMPLATE = r'''#!%(python)s
import hashlib, json, os, sys, time
CFG = %(cfg)s
LOG = CFG['log']
argv = sys.argv[1:]
flavor = CFG['flavor']

def flag(name):
    return name in argv

def value(name):
    return argv[argv.index(name) + 1] if name in argv and argv.index(name) + 1 < len(argv) else None

# --- argument checks that mimic the real parsers (unknown option -> exit 1)
known = CFG['known']
for a in argv:
    if a.startswith('-') and a not in known and not a.startswith('--image='):
        sys.stderr.write("error: unknown option '%%s'\n" %% a); sys.exit(1)

n = 0
try:
    n = int(open(LOG + '.count').read())
except Exception:
    pass
open(LOG + '.count', 'w').write(str(n + 1))

stdin = '' if sys.stdin.isatty() else sys.stdin.read()
prompt = stdin
if flavor == 'agy':
    prompt = value('--print') or value('-p') or ''
elif flavor == 'gemini':
    prompt = stdin + '\n' + (value('-p') or '')
elif flavor == 'grok':
    pf = value('--prompt-file')
    prompt = open(pf).read() if pf else ''
cwd_files = sorted(os.listdir(os.getcwd()))
rec = {'argv': argv, 'cwd': os.getcwd(), 'cwd_files': cwd_files, 'env': sorted(os.environ), 'prompt': prompt,
       'stdin_len': len(stdin), 'n': n + 1,
       'secret_seen': os.environ.get('FAKE_SECRET'), 'stdin_tty': sys.stdin.isatty()}
with open(LOG, 'a') as f:
    f.write(json.dumps(rec) + '\n')

beh = CFG['behaviour']
mode = beh.get('mode', 'ok')
if beh.get('fail_first') and n < beh['fail_first']:
    mode = beh.get('first_mode', 'prose')

def objects():
    a = prompt.index('OBJECTS (JSON, begin):') + len('OBJECTS (JSON, begin):')
    b = prompt.index('OBJECTS (end)')
    return json.loads(prompt[a:b].strip())

def cols():
    # "  - NAME (type[ in [0, 1]][, optional]): doc"
    out = []
    for ln in prompt.splitlines():
        if ln.startswith('  - '):
            nm, rest = ln[4:].split(' (', 1)
            ty = rest.split(')')[0].split(' ')[0].split(',')[0]
            out.append((nm, ty, 'optional' in rest.split('):')[0]))
    return out

def answer():
    res = []
    for i, r in enumerate(objects()):
        it = {'id': r['id']}
        for nm, ty, opt in cols():
            if opt:
                continue
            if ty == 'float':
                it[nm] = round(0.1 + 0.01 * ((int(r['id']) if r['id'].isdigit() else i) %% 50), 4)
            elif ty == 'int':
                it[nm] = i
            elif ty == 'bool':
                it[nm] = True
            else:
                it[nm] = 'lbl%%d' %% (i %% 3)
        res.append(it)
    return {'results': res, 'model': beh.get('model', 'fake-model-1')}

def text_for(mode):
    a = answer()
    if mode == 'ok':
        return json.dumps(a)
    if mode == 'fence':
        return '```json\n' + json.dumps(a, indent=1) + '\n```'
    if mode == 'prose':
        return 'Sure! Here is the JSON you asked for:\n' + json.dumps(a) + '\nHope that helps.'
    if mode == 'garbage':
        return 'I cannot do that.'
    if mode == 'bad_range':
        for it in a['results']:
            for k in list(it):
                if k in ('STAR_PROB', 'REALBOGUS_SCORE', 'MORPH_CONF', 'CLASS_PROB'):
                    it[k] = 1.7
        return json.dumps(a)
    if mode == 'missing_id':
        a['results'] = a['results'][:-1]
        return json.dumps(a)
    if mode == 'wrong_type':
        for it in a['results']:
            for k in list(it):
                if isinstance(it[k], float):
                    it[k] = 'high'
        return json.dumps(a)
    if mode == 'extra_id':
        a['results'].append(dict(a['results'][0], id='999999'))
        return json.dumps(a)
    if mode == 'per_object_error':
        a['results'][0] = {'id': a['results'][0]['id'], 'error': 'cannot judge this one'}
        return json.dumps(a)
    return json.dumps(a)

if mode == 'sleep':
    time.sleep(beh.get('seconds', 30))
if mode == 'exit1':
    sys.stderr.write('internal error: boom\n'); sys.exit(1)
if mode == 'auth':
    msg = 'Not logged in - please run the login command'
    if flavor == 'claude':
        print(json.dumps({'type': 'result', 'subtype': 'success', 'is_error': True, 'result': msg}))
        sys.exit(1)
    sys.stderr.write(msg + '\n'); sys.exit(1)
if mode == 'rate':
    sys.stderr.write('429 rate limit exceeded, try again later\n'); sys.exit(1)
if mode == 'secret_echo':
    sys.stderr.write('failed with key ' + os.environ.get('FAKE_SECRET', '') + '\n'); sys.exit(1)
if mode == 'inject_probe':
    # prove the bridge never executed the catalog text: the marker file must not exist
    pass

t = text_for(mode)
if flavor == 'codex':
    sys.stderr.write('OpenAI Codex vX (fake)\n')
    sys.stdout.write(t + '\n')
elif flavor == 'claude':
    env = {'type': 'result', 'subtype': 'success', 'is_error': False, 'duration_ms': 5, 'num_turns': 1, 'result': t,
           'session_id': 'sess-' + hashlib.md5(prompt.encode()).hexdigest()[:8], 'total_cost_usd': 0.001, 'modelUsage': {beh.get('model', 'claude-fake-1'): {}}}
    if beh.get('structured'):
        env['structured_output'] = json.loads(t)
    sys.stdout.write(json.dumps(env) + '\n')
elif flavor == 'agy':
    sys.stdout.write(json.dumps({'conversation_id': 'conv-' + hashlib.md5(prompt.encode()).hexdigest()[:8], 'status': 'SUCCESS', 'response': t + '\n',
                                 'duration_seconds': 0.1, 'num_turns': 1,
                                 'usage': {'input_tokens': 1, 'output_tokens': 1}}) + '\n')
elif flavor == 'gemini':
    sys.stdout.write(json.dumps({'response': t, 'stats': {'models': {beh.get('model', 'gemini-fake-2.5'): {}}}}) + '\n')
elif flavor == 'grok':
    sys.stdout.write(json.dumps({'result': t}) + '\n')
'''

KNOWN = {
    'codex': ['exec', '--skip-git-repo-check', '--ephemeral', '--ignore-user-config', '--ignore-rules', '--sandbox',
              '--color', '-C', '-', '--image', '-m', '--model', '--output-schema', '-o'],
    'claude': ['-p', '--print', '--output-format', '--no-session-persistence', '--disable-slash-commands',
               '--strict-mcp-config', '--permission-mode', '--tools', '--model', '--max-turns', '--bare'],
    'agy': ['--print', '-p', '--output-format', '--print-timeout', '--disable-slash-commands', '--model', '--effort',
            '--sandbox'],
    'gemini': ['-p', '--prompt', '--output-format', '-o', '--approval-mode', '--skip-trust', '-m', '--model'],
    'grok': ['--prompt-file', '--output-format', '--permission-mode', '--tools', '--no-auto-update', '--cwd', '-m',
             '--model', '--max-turns'],
}


def make_fake(directory, name, flavor=None, behaviour=None, log=None, known_extra=()):
    """Write executable `directory/name`; returns (path, log_path).  `log` receives one JSON line per call."""
    flavor = flavor or name
    log = log or os.path.join(directory, name + '.calls.jsonl')
    cfg = {'flavor': flavor, 'behaviour': behaviour or {}, 'log': log, 'known': KNOWN[flavor] + list(known_extra)}
    path = os.path.join(directory, name)
    with open(path, 'w') as f:
        f.write(TEMPLATE % {'python': sys.executable, 'cfg': repr(cfg)})
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path, log


def calls(log):
    try:
        return [json.loads(l) for l in open(log) if l.strip()]
    except OSError:
        return []
