#!/usr/bin/env python3
"""Agent-CLI backends: GUI test under Xvfb + exported-script replay + pipeline mode, all with FAKE agent CLIs on PATH.

  verify_agent_gui.py [--python PY] [--display :77] [--workdir DIR]

NO real Codex / Claude Code / agy / Gemini / Grok binary is run and no model is contacted: the executables on PATH
are small scripts (ai_bridge/tests/fake_agents.py) that mimic each CLI's documented output format.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from ai_bridge.tests import fake_agents as fa  # noqa: E402

results = []


def rec(check, ok, detail=''):
    results.append((check, bool(ok), detail))
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', check, ('  -- ' + detail) if detail else ''), flush=True)
    return ok


def ncalls(log):
    return len(fa.calls(log))


def col(path, name):
    L = open(path).read().split('\n')
    h = L[0].split('\t')
    i = h.index(name)
    return [l.split('\t')[i] for l in L[1:] if l]


def drop_time(path):
    """catalogue text without the *_TIME provenance columns (a fresh request has a new wall-clock time)"""
    L = open(path).read().rstrip('\n').split('\n')
    h = L[0].split('\t')
    keep = [i for i, c in enumerate(h) if not c.endswith('_TIME')]
    return '\n'.join('\t'.join(r.split('\t')[i] for i in keep) for r in L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--python', default=os.environ.get('OGFINDER_PYTHON', '/workspace/ogf_venv/bin/python'))
    ap.add_argument('--display', default=':77')
    ap.add_argument('--workdir', default=None)
    a = ap.parse_args()
    base = a.workdir or tempfile.mkdtemp(prefix='ogf_agent_verify_')
    shutil.rmtree(base, ignore_errors=True)
    os.makedirs(base)
    print('workdir', base)
    binp = os.path.join(base, 'bin')
    os.makedirs(binp)
    logs = {}
    for name, flavor in (('claude', 'claude'), ('codex', 'codex'), ('gemini', 'gemini')):      # agy and grok: NOT installed
        _p, logs[name] = fa.make_fake(binp, name, flavor)
    sf = os.path.join(base, 'ai_services.json')
    json.dump({'services': []}, open(sf, 'w'))
    if subprocess.run(['pgrep', '-f', 'Xvfb %s' % a.display], stdout=subprocess.DEVNULL).returncode != 0:
        subprocess.Popen(['Xvfb', a.display, '-screen', '0', '1400x1000x24'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        time.sleep(1.5)
    out = os.path.join(base, 'gui')
    os.makedirs(out)
    for f in ('ds9.auto', 'ds9.auto.dir'):
        shutil.rmtree(os.path.expanduser('~/' + f), ignore_errors=True)
    home = os.path.join(base, 'home')
    os.makedirs(home)
    path = binp + os.pathsep + '/usr/bin:/bin'
    env = dict(os.environ, DISPLAY=a.display, OGF_AGENT_OUT=out, OGF_AGENT_BIN=binp, OGFINDER_PYTHON=a.python,
               OGF_AI_SERVICES_FILE=sf, OGF_AI_CACHE_DIR=os.path.join(base, 'cache'), PATH=path)
    print('\n=== GUI session (Xvfb %s, fake agent CLIs: claude codex gemini on PATH; agy, grok missing) ===' % a.display)
    t0 = time.time()
    p = subprocess.run([os.path.join(ROOT, 'bin', 'ds9'), '/workspace/fits/m51.fits', '-geometry', '1300x950', '-source',
                        os.path.join(HERE, 'verify_agent_gui.tcl')], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)
    dt = time.time() - t0
    se = p.stderr.decode()
    open(os.path.join(out, 'gui_stderr.txt'), 'w').write(se)
    open(os.path.join(out, 'gui_stdout.txt'), 'w').write(p.stdout.decode())
    log = open(os.path.join(out, 'agent_gui_steps.log')).read()
    rec('ds9 session exit 0 and DONE logged (%.1f s)' % dt, p.returncode == 0 and 'DONE' in log)
    bad_err = [l for l in se.splitlines() if l.strip() and not l.startswith('Filtered:')]
    rec('no Tcl errors on stderr', not bad_err, '; '.join(bad_err[:3]))

    def has(pat, flags=0):
        return re.search(pat, log, flags) is not None

    rec('detection: codex / claude / agy (via gemini) installed, grok not found',
        has(r'detect agent_codex: installed=1') and has(r'detect agent_claude: installed=1') and has(r'detect agent_agy: installed=1 path=\S*/gemini')
        and has(r'detect agent_grok: installed=0'))
    rec('registry lists the four ready-made agent profiles + mock, with "AGENT CLI installed / not found"',
        has(r'registry rows: .*agent_codex') and has(r'registry note agent_claude: AGENT CLI installed / agent_grok: AGENT CLI not found'))
    rec('"Detect Agent CLIs" prints installed / NOT FOUND per CLI', has(r'registry detect text: .*Codex CLI +installed') and has(r'Grok \(xAI CLI\) +NOT FOUND'))
    rec('registry "Test Connection" on an agent runs no model and opens no dialog',
        has(r'registry test claude: OK: .*no model call made') and has(r'confirmation dialogs after detect/test: 0'))
    rec('Run Task dialog Backend dropdown lists Codex CLI, Claude Code, agy / Gemini CLI, Grok, Generic services with installed / NOT FOUND',
        has(r'dialog backends: \{Codex CLI \(codex exec\)  -  installed\} \{Claude Code \(claude -p\)  -  installed\} '
            r'\{agy / Gemini CLI \(agy, gemini -p\)  -  installed\} \{Grok \(xAI CLI\)  -  NOT FOUND on PATH\} \{Generic services'))
    rec('choosing Claude Code shows the agent frame with the found path', has(r'dialog claude: services=agent_claude agframe=1 st=Found: \S*/claude'))
    rec('choosing Generic hides the agent frame and lists only non-agent services', has(r'dialog generic: services=mock agframe=0'))
    rec('choosing Grok shows NOT FOUND with the install hint', has(r'dialog grok: services=agent_grok st=NOT FOUND on PATH'))
    rec('Grok run (not installed): clear status, no confirmation, no columns', has(r'grok run \(not installed\): status=AI Services: the agent CLI of .agent_grok. was not found on PATH.* confirmations=0'))
    rec('dry run: prompt + argv in the log, no confirmation, no columns, CLI not executed',
        has(r'dry run: status=AI Services: dry run done .*cols=\d+ confirmations=0') and has(r'dry log has prompt: 1 argv: 1'))
    rec('declined confirmation: nothing runs, no columns', has(r'declined: status=AI Services: cancelled confirmations=1'))
    rec('dry run / declined / not-installed / registry test: no fake CLI received a prompt',
        has(r'prompt-carrying calls of the fakes before any confirmed run: claude=0 codex=0 gemini=0'))
    m = re.search(r'claude run: status=AI Services: (\d+) column\(s\) added \(service agent_claude\).* new=(.*)', log)
    rec('Claude run added STAR_PROB, CLASS_LABEL and AI_SG_* provenance', bool(m) and 'STAR_PROB CLASS_LABEL AI_SG_SERVICE AI_SG_MODEL' in m.group(2), m.group(0)[:160] if m else '')
    conf = re.search(r'CONFIRMTEXT: (.*)', log).group(1)
    rec('confirmation names the CLI path, the provider, columns sent, "NONE" for images and the login note',
        'agent CLI claude' in conf and "Anthropic's cloud model" in conf and 'LEAVES YOUR COMPUTER' in conf and 'image data: NONE' in conf
        and 'mags' in conf and 'No credentials are stored' in conf and 'catalog row text: no' in conf, conf[:200])
    rec('second Claude run: no second confirmation, columns updated in place', has(r'claude 2nd run: confirmations=0 cols=(\d+) \(was \1\)'))
    cl = [x for x in fa.calls(logs['claude']) if '-p' in x['argv']]
    rec('Claude fake got: prompt on stdin, tools disabled, JSON output, empty temp dir, no env secrets', len(cl) >= 2 and cl[0]['stdin_len'] > 500 and
        cl[0]['argv'][cl[0]['argv'].index('--tools') + 1] == '' and cl[0]['cwd_files'] == [] and 'ANTHROPIC_API_KEY' not in cl[0]['env'])
    conf2 = re.search(r'CONFIRMTEXT2: (.*)', log).group(1)
    rec('codex morphology without images: confirmation says image data NONE; the CLI got no image file',
        'image data: NONE' in conf2 and has(r'codex morphology no images: status=AI Services: \d+ column') and
        (lambda cc: len(cc) == 12 and all(x['cwd_files'] == [] and not any(a.startswith('--image') for a in x['argv']) for x in cc))(
            [x for x in fa.calls(logs['codex']) if 'exec' in x['argv']][:12]))
    conf3 = re.search(r'CONFIRMTEXT3: (.*)', log).group(1)
    rec('ticking "Send image cutouts" asks again and names the image files; the CLI then gets them',
        has(r'codex morphology WITH images: confirmations=1') and 'IMAGE DATA: ' in conf3 and 'cutout file(s)' in conf3 and
        any(any(f.startswith('img_') for f in x['cwd_files']) and any(a.startswith('--image=') for a in x['argv'])
            for x in fa.calls(logs['codex']) if 'exec' in x['argv']))
    rec('agy service ran through the gemini executable (agy missing)', has(r'agy\(gemini\) real_bogus: status=AI Services: \d+ column') and
        any(x['argv'][:1] == ['-p'] for x in fa.calls(logs['gemini'])))
    rec('"Save to profile" wrote executable + extra_args only (no secrets); the next run used them',
        has(r'save to profile: Saved: agent_claude: executable=\S+ extra_args=\["--model", "haiku"\]') and
        any(x['argv'][-2:] == ['--model', 'haiku'] for x in fa.calls(logs['claude'])))
    prof = open(sf).read()
    rec('profile file holds executable/extra_args and no key / token text', '"extra_args"' in prof and not re.search(r'sk-|Bearer|api_key|token', prof, re.I))
    # ---- recorder + export
    script = os.path.join(out, 'ogfinder_session.py')
    rec('session script exported', os.path.exists(script))
    txt = open(script).read() if os.path.exists(script) else ''
    steps = re.findall(r'^#\s+\d+\s+(\S+)\s+(\S+)\s+(.*)$', txt, re.M)
    ai = [s for s in steps if s[1] == 'ai.run']
    print('  recorded steps:')
    for s in steps:
        print('     %s %s %s' % s)
    rec('ai.run steps recorded for every confirmed agent run (claude x2 + codex x2 + agy + claude photoz = 6), all AUTO', len(ai) == 6 and all(s[0] == 'AUTO' for s in ai), str(len(ai)))
    rec('recorded argv keeps --allow-agent-cli and (only for the ticked run) --send-images',
        len([l for l in txt.splitlines() if l.strip().startswith('"argv_t"') and '"--allow-agent-cli"' in l]) == 6 and
        len([l for l in txt.splitlines() if l.strip().startswith('"argv_t"') and '"--send-images"' in l]) == 1)
    rec('script holds no key / token text', not re.search(r'sk-[A-Za-z0-9]{10}|Bearer ', txt))
    # ---- replay in a plain shell
    print('\n=== replay (plain shell, empty HOME, fake CLIs on PATH) ===')
    for l in logs.values():
        for f in (l, l + '.count'):
            if os.path.exists(f):
                os.remove(f)
    renv = {'PATH': path, 'HOME': home, 'LANG': 'C.UTF-8', 'OGF_AI_SERVICES_FILE': sf, 'OGF_AI_CACHE_DIR': os.path.join(base, 'cache_replay')}
    renv0 = dict(renv)
    rep = os.path.join(base, 'replay')
    r = subprocess.run([a.python, script, '--python', a.python, '--mode', 'replay', '--outdir', rep, '/workspace/fits/m51.fits'],
                       env=renv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=base)
    t = r.stdout.decode()
    open(rep + '.log', 'w').write(t)
    rec('replay exit code 0', r.returncode == 0, t.strip().splitlines()[-1][:120] if t.strip() else '')
    gui_cat = os.path.join(out, 'catalog_gui.tsv')
    rc = os.path.join(rep, 'm51', 'catalog_final.tsv')
    rec('replay catalogue identical to the GUI catalogue (all columns, request time column excluded)',
        os.path.exists(rc) and drop_time(rc) == drop_time(gui_cat))
    # 111 objects / batch 10 = 12 requests per step; the 2nd identical claude star_galaxy step is a cache hit
    nc = sum(1 for x in fa.calls(logs['claude']) if '-p' in x['argv'])
    nx = sum(1 for x in fa.calls(logs['codex']) if 'exec' in x['argv'])
    ng = sum(1 for x in fa.calls(logs['gemini']) if '-p' in x['argv'])
    rec('replay called the fake CLIs again: claude 2 distinct steps x 12 requests, codex 2 x 12 (with / without images), gemini 1 x 12',
        (nc, nx, ng) == (24, 24, 12), 'claude %d codex %d gemini %d' % (nc, nx, ng))
    # ---- pipeline mode
    print('\n=== pipeline mode (new data = m51.fits again) ===')
    outs = {}
    for flag, label in (([], 'plain'), (['--allow-network'], 'allow-network only'), (['--allow-agent-cli'], 'allow-agent-cli'),
                        (['--allow-agent-cli', '--allow-agent-images'], 'allow-agent-cli + images')):
        for l in logs.values():
            for f in (l, l + '.count'):
                if os.path.exists(f):
                    os.remove(f)
        pdir = os.path.join(base, 'pipe_' + label.replace(' ', '_').replace('+', 'p'))
        renv = dict(renv0, OGF_AI_CACHE_DIR=os.path.join(base, 'cache_' + label.replace(' ', '_').replace('+', 'p')))   # fresh cache: count real calls
        r = subprocess.run([a.python, script, '--python', a.python, '--outdir', pdir, '/workspace/fits/m51.fits'] + flag,
                           env=renv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=base)
        t = r.stdout.decode()
        open(pdir + '.log', 'w').write(t)
        cat = os.path.join(pdir, 'm51', 'catalog_final.tsv')
        hdr = open(cat).readline().rstrip('\n').split('\t') if os.path.exists(cat) else []
        outs[label] = (r.returncode, t, hdr, sum(1 for x in fa.calls(logs['claude']) + fa.calls(logs['codex']) + fa.calls(logs['gemini']) if x['argv'] != ['--version']) // 12)
    rc0, t0_, h0, n0 = outs['plain']
    rec('pipeline without flags: every agent step skipped with the --allow-agent-cli message, no CLI executed, no AI columns',
        rc0 == 0 and t0_.count('use --allow-agent-cli') == 6 and n0 == 0 and 'STAR_PROB' not in h0, 'skipped=%d CLI steps run=%d' % (t0_.count('use --allow-agent-cli'), n0))
    rc1, t1, h1, n1 = outs['allow-network only']
    rec('--allow-network alone does not enable agent CLIs', rc1 == 0 and n1 == 0 and 'STAR_PROB' not in h1)
    rc2, t2, h2, n2 = outs['allow-agent-cli']
    rec('--allow-agent-cli: the 5 steps without images ran (4 CLI call batches + 1 cache hit) (STAR_PROB, PHOTOZ, REALBOGUS_SCORE present); the image step is skipped with its own message',
        rc2 == 0 and 'STAR_PROB' in h2 and 'PHOTOZ' in h2 and 'REALBOGUS_SCORE' in h2 and t2.count('use --allow-agent-images as well') == 1 and n2 == 4, 'steps that really called a CLI=%d (the 2nd identical claude step is a cache hit)' % n2)
    rc3, t3, h3, n3 = outs['allow-agent-cli + images']
    rec('--allow-agent-cli --allow-agent-images: all 6 steps ran (5 CLI call batches + 1 cache hit)', rc3 == 0 and n3 == 5 and 'MORPH_TYPE' in h3, 'CLI steps run=%d' % n3)
    npass = sum(1 for _, ok, _ in results if ok)
    print('\nSUMMARY: %d/%d checks passed' % (npass, len(results)))
    json.dump([{'check': c, 'ok': o, 'detail': d} for c, o, d in results], open(os.path.join(base, 'verify_agent_gui_results.json'), 'w'), indent=1)
    return 0 if npass == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())
