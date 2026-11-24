import json
from pathlib import Path
import subprocess
import sys

import pytest

from receive_s16a_cosmos_companion import APPROVAL_HASH, DIRECTORY, credential, load_plan

HERE = Path(__file__).resolve().parent
PROBE = Path('/home/kjhan/BACKUP/Astrafek/s16a_cosmos_companion_metadata_20261124_02.json')


@pytest.mark.skipif(not (DIRECTORY/'plan.json').is_file(), reason='pinned plan is local data')
def test_pinned_plan_matches_authenticated_head_evidence():
    plan = load_plan(DIRECTORY/'plan.json')
    assert [f['filename'] for f in plan['files']] == ['target_wide_s17a_9812.fits', 'target_wide_s17a_9813.fits']
    assert plan['maximum_total_payload_bytes'] == 552268800
    if PROBE.is_file():
        probe = {f['filename']: f for f in json.loads(PROBE.read_text())['files']}
        for entry in plan['files']:
            head = probe[entry['filename']]
            assert (entry['url'], entry['content_length'], entry['etag']) == (head['url'], head['file_size_bytes'], head['etag'])


def test_edited_plan_is_refused(tmp_path):
    plan = dict(schema='x', files=[dict(filename='a.fits', content_length=3)], maximum_total_payload_bytes=3)
    plan['plan_sha256'] = APPROVAL_HASH
    (tmp_path/'plan.json').write_text(json.dumps(plan))
    with pytest.raises(ValueError, match='pinned two-file'):
        load_plan(tmp_path/'plan.json')


def test_environment_credential_requires_password(monkeypatch):
    class Args:
        use_env_auth, user = True, 'someone'
    monkeypatch.delenv('HSC_PASSWORD', raising=False)
    with pytest.raises(ValueError, match='not set'):
        credential(Args)
    monkeypatch.setenv('HSC_PASSWORD', 'fixture')
    monkeypatch.delenv('HSC_USER', raising=False)
    assert credential(Args) == 'Basic c29tZW9uZTpmaXh0dXJl'


def test_prompt_credential_refuses_without_terminal():
    code = ('import sys; sys.path.insert(0, %r)\n'
            'from receive_s16a_cosmos_companion import credential\n'
            'class A: use_env_auth=False; user="u"\n'
            'credential(A)') % str(HERE)
    result = subprocess.run([sys.executable, '-c', code], stdin=subprocess.DEVNULL, capture_output=True, text=True)
    assert result.returncode != 0 and 'interactive terminal' in result.stderr


@pytest.mark.skipif(not (DIRECTORY/'plan.json').is_file(), reason='pinned plan is local data')
def test_receive_requires_explicit_plan_hash():
    result = subprocess.run([sys.executable, str(HERE/'receive_s16a_cosmos_companion.py'), '--receive'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 2 and 'explicit user approval' in result.stderr
