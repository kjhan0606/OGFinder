import importlib.util
import io
import json
from pathlib import Path
import urllib.error
import pytest

SOURCE = Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug/fetch_ancillary_catalogues.py'
spec = importlib.util.spec_from_file_location('fetch_ancillary_catalogues', SOURCE)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class Fake:
    credential = {'account_name':'someone', 'password':'literal/slash$secret'}
    def __init__(self, failure=None): self.calls=[]; self.failure=failure
    def json(self, action, **fields):
        self.calls.append(action)
        if action=='submit':
            if self.failure: raise urllib.error.HTTPError('https://example',401,'Unauthorized',{},None)
            return {'id':7}
        return {'status':'done'}
    def post(self, action, **fields):
        self.calls.append(action)
        return io.BytesIO(b'# target_cluster_id,object_id,ra,dec,a_i\n1,2,3,4,.1\n')


def test_download_verify_and_resume_without_resubmit(tmp_path):
    c=Fake();s={'jobs':{}};p=tmp_path/'status.json'
    m.fetch_one(c,'extinction_000','SELECT 1',tmp_path,s,p)
    assert s['jobs']['extinction_000']['rows']==1
    m.fetch_one(c,'extinction_000','SELECT 1',tmp_path,s,p)
    assert c.calls==['submit','status','download']
    assert all(secret not in p.read_text() for secret in c.credential.values())
    (tmp_path/'extinction_000.csv').write_text('changed')
    with pytest.raises(ValueError,match='changed'):m.fetch_one(c,'extinction_000','SELECT 1',tmp_path,s,p)


def test_existing_remote_job_is_resumed(tmp_path):
    c=Fake();s={'jobs':{'extinction_000':{'job_id':7,'status':'submitted'}}}
    m.fetch_one(c,'extinction_000','SELECT 1',tmp_path,s,tmp_path/'status.json')
    assert 'submit' not in c.calls


def test_auth_failure_is_retryable_without_saving_password(tmp_path):
    c=Fake(True);s={'jobs':{}};p=tmp_path/'status.json'
    with pytest.raises(urllib.error.HTTPError):m.fetch_one(c,'extinction_000','SELECT 1',tmp_path,s,p)
    assert s['jobs']['extinction_000']['status']=='pending'
    assert c.credential['password'] not in p.read_text()


def test_ambiguous_submission_is_not_duplicated(tmp_path):
    c=Fake();s={'jobs':{'extinction_000':{'status':'submitting'}}}
    with pytest.raises(ValueError,match='unknown'):m.fetch_one(c,'extinction_000','SELECT 1',tmp_path,s,tmp_path/'status.json')
    assert not c.calls


def test_bad_csv_and_non_select_are_rejected(tmp_path):
    p=tmp_path/'bad.csv';p.write_text('<html>error</html>')
    with pytest.raises(ValueError):m.validate_csv(p,'extinction_000')
    with pytest.raises(ValueError):m.plan([{'cluster_id':'1','extinction_sql':'DELETE FROM x;'}],'SELECT 1;',10)


@pytest.mark.parametrize('prefix', ['', '# ', '\ufeff# '])
def test_actual_hsc_header_convention(tmp_path, prefix):
    p=tmp_path/'table.csv'
    p.write_text(prefix+'target_cluster_id,object_id,ra,dec,a_i\n1439,40973322283929968,132.65001852423325,-0.5704721319393459,0.04622873\n')
    assert m.validate_csv(p,'extinction_000')==1
    with p.open(encoding='utf-8-sig') as f:
        assert next(m.catalogue_reader(f))['target_cluster_id']=='1439'


def test_commented_patch_header_and_empty_extinction(tmp_path):
    p=tmp_path/'table.csv';p.write_text('# skymap_id,ra,dec\n82780004,29.13,-6.69\n')
    assert m.validate_csv(p,'patch_qa')==1
    p.write_text('# target_cluster_id,object_id,ra,dec,a_i\n')
    assert m.validate_csv(p,'extinction_000')==0


def test_header_cleanup_does_not_hide_duplicates_or_ragged_rows(tmp_path):
    p=tmp_path/'table.csv';p.write_text('# target_cluster_id,target_cluster_id\n1,2\n')
    with pytest.raises(ValueError,match='duplicate'):m.validate_csv(p,'extinction_000')
    p.write_text('# target_cluster_id,object_id,ra,dec,a_i\n1,2,3,4\n')
    with pytest.raises(ValueError,match='Malformed'):m.validate_csv(p,'extinction_000')
