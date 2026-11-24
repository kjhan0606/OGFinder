import hashlib
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from build_s16a_calibration_evidence import SELECTION_PRODUCTS, selection_evidence
from ogfkit import cliexpand

PILOT = dict(plan_sha256='p'*64, receipt_sha256='r'*64, object_qa_sha256='q'*64)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(directory, payload, artifact):
    directory.mkdir(parents=True)
    (directory/artifact).write_text(artifact+'\n')
    payload = dict(payload, science_measurement=False, calibration_applied=False,
                   output_sha256={artifact: digest(directory/artifact)})
    (directory/'summary.json').write_text(json.dumps(payload))
    return directory/'summary.json'


def selection_tree(root):
    intake = root/'intake_summary.json'
    intake.write_text('{}')
    missing = write(root/SELECTION_PRODUCTS['missingness'],
                    dict(provenance=PILOT, selected_objects=10, source_object_counts=dict(missing=dict(n=2)),
                         limitations=['a']), 'center_missingness.tsv')
    table = json.loads(missing.read_text())['output_sha256']['center_missingness.tsv']
    write(root/SELECTION_PRODUCTS['matched_missingness'],
          dict(input_summary_sha256=digest(missing), input_table_sha256=table, n_lenses=1, all_five_supported=0,
               limitations=['b']), 'matched_missingness.tsv')
    write(root/SELECTION_PRODUCTS['selection_fields'],
          dict(intake_summary_sha256=digest(intake), n_literal_names_present=0, n_required_fields=30,
               limitations=['c']), 'required_selection_fields.tsv')
    return intake


def test_selection_evidence_verifies_lineage(tmp_path):
    intake = selection_tree(tmp_path)
    result = selection_evidence(tmp_path, intake, PILOT)
    assert result['status'] == 'verified'
    assert (result['selected_objects'], result['missing_objects'], result['n_required_quality_fields']) == (10, 2, 30)
    assert result['limitations'] == ['a', 'b', 'c']


def test_changed_selection_artifact_is_refused(tmp_path):
    intake = selection_tree(tmp_path)
    (tmp_path/SELECTION_PRODUCTS['matched_missingness']/'matched_missingness.tsv').write_text('edited\n')
    with pytest.raises(ValueError, match='Changed selection artifact'):
        selection_evidence(tmp_path, intake, PILOT)


def test_other_pilot_lineage_is_refused(tmp_path):
    intake = selection_tree(tmp_path)
    with pytest.raises(ValueError, match='Selection pilot lineage mismatch'):
        selection_evidence(tmp_path, intake, dict(PILOT, receipt_sha256='x'*64))


def test_other_intake_is_refused(tmp_path):
    selection_tree(tmp_path)
    other = tmp_path/'other_intake.json'
    other.write_text('{"different": true}')
    with pytest.raises(ValueError, match='intake lineage mismatch'):
        selection_evidence(tmp_path, other, PILOT)


def test_missing_selection_product_is_refused(tmp_path):
    intake = selection_tree(tmp_path)
    (tmp_path/SELECTION_PRODUCTS['selection_fields']/'summary.json').unlink()
    with pytest.raises(FileNotFoundError):
        selection_evidence(tmp_path, intake, PILOT)


@pytest.mark.parametrize('root,expected', [('', None), ('/saved/selection', '/saved/selection')])
def test_desktop_step_passes_selection_root_only_when_set(root, expected):
    repo = Path(__file__).resolve().parents[3]
    manifests = cliexpand.load_manifests(str(repo))
    argv = cliexpand.build_argv(manifests, 'icl', 's16a-calibration',
                                dict(python='python3', work='/w', image='', catalog='', root=str(repo)),
                                {'s16a-selection-root': root})
    if expected is None:
        assert '--selection-root' not in argv
    else:
        assert argv[argv.index('--selection-root')+1] == expected
    assert argv[argv.index('--output-dir')+1] == '/w/s16a_calibration'
