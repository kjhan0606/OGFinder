import csv
import json

from astropy.table import MaskedColumn, Table
import pytest

import intake_s16a_ephor_ab_cosmos as intake


def fixture_tables():
    base = 2**53
    points = [Table(dict(ID=[base, base+1, base+3], PHOTOZ_BEST=[.8, 1., 1.2])),
              Table(dict(ID=[base+2], PHOTOZ_BEST=[1.2], PHOTOZ_RISK_BEST=[.2]))]
    reference = Table(dict(S17a_objid=[base+1, base+2, base+4], COSMOS_photoz=[.9, 1.2, 1.4],
                           weight_source=[2., 1., 1.], SOM_weight=[3., 1., 1.]))
    reference['COSMOS_photoz'] = MaskedColumn(reference['COSMOS_photoz'], mask=[False, True, False])
    return points, reference


def test_fits_intake_preserves_exact_ids_and_records_unmatched_and_masked_rows(tmp_path, monkeypatch):
    points, reference = fixture_tables()
    folder = tmp_path/'received'
    folder.mkdir()
    entries, files = [], {}
    for i, table in enumerate([*points, reference]):
        name = f'table{i}.fits'
        file = folder/name
        table.write(file)
        size = file.stat().st_size
        entries.append(dict(filename=name, content_length=size))
        files[name] = dict(status='complete', sha256=intake.sha(file))
    total = sum(e['content_length'] for e in entries)
    receipt = dict(plan_sha256='test', reserved_bytes=total, files=files)
    (folder/intake.CHECKPOINT).write_text(json.dumps(receipt))
    monkeypatch.setattr(intake, 'load_plan', lambda _: dict(plan_sha256='test', destination=str(folder),
                                                          files=entries, maximum_total_payload_bytes=total))
    output = tmp_path/'output'
    intake.run(output)
    with (output/'cosmos_join.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert [row['s17a_objid'] for row in rows] == [str(2**53+1), str(2**53+2)]
    assert [row['eligible_for_selection_review'] for row in rows] == ['True', 'False']
    assert (output/'reference_only_ids.tsv').read_text().splitlines()[-1] == str(2**53+4)
    result = json.loads((output/'summary.json').read_text())
    assert result['point_only_rows'] == 2 and result['eligible_rows'] == 1
    assert result['invalid_by_column']['reference_z'] == 1
    assert result['effective_count_shape_som'] == 1
    assert result['joined_reference_weight_fraction'] == pytest.approx(7/8)
    assert result['eligible_reference_weight_fraction'] == pytest.approx(6/8)


def test_join_rejects_duplicate_int_ids_and_floating_id_columns():
    points, reference = fixture_tables()
    points[1]['ID'][0] = points[0]['ID'][0]
    with pytest.raises(ValueError, match='duplicate IDs'):
        intake.join_tables(points, reference)
    points, reference = fixture_tables()
    reference['S17a_objid'] = reference['S17a_objid'].astype(float)
    with pytest.raises(ValueError, match='exact unmasked integer'):
        intake.join_tables(points, reference)
