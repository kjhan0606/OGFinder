import numpy as np
import pytest

from audit_s16a_cosmos_afterburner_photoz import unavailable
from audit_s16a_cosmos_companion_schema import SUPPLEMENT_SQL, match_fields, source_columns
from s16a_y1_source_contract import REQUIRED


def table(*names):
    return {'t.fits': dict(rows=1, columns=[dict(name=n, format='D', unit=None) for n in names])}


@pytest.mark.skipif(not SUPPLEMENT_SQL.is_file(), reason='executed S16A supplement SQL is local data')
def test_every_required_alias_has_an_s16a_source_column():
    mapping = source_columns()
    assert set(REQUIRED) <= set(mapping)
    assert mapping['pixel_edge'] == dict(table='s16a_wide.meas', column='iflags_pixel_edge')
    assert mapping['e1']['column'] == 'ishape_hsm_regauss_e1'


def test_match_prefers_source_name_and_reports_absence():
    mapping = {k: dict(table='x', column='src_'+k) for k in REQUIRED}
    rows = {r['required_field']: r for r in match_fields(mapping, table('SRC_I_FLUX', 'pixel_edge', 'src_e1_extra'))}
    assert rows['i_flux']['status'] == 's16a_column_name_present' and rows['i_flux']['matched_column'] == 'SRC_I_FLUX'
    assert rows['pixel_edge']['status'] == 'alias_name_present'
    assert rows['e1']['status'] == 'substring_candidates_only' and rows['e1']['matched_column'] == ''
    assert rows['blendedness']['status'] == 'absent'


def test_afterburner_unavailable_uses_null_flag_and_nonfinite_values():
    data = {'f': np.array([1.0, np.nan, 2.0, np.inf]), 'f_isnull': np.array([False, False, True, False])}
    assert unavailable(data, 'f').tolist() == [False, True, True, True]
