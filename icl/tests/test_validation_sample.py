import importlib.util
from pathlib import Path
import pytest

source=Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug/prepare_validation_sample.py'
spec=importlib.util.spec_from_file_location('prepare_validation_sample',source)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_stratification_is_order_independent_and_preserves_small_strata():
    rows=[dict(id=str(i),z_cl=.31 if i<19 else .51,N_mem=18) for i in range(20)]
    a,table=m.select(rows,6,123)
    b,_=m.select(rows[::-1],6,123)
    assert a==b and len({r['id'] for r in a})==6
    assert any(r['id']=='19' for r in a)
    assert sum(r['selected'] for r in table)==6
    assert all(0<r['conditional_inclusion_probability']<=1 for r in table)


def test_exact_bin_boundaries_and_insufficient_pool():
    assert m.strata(dict(z_cl=.35,N_mem=20))==(1,1)
    assert m.strata(dict(z_cl=.50,N_mem=50))==(4,3)
    with pytest.raises(ValueError):m.select([],1,1)
