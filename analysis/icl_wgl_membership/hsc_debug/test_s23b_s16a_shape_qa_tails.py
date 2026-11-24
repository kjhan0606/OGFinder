from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_s23b_s16a_shape_qa_tails import flags


def test_tail_boundaries_are_flags_and_nan_is_not_an_extreme_value():
    assert flags(-1,2,0,.4)==['m_le_minus_one','ellipticity_norm_ge_two']
    assert flags(-.99,1.99,0,.4001)==['derived_sigma_e_gt_0p4']
    assert flags(float('nan'),float('nan'),0,float('nan'))==[]
