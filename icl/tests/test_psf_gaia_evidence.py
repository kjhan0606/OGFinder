from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from summarize_psf_gaia import classify


def source():
    return dict(separation_arcsec='.2',ruwe='1.1',classprob_dsc_combmod_star='.999',parallax_over_error='6',
        pmra='1',pmra_error='1',pmdec='1',pmdec_error='1',duplicated_source='false',in_qso_candidates='false',in_galaxy_candidates='false')


def test_match_evidence_does_not_certify_psf():
    row=source();label,evidence=classify([row])
    assert label=='stellar_evidence_requires_epoch_and_psf_review' and evidence==['positive_parallax_snr_ge5']
    row['parallax_over_error']='-9';row['pmdec']='-7'
    assert classify([row])[1]==['pmdec_component_snr_ge5']


def test_missing_ambiguous_and_large_offset_are_not_stars():
    row=source();assert classify([])[0]=='no_gaia_match'
    assert classify([row,row])[0]=='ambiguous_position'
    row['separation_arcsec']='1.5';assert classify([row])[0]=='offset_requires_epoch_review'


def test_bad_quality_or_missing_data_blocks_support():
    for key,val in [('ruwe','nan'),('ruwe','2'),('classprob_dsc_combmod_star',''),('duplicated_source','true'),('in_qso_candidates','true'),('in_galaxy_candidates','')]:
        row=source();row[key]=val
        assert classify([row])[0]=='matched_but_stellar_evidence_insufficient'
    row=source();row['parallax_over_error']='-10';row['pmra_error']='0'
    assert classify([row])[0]=='matched_but_stellar_evidence_insufficient'
