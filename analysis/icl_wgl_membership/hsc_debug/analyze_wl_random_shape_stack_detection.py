#!/usr/bin/env python3
"""Global tangential-signal and cross-shear null checks for saved WL stacks."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import chi2

DEFAULT_ROOT=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
DEFAULT_RUNS={
 '5deg':DEFAULT_ROOT/'wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json',
 '2p5deg':DEFAULT_ROOT/'wl_random_shape_stack_lens_reuse_2026-11-04_02p5/summary.json',
}


def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()


def compute(path):
 s=json.loads(path.read_text()); j=s['paired_spatial_jackknife']
 if s.get('science_measurement') is not False or j.get('science_measurement') is not False:
  raise ValueError('science gate must remain false')
 if s.get('random_lens_subtraction_applied') is not True or j.get('n_occupied_blocks',0)<3:
  raise ValueError('expected random-subtracted stack with paired spatial covariance')
 profile=j['profile']; C=np.asarray(j['covariance'],float)
 if C.shape!=(2*len(profile),2*len(profile)) or not np.all(np.isfinite(C)):
  raise ValueError('unexpected or non-finite covariance matrix')
 if not np.allclose(C,C.T,rtol=1e-12,atol=0): raise ValueError('covariance is not symmetric')
 eig=np.linalg.eigvalsh(C)
 if eig[0]<=0: raise ValueError('covariance is not positive definite')
 ds=np.asarray([r['delta_sigma_random_subtracted_msun_mpc2'] for r in profile],float)
 cross=np.asarray([r['cross_delta_sigma_random_subtracted_msun_mpc2'] for r in profile],float)
 ds_idx=np.arange(0,2*len(profile),2); x_idx=ds_idx+1
 band_results={}
 for key,select in [('0.1-3.0_mpc',np.ones(len(profile),bool)),
                    ('0.2-3.0_mpc',np.asarray([r['r_min_mpc']>=.2 for r in profile])),
                    ('0.4-3.0_mpc',np.asarray([r['r_min_mpc']>=.4 for r in profile])),
                    ('0.7-3.0_mpc',np.asarray([r['r_min_mpc']>=.7 for r in profile]))]:
  ids=ds_idx[select]; sub=C[np.ix_(ids,ids)]; values=ds[select]
  stat=float(values@np.linalg.solve(sub,values)); dof=len(values)
  band_results[key]={'n_bins':dof,'chi_square_vs_zero':stat,'p_value_chi_square_approx':float(chi2.sf(stat,dof)),
                     'diagonal_signal_to_error':[float(a/b) for a,b in zip(values,np.sqrt(np.diag(sub)))],
                     'minimum_covariance_eigenvalue':float(np.linalg.eigvalsh(sub)[0]),
                     'condition_number':float(np.linalg.cond(sub))}
 xsub=C[np.ix_(x_idx,x_idx)]; xstat=float(cross@np.linalg.solve(xsub,cross))
 return {'summary_path':str(path),'summary_sha256':sha(path),'block_method':j['block_method'],
         'n_occupied_blocks':j['n_occupied_blocks'],'n_common_lenses':s['input']['n_matched_lenses_in_cluster_and_random_stacks'],
         'n_random_centers':s['input']['n_random_centers'],'n_random_source_pairs':s['input']['matched_random_source_pairs'],
         'science_measurement':False,'covariance_shape':list(C.shape),'covariance_rank':int(np.linalg.matrix_rank(C)),
         'minimum_full_covariance_eigenvalue':float(eig[0]),'maximum_covariance_condition_number':float(np.linalg.cond(C)),
         'tangential_zero_signal_tests':band_results,
         'cross_component_zero_tests':{'n_bins':len(cross),'chi_square':xstat,'p_value_chi_square_approx':float(chi2.sf(xstat,len(cross)))},
         'point_profile':profile}


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--run-5deg',type=Path,default=DEFAULT_RUNS['5deg'])
 p.add_argument('--run-2p5deg',type=Path,default=DEFAULT_RUNS['2p5deg'])
 p.add_argument('--output-dir',type=Path,required=True)
 a=p.parse_args()
 if a.output_dir.exists():raise FileExistsError(a.output_dir)
 data={'schema':'astrafek-wl-random-subtracted-detection-audit-v1','status':'complete','science_measurement':False,
       'interpretation':'conditional global profile-null diagnostics using stored paired spatial jackknife covariance',
       'runs':{'5deg':compute(a.run_5deg),'2p5deg':compute(a.run_2p5deg)}}
 p5=data['runs']['5deg']['point_profile'];p2=data['runs']['2p5deg']['point_profile']
 if len(p5)!=len(p2) or any(a['delta_sigma_random_subtracted_msun_mpc2']!=b['delta_sigma_random_subtracted_msun_mpc2'] for a,b in zip(p5,p2)):
  raise ValueError('block-width runs do not share identical point estimates')
 data['point_estimates_identical_across_partitions']=True
 a.output_dir.mkdir(parents=True)
 out=a.output_dir/'summary.json';out.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
 lines=['# Random-subtracted weak-lensing profile null audit','',
        'Uses the saved 5° and 2.5° paired endpoint-jackknife covariance matrices. It evaluates the full tangential ΔΣ vector against zero with the covariance submatrix, and repeats after dropping inner radial bins. The cross component is tested as a null. These are covariance-conditional chi-square diagnostics; the jackknife partitions are not yet validated contiguous HSC survey regions and calibration/systematic uncertainties are incomplete.', '',
        '| Block width | Bins included | χ² vs zero | dof | nominal p | diagonal S/N range |', '|---|---|---:|---:|---:|---:|']
 for run in ('5deg','2p5deg'):
  for bounds,item in data['runs'][run]['tangential_zero_signal_tests'].items():
   sn=item['diagonal_signal_to_error']
   lines.append(f"| {run} | {bounds} | {item['chi_square_vs_zero']:.2f} | {item['n_bins']} | {item['p_value_chi_square_approx']:.3e} | {min(sn):.2f}–{max(sn):.2f} |")
 lines += ['', 'The 2.5° point estimates exactly match the 5° values, as expected; covariance changes with the block partition. The nominal p-values are not a science detection claim: block geometry, correlated coadd noise, n(z) calibration, reduced shear, B-mode polarity, and random-center matching remain systematic gates. The stack and this audit retain `science_measurement=false`.', '']
 (a.output_dir/'assessment.md').write_text('\n'.join(lines))
 print(json.dumps({'output_dir':str(a.output_dir),'summary_sha256':sha(out),'assessment_sha256':sha(a.output_dir/'assessment.md'),
                   'tangential_chi2_2p5_all':data['runs']['2p5deg']['tangential_zero_signal_tests']['0.1-3.0_mpc']['chi_square_vs_zero'],
                   'p_2p5_all':data['runs']['2p5deg']['tangential_zero_signal_tests']['0.1-3.0_mpc']['p_value_chi_square_approx']},indent=2))

if __name__=='__main__':main()
