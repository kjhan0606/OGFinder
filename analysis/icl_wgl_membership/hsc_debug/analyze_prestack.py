#!/usr/bin/env python3
"""Summarize re-detected injections and complete-case empirical stack covariance."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import platform

import numpy as np
import prepare_detection_psf as prep
from download_pilot import file_sha256 as sha
from prestack_core import EDGES, FACTORS, covariance_diagnostics, spatial_labels


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extraction', type=Path, required=True); p.add_argument('--work', type=Path, required=True)
    p.add_argument('--clusters', type=Path, help='existing 50-cluster profile output for paired covariance')
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False)
    src = a.extraction/'extraction_summary.json'; doc = json.loads(src.read_text())
    if doc['status'] != 'complete':
        raise ValueError('Extraction did not complete')
    protocol = a.extraction/'protocol.json'
    if sha(protocol) != doc['protocol_sha256']:
        raise ValueError('Protocol changed')
    for path, digest in doc['script_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Extraction code changed: '+path)
    records = doc['records']; hashes = {str(src):sha(src), str(protocol):sha(protocol)}
    fields = []; recovery = []
    for record in records:
        path = a.extraction/record['control_id']/'baseline.npz'
        if sha(path) != record['baseline_sha256']:
            raise ValueError('Profile file changed')
        hashes[str(path)] = record['baseline_sha256']
        with np.load(path) as z:
            fields.append({key:z[key].copy() for key in ('raw','corrected','coverage','nominal_raw','nominal_corrected')})
        if record['n_injection_models']:
            path = a.extraction/record['control_id']/'recovery.tsv'
            hashes[str(path)] = sha(path)
            with path.open() as f:
                local=list(csv.DictReader(f, delimiter='\t'))
            open_truth={}
            for mi in range(record['n_injection_models']):
                case_path=a.extraction/record['control_id']/f'injection_{mi:02d}.npz'
                if sha(case_path)!=record['injection_sha256'][str(case_path)]:
                    raise ValueError('Injection artifact changed')
                hashes[str(case_path)]=sha(case_path)
                with np.load(case_path) as z:
                    open_truth[mi]=np.array([z[f'fixed_signal_raw_{j}'] for j in range(3)])
            for row in local:
                mi=int(row['model_index']);fi=FACTORS.index(float(row['kron_factor']))
                bi=list(EDGES[:-1]).index(float(row['r_min_kpc']))
                expected=float(open_truth[mi][fi,bi])
                row['expected_open_area_signal']=expected
                row['mask_aware_recovery_fraction']=(float(row['recovered_increment'])/expected
                    if row['recovered_increment'] not in ('',None) and np.isfinite(expected) and expected>0 else None)
            recovery.extend(local)
    grouped = defaultdict(list)
    for row in recovery:
        grouped[(row['cohort'],row['model_index'],row['kron_factor'],row['extraction'],row['r_min_kpc'])].append(row)
    recovery_groups = []
    for key, rows in grouped.items():
        eligible = [r for r in rows if r['area_50pct_eligible']=='True' and np.isfinite(float(r['recovery_fraction']))]
        bias = np.array([float(r['recovery_fraction'])-1 for r in eligible])
        def finite_values(name):
            return [float(r[name]) for r in eligible if r[name] not in ('',None) and np.isfinite(float(r[name]))]
        fixed=finite_values('fixed_support_fraction');common=finite_values('common_support_fraction')
        mask_aware=finite_values('mask_aware_recovery_fraction')
        item = dict(cohort=key[0], model_index=int(key[1]), kron_factor=float(key[2]), extraction=key[3],
                    r_min_kpc=float(key[4]), r_max_kpc=float(rows[0]['r_max_kpc']),
                    scale_kpc=float(rows[0]['scale_kpc']), axis_ratio=float(rows[0]['axis_ratio']),
                    amplitude_at_100kpc=float(rows[0]['amplitude_at_100kpc']), n_total=len(rows), n_eligible=len(eligible),
                    median_fractional_bias=float(np.median(bias)) if len(bias) else None,
                    p16_fractional_bias=float(np.quantile(bias,.16)) if len(bias) else None,
                    p84_fractional_bias=float(np.quantile(bias,.84)) if len(bias) else None,
                    median_fixed_support_fraction=float(np.median(fixed)) if fixed else None,
                    median_common_support_fraction=float(np.median(common)) if common else None,
                    mask_aware_median_fractional_bias=float(np.median(mask_aware))-1 if mask_aware else None,
                    n_common_support_eligible=len(common))
        item['status'] = ('insufficient_area_support' if len(eligible)<10 else
                          'pass' if abs(item['median_fractional_bias']) <= .10 else 'bias_exceeds_10pct')
        recovery_groups.append(item)
    if recovery:
        prep.table(a.work/'recovery.tsv', recovery, list(recovery[0]))
        prep.table(a.work/'recovery_groups.tsv', recovery_groups, list(recovery_groups[0]))
    profiles = {key:np.array([r[key] for r in fields]) for key in fields[0]}
    ra = np.array([r['ra_deg'] for r in records]); dec = np.array([r['dec_deg'] for r in records])
    validation = np.array([r['n_injection_models'] > 0 for r in records])
    covariance_rows = []; archive = {}

    def run_cov(name, factor_index, mode, population, bins, width=5., offset=0.):
        coverage = profiles['coverage'][:,factor_index][:,bins]
        vals = profiles[mode][:,factor_index][:,bins]
        nominal = profiles['nominal_'+mode][:,factor_index][:,bins][:,:,bins]
        selected = population & np.all(coverage >= .5, axis=1) & np.isfinite(vals).all(axis=1)
        selected &= np.isfinite(nominal).all(axis=(1,2))
        info, arrays = covariance_diagnostics(vals[selected], spatial_labels((ra[selected]+offset)%360,dec[selected]+offset,width), nominal[selected])
        for key, matrix in arrays.items():
            archive[name+'__'+key] = matrix
        archive[name+'__selected_indices'] = np.flatnonzero(selected)
        item = dict(name=name,sample_kind='control_only',kron_factor=FACTORS[factor_index],extraction=mode,bin_indices=bins,
                    sky_cell_width_deg=width,sky_grid_offset_deg=offset,**info)
        covariance_rows.append(item)

    for fi, factor in enumerate(FACTORS):
        for mode in ('raw','corrected'):
            for population, selected in [('all',np.ones(len(records),bool)),('validation50',validation)]:
                run_cov(f'{population}_k{factor}_{mode}',fi,mode,selected,list(range(6)))
            for b in range(6):
                run_cov(f'bin{b}_k{factor}_{mode}',fi,mode,np.ones(len(records),bool),[b])
    for width in (2.5,10.):
        for mode in ('raw','corrected'):
            run_cov(f'all_k3.5_{mode}_cell{width}',1,mode,np.ones(len(records),bool),list(range(6)),width)
    for mode in ('raw','corrected'):
        run_cov(f'all_k3.5_{mode}_shifted',1,mode,np.ones(len(records),bool),list(range(6)),5.,2.5)
    for header_filter in ('i','i2'):
        population=np.array([r['header_filter']==header_filter for r in records])
        for mode in ('raw','corrected'):
            run_cov(f'filter_{header_filter}_k3.5_{mode}',1,mode,population,list(range(6)))
    for lo,hi in ((.3,.4),(.4,.55)):
        population=np.array([lo<=r['assigned_z']<hi for r in records])
        for mode in ('raw','corrected'):
            run_cov(f'assigned_z_{lo}_{hi}_k3.5_{mode}',1,mode,population,list(range(6)))
    if a.clusters:
        cp=a.clusters/'summary.json';cd=json.loads(cp.read_text());pp=a.clusters/'profiles.npz'
        if cd['status']!='complete' or sha(pp)!=cd['output_sha256']:
            raise ValueError('Cluster profile input is incomplete or changed')
        hashes[str(cp)]=sha(cp);hashes[str(pp)]=sha(pp)
        by_target={r['target_cluster_id']:i for i,r in enumerate(records)}
        linked=np.array([by_target[r['id']] for r in cd['records']])
        archive['paired_control_indices']=linked
        labels=spatial_labels([r['ra_deg'] for r in cd['records']],[r['dec_deg'] for r in cd['records']])
        with np.load(pp) as c:
            for fi,factor in enumerate(FACTORS):
                for mode in ('raw','corrected'):
                    for bins in [list(range(6))]+[[b] for b in range(6)]:
                        vals=(c[mode][:,fi]-profiles[mode][linked,fi])[:,bins]
                        coverage=np.minimum(c['coverage'][:,fi],profiles['coverage'][linked,fi])[:,bins]
                        nominal=(c['nominal_'+mode][:,fi]+profiles['nominal_'+mode][linked,fi])[:,bins][:,:,bins]
                        selected=np.all(coverage>=.5,axis=1)&np.isfinite(vals).all(axis=1)&np.isfinite(nominal).all(axis=(1,2))
                        info,matrices=covariance_diagnostics(vals[selected],labels[selected],nominal[selected])
                        name=f'paired50_k{factor}_{mode}_bins'+''.join(map(str,bins))
                        for key,value in matrices.items():archive[name+'__'+key]=value
                        archive[name+'__selected_indices']=np.flatnonzero(selected)
                        covariance_rows.append(dict(name=name,sample_kind='paired_cluster_minus_control',kron_factor=factor,
                            extraction=mode,bin_indices=bins,sky_cell_width_deg=5.,**info))
    np.savez_compressed(a.work/'covariances.npz', **archive)
    prep.write_json(a.work/'covariance_report.json',covariance_rows)
    # Existing native-profile/paired-random pipeline uses raw annular means.
    # The fitted outer plane is a sensitivity branch, never an activated correction.
    primary = next(r for r in covariance_rows if r['name']=='all_k3.5_raw')
    plane_covariance = next(r for r in covariance_rows if r['name']=='all_k3.5_corrected')
    boot = primary.get('iid_bootstrap_diagonal_ratio',[])
    cov_gate = bool(primary.get('covariance_sample_size_gate') and primary.get('rank')==6
                    and boot and max(abs(np.array(boot)-1))<=.1
                    and len(primary.get('folds',[]))==5
                    and primary.get('cv_fraction_within_1p96',0)>=.8)
    primary_recovery = [r for r in recovery_groups if r['kron_factor']==3.5 and r['extraction']=='raw']
    status_counts = {key:sum(r['status']==key for r in primary_recovery)
                     for key in ('pass','bias_exceeds_10pct','insufficient_area_support')}
    summary = dict(status='complete',n_control_fields=len(records),n_injection_fields=doc['n_injection_fields'],
        n_injection_images=doc['n_injection_images'], n_recovery_rows=len(recovery),
        n_recovery_groups=len(recovery_groups),primary_recovery_status_counts=status_counts,
        primary_covariance=primary, outer_plane_covariance_sensitivity=plane_covariance,
        primary_extraction='raw_native_profile; ensemble matched-control subtraction remains a separate stack operation',
        paired_cluster_input=str(a.clusters) if a.clusters else None,
        covariance_pilot_engineering_gate=cov_gate,
        coadd_injection_engineering_gate=bool(primary_recovery and all(r['status']=='pass' for r in primary_recovery)),
        source_sha256=hashes,script_sha256={str(Path(__file__)):sha(Path(__file__)),
            str(Path(__file__).with_name('prestack_core.py')):sha(Path(__file__).with_name('prestack_core.py'))},
        science_measurement=False,sky_correction_approved=False,controls_ready_for_science=False,
        limitations=[
            'Injections enter existing coadds: upstream HSC sky subtraction and PSF-wing transfer are not validated.',
            'Smooth observed-SB templates have no additional source photon noise; no total-flux or PSF deconvolution claim.',
            'Re-detected masks and outer planes are recomputed. Fixed/common support outputs are linear operator diagnostics only.',
            'The predeclared geometric-ring bias gate is retained. Open-area-normalized recovery is a secondary diagnostic separating footprint selection from pipeline transfer.',
            'All 798 controls feed covariance; complete-case selection at the unchanged 50% area floor changes the selected population.',
            'Control covariance is the control/background contribution. Paired50 diagnostics separately resample complete cluster-minus-control profiles.',
            'Paired50 scatter includes real cluster diversity; its nominal comparator assumes zero cross-image pixel covariance and is not a calibrated total science error.',
            'Empirical control scatter includes residual sources, cirrus, selection and survey systematics; not a pure pixel-noise kernel.',
            'Five-degree cells and 2.5/10-degree sensitivity cells are engineering blocks; residual inter-cell dependence remains unproven.',
            'Delete-cell jackknife uses unequal field counts and is reported as a sensitivity diagnostic.',
            'Spatial hold-out tests predict a held-out mean from the training mean; they test exchangeability, not an absolute zero sky.',
            'No inverse covariance is used to claim a science detection; satellite/ICL joint leakage, PSF and extinction gates remain.'])
    prep.write_json(a.work/'summary.json',summary)
    plot(a.work,recovery_groups,covariance_rows,archive)
    print(json.dumps(prep.plain({k:v for k,v in summary.items() if k not in ('source_sha256','script_sha256','limitations')})))


def plot(work,recovery,covariances,arrays):
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for h,color in [(100,'#286a9d'),(250,'#c06b20')]:
        for mode,style in [('raw','o-'),('corrected','s--')]:
            rows=sorted([r for r in recovery if r['cohort']=='area_qualified' and r['kron_factor']==3.5
                         and r['extraction']==mode and r['scale_kpc']==h and r['amplitude_at_100kpc']==3
                         and r['n_eligible']>0],key=lambda r:r['r_min_kpc'])
            axes[0].plot([(r['r_min_kpc']+r['r_max_kpc'])/2 for r in rows],
                         [1+r['median_fractional_bias'] for r in rows],style,color=color,label=f'h={h}, {mode}')
    axes[0].axhline(1,color='k',lw=.7);axes[0].axhspan(.9,1.1,alpha=.1,color='gray')
    axes[0].set_xlabel('Radius (kpc)');axes[0].set_ylabel('Re-detected increment / geometric truth');axes[0].legend(fontsize=8)
    axes[0].set_title('Validation40; eligible counts vary by bin')
    for ax,mode in zip(axes[1:],('raw','corrected')):
        prefix='all_k3.5_'+mode+'__'; key=prefix+'sample_covariance'
        if key in arrays:
            matrix=arrays[key]; diagonal=np.sqrt(np.diag(matrix)); corr=matrix/np.outer(diagonal,diagonal)
            im=ax.imshow(corr,vmin=-1,vmax=1,cmap='coolwarm');fig.colorbar(im,ax=ax,shrink=.8)
            primary=next(r for r in covariances if r['name']=='all_k3.5_'+mode)
            ax.set_title(f'{mode}: n={primary["n"]}, cells={primary["n_sky_cells"]}')
        ax.set_xlabel('Radial bin');ax.set_ylabel('Radial bin')
    fig.suptitle('Coadd-level pre-stack validation; science acceptance remains pending',fontsize=10)
    fig.tight_layout();fig.savefig(work/'prestack_validation.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    main()
