#!/usr/bin/env python3
"""S16A primary-definition stacked Delta Sigma of the 11 pilot lenses minus the 76 new randoms. Local only.

Pre-registered in HANDOVER 2026-11-27 (user approval "승인"). Fixed before the run:
- Sources: primary definition (standard selection, photoz_best - z_l > 0.2, Ephor_AB err95_min > z_l, dereddened i > 23).
  Main result all colours; comparison blue g-i <= random primary tertile; sensitivity without S23B_9307.
- Estimator per proper-radius bin (RADIAL_EDGES): Sigma_crit^-1 from photoz_best (flat LCDM H0=70, Om=0.3, proper),
  w_ls = w Sigma_crit^-2, DS = sum w_ls (e_t - c_t) Sigma_crit / [2 R (1 + K) sum w_ls], R = 1 - <rms_e^2>, K = <m>,
  both w_ls-weighted over the lens sample itself; the same for randoms (pooled, each with its matched lens z) and DS_lens -
  DS_random. Boost = 1. Cross component DS_x the same way.
- Photo-z multiplier: COSMOS (S17A Wide-COSMOS intake, shape x SOM weights) selected with the primary cuts per lens
  (photo_z > z_l + 0.2, PHOTOZ_ERR95_MIN > z_l, i_cmodel - a_i > 23 from the companion photometry), response of
  ogfkit point_redshift_response, stacked per bin with the observed lens w_ls sums; DS multiplied by 1/response.
  COSMOS object bootstrap (200 draws, seed 2026112701) for its uncertainty; the standard-cut baseline is reproduced as a check.
- Errors: group jackknife over the 11 groups (lens with its randoms), analytic shape noise, multiplier bootstrap.
- Fit range 0.4-3 Mpc (S/N reported); 0-0.4 Mpc reported but flagged (systematics not settled, HANDOVER limitations 1-4).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_outer_excess_split import colour, new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, load_new, load_new_centres, receipt, rows
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from select_s16a_random_centres import sha

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wl_photoz import point_redshift_response, stack_photoz_response  # noqa: E402
from ogfkit.wlcluster import C_KMS, G_MPC, Cosmology, tangential  # noqa: E402

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
INTAKE = DATA/'s16a_ephor_ab_cosmos_intake_2026-11-19_01'
EPHOR = DATA/'s16a_y1_ephor_ab_cosmos_plan_2026-11-18_01/results'
COMPANION = DATA/'s16a_cosmos_companion_plan_2026-11-24_01/results'
OUTLIER = 'S23B_9307'
SEED = 2026112701
BOOT = 200
FIT = (0.4, 3.0)


def cosmos():
    intake = json.loads((INTAKE/'summary.json').read_text())
    table = INTAKE/'cosmos_join.tsv'
    if sha(table) != intake['output_sha256'][table.name]:
        raise ValueError('COSMOS table changed')
    with table.open() as stream:
        rows_ = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review'] == 'True']
    if len(rows_) != intake['eligible_rows']:
        raise ValueError('COSMOS eligibility changed')
    extra = {}
    for tract in ('9812', '9813'):
        ename = f'pdf-s17a_wide-{tract}.cat.fits'
        if sha(EPHOR/ename) != intake['input_sha256'][ename]:
            raise ValueError('Ephor table changed '+ename)
        with fits.open(EPHOR/ename, memmap=True) as h:
            e = h[1].data
            err95 = dict(zip(np.asarray(e['ID'], np.int64).tolist(), np.asarray(e['PHOTOZ_ERR95_MIN'], float).tolist()))
        with fits.open(COMPANION/f'target_wide_s17a_{tract}.fits', memmap=True) as h:
            c = h[1].data
            f, a = np.asarray(c['icmodel_flux'], float), np.asarray(c['a_i'], float)
            bad = np.asarray(c['icmodel_flux_isnull'], bool) | np.asarray(c['a_i_isnull'], bool) | ~(f > 0)
            imag = np.where(bad, np.nan, -2.5*np.log10(np.where(bad, 1.0, f))-48.6-a)
            idmag = dict(zip(np.asarray(c['object_id'], np.int64).tolist(), imag.tolist()))
        extra[tract] = (err95, idmag)
    out = {k: [] for k in ('photo', 'ref', 'shape', 'som', 'err95', 'imag')}
    for r in rows_:
        err95, idmag = extra[r['tract']]
        oid = int(r['s17a_objid'])
        out['photo'].append(float(r['photo_z'])); out['ref'].append(float(r['reference_z']))
        out['shape'].append(float(r['shape_weight'])); out['som'].append(float(r['som_weight']))
        out['err95'].append(err95.get(oid, np.nan)); out['imag'].append(idmag.get(oid, np.nan))
    out = {k: np.array(v) for k, v in out.items()}
    return out, dict(rows=len(rows_), missing_err95=int(np.isnan(out['err95']).sum()), missing_imag=int(np.isnan(out['imag']).sum()))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    cosmo = Cosmology()
    centres, selection_sha = load_new_centres()
    new, _, provenance, _ = load_new(centres)
    pilot_centres, pilot, _, _ = bt.load_centres()
    lens_row = {c['matched_lens_id']: c for c in pilot_centres if c['role'] == 'lens'}
    groups = sorted(lens_row)
    new_rand = defaultdict(list)
    for cid, c in centres.items():
        new_rand[c['matched_lens_id']].append(cid)
    _, ppairs, _, _, _, _ = load_inputs(PLAN, OBJECT_QA)
    res, summ = receipt(SOURCE_PLAN, SOURCE_SHA)
    npairs = [{k: r[k] for k in PAIR_FIELDS} for r in rows(res, sorted(n for n in summ['output_sha256'] if 'pairs' in n))]
    pair = {(p['center_id'], p['object_id']): p for p in ppairs+npairs}
    pq, nq = pilot_quality(), new_quality()
    pz, pz_sha = load_pz()
    lens_ids = {g: lens_row[g]['center_id'] for g in groups}
    entries = {**{lens_ids[g]: pilot[lens_ids[g]] for g in groups}, **{r: new[r] for g in groups for r in new_rand[g]}}
    data = {}
    for cid, e in entries.items():
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        zs = g('photoz_best', q)
        keep = ~e['missing'] & (zs-zl > .2) & (g('ephor_ab_photoz_err95_min', pz)-zl > 0) & (g('i_mag', q)-g('a_i', q) > 23)
        ids = np.array(e['ids'])[keep]
        P = [pair[(cid, o)] for o in ids]
        f = {k: np.array([num(p[k]) for p in P]) for k in ('e1', 'e2', 'c1', 'c2', 'm', 'rms_e', 'sigma_e', 'weight')}
        xy = e['xy'][keep]
        chil, chis = float(cosmo.comoving(zl)), cosmo.comoving(zs[keep])
        sinv = 4*np.pi*G_MPC/C_KMS**2*chil/(1+zl)*(1-chil/chis)  # Mpc^2/Msun, proper
        et, ex = tangential(xy[:, 0], xy[:, 1], f['e1']-f['c1'], f['e2']-f['c2'])
        data[cid] = dict(zl=zl, rprop=np.hypot(*xy.T)*3.0/float(e['radius']), et=et, ex=ex, sinv=sinv, wls=f['weight']*sinv**2,
                         m=f['m'], rms=f['rms_e'], sig=f['sigma_e'], gi=np.array([colour(q[o]) for o in ids]))
    rand_cids = [r for g in groups for r in new_rand[g]]
    t1 = float(np.nanpercentile(np.concatenate([data[c]['gi'] for c in rand_cids]), 100/3))
    nb = len(RADIAL_EDGES)-1
    for d in data.values():
        d['bin'] = np.digitize(d['rprop'], RADIAL_EDGES)-1
        d['blue'] = d['gi'] <= t1

    # photo-z multiplier
    C, cinfo = cosmos()
    def responses(zl, idx=None):
        p, r, sh, so, e95, im = (C[k] if idx is None else C[k][idx] for k in ('photo', 'ref', 'shape', 'som', 'err95', 'imag'))
        std = p > zl+0.2
        prim = std & (e95 > zl) & (im > 23)
        return (point_redshift_response(zl, p[prim], r[prim], sh[prim], so[prim], cosmo)['response'],
                point_redshift_response(zl, p[std], r[std], sh[std], so[std], cosmo)['response'], int(prim.sum()), int(std.sum()))
    lens_resp = {}
    for g in groups:
        rp, rs, npm, nst = responses(float(lens_row[g]['z_lens']))
        lens_resp[g] = dict(lens_z=float(lens_row[g]['z_lens']), primary_response=rp, primary_multiplier=1/rp, standard_multiplier=1/rs,
                            cosmos_primary_rows=npm, cosmos_standard_rows=nst)
    rng = np.random.default_rng(SEED)
    nC = len(C['photo'])
    boot = [{g: responses(lens_resp[g]['lens_z'], idx)[0] for g in groups} for idx in (rng.integers(0, nC, nC) for _ in range(BOOT))]

    def sums(cids, k, cname):
        s = defaultdict(float)
        for c in cids:
            d = data[c]
            sel = (d['bin'] == k) & (d['blue'] if cname == 'blue' else True)
            w = d['wls'][sel]
            s['w'] += w.sum(); s['wt'] += (w*d['et'][sel]/d['sinv'][sel]).sum(); s['wx'] += (w*d['ex'][sel]/d['sinv'][sel]).sum()
            s['wm'] += (w*d['m'][sel]).sum(); s['wr'] += (w*d['rms'][sel]**2).sum(); s['n'] += sel.sum()
            s['var'] += (w**2/d['sinv'][sel]**2*(d['rms'][sel]**2+d['sig'][sel]**2)).sum()
        return s

    def ds(s):
        if s['w'] <= 0:
            return np.nan, np.nan, np.nan
        norm = 2*(1-s['wr']/s['w'])*(1+s['wm']/s['w'])*s['w']
        return s['wt']/norm, s['wx']/norm, np.sqrt(s['var'])/norm

    def profile(gs, cname, resp=None):
        resp = resp or {g: lens_resp[g]['primary_response'] for g in gs}
        out = []
        for k in range(nb):
            per = {g: sums([lens_ids[g]], k, cname)['w'] for g in gs}
            mult = 1/stack_photoz_response([resp[g] for g in gs], [per[g] for g in gs])['response'] if sum(per.values()) > 0 else np.nan
            L, R = sums([lens_ids[g] for g in gs], k, cname), sums([r for g in gs for r in new_rand[g]], k, cname)
            (lt, lx, le), (rt, rx, re) = ds(L), ds(R)
            out.append(dict(k=k, mult=mult, lens_t=lt*mult, lens_x=lx*mult, rand_t=rt*mult, rand_x=rx*mult,
                            dst=(lt-rt)*mult, dsx=(lx-rx)*mult, shape_se=np.hypot(le, re)*mult, n_lens=int(L['n']), n_rand=int(R['n'])))
        return out

    def jk(gs, cname, key):
        reps = np.array([[r[key] for r in profile([g for g in gs if g != d], cname)] for d in gs])
        n = len(gs)
        return np.sqrt((n-1)/n*np.nansum((reps-np.nanmean(reps, 0))**2, 0))
    results, cov_out = [], {}
    for label, gs, cname in (('all colours', groups, 'all'), ('blue', groups, 'blue'),
                             (f'all colours, no {OUTLIER}', [g for g in groups if g != OUTLIER], 'all')):
        base = profile(gs, cname)
        jt, jx = jk(gs, cname, 'dst'), jk(gs, cname, 'dsx')
        reps = np.array([[r['dst'] for r in profile(gs, cname, resp=b)] for b in boot])
        bse = np.nanstd(reps, 0, ddof=1)
        fit = [r['k'] for r in base if RADIAL_EDGES[r['k']] >= FIT[0]-1e-9 and RADIAL_EDGES[r['k']+1] <= FIT[1]+1e-9]
        tot = np.sqrt(jt**2+bse**2)
        snr = float(np.sqrt(sum((base[k]['dst']/tot[k])**2 for k in fit)))
        for r in base:
            k = r['k']
            results.append(dict(sample=label, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], in_fit_range=k in fit,
                                delta_sigma_t=r['dst'], jackknife_se=jt[k], shape_noise_se=r['shape_se'], multiplier_se=bse[k],
                                total_se=tot[k], delta_sigma_x=r['dsx'], jackknife_se_x=jx[k], lens_t=r['lens_t'], random_t=r['rand_t'],
                                multiplier=r['mult'], n_lens=r['n_lens'], n_random=r['n_rand']))
        cov_out[label] = dict(snr_fit_range_diag=snr, chi2_cross_fit_range=float(sum((base[k]['dsx']/jx[k])**2 for k in fit)), dof=len(fit))
    output.mkdir(parents=True)
    with (output/'delta_sigma.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(results)
    summary = dict(schema='astrafek-s16a-primary-delta-sigma-v1', science_measurement=True, calibration_applied=True, new_query=False,
                   units='Msun / pc^2 after division by 1e12 in the report; table in Msun/Mpc^2', radial_edges_mpc=list(RADIAL_EDGES),
                   fit_range_mpc=FIT, inner_flag='0-0.4 Mpc: systematics not settled (HANDOVER 2026-11-27 limitations 1-4)',
                   cosmology=dict(H0=70.0, Om=0.3), blue_gi_max=t1, lens_photoz=lens_resp, cosmos=cinfo, bootstrap=dict(n=BOOT, seed=SEED),
                   fit_range_statistics=cov_out, groups=groups, pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha,
                   receipt_summary_sha256=provenance,
                   limitations=['Boost = 1 (primary definition: loss, not dilution); crowding loss uncorrected.',
                                'Extra-light tangential bias upper bound ~10% at 0.4-1.7 Mpc for all colours (0 for blue), not corrected.',
                                'COSMOS calibration seeing transfer untested; COSMOS object bootstrap ignores field variance.',
                                'Diagonal S/N only; no covariance between bins; randoms pooled.',
                                'COSMOS i from S17A companion icmodel_flux (S16A uses meas icmodel_mag).'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'delta_sigma.tsv': sha(output/'delta_sigma.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return results, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_primary_delta_sigma_20261127_01'))
    results, summary = run(parser.parse_args().output)
    u = 1e-12
    for r in results:
        print(f"{r['sample']:24s} {r['r_lo']:>4}-{r['r_hi']:<4} {'fit' if r['in_fit_range'] else '   '} DS={r['delta_sigma_t']*u:7.2f}"
              f"±{r['total_se']*u:5.2f} (jk {r['jackknife_se']*u:5.2f} sn {r['shape_noise_se']*u:5.2f} pz {r['multiplier_se']*u:5.2f})"
              f" DSx={r['delta_sigma_x']*u:6.2f}±{r['jackknife_se_x']*u:5.2f} L={r['lens_t']*u:6.2f} R={r['random_t']*u:6.2f}"
              f" mult={r['multiplier']:.4f} n={r['n_lens']}/{r['n_random']}")
    print(json.dumps(dict(fit=summary['fit_range_statistics'], cosmos=summary['cosmos'], lens_photoz=summary['lens_photoz']), indent=1, default=float))
