#!/usr/bin/env python3
"""Reduced-shear sensitivity for the redshift-matched richness mass contrast.

Uses the saved 204-pair high-minus-low random-subtracted DeltaSigma profile,
its shared covariance, exact HSC source-pair weights, and calibrated source
n(z). It changes only the NFW forward model from weak shear to reduced shear.
This remains an exploratory Astrafek diagnostic, not a science measurement.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2

from fit_wl_richness_split_nfw_mass_contrast import (
    CONCENTRATIONS, EDGES_MPC, FIT_BINS, NZ_SACC, SPLIT_SUMMARY,
    load_pairs_and_provenance,
)
from fit_wl_random_stack_reduced_shear_sensitivity import beta_moments_by_pair
from fit_wl_random_stack_nfw_mass import (
    EXPECTED_FULL_LENS_SUMMARY, RANDOM_EXPORT_AUDIT, RANDOM_SUMMARY, TOMO_SUMMARY,
    LENS_PLAN,
)
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology, nfw_sigma
from s19a_shape_stack import (
    C_KMS, G_MPC, hsc_selection_bias, load_hsc_y3_nz_sacc,
    pair_geometry, tomographic_nz_lensing_efficiency,
)

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCEBIN1_SUMMARY = DATA / "wl_random_subtracted_richness_split_zbin1_excluded_2026-11-05_01/summary.json"
EXPECTED_SPLIT_SHA = "b38a78f6005dc17d8a12ddd6117debd63a005d0b21ff8dbda05961bd31eb65b7"
EXPECTED_SOURCEBIN1_SHA = "814e9162aad25cc6c80c65a19ddfc52b88f1acc6f7d324723253395cea456df6"
OUTPUT = DATA / "wl_richness_split_reduced_shear_2026-11-05_01"
SERIES_ORDER = 16


def beta_support_max_by_pair(lens_z, source_bin, z_grid, nz, cosmology):
    """Maximum geometric efficiency over nonzero tomographic n(z) support."""
    lens_z=np.asarray(lens_z,float); bins=np.asarray(source_bin,int)
    z=np.asarray(z_grid,float); nz=np.asarray(nz,float); chi_s=np.asarray(cosmology.comoving(z),float)
    unique,inverse=np.unique(lens_z,return_inverse=True); chi_l=np.asarray(cosmology.comoving(unique),float)
    table=np.zeros((len(unique),nz.shape[0]),float)
    for li,(zl,cl) in enumerate(zip(unique,chi_l)):
        beta=np.divide(chi_s-cl,chi_s,out=np.zeros_like(chi_s),where=(z>zl)&(chi_s>0))
        beta=np.clip(beta,0,1)
        for bi in range(nz.shape[0]):
            support=nz[bi]>0
            if np.any(support): table[li,bi]=float(np.max(beta[support]))
    out=np.zeros(len(lens_z),float); valid=(bins>=1)&(bins<=nz.shape[0])
    out[valid]=table[inverse[valid],bins[valid]-1]
    return out


def prepare_group_rows(pairs, lens_ids, z_grid, nz, cosmology, source_bin_min):
    lens_id = np.asarray(pairs["lens_id"]).astype(str)
    source_bin = np.asarray(pairs["source_zbin"], dtype=int)
    take = np.isin(lens_id, list(lens_ids)) & (source_bin >= source_bin_min) & (source_bin <= 4)
    sub = {key: np.asarray(value)[take] for key, value in pairs.items()}
    bins = source_bin[take]
    beta = tomographic_nz_lensing_efficiency(sub["lens_z"], bins, z_grid, nz, cosmology)
    moments = beta_moments_by_pair(sub["lens_z"], bins, z_grid, nz, cosmology, SERIES_ORDER + 1)
    beta_max = beta_support_max_by_pair(sub["lens_z"], bins, z_grid, nz, cosmology)
    radius, _ = pair_geometry(sub, cosmology)
    amp = np.hypot(sub["e1"], sub["e2"])
    shape_ok = ((sub["sigma_e"] >= 0) & (sub["sigma_e"] <= .4) & (amp < 2))
    source_ok = shape_ok & (sub["mag_a10"] <= 25.5) & (sub["resolution"] >= .3)
    q_geom = 4 * math.pi * G_MPC / C_KMS**2 * np.asarray(cosmology.angular_diameter(sub["lens_z"]))
    usable = (np.asarray(sub["b_mode_pass"], dtype=bool) & source_ok & (beta > 0)
              & (moments[:, 1] > 0) & (radius >= EDGES_MPC[0]) & (radius <= EDGES_MPC[-1]))
    rows = []
    used_lenses = set()
    for b, (lo, hi) in enumerate(zip(EDGES_MPC[:-1], EDGES_MPC[1:])):
        in_bin = usable & (radius >= lo) & ((radius < hi) if b < len(EDGES_MPC)-2 else (radius <= hi))
        if not np.any(in_bin):
            raise ValueError(f"no usable source rows in radial bin {b}")
        weight = np.asarray(sub["weight"][in_bin], float)
        q = q_geom[in_bin] * beta[in_bin]
        lw = weight * q**2
        rms = np.asarray(sub["rms_e"][in_bin], float)
        response = 1 - np.sum(lw * rms**2) / np.sum(lw)
        if not np.isfinite(response) or response <= 0:
            raise ValueError(f"invalid shape response in radial bin {b}")
        sel = hsc_selection_bias(lw, sub["mag_a10"][in_bin], sub["resolution"][in_bin],
                                 sub["psf_shape11"][in_bin], sub["psf_shape22"][in_bin],
                                 sub["psf_shape12"][in_bin])
        m = np.asarray(sub["m"][in_bin], float)
        w = weight * (2 * response * (1 + m) * (1 + float(sel["m_sel"])))**2
        rows.append({"bin": b, "radius": radius[in_bin], "lens_z": np.asarray(sub["lens_z"][in_bin], float),
                     "geometry": q_geom[in_bin], "beta_ratios": moments[in_bin] / beta[in_bin, None],
                     "q_eff": q, "beta_support_max": beta_max[in_bin],
                     "weight": w, "weight_sum": float(np.sum(w * q**2))})
        used_lenses.update(np.asarray(sub["lens_id"])[in_bin].astype(str))
    return rows, used_lenses


def model_profile(rows, cosmology, mass, concentration):
    model = np.full(len(rows), np.nan)
    max_kappa = 0.0
    for item in rows:
        z = item["lens_z"]
        hz = cosmology.H0 * np.sqrt(cosmology.Om * (1 + z)**3 + (1 - cosmology.Om))
        rho = 3 * hz**2 / (8 * math.pi * G_MPC)
        r200 = (3 * mass / (4 * math.pi * 200 * rho))**(1/3)
        rs = r200 / concentration
        dc = (200/3) * concentration**3 / (np.log1p(concentration) - concentration/(1+concentration))
        sigma, ds = nfw_sigma(item["radius"], dc * rho, rs)
        kappa_inf = np.asarray(sigma) * item["geometry"]
        if item["bin"] in FIT_BINS:
            max_kappa = max(max_kappa, float(np.max(kappa_inf * item["beta_support_max"])))
        factor = np.zeros_like(kappa_inf)
        for order in range(SERIES_ORDER + 1):
            factor += kappa_inf**order * item["beta_ratios"][:, order]
        reduced_ds = np.asarray(ds) * factor
        model[item["bin"]] = float(np.sum(item["weight"] * item["q_eff"]**2 * reduced_ds) / item["weight_sum"])
    return model, max_kappa


def fit_contrast(low, high, observed, cov, cosmology, concentration):
    idx = FIT_BINS
    y = np.asarray(observed)[idx]
    c = np.asarray(cov)[np.ix_(idx, idx)]
    if np.min(np.linalg.eigvalsh(c)) <= 0:
        raise ValueError("contrast covariance is not positive definite")
    inv = np.linalg.inv(c)
    cache = {}
    def prof(group, rows, lm):
        k=(group,float(lm))
        if k not in cache: cache[k]=model_profile(rows,cosmology,10**float(lm),concentration)
        return cache[k][0][idx]
    def q(logm):
        res=y-(prof("high",high,logm[1])-prof("low",low,logm[0]))
        return float(res@inv@res)
    from scipy.optimize import least_squares
    chol=np.linalg.cholesky(c)
    def residual(lm): return np.linalg.solve(chol,y-(prof("high",high,lm[1])-prof("low",low,lm[0])))
    fit=least_squares(residual,[14.0,14.3],bounds=([12.5,12.5],[16,16]),max_nfev=140,
                      ftol=2e-6,xtol=2e-6,gtol=2e-6,x_scale="jac")
    if not fit.success: raise RuntimeError(f"reduced-shear fit failed: {fit.message}")
    best=fit.x; bestq=float(residual(best)@residual(best)); target=bestq+1
    lr=float(best[1]-best[0])
    def ratio_q(d):
        lo=max(12.5,12.5-d); hi=min(16,16-d)
        if lo>=hi:return float("inf")
        def obj(x):
            v=residual([x,x+d]); return float(v@v)
        r=minimize_scalar(obj,bounds=(lo,hi),method="bounded",options={"xatol":2e-4})
        return float(r.fun)
    ends=[]
    for sign in [-1,1]:
        inner=lr; step=.02; bound=-3.5 if sign<0 else 3.5
        while True:
            outer=lr+sign*step
            if (sign<0 and outer<=bound) or (sign>0 and outer>=bound): outer=bound
            if ratio_q(outer)>=target:
                ends.append(float(brentq(lambda d:ratio_q(d)-target,min(inner,outer),max(inner,outer),xtol=2e-4))); break
            if outer==bound: ends.append(bound); break
            inner=outer; step*=2
    max_k=max(model_profile(rows,cosmology,10**best[j],concentration)[1] for j,rows in enumerate([low,high]))
    return {"concentration_fixed":concentration,"low_richness_M200c_Msun":10**float(best[0]),
            "high_richness_M200c_Msun":10**float(best[1]),"mass_ratio_high_over_low":10**lr,
            "mass_ratio_conditional_delta_chi2_1_interval_log10":[*ends],"chi2":bestq,"dof":len(idx)-2,
            "max_kappa_times_beta_on_nz_support_in_fit_range":max_k,
            "relative_truncation_bound_on_nz_support":max_k**(SERIES_ORDER+1)/(1-max_k) if max_k < 1 else None,
            "reduced_shear_series_order":SERIES_ORDER,
            "observed_fit_profile_msun_mpc2":y.tolist(),"model_fit_profile_msun_mpc2":(prof("high",high,best[1])-prof("low",low,best[0])).tolist()}


def run(output: Path, exclude_source_bin_1=False):
    output=Path(output).resolve()
    if output.exists(): raise FileExistsError(output)
    split=json.loads(SPLIT_SUMMARY.read_text())
    if sha(SPLIT_SUMMARY)!=EXPECTED_SPLIT_SHA or split.get("science_measurement") is not False or len(split.get("matched_lens_pairs",[]))!=204:
        raise ValueError("richness split failed exact audited-input gate")
    if exclude_source_bin_1:
        sb=json.loads(SOURCEBIN1_SUMMARY.read_text())
        if sha(SOURCEBIN1_SUMMARY)!=EXPECTED_SOURCEBIN1_SHA or sb.get("science_measurement") is not False:
            raise ValueError("source-bin-1-excluded profile failed input gate")
        split=sb; source_min=2
    else: source_min=1
    pairs,_,_=load_pairs_and_provenance()
    low_ids={str(r["lower_lens_id"]) for r in json.loads(SPLIT_SUMMARY.read_text())["matched_lens_pairs"]}
    high_ids={str(r["upper_lens_id"]) for r in json.loads(SPLIT_SUMMARY.read_text())["matched_lens_pairs"]}
    z,nz=load_hsc_y3_nz_sacc(NZ_SACC); cosmo=Cosmology(H0=70,Om=.3)
    low,nlow=prepare_group_rows(pairs,low_ids,z,nz,cosmo,source_min)
    high,nhigh=prepare_group_rows(pairs,high_ids,z,nz,cosmo,source_min)
    if len(nlow)!=204 or len(nhigh)!=204: raise ValueError("not all matched lenses contribute estimator rows")
    observed=np.array([r["upper_minus_lower_delta_sigma_msun_mpc2"] for r in split["profile"]])
    cov=np.asarray(split["upper_minus_lower_covariance"],float)[np.ix_(np.arange(0,14,2),np.arange(0,14,2))]
    fits={str(c):fit_contrast(low,high,observed,cov,cosmo,c) for c in CONCENTRATIONS}
    out={"schema":"astrafek-wl-richness-split-reduced-shear-sensitivity-v1","status":"complete_exploratory_sensitivity",
         "science_measurement":False,"hsc_network_used":False,"images_read":False,"new_data_requested":False,
         "source_bins_used":list(range(source_min,5)),"matched_lenses_per_group":204,"radial_fit_range_mpc_proper":[.2,1.5],
         "model":f"two single-mass spherical NFW profiles jointly fit to high-minus-low random-subtracted DeltaSigma; reduced shear integrated over calibrated HSC Y3 n(z) beta moments through power {SERIES_ORDER+1} (series order {SERIES_ORDER})",
         "covariance":"saved shared 50-block spatial jackknife covariance; held fixed from the corresponding richness-split profile",
         "fits_by_fixed_concentration":fits,"limitations":["fixed single concentration and one effective mass per group do not represent the richness-selected mass distributions","source n(z) and shear calibration uncertainty are not propagated","the saved spatial covariance and matched-random depth/mask/PSF/sky remain unvalidated","conditional profile intervals are sensitivity diagnostics, not detection significance"],
         "input_sha256":{str(SPLIT_SUMMARY):sha(SPLIT_SUMMARY),str(NZ_SACC):sha(NZ_SACC),str(Path(__file__).resolve()):sha(Path(__file__).resolve())}}
    out["input_sha256"].update({str(path):sha(path) for path in
        (EXPECTED_FULL_LENS_SUMMARY,RANDOM_EXPORT_AUDIT,RANDOM_SUMMARY,TOMO_SUMMARY,LENS_PLAN/"row_export_plan.json")})
    if exclude_source_bin_1: out["input_sha256"][str(SOURCEBIN1_SUMMARY)]=sha(SOURCEBIN1_SUMMARY)
    output.mkdir(parents=True)
    (output/"summary.json").write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"source_bins":out["source_bins_used"],"fits":{c:{"ratio":v["mass_ratio_high_over_low"],"interval_log10":v["mass_ratio_conditional_delta_chi2_1_interval_log10"]} for c,v in fits.items()}},indent=2))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--output",type=Path,default=OUTPUT); ap.add_argument("--exclude-source-bin-1",action="store_true"); a=ap.parse_args();run(a.output,a.exclude_source_bin_1)

if __name__=="__main__": main()
