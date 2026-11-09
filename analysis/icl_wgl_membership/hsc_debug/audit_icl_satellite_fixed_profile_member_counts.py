#!/usr/bin/env python3
"""Check whether per-radius satellite-profile counts differ by ICL strength.

The audit uses the fixed complete-profile sample and median ICL groups from
the existing saved 94-pair table. At every radial bin it models the number of
satellite profiles contributing to the cluster-level mean as a function of
ICL group, redshift, and richness. Inference is spatial-block wild-bootstrap
and familywise-adjusted over 14 filter/radius cells. Counts are a support and
composition diagnostic, not a physical satellite abundance measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

DATA=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ASSOC=DATA/"icl_satellite_supported_94_association_2026-11-05_01"
SAT=DATA/"satellite_outer_profile_perbin_full_2026-10-18_01"
OUT=DATA/"icl_satellite_fixed_profile_member_counts_2026-11-06_01"
BANDS=("i","i2"); RADIAL=("5-10","10-15","15-20","20-30","30-40","40-50","50-75")
WEBB=np.asarray([-math.sqrt(1.5),-1.,-math.sqrt(.5),math.sqrt(.5),1.,math.sqrt(1.5)])

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def tsv(p):
    with Path(p).open(encoding="utf-8",newline="") as f: return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows):
    with Path(p).open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",extrasaction="raise");w.writeheader();w.writerows(rows)
def cluster_t(beta,resid,x,bread,blocks):
    groups=np.unique(blocks); n,k=x.shape
    if n<=k or len(groups)<2:return math.nan
    scores=np.vstack([x[blocks==g].T@resid[blocks==g] for g in groups])
    factor=(len(groups)/(len(groups)-1))*((n-1)/(n-k))
    v=factor*(bread@(scores.T@scores)@bread)[-1,-1]
    return float(beta[-1]/math.sqrt(v)) if v>0 else math.nan
def prepare(rows,band,rad,minblocks):
    y=np.asarray([r["counts"][rad] for r in rows],float)
    z=np.asarray([float(r["redshift"]) for r in rows]); lr=np.log([float(r["richness"]) for r in rows])
    if min(float(z.std()),float(lr.std()))<=0:raise ValueError("zero-variance covariate")
    high=np.asarray([r["high"] for r in rows],float)
    nuisance=np.column_stack((np.ones(len(rows)),(z-z.mean())/z.std(),(lr-np.mean(lr))/np.std(lr)))
    x=np.column_stack((nuisance,high))
    if np.linalg.matrix_rank(x)<4:raise ValueError(f"rank-deficient model {band}/{rad}")
    blocks=np.asarray([r["block"] for r in rows],str); unique=np.unique(blocks)
    if len(unique)<minblocks:raise ValueError(f"insufficient blocks {band}/{rad}")
    bread=np.linalg.inv(x.T@x); nb=np.linalg.inv(nuisance.T@nuisance)
    null=nuisance@(nb@nuisance.T@y); residual0=y-null
    leverage=np.sum((nuisance@nb)*nuisance,axis=1)
    adjusted=residual0/np.sqrt(np.maximum(1.-leverage,.05))
    beta=bread@x.T@y; resid=y-x@beta; t=cluster_t(beta,resid,x,bread,blocks)
    if not math.isfinite(t):raise ValueError(f"invalid CR1 statistic {band}/{rad}")
    return {"band":band,"rad":rad,"x":x,"bread":bread,"blocks":blocks,"null":null,
            "adjusted":adjusted,"beta":beta,"resid":resid,"t":t,"effect":float(beta[-1]),
            "se":abs(float(beta[-1])/t) if t else math.nan,"n":len(rows),"nblocks":len(unique),
            "mean_low":float(np.mean(y[high==0])) if False else None}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output",type=Path,default=OUT)
    ap.add_argument("--replicas",type=int,default=10000);ap.add_argument("--seed",type=int,default=20261112)
    a=ap.parse_args()
    if a.replicas<1000:ap.error("replicas must be >=1000")
    out=a.output.resolve()
    if out.exists():raise FileExistsError(f"refusing to overwrite {out}")
    assocpath=ASSOC/"supported_pairs.tsv"; aspath=ASSOC/"summary.json"; satpath=SAT/"cluster_profiles.tsv"; sspath=SAT/"summary.json"
    ass=json.loads(aspath.read_text()); ss=json.loads(sspath.read_text())
    if (ass.get("status")!="complete_diagnostic" or ass.get("science_measurement") is not False or ass.get("n_pairs")!=94
        or ass.get("output_sha256",{}).get(str(assocpath))!=sha(assocpath)
        or ss.get("status")!="complete" or ss.get("science_measurement") is not False
        or ss.get("output_sha256",{}).get(str(satpath))!=sha(satpath)):
        raise ValueError("pinned exploratory input status/hash mismatch")
    assoc=tsv(assocpath); profiles={}
    for r in tsv(satpath):
        if r.get("cohort")!="full_sample" or r.get("radial_bin_kpc") not in RADIAL:continue
        key=(str(r["cluster_id"]),r["filter"],r["radial_bin_kpc"])
        if key in profiles:raise ValueError(f"duplicate profile {key}")
        profiles[key]={"count":int(r["n_satellites"]),"sb":float(r["mean_delta_sb_njy_arcsec2"])}
    from evaluate_prestack_pair_support import pair_block_labels
    tar=[{"ra_deg":float(r["target_ra_deg"]),"dec_deg":float(r["target_dec_deg"])} for r in assoc]
    ctl=[{"ra_deg":float(r["control_ra_deg"]),"dec_deg":float(r["control_dec_deg"])} for r in assoc]
    ix,labels,*_=pair_block_labels(tar,ctl,np.ones(len(assoc),dtype=bool))
    if len(ix)!=len(assoc):raise ValueError("block helper changed row order")
    blocks={(r["cluster_id"],r["filter"]):str(b) for r,b in zip(assoc,labels)}
    samples={}; medians={}
    for band in BANDS:
        aa=[r for r in assoc if r["filter"]==band]
        sets=[{(r["cluster_id"],band) for r in aa if (r["cluster_id"],band,rad) in profiles and profiles[(r["cluster_id"],band,rad)]["count"]>0} for rad in RADIAL]
        common=set.intersection(*sets); sample=[]
        for r in aa:
            key=(r["cluster_id"],band)
            if key not in common:continue
            sample.append({**r,"block":blocks[key],"icl":float(r["icl_delta_200_300"]),
                           "counts":{rad:profiles[(key[0],band,rad)]["count"] for rad in RADIAL}})
        med=float(np.median([r["icl"] for r in sample])); medians[band]=med
        for r in sample:r["high"]=int(r["icl"]>=med)
        samples[band]=sample
    groups={(b,r):prepare(samples[b],b,r,6) for b in BANDS for r in RADIAL}
    allblocks=sorted(set(blocks.values()));rng=np.random.default_rng(a.seed)
    boot={k:np.empty(a.replicas) for k in groups}; max_t=np.empty(a.replicas)
    for j in range(a.replicas):
        weights=dict(zip(allblocks,rng.choice(WEBB,len(allblocks),replace=True)));mx=0.
        for k,g in groups.items():
            w=np.asarray([weights[b] for b in g["blocks"]]);ys=g["null"]+g["adjusted"]*w
            beta=g["bread"]@g["x"].T@ys;res=ys-g["x"]@beta
            value=cluster_t(beta,res,g["x"],g["bread"],g["blocks"])
            if not math.isfinite(value):raise ValueError(f"invalid bootstrap {k}")
            boot[k][j]=value;mx=max(mx,abs(value))
        max_t[j]=mx
    raw={k:(1+np.count_nonzero(np.abs(boot[k])>=abs(g["t"])))/(a.replicas+1) for k,g in groups.items()}
    ordered=sorted(raw,key=raw.get);holm={};run=0.
    for rank,k in enumerate(ordered):run=max(run,(len(ordered)-rank)*raw[k]);holm[k]=min(1.,run)
    rows=[]
    for key,g in groups.items():
        sample=samples[key[0]];vals=np.asarray([r["counts"][key[1]] for r in sample]);hi=np.asarray([r["high"] for r in sample],bool)
        crit=float(student_t.ppf(.975,g["nblocks"]-1))
        rows.append({"filter":key[0],"radial_bin_kpc":key[1],"n_fixed_clusters":len(sample),
                     "n_low_icl":int((~hi).sum()),"n_high_icl":int(hi.sum()),"n_connected_blocks":g["nblocks"],
                     "mean_n_satellite_profiles_low_icl":float(vals[~hi].mean()),"mean_n_satellite_profiles_high_icl":float(vals[hi].mean()),
                     "adjusted_high_minus_low_count":g["effect"],"cr1_se":g["se"],
                     "cr1_ci95_low":g["effect"]-crit*g["se"],"cr1_ci95_high":g["effect"]+crit*g["se"],
                     "t_cr1":g["t"],"p_wild_raw":raw[key],"p_holm_14":holm[key],
                     "p_maxT_14":(1+np.count_nonzero(max_t>=abs(g["t"])))/(a.replicas+1)})
    out.mkdir(parents=True,exist_ok=False);table=out/"fixed_profile_member_counts.tsv";write(table,rows)
    null=out/"maxT_null.tsv";write(null,[{"replica":i+1,"max_abs_t_14":float(v)} for i,v in enumerate(max_t)])
    products=[table,null]
    summary={"schema":"astrafek-icl-satellite-fixed-profile-member-counts-v1","status":"complete_exploratory",
             "science_measurement":False,"hsc_network_used":False,"new_image_download":False,
             "method":"OLS contributing satellite-profile count ~ fixed high/low ICL group + standardized redshift + standardized log richness; complete seven-bin sample fixed within each filter; CR1 connected 5-degree blocks; restricted HC2 Webb bootstrap-t; Holm/maxT across 14 tests",
             "fixed_sample":{b:{"n_clusters":len(samples[b]),"n_blocks":len({r['block'] for r in samples[b]}),"n_low":sum(not r['high'] for r in samples[b]),"n_high":sum(r['high'] for r in samples[b]),"icl_median":medians[b]} for b in BANDS},
             "replicas":a.replicas,"seed":a.seed,"minimum_maxT_p_14":float(min(r["p_maxT_14"] for r in rows)),
             "input_sha256":{str(p):sha(p) for p in (assocpath,aspath,satpath,sspath)},"runner_sha256":sha(Path(__file__).resolve()),
             "output_sha256":{p.name:sha(p) for p in products},
             "limitations":["counts are the number of satellite profiles contributing at each radius, not total member richness or complete satellite abundance",
                            "fixed sample is a selected complete-profile subset of 94 supported pairs, not representative of the 800-cluster parent",
                            "five-degree blocks approximate HSC survey fields",
                            "counts inherit mask, detection, profile support, and satellite-selection changes",
                            "the 14-test adjustment does not include prior exploratory research"]}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"summary":str(out/"summary.json"),"minimum_maxT_p_14":summary["minimum_maxT_p_14"],"results":rows},indent=2))
if __name__=="__main__":main()
