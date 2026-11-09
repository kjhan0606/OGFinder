#!/usr/bin/env python3
"""Compare fixed-sample satellite radial profiles above/below median ICL.

Uses the same seven-bin complete-profile sample as
``analyze_icl_satellite_fixed_profile.py``. Within-filter median ICL groups are
fixed once, and adjusted high-minus-low satellite excess is estimated at each
radius with redshift/richness covariates. Restricted Webb wild-cluster
bootstrap-t and max-T inference cover all 14 filter/radius comparisons.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as student_t

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ASSOC = DATA / "icl_satellite_supported_94_association_2026-11-05_01"
SAT = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
OUT = DATA / "icl_satellite_fixed_profile_split_2026-11-06_01"
BANDS = ("i", "i2")
RADIAL = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
WEBB = np.asarray([-math.sqrt(1.5), -1., -math.sqrt(.5), math.sqrt(.5), 1., math.sqrt(1.5)])


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f: return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        w.writeheader(); w.writerows(rows)


def cluster_t(beta: np.ndarray, resid: np.ndarray, x: np.ndarray,
              bread: np.ndarray, labels: np.ndarray) -> float:
    groups = np.unique(labels); n,k=x.shape
    if n<=k or len(groups)<2: return math.nan
    scores=np.vstack([x[labels==g].T@resid[labels==g] for g in groups])
    corr=(len(groups)/(len(groups)-1))*((n-1)/(n-k))
    var=corr*(bread@(scores.T@scores)@bread)[-1,-1]
    return float(beta[-1]/math.sqrt(var)) if var>0 else math.nan


def prepare(rows: list[dict], band: str, radial: str, min_blocks: int) -> dict:
    y=np.asarray([r["profile"][radial] for r in rows],float)
    z=np.asarray([float(r["redshift"]) for r in rows]); lr=np.log([float(r["richness"]) for r in rows])
    zsd,rsd=float(z.std()),float(np.std(lr))
    if min(zsd,rsd)<=0: raise ValueError(f"zero-variance covariate: {band}/{radial}")
    hi=np.asarray([r["icl_group"] for r in rows],float)
    nuisance=np.column_stack((np.ones(len(rows)),(z-z.mean())/zsd,(lr-np.mean(lr))/rsd))
    x=np.column_stack((nuisance,hi))
    if np.linalg.matrix_rank(x)!=x.shape[1]: raise ValueError(f"rank deficient: {band}/{radial}")
    blocks=np.asarray([r["block"] for r in rows],str); unique=np.unique(blocks)
    if len(unique)<min_blocks: raise ValueError(f"insufficient blocks {band}/{radial}: {len(unique)}")
    bread=np.linalg.inv(x.T@x); nb=np.linalg.inv(nuisance.T@nuisance)
    nullfit=nuisance@(nb@nuisance.T@y); nullresid=y-nullfit
    lev=np.sum((nuisance@nb)*nuisance,axis=1)
    adjusted=nullresid/np.sqrt(np.maximum(1.-lev,.05))
    beta=bread@x.T@y; resid=y-x@beta; tval=cluster_t(beta,resid,x,bread,blocks)
    if not math.isfinite(tval): raise ValueError(f"invalid CR1 variance: {band}/{radial}")
    se=abs(float(beta[-1])/tval) if tval else math.nan
    return {"band":band,"radial":radial,"x":x,"bread":bread,"blocks":blocks,
            "nullfit":nullfit,"adjusted":adjusted,"beta":beta,"resid":resid,
            "t":tval,"effect":float(beta[-1]),"se":se,"n":len(rows),
            "nblocks":len(unique),"largest_block":max(np.count_nonzero(blocks==b) for b in unique),
            "n_high":int(hi.sum()),"n_low":int(len(hi)-hi.sum())}


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=OUT)
    p.add_argument("--replicas",type=int,default=10000)
    p.add_argument("--seed",type=int,default=20261111)
    p.add_argument("--minimum-blocks",type=int,default=6)
    a=p.parse_args()
    if a.replicas<1000 or a.minimum_blocks<6: p.error("replicas >=1000 and minimum-blocks >=6 required")
    out=a.output.resolve()
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    assoc_path=ASSOC/"supported_pairs.tsv"; assoc_sum_path=ASSOC/"summary.json"
    sat_path=SAT/"cluster_profiles.tsv"; sat_sum_path=SAT/"summary.json"
    ass=json.loads(assoc_sum_path.read_text()); ss=json.loads(sat_sum_path.read_text())
    if (ass.get("status")!="complete_diagnostic" or ass.get("science_measurement") is not False
            or ass.get("n_pairs")!=94 or ass.get("output_sha256",{}).get(str(assoc_path))!=sha(assoc_path)
            or ss.get("status")!="complete" or ss.get("science_measurement") is not False
            or ss.get("output_sha256",{}).get(str(sat_path))!=sha(sat_path)):
        raise ValueError("pinned input summaries/status/hash mismatch")
    assoc=read_tsv(assoc_path)
    if len(assoc)!=94 or len({(r["cluster_id"],r["filter"]) for r in assoc})!=94: raise ValueError("expected 94 unique pairs")
    prof={}
    for r in read_tsv(sat_path):
        if r.get("cohort")!="full_sample" or r.get("radial_bin_kpc") not in RADIAL: continue
        key=(str(r["cluster_id"]),r["filter"],r["radial_bin_kpc"])
        if key in prof: raise ValueError(f"duplicate profile {key}")
        if int(r["n_satellites"])>0: prof[key]=float(r["mean_delta_sb_njy_arcsec2"])

    from evaluate_prestack_pair_support import pair_block_labels
    t=[{"ra_deg":float(r["target_ra_deg"]),"dec_deg":float(r["target_dec_deg"])} for r in assoc]
    c=[{"ra_deg":float(r["control_ra_deg"]),"dec_deg":float(r["control_dec_deg"])} for r in assoc]
    ix,labels,*_=pair_block_labels(t,c,np.ones(len(assoc),dtype=bool))
    if len(ix)!=len(assoc): raise ValueError("block helper changed row order")
    block={(r["cluster_id"],r["filter"]):str(b) for r,b in zip(assoc,labels)}
    samples={}; medians={}
    for band in BANDS:
        rows=[r for r in assoc if r["filter"]==band]
        supports=[{(r["cluster_id"],band) for r in rows if (r["cluster_id"],band,rad) in prof} for rad in RADIAL]
        common=set.intersection(*supports)
        sample=[]
        for r in rows:
            key=(r["cluster_id"],band)
            if key not in common: continue
            sample.append({**r,"block":block[key],"icl":float(r["icl_delta_200_300"]),
                           "profile":{rad:prof[(key[0],band,rad)] for rad in RADIAL}})
        if len(sample)<10: raise ValueError(f"fixed-profile support too small in {band}")
        median=float(np.median([r["icl"] for r in sample])); medians[band]=median
        for r in sample: r["icl_group"]=int(r["icl"]>=median)
        if min(sum(r["icl_group"]==0 for r in sample),sum(r["icl_group"]==1 for r in sample))<4:
            raise ValueError(f"median split too imbalanced in {band}")
        samples[band]=sample
    groups={(b,r):prepare(samples[b],b,r,a.minimum_blocks) for b in BANDS for r in RADIAL}
    global_blocks=sorted(set(block.values())); rng=np.random.default_rng(a.seed)
    draws={k:np.empty(a.replicas) for k in groups}; max_t=np.empty(a.replicas)
    for j in range(a.replicas):
        weights=dict(zip(global_blocks,rng.choice(WEBB,len(global_blocks),replace=True))); mx=0.
        for k,g in groups.items():
            w=np.asarray([weights[b] for b in g["blocks"]]); ys=g["nullfit"]+g["adjusted"]*w
            beta=g["bread"]@g["x"].T@ys; residual=ys-g["x"]@beta
            ts=cluster_t(beta,residual,g["x"],g["bread"],g["blocks"])
            if not math.isfinite(ts): raise ValueError(f"invalid bootstrap t {k}")
            draws[k][j]=ts; mx=max(mx,abs(ts))
        max_t[j]=mx
    raw={k:(1+np.count_nonzero(np.abs(draws[k])>=abs(g["t"])))/(a.replicas+1) for k,g in groups.items()}
    order=sorted(raw,key=raw.get); holm={}; running=0.
    for rank,k in enumerate(order): running=max(running,(len(order)-rank)*raw[k]); holm[k]=min(1.,running)
    rows=[]
    for k,g in groups.items():
        crit=float(student_t.ppf(.975,g["nblocks"]-1))
        rows.append({"filter":k[0],"radial_bin_kpc":k[1],"n_fixed_clusters":g["n"],
                     "n_low_icl":g["n_low"],"n_high_icl":g["n_high"],"n_connected_blocks":g["nblocks"],
                     "icl_split_median_njy_arcsec2":medians[k[0]],"adjusted_high_minus_low_njy_arcsec2":g["effect"],
                     "cr1_se":g["se"],"cr1_ci95_low":g["effect"]-crit*g["se"],"cr1_ci95_high":g["effect"]+crit*g["se"],
                     "t_cr1":g["t"],"p_wild_raw":raw[k],"p_holm_14":holm[k],
                     "p_maxT_14":(1+np.count_nonzero(max_t>=abs(g["t"])))/(a.replicas+1)})
    out.mkdir(parents=True,exist_ok=False); table=out/"fixed_sample_icl_split_profiles.tsv"; write_tsv(table,rows)
    null=out/"maxT_null.tsv"; write_tsv(null,[{"replica":i+1,"max_abs_t_14":float(v)} for i,v in enumerate(max_t)])
    fig,axes=plt.subplots(2,1,figsize=(7.4,7.0),sharex=True,sharey=False,constrained_layout=True)
    center=np.asarray([7.5,12.5,17.5,25,35,45,62.5]); half=np.asarray([2.5,2.5,2.5,5,5,5,12.5])
    support={b:len(samples[b]) for b in BANDS}
    for ax,b in zip(axes,BANDS):
        rr=[r for r in rows if r["filter"]==b]; y=np.asarray([r["adjusted_high_minus_low_njy_arcsec2"] for r in rr])
        lo=np.asarray([r["cr1_ci95_low"] for r in rr]); hi=np.asarray([r["cr1_ci95_high"] for r in rr])
        ax.axhline(0,color="0.35",lw=1,zorder=0)
        ax.errorbar(center,y,yerr=np.vstack((y-lo,hi-y)),xerr=half,fmt="o",color="#235789",ecolor="#607d9b",capsize=3,ms=5)
        ax.text(.02,.94,f"{b}: fixed N={support[b]}",transform=ax.transAxes,ha="left",va="top",fontsize=10)
        ax.set_ylabel("adjusted high−low satellite excess\n[nJy arcsec⁻²]")
        extent=float(np.nanmax(np.abs(np.r_[lo,hi]))); ax.set_ylim(-1.08*extent,1.08*extent); ax.grid(axis="y",alpha=.2)
    axes[-1].set_xlabel("Satellite-centered radius [proper kpc]")
    fig.savefig(out/"icl_satellite_fixed_sample_icl_split.png",dpi=180); fig.savefig(out/"icl_satellite_fixed_sample_icl_split.pdf"); plt.close(fig)
    products=[table,null,out/"icl_satellite_fixed_sample_icl_split.png",out/"icl_satellite_fixed_sample_icl_split.pdf"]
    summary={"schema":"astrafek-icl-satellite-fixed-profile-split-v1","status":"complete_exploratory",
             "science_measurement":False,"hsc_network_used":False,"new_image_download":False,
             "method":"within-filter median ICL split fixed across all seven satellite radial bins; adjusted high-minus-low means from OLS with standardized z/log-richness; CR1 connected target/control 5-degree blocks; restricted HC2 Webb wild-cluster bootstrap-t; Holm/maxT across 14 tests",
             "fixed_support":{b:{"n_clusters":len(samples[b]),"n_blocks":len({r['block'] for r in samples[b]}),
                                  "n_low_icl":sum(r['icl_group']==0 for r in samples[b]),
                                  "n_high_icl":sum(r['icl_group']==1 for r in samples[b]),"icl_median":medians[b]} for b in BANDS},
             "replicas":a.replicas,"seed":a.seed,
             "minimum_maxT_p_14":float(min(r["p_maxT_14"] for r in rows)),
             "input_sha256":{str(p):sha(p) for p in (ASSOC/"supported_pairs.tsv",ASSOC/"summary.json",SAT/"cluster_profiles.tsv",SAT/"summary.json")},
             "runner_sha256":sha(Path(__file__).resolve()),"output_sha256":{p.name:sha(p) for p in products},
             "limitations":["selected complete-profile intersection of prior 94-pair sample; not representative of all 800 images",
                             "fixed clusters do not imply constant member detection or masks across radius",
                             "equal-area 5-degree blocks approximate HSC survey field dependence",
                             "intervals omit source Poisson, coadd weights, sky transfer, extended PSF, completeness, and injection-response uncertainties",
                             "familywise adjustment covers the declared 14 tests, not prior exploratory research"]}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"summary":str(out/"summary.json"),"fixed_support":summary["fixed_support"],
                      "minimum_maxT_p":min(r["p_maxT_14"] for r in rows)},indent=2))


if __name__=="__main__": main()
