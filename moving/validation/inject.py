import sys,time,pickle,os,json; sys.path.insert(0,'/workspace/OGFinder')
import numpy as np
from scipy import ndimage as ndi
from moving import imaging as I, align as A, pipeline as P, util
C0=os.path.expanduser('~/.ds9/mast_cache/')
names=['j8pu38c7q','j8pu38caq','j8pu38ceq','j8pu38ciq']
seed=int(sys.argv[1]); nobj=int(sys.argv[2]); out=sys.argv[3]
rng=np.random.default_rng(seed)
chips=[]
for n in names: chips+=I.load_chips(C0+n+'_flc.fits')
A.align_field(chips)
psf=I.estimate_psf(chips,snr_min=6.0)
chips=[c for c in chips if c.name.endswith('[SCI,1]')]
offs,tm=P.hst_offsets(chips)
print('offs',offs)
center=(150.1375,2.3610)
# reference chip pixel region
c0=chips[0]; x,y=c0.wcs.all_world2pix([center[0]],[center[1]],0)
ra0,de0=center
a0,d0=np.radians(ra0),np.radians(de0)
ex_=np.array([-np.sin(a0),np.cos(a0),0]);ey_=np.array([-np.sin(d0)*np.cos(a0),-np.sin(d0)*np.sin(a0),np.cos(d0)])
fw=max(psf[1] if psf[0] is not None else 1.7, 0.085/chips[0].pixscale)
def gauss_psf(fwhm,size=41):
    r=size//2; yy,xx=np.mgrid[-r:r+1,-r:r+1]; s=fwhm/2.355
    g=np.exp(-(xx**2+yy**2)/(2*s*s)); return g/g.sum()
truth=[]
exp_t=chips[0].texp
for k in range(nobj):
    mag=rng.uniform(21.0,26.5)
    rate=np.exp(rng.uniform(np.log(1.5),np.log(40)))      # arcsec/h geocentric
    pa=rng.uniform(0,360); delta=rng.uniform(1.5,4.5)
    xi0=rng.uniform(-18,18); eta0=rng.uniform(-18,18)   # arcsec offsets from center, at t_ref
    truth.append(dict(mag=mag,rate=rate,pa=pa,delta=delta,xi0=xi0,eta0=eta0))
tref=np.mean(tm)
from moving.align import tangent_inverse
inj=[]
for ci,c in enumerate(chips):
    img=np.zeros(c.shape,np.float32)
    ex_i=ci
    for k,tk in enumerate(truth):
        # position at exposure mid (and trail endpoints) in arcsec tangent coords
        vx=tk['rate']*np.sin(np.radians(tk['pa']))*24.0; vy=tk['rate']*np.cos(np.radians(tk['pa']))*24.0      # arcsec/hour -> arcsec/day
        o=offs[ex_i]*206264.806
        nsub=40
        tt=np.linspace(-0.5,0.5,nsub)*c.texp/86400.0
        xs=[];ys=[]
        for dt in tt:
            tday=tm[ex_i]-tref+dt
            # geocentric motion minus parallax of the observer (k = 1/Delta)
            ox=o@ex_; oy=o@ey_
            # HST offset at sub-time: linearise using offset at mid exposure (HST moves ~ 7km/s: ok within 8 min? 3000 km -> parallax 0.2" at 1.5AU... include linear)
            xs.append(tk['xi0']+vx*tday-ox/tk['delta']); ys.append(tk['eta0']+vy*tday-oy/tk['delta'])
        # account HST motion during exposure with the offsets at neighbouring exposures (small): ignored, documented
        ra_s,de_s=tangent_inverse(np.array(xs)/3600.0*0+np.array(xs),np.array(ys),ra0,de0)
        px,py=c.wcs.all_world2pix(ra_s,de_s,0)
        flux=10**(-0.4*(tk['mag']-c.zp_ab))   # e-/s total
        # draw: accumulate sub-pixel shifted PSFs
        g=gauss_psf(fw if True else 1.7)
        tr=np.zeros(c.shape,np.float32)
        for xp,yp in zip(px,py):
            xi_,yi_=int(round(xp)),int(round(yp)); r=20
            if xi_<r or yi_<r or xi_>=c.shape[1]-r or yi_>=c.shape[0]-r: continue
            sh=ndi.shift(g,(yp-yi_,xp-xi_),order=1)
            tr[yi_-r:yi_+r+1,xi_-r:xi_+r+1]+=sh*(flux/nsub)
        # shot noise of the source (e-)
        n_e=tr*c.texp
        noise=rng.normal(0,np.sqrt(np.clip(n_e,0,None)))/c.texp
        img+=tr+noise.astype(np.float32)
        if ci==0:
            tk['x0']=float(px[nsub//2]);tk['y0']=float(py[nsub//2])
        tk.setdefault('pos',{})[ci]=(float(np.mean(px)),float(np.mean(py)),float(np.mean(ra_s)),float(np.mean(de_s)))
    c.data=c.data+img
res=P.detect_in_region(chips,center=center,half_pix=620,psf=psf)
dets,infos=res
pickle.dump((dets,truth,offs),open(out+'.pkl','wb'))
print(len(dets),'dets')
trs=P.link_detections(dets,offs,snr_min=8.0,tol_arcsec=1.0,min_exposures=3,max_per_exposure=600,max_tracklets=400)
print(len(trs),'tracklets')
# recovery
rec=[]
for k,tk in enumerate(truth):
    pos=tk['pos']
    det_hit=[]
    for ci in range(4):
        px,py,ra,de=pos[ci]
        best=None
        for d in dets:
            if d['ex']!=ci or d['sign']<0: continue
            s=np.hypot((d['ra']-ra)*np.cos(np.radians(de)),d['dec']-de)*3600
            if s<1.0 and (best is None or s<best[0]): best=(s,d)
        det_hit.append(best)
    nd=sum(b is not None for b in det_hit)
    # tracklet match
    tr_hit=None
    for t in trs:
        m=0
        for e,r,d in zip(t['ex'],t['ra'],t['dec']):
            px,py,ra,de=pos[e]
            if np.hypot((r-ra)*np.cos(np.radians(de)),d-de)*3600<1.0: m+=1
        if m>=3: tr_hit=t;break
    rec.append(dict(mag=tk['mag'],rate=tk['rate'],pa=tk['pa'],delta=tk['delta'],n_det=nd,
        det_astrom=[round(b[0],3) if b else None for b in det_hit],
        det_cls=[b[1]['cls'] if b else None for b in det_hit],det_snr=[round(b[1]['snr'],1) if b else None for b in det_hit],
        linked=tr_hit is not None,
        rate_fit=tr_hit['rate_ash'] if tr_hit else None, pa_fit=tr_hit['pa_deg'] if tr_hit else None))
json.dump(dict(rec=rec,n_tracklets=len(trs),n_false_tracklets=len(trs)-sum(r['linked'] for r in rec),n_dets=len(dets)),open(out,'w'),indent=1)
for r in rec: print(r)
