import sys,pickle,json; sys.path.insert(0,'/workspace/OGFinder')
import numpy as np
from moving import util, horizons as H, orbit as O, orbitfit as OF, identify, kepler as K
dets=pickle.load(open('dets_real.pkl','rb'))
ts,e=pickle.load(open('truth.pkl','rb'))
sel=[]
for k in range(4):
    c=[(np.hypot((d['ra']-e[k,0])*np.cos(np.radians(e[k,1])),d['dec']-e[k,1])*3600,d) for d in dets if d['ex']==k and d['sign']>0]
    c.sort(key=lambda x:x[0]); s,d=c[0]; sel.append((s,d))
print('astrometric offsets vs Horizons(-48) for 2015 BB89 (arcsec):',[round(s,3) for s,_ in sel])
tr=dict(t=[d['t'] for _,d in sel],ra=[d['ra'] for _,d in sel],dec=[d['dec'] for _,d in sel],sig=[max(d['sig_pos_arcsec'],0.1) for _,d in sel])
tr['sig']=[0.15]*4
idr=identify.identify_tracklet(tr,radius_arcmin=3.0,tol_arcsec=1.5,max_cand=12)
print('IDENT',idr['status'],[(m['name'],m['num'],round(m['max_sep'],2),m['matched']) for m in idr['matches'][:3]])
jd=util.utc_mjd_to_tdb_jd(np.array(tr['t'])); sat=H.vectors('-48',jd,center='500@399')[:,:3]*util.AU_KM
obs=O.Obs(tr['t'],tr['ra'],tr['dec'],tr['sig'],tr['sig'],['250']*4,sat_xyz_km=list(sat))
p=O.Propagator()
res=OF.fit_orbit(obs,p,ranging_samples=250,n_starts=4)
el=res['elements']; sg=res['sigmas']
js=H.sbdb('558706'); jel,_=H.sbdb_elements(js)
print('fit rms',res['rms_arcsec'],'chi2red',res['chi2_red'],res['orbit_class'],res['info'].get('ranging',{}).get('class_probs'))
out={}
for k in('a','e','i','om','w','ma','q'):
    print(k,round(el[k],4),'+-',round(sg[k],4),'JPL',round(jel[k],4)); out[k]=dict(fit=el[k],sig=sg[k],jpl=jel[k])
json.dump(dict(offsets=[s for s,_ in sel],ident=idr,elements=out,rms=res['rms_arcsec'],cls=res['orbit_class'],ranging=res['info'].get('ranging')),open('real_bb89.json','w'),indent=1,default=str)
