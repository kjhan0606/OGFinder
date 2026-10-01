import sys,time,json; sys.path.insert(0,'/workspace/OGFinder')
import numpy as np
from moving import obs as OBS, orbit as O, horizons as H, kepler as K, orbitfit as OF
des=sys.argv[1]; t0m=float(sys.argv[2]); t1m=float(sys.argv[3]); tag=sys.argv[4]
j=json.load(open('/tmp/obs%s.json'%des))[0]
ob=[o for o in OBS.parse_obs80(j['OBS80']) if o['stn'] not in ('500',)]
sel=[o for o in ob if t0m<o['mjd_utc']<t1m and (OBS.obscodes().get(o['stn'],(None,))[0] is not None or o['sat_xyz_km'] is not None)]
import os
thin=int(os.environ.get('THIN','1'))
sel=sel[::thin]
print(len(sel),'obs; stations',sorted(set(o['stn'] for o in sel)))
ra=[];de=[]
for o in sel:
    dr,dd=OBS.debias(o['ra'],o['dec'],o['mjd_utc'],o['cat'])
    ra.append(o['ra']-dr/np.cos(np.radians(o['dec']))/3600); de.append(o['dec']-dd/3600)
yr=int(2000+(t0m-51544)/365.25)
sig=[max(OBS.station_sigma(o['stn'],yr),0.1) for o in sel]
obs=O.Obs([o['mjd_utc'] for o in sel],ra,de,sig,sig,[o['stn'] for o in sel],sat_xyz_km=[o['sat_xyz_km'] for o in sel])
p=O.Propagator()
js=H.sbdb(des); el2,sg2=H.sbdb_elements(js)
jd0=el2['epoch_jd']
t=time.time()
res=OF.fit_orbit(obs,p,epoch_jd=jd0,ranging_samples=int(os.environ.get('NS','200')),n_starts=3,verbose=True)
print('fit time',time.time()-t,'rms',res['rms_arcsec'],'chi2red',res['chi2_red'],res['n_used'],res['n_total'],res['orbit_class'])
print('ranging',res['info'].get('ranging'))
e,sg=res['elements'],res['sigmas']
out={}
for k in ('a','e','i','om','w','ma'):
    d=e[k]-el2[k]
    if k in('om','w','ma'): d=(d+180)%360-180
    print(k,e[k],sg[k],el2[k],sg2.get(k),'dev/sigma_fit',d/sg[k] if sg[k]>0 else None)
    out[k]=dict(fit=e[k],sig=sg[k],jpl=el2[k],jpl_sig=sg2.get(k),dev_sigma=d/sg[k])
json.dump(dict(des=des,n_obs=len(sel),n_used=res['n_used'],rms=res['rms_arcsec'],chi2_red=res['chi2_red'],elements=out,arc_days=t1m-t0m,epoch_jd=jd0,cls=res['orbit_class']),open('val_%s.json'%tag,'w'),indent=1)
