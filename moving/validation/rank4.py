import sys,pickle,json; sys.path.insert(0,'/workspace/OGFinder')
import numpy as np
from moving import pipeline as P
dets,truth,offs=pickle.load(open('inj4.json.pkl','rb'))
trs=P.link_detections(dets,offs,snr_min=8.0,tol_arcsec=1.0,min_exposures=3,max_per_exposure=600,max_tracklets=400)
rk=[]
for i,t in enumerate(trs):
    for k,tk in enumerate(truth):
        m=sum(np.hypot((r-tk['pos'][e][2])*np.cos(np.radians(tk['pos'][e][3])),d-tk['pos'][e][3])*3600<1.0 for e,r,d in zip(t['ex'],t['ra'],t['dec']))
        if m>=3: rk.append((i,k,round(truth[k]['mag'],1),t['n'],round(t['score'],1)))
print(len(trs),rk)
