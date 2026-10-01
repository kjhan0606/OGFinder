#!/usr/bin/env python3
"""Move Tcl procs out of ds9/library/layout.tcl into plugins/<id>/<file>.tcl.
usage: migrate.py MAPFILE   (lines: "<plugin>/<file.tcl> ProcName ProcName ..."; '#' comments)
Each proc block = the proc plus the immediately preceding '#' comment lines; blank line after it is dropped.
Blocks are appended verbatim (same order as in layout.tcl) to the target; manifest 'tcl' is updated by hand/by --manifest."""
import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.environ.get('OGF_ROOT',_os.path.dirname(_HERE))
import subprocess,sys,os,re,json
ROOT=_ROOT
LAY=ROOT+'/ds9/library/layout.tcl'
def procs(path):
    out=subprocess.check_output(['tclsh',_HERE+'/tclprocs.tcl',path]).decode().split('\n')
    return {l.split()[0]:(int(l.split()[1]),int(l.split()[2])) for l in out if l.strip()}
def main(mapfile):
    P=procs(LAY); L=open(LAY).read().split('\n')
    groups=[]
    for ln in open(mapfile):
        ln=ln.split('#')[0].strip()
        if not ln: continue
        t=ln.split(); groups.append((t[0],t[1:]))
    drop=set(); moved=0
    for tgt,names in groups:
        miss=[n for n in names if n not in P]
        if miss: sys.exit('missing in layout.tcl: %s'%miss)
        names=sorted(names,key=lambda n:P[n][0])
        path=os.path.join(ROOT,'plugins',tgt)
        os.makedirs(os.path.dirname(path),exist_ok=True)
        blocks=[]
        for n in names:
            s,e=P[n]
            blocks.append('\n'.join(L[s-1:e]))
            drop.update(range(s,e+1))
            if e<len(L) and L[e].strip()=='' : drop.add(e+1)
            moved+=e-s+1
        new=not os.path.exists(path)
        with open(path,'a') as f:
            if new:
                f.write('# %s -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).\n# Loaded through the "tcl" field of plugins/%s/plugin.json.\n\n'%(tgt,tgt.split('/')[0]))
            f.write('\n\n'.join(blocks)+'\n\n')
        print(tgt,len(names),'procs')
        pid,fn=tgt.split('/')
        mp=os.path.join(ROOT,'plugins',pid,'plugin.json'); m=json.load(open(mp))
        cur=m.get('tcl'); cur=[] if cur is None else ([cur] if isinstance(cur,str) else cur)
        if fn not in cur: cur.append(fn)
        m['tcl']=cur[0] if len(cur)==1 else cur
        if not m.get('required'): m['tcl_always']=1
        open(mp,'w').write(json.dumps(m,indent=1,ensure_ascii=False)+'\n')
    keep=[l for i,l in enumerate(L,1) if i not in drop]
    open(LAY,'w').write('\n'.join(keep))
    print('layout.tcl lines',len(L),'->',len(keep),'moved',moved)
main(sys.argv[1])
