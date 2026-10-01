#!/usr/bin/env python3
"""rmprocs.py FILE name...  : delete the named top-level procs (with their preceding comment block) from FILE"""
import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.environ.get('OGF_ROOT',_os.path.dirname(_HERE))
import subprocess,sys
f=sys.argv[1]; names=set(sys.argv[2:])
P={l.split()[0]:(int(l.split()[1]),int(l.split()[2])) for l in subprocess.check_output(['tclsh',_HERE+'/tclprocs.tcl',f]).decode().split('\n') if l.strip()}
L=open(f).read().split('\n'); drop=set()
for n in names:
    if n not in P: sys.exit('no proc %s in %s'%(n,f))
    s,e=P[n]; drop.update(range(s,e+1))
    if e<len(L) and L[e].strip()=='': drop.add(e+1)
open(f,'w').write('\n'.join(l for i,l in enumerate(L,1) if i not in drop))
print(f,'removed',len(names),'procs')
