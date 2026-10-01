import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.environ.get('OGF_ROOT',_os.path.dirname(_HERE))
import subprocess,re,sys,glob,collections
ROOT=_ROOT
def procs(path):
    out=subprocess.check_output(['tclsh',_HERE+'/tclprocs.tcl',path]).decode().split('\n')
    return {l.split()[0]:(int(l.split()[1]),int(l.split()[2])) for l in out if l.strip()}
files=[ROOT+'/ds9/library/layout.tcl']+sorted(glob.glob(ROOT+'/ds9/library/ogf_*.tcl'))+sorted(glob.glob(ROOT+'/plugins/*/*.tcl'))+[ROOT+'/ds9/library/frame.tcl',ROOT+'/ds9/library/ds9.tcl']
owner={}   # proc -> (file)
bodies={}
for f in files:
    P=procs(f); L=open(f).read().split('\n')
    for n,(s,e) in P.items():
        owner[n]=f; bodies[n]='\n'.join(L[s-1:e])
tok=re.compile(r'[A-Za-z_][A-Za-z0-9_:]*')
refs=collections.defaultdict(set)
for n,b in bodies.items():
    for t in set(tok.findall(b)):
        if t in owner and t!=n: refs[n].add(t)
import json
json.dump({k:sorted(v) for k,v in refs.items()},open('/tmp/refs.json','w'))
json.dump(owner,open('/tmp/owner.json','w'))
