#!/usr/bin/env python3
"""Derive feature inventory from layout.tcl / ogf_*.tcl (used to write docs/architecture.md)."""
import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.environ.get('OGF_ROOT',_os.path.dirname(_HERE))
import re, sys, collections, json, os
LIB=_ROOT+'/ds9/library'
FEATURES=[ # (feature, regex on proc name)
 ('sessionlog', r'^(CatalogPanelSession|OGFSess)'),
 ('ai-merge', r'^CatalogPanelAI(Merge|Bind|Unbind|Show|Next|Prev|Accept|Reject|Done)'),
 ('ai-services', r'^OGFAI'),
 ('moving', r'^OGFMov'),
 ('bands', r'^(OGFBands|CatalogPanelBands)'),
 ('mask', r'^(OGFMask|CatalogPanelMask)'),
 ('link', r'^OGFLink'),
 ('icl', r'^(CatalogPanelICL|CatalogPanelIcl)'),
 ('lsbg', r'^(CatalogPanelLSBG|CatalogPanelLsbg)'),
 ('cli-script', r'^(CatalogPanelExportCLIScript|CatalogPanelImportCLIScript|CatalogPanelCmdLog|ShellQuote)'),
 ('deconv', r'^CatalogPanel(Deconv|QuickDeconv|PSFDeconv|PSFSettings)'),
 ('star-psf', r'^CatalogPanel(FindStars|ShowStars|ClearStars|StarFinder|PSF|BuildPSF|BuildExtendedPSF|CheckSimAvail|SimPSF|ViewPSF|SavePSF|LoadPSF|StarPSF)'),
 ('galaxy-model', r'^CatalogPanel(Galaxy|Morph(Parse|Add|Color))'),
 ('separate/edit', r'^CatalogPanel(Separate|DeleteSelected|AddObject|GetSelectedSource)'),
 ('trim', r'^CatalogPanelTrim'),
 ('merge', r'^CatalogPanel(Merge|SetLogScale)'),
 ('photo-z', r'^CatalogPanel(PhotoZ)'),
 ('sed', r'^CatalogPanel(SED)'),
 ('bulge-disk', r'^CatalogPanel(BD|BulgeDisk)'),
 ('plot/viewer', r'^CatalogPanel(Plot|AnalysisViewer)'),
 ('measure-misc', r'^CatalogPanel(Morphometry|SersicFit|PSFPhotometry|MultiBand|CrowdedPhot|CrossMatch|SegmentationMap|Completeness)'),
 ('extract', r'^CatalogPanel(Extract|DualExtract|Param|SettingsDialog|AutoExtract)'),
 ('catalog-core', r'^CatalogPanel'),
 ('ds9-layout', r'.*'),
]
def procs(path):
    L=open(path,encoding='utf-8',errors='replace').read().split('\n')
    idx=[(i+1,m.group(1)) for i,l in enumerate(L) for m in [re.match(r'^proc (\S+)',l)] if m]
    out=[]
    for k,(ln,name) in enumerate(idx):
        end=(idx[k+1][0]-1) if k+1<len(idx) else len(L)
        out.append((name,ln,end,'\n'.join(L[ln-1:end])))
    return out
def feat(name):
    for f,r in FEATURES:
        if re.search(r,name): return f
allp=[]
for fn in ['layout.tcl']+sorted(x for x in os.listdir(LIB) if x.startswith('ogf_') and x.endswith('.tcl')):
    for name,a,b,body in procs(os.path.join(LIB,fn)):
        allp.append((fn,name,a,b,body,feat(name) if fn=='layout.tcl' or True else ''))
procnames={p[1]:p for p in allp}
keyre=re.compile(r'catpanel\(([A-Za-z0-9_,\-]+)\)')
def prefix(k):
    return k.split(',')[0]
res=collections.OrderedDict()
for fn,name,a,b,body,f in allp:
    if f=='ds9-layout' and fn=='layout.tcl' : pass
    d=res.setdefault(f,{'procs':[],'files':set(),'read':collections.Counter(),'write':collections.Counter(),'calls':collections.Counter(),'lines':0})
    d['procs'].append((fn,name,a,b)); d['files'].add(fn); d['lines']+=b-a+1
    for m in re.finditer(r'set\s+catpanel\(([A-Za-z0-9_,\-]+)\)',body): d['write'][prefix(m.group(1))]+=1
    for m in keyre.finditer(body): d['read'][prefix(m.group(1))]+=1
    for m in re.finditer(r'\b(CatalogPanel\w+|OGF\w+)\b',body):
        t=m.group(1)
        if t in procnames and feat(t)!=f and t!=name: d['calls'][feat(t)+':'+t]+=1
json.dump({f:{'lines':d['lines'],'n':len(d['procs']),'procs':d['procs'],'write':d['write'],'read':d['read'],
  'calls':d['calls']} for f,d in res.items()},open(_HERE+'/audit.json','w'),default=list,indent=1)
for f,d in res.items():
    rs=sorted(d['procs'],key=lambda p:(p[0],p[2]))
    print(f, d['lines'], 'lines', len(rs),'procs', rs[0][0],rs[0][2],'..',rs[-1][3])
