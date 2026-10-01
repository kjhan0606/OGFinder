import sys
def load(p):
    rows=[]
    for l in open(p,encoding='utf-8'):
        f=l.rstrip('\n').split('\t')
        while len(f)<5: f.append('')
        rows.append(f)
    return rows
base=load(sys.argv[1]); new=load(sys.argv[2])
REPLACED={'OGFMovTable':'shared time-domain table + kind filter (intentional, requested)',
 'CatalogPanelSettingsDialog':'OGFParamDialog extract (declarative dialog)',
 'CatalogPanelSeparateSettings':'OGFUISettings objects (chip menu Settings...)',
 'CatalogPanelDeconvSettings':'OGFParamDialog deconv (declarative dialog)'}
newcmds={}
for path,t,lab,cmd,var in new:
    if t=='cascade': continue
    newcmds.setdefault((cmd or var).strip(),[]).append(path+' | '+lab)
miss=[]; ok=0; skipped=0
for path,t,lab,cmd,var in base:
    if t=='cascade': skipped+=1; continue
    key=(cmd or var).strip()
    if not key and lab=='(none)': skipped+=1; continue   # disabled placeholder of a dynamic band menu
    if not key: miss.append((path,lab,'(no command)')); continue
    if key in REPLACED: ok+=1; continue
    if key in newcmds: ok+=1
    else: miss.append((path,lab,key))
print('baseline leaf entries:',len(base)-skipped,' reachable in new UI:',ok,' missing:',len(miss))
for m in miss: print('  MISSING',m)
# entries only in new
bk={(c or v).strip() for p,t,l,c,v in base if t!='cascade'}
extra=[(p,l) for p,t,l,c,v in new if t!='cascade' and (c or v).strip() not in bk]
print('entries only in new UI:',len(extra))
for e in extra: print('  NEW',e)
