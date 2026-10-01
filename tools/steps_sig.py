import sys,re,json,ast
def load(p):
    s=open(p).read()
    i=s.index('STEPS = [')
    # find matching: exec the DATA part
    ns={}
    j=s.index('\n]\n',i)+3
    exec('false=False;true=True;null=None\n'+s[s.index('SESSION_META'):j],ns)
    return ns['STEPS']
for p in sys.argv[1:]:
    for r in load(p):
        print(r['seq'],r['step'],r['class'],' '.join(map(str,r['argv_t'])), json.dumps(r.get('post'),sort_keys=True))
