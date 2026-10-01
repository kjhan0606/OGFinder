import re,sys
def parse(path):
    """return list of (name,start,end) 1-based inclusive line ranges of top-level procs, start includes preceding comment block"""
    L=open(path).read().split('\n')
    out=[];i=0;n=len(L)
    while i<n:
        m=re.match(r'proc\s+(\S+)\s',L[i])
        if m:
            j=i
            while j<n and L[j]!='}': j+=1
            s=i
            k=i-1
            while k>=0 and L[k].startswith('#'): k-=1
            s=k+1
            out.append((m.group(1),s+1,j+1))
            i=j+1
        else: i+=1
    return out
if __name__=='__main__':
    p=sys.argv[1]
    L=open(p).read().split('\n')
    pr=parse(p)
    cov=set()
    for nm,s,e in pr: cov.update(range(s,e+1))
    # report non-proc, non-blank, non-comment top-level lines
    for i,l in enumerate(L,1):
        if i not in cov and l.strip() and not l.startswith('#'):
            print(i,l[:100])
