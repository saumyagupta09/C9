import os,subprocess,re,itertools
GC=os.path.dirname(os.path.abspath(__file__))
fa=os.path.join(GC,"all_finalized_cds.fa")
def read_fa(p):
    d={};n=None
    for ln in open(p):
        ln=ln.rstrip()
        if ln.startswith(">"):n=ln[1:].split()[0];d[n]=""
        elif n:d[n]+=ln.upper()
    return d
seqs=read_fa(fa)
# group by species
sp={}
for k in seqs:
    s=re.sub(r"_C9[ABC].*$","",k)
    sp.setdefault(s,[]).append(k)
multi={s:v for s,v in sp.items() if len(v)>=2}

def mafft(d):
    inp=os.path.join(GC,"tmp_in.fa"); out=os.path.join(GC,"tmp_out.fa")
    with open(inp,"w") as f:
        for k,v in d.items(): f.write(f">{k}\n{v}\n")
    subprocess.run(f"mafft --quiet --auto {inp} > {out}",shell=True,check=True)
    return read_fa(out)

def pid(a,b):
    m=t=0
    for x,y in zip(a,b):
        if x=="-" or y=="-": continue
        t+=1; m+= (x==y)
    return (100*m/t if t else 0), t

def windows(a,b,win=90,step=30):
    res=[]
    L=len(a)
    for i in range(0,L-win+1,step):
        wa=a[i:i+win]; wb=b[i:i+win]
        p,t=pid(wa,wb)
        res.append((i,p,t))
    return res

rows=[]
tracts=[]
for s,copies in sorted(multi.items()):
    aln=mafft({k:seqs[k] for k in copies})
    keys=list(aln)
    for k1,k2 in itertools.combinations(keys,2):
        a,b=aln[k1],aln[k2]
        bg,_=pid(a,b)
        ws=windows(a,b)
        if not ws: continue
        wvals=[w[1] for w in ws if w[2]>=60]
        # conversion tract: contiguous windows >=99% where background <96%
        hi=[(i,p) for i,p,t in ws if t>=60 and p>=99.0]
        # collapse contiguous hi windows into tracts
        tract_spans=[]
        if bg<96.0 and hi:
            cur=[hi[0][0],hi[0][0]]
            for pos,_ in hi[1:]:
                if pos-cur[1]<=60: cur[1]=pos
                else: tract_spans.append(tuple(cur)); cur=[pos,pos]
            tract_spans.append(tuple(cur))
        maxw=max(wvals) if wvals else 0; minw=min(wvals) if wvals else 0
        rows.append((s,k1,k2,f"{bg:.1f}",f"{minw:.0f}",f"{maxw:.0f}",
                     str(len(tract_spans)), ";".join(f"{x}-{y+90}" for x,y in tract_spans) or "-"))
        for x,y in tract_spans:
            tracts.append((s,k1,k2,x,y+90,f"bg={bg:.1f}"))

with open(os.path.join(GC,"paralog_divergence.tsv"),"w") as f:
    f.write("species\tcopy1\tcopy2\tCDS_bg_pid\tmin_win_pid\tmax_win_pid\tn_conv_tracts\ttract_spans(aln_coord,>=99%&bg<96%)\n")
    for r in rows: f.write("\t".join(r)+"\n")
with open(os.path.join(GC,"conversion_tracts.tsv"),"w") as f:
    f.write("species\tcopy1\tcopy2\ttract_start\ttract_end\tnote\n")
    for r in tracts: f.write("\t".join(map(str,r))+"\n")
print("multi-copy species:",len(multi))
print("paralog pairs:",len(rows))
print("pairs with >=1 conversion tract:",sum(1 for r in rows if r[6]!="0"))
print()
print("species\tpair\tbg_pid\tmin/max_win\tn_tracts")
for r in rows:
    flag=" <== TRACT" if r[6]!="0" else ""
    print(f"{r[0]:24} {r[1].split('_C9')[1]+'/'+r[2].split('_C9')[1]:8} bg={r[3]:5} win[{r[4]}-{r[5]}] tracts={r[6]}{flag}")
