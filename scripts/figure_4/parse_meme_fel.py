#!/usr/bin/env python3
import json, os, csv, sys
SR=sys.argv[1]
REF="Pituophis_catenifer_C9A"
# approximate C9 domain intervals (reference codon positions, ~583-aa squamate C9)
DOMAINS=[("signal_peptide",1,22),("TSP1_N",23,83),("LDLRA",84,125),("MACPF",126,430),
         ("EGF",431,520),("TSP1_C",521,999)]
def domain(pos):
    for n,a,b in DOMAINS:
        if a<=pos<=b: return n
    return "?"
def readaln(p):
    d={};cur=None
    for l in open(p):
        if l.startswith(">"): cur=l[1:].strip(); d[cur]=""
        else: d[cur]+=l.strip().upper()
    return d
def ref_map(alnfa):
    d=readaln(alnfa); ref=d.get(REF)
    if ref is None: return None
    m={}; rp=0
    ncod=len(ref)//3
    for c in range(ncod):
        cod=ref[c*3:c*3+3]
        if cod!="---" and "-" not in cod: rp+=1; m[c+1]=rp   # aln codon (1-based) -> ref codon pos
        elif cod!="---": rp+=1; m[c+1]=rp
        else: m[c+1]=None   # gap in ref -> insertion
    return m
DS={"snakes_C9A":"snakes_C9A_relaxready.fasta","duplicated_snakes_C9A_C9B":"dupsnakes_C9A_C9B_relaxready.fasta",
    "colubrid_C9A_C9B_C9C":"colubrid_ABC_relaxready.fasta"}
memerows=[]; felrows=[]; dom_counter={}
for ds,fa in DS.items():
    rmap=ref_map(os.path.join(SR,"RELAX",ds,fa))
    # MEME
    mj=json.load(open(os.path.join(SR,"MEME",ds+".MEME.json")))
    for i,row in enumerate(mj["MLE"]["content"]["0"],start=1):
        p=row[6]
        if p is not None and p<=0.05:
            rp=rmap.get(i); dom=domain(rp) if rp else "insertion"
            memerows.append({"dataset":ds,"aln_codon":i,"ref_codon":rp if rp else "","domain":dom,
                             "beta_plus":round(row[3],3),"p_value":f"{p:.3g}","branches_sel":row[7]})
            if rp: dom_counter[(ds,"MEME",dom)]=dom_counter.get((ds,"MEME",dom),0)+1
    # FEL
    fj=json.load(open(os.path.join(SR,"FEL",ds+".FEL.json")))
    for i,row in enumerate(fj["MLE"]["content"]["0"],start=1):
        alpha,beta,p=row[0],row[1],row[4]
        if p is not None and p<=0.05:
            kind="positive" if beta>alpha else "purifying"
            rp=rmap.get(i); dom=domain(rp) if rp else "insertion"
            felrows.append({"dataset":ds,"aln_codon":i,"ref_codon":rp if rp else "","domain":dom,
                            "alpha":round(alpha,3),"beta":round(beta,3),"p_value":f"{p:.3g}","selection":kind})
            if rp and kind=="positive": dom_counter[(ds,"FEL_pos",dom)]=dom_counter.get((ds,"FEL_pos",dom),0)+1
def dump(path,rows,fields):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
dump(SR+"/MEME_sites.csv",memerows,["dataset","aln_codon","ref_codon","domain","beta_plus","p_value","branches_sel"])
dump(SR+"/FEL_sites.csv",felrows,["dataset","aln_codon","ref_codon","domain","alpha","beta","p_value","selection"])
# domain summary
drows=[]
for (ds,test,dom),n in sorted(dom_counter.items()):
    drows.append({"dataset":ds,"test":test,"domain":dom,"n_sites":n})
dump(SR+"/sites_by_domain.csv",drows,["dataset","test","domain","n_sites"])
print(f"MEME sig sites: {len(memerows)}  FEL sig sites: {len(felrows)} (positive: {sum(1 for r in felrows if r['selection']=='positive')})")
print("=== sites by domain (positive/episodic) ===")
for r in drows: print(f"  {r['dataset']:26} {r['test']:8} {r['domain']:14} {r['n_sites']}")
