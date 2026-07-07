#!/usr/bin/env python3
"""
Build the per-residue selection table for the C9 domain + MEME/FEL figure.
Maps every FEL/MEME site (alignment codon) onto the Pituophis catenifer C9A reference
residue, tags it by InterProScan domain, and records FEL dN-dS + class and MEME p + class.

Inputs (all in figure_4/):
  snakes_C9A_codon.aln            codon alignment (contains Pituophis_catenifer_C9A)
  snakes_C9A.FEL.sites.tsv        FEL per-site table  (site, alpha, beta, ..., p-value)
  snakes_C9A.MEME.sites.tsv       MEME per-site table (site, ..., beta+, ..., p-value)
  Pituophis_C9A_interproscan.tsv  (domain intervals hard-coded below, from this file)
Output:
  c9_selection_sites.csv          one row per reference residue
"""
import csv, sys

REF = "Pituophis_catenifer_C9A"
# dataset presets: tag -> input tsvs, codon alignment, output CSV  (default = snakes)
DATASETS = {
  "snakes":   dict(fel="figure_4/snakes_C9A.FEL.sites.tsv",
                   meme="figure_4/snakes_C9A.MEME.sites.tsv",
                   aln="figure_4/snakes_C9A_codon.aln",
                   out="figure_4/c9_selection_sites.csv"),
  "colubrid": dict(fel="figure_4/colubrid_C9A_C9B_C9C.FEL.sites.tsv",
                   meme="figure_4/colubrid_C9A_C9B_C9C.MEME.sites.tsv",
                   aln="figure_4/colubrid_ABC_codon.aln",
                   out="figure_4/c9_selection_sites_colubrid.csv"),
}
tag = sys.argv[1] if len(sys.argv) > 1 else "snakes"
cfg = DATASETS[tag]

# InterProScan domain intervals on the Pituophis C9A reference (see Pituophis_C9A_interproscan.tsv)
DOM = [("signal_peptide",1,20),("TSP1_N",38,91),("LDLRA",96,133),
       ("MACPF",135,495),("EGF",496,536),("TSP1_C",537,583)]
def domain_of(pos):
    for n,a,b in DOM:
        if a <= pos <= b: return n
    return "linker"

# --- read codon alignment, map alignment-codon index -> reference residue ---
seq={}; name=None
for l in open(cfg["aln"]):
    l=l.rstrip("\n")
    if l.startswith(">"): name=l[1:]; seq[name]=""
    elif name: seq[name]+=l
ref=seq[REF]
cods=[ref[i:i+3] for i in range(0,len(ref),3)]
aln2ref={}; r=0
for i,c in enumerate(cods, start=1):
    if c!="---":
        r+=1; aln2ref[i]=r
    else:
        aln2ref[i]=None

def load(path, cols):
    out={}
    for row in csv.DictReader(open(path), delimiter="\t"):
        out[int(row["site"])] = tuple(float(row[c]) for c in cols)
    return out
fel  = load(cfg["fel"],  ["alpha","beta","p-value"])
meme = load(cfg["meme"], ["beta+","p-value"])

rows=[]
for i in sorted(aln2ref):
    ref_pos=aln2ref[i]
    if ref_pos is None: continue
    a,b,fp = fel.get(i,(0.0,0.0,1.0))
    bp,mp  = meme.get(i,(0.0,1.0))
    fel_class = "positive" if (b>a and fp<=0.05) else ("purifying" if (b<a and fp<=0.05) else "ns")
    rows.append({
        "ref_pos":ref_pos, "domain":domain_of(ref_pos),
        "fel_alpha":round(a,4), "fel_beta":round(b,4), "fel_dNdS":round(b-a,4),
        "fel_p":fp, "fel_class":fel_class,
        "meme_betaplus":round(bp,4), "meme_p":mp,
        "meme_class":("episodic" if mp<=0.05 else "ns"),
    })
with open(cfg["out"],"w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
npos = sum(1 for r in rows if r["fel_class"]=="positive")
nmeme= sum(1 for r in rows if r["meme_class"]=="episodic")
print(f"wrote {cfg['out']}: {len(rows)} reference residues; "
      f"FEL positive={npos}, purifying={sum(1 for r in rows if r['fel_class']=='purifying')}, "
      f"MEME episodic={nmeme}")
