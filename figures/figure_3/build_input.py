#!/usr/bin/env python3
# Assemble the Fig-5-style C8/C9 nucleotide (CDS) input:
#   our 81 intact squamate C9 CDS + Sharma-2022 C8a/C8b + non-squamate C9 outgroups.
import glob, os, re
def readfa(p):
    d={}; n=None
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=""
        elif n: d[n]+=l.strip()
    return d
def strip_stop(s):
    s=s.upper()
    return s[:-3] if len(s)>=3 and s[-3:] in ("TAA","TAG","TGA") else s
out=[]; comp={"C9A":0,"C9B":0,"C9C":0}; outg=[]
for f in sorted(glob.glob("exonwise_sequence_final/*/*.cds.fa")):
    b=os.path.basename(f).replace(".cds.fa","")
    if any(x in b.lower() for x in ("pseudogene","remnant","readthrough","alt")): continue
    for n,s in readfa(f).items():
        cp="C9C" if n.endswith("C9C") else "C9B" if n.endswith("C9B") else "C9A" if n.endswith("C9A") else None
        if cp is None: continue
        comp[cp]+=1; out.append((n, strip_stop(s)))
nsq=len(out)
sh=readfa("data/Squmata_C9_functional_duplicated_C8A_C8B.fa")
nc8=nog=0
for n,cds in sh.items():
    if re.search(r'_C8[AB]$', n): out.append((n, strip_stop(cds))); nc8+=1; outg.append(n)
    elif n.endswith("_C9"):       out.append((n, strip_stop(cds))); nog+=1
with open("figure_3/c8_c9_cds.fna","w") as fh:
    for n,s in out: fh.write(f">{n}\n{s}\n")
open("figure_3/outgroup_C8.txt","w").write(",".join(outg)+"\n")
print(f"squamate C9 (ours): {nsq}  {comp}")
print(f"C8 (Sharma): {nc8}  |  C9 outgroup (Sharma): {nog}  |  TOTAL: {len(out)}")
