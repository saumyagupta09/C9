#!/usr/bin/env bash
# Fig-5-style C8/C9 nucleotide ML tree: codon-aware alignment (protein-guided) + IQ-TREE.
set -euo pipefail
WORK=/media/ashutosh/disk3/Open_seminar/C9_related/Squamata/Squamata_genomes/C9_copy_number_2026/runs/v2_ashu_C9A/figure3_c8c9_tree
cd "$WORK"
IN=c8_c9_cds.fna
echo "[$(date)] start; iqtree=$(iqtree --version 2>&1 | head -1)" > run.log

# 1) translate CDS -> protein (in-frame)
python3 - "$IN" > c8_c9_prot.faa <<'PY'
import sys
b="TCAG"; aa="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
cod={"".join(b[(i>>s)&3] for s in (4,2,0)):a for i,a in enumerate(aa)}
def rd(p):
    d={}; n=None
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=""
        elif n: d[n]+=l.strip()
    return d
for n,s in rd(sys.argv[1]).items():
    s=s.upper(); print(">"+n); print("".join(cod.get(s[i:i+3],"X") for i in range(0,len(s)-2,3)))
PY

# 2) MAFFT align protein
mafft --auto --anysymbol --thread 8 c8_c9_prot.faa > c8_c9_prot.aln 2> mafft.log
echo "[$(date)] mafft done" >> run.log

# 3) back-translate protein alignment -> codon (nucleotide) alignment
python3 - c8_c9_prot.aln "$IN" > c8_c9_codon.aln <<'PY'
import sys
def rd(p):
    d={}; n=None
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=""
        elif n: d[n]+=l.strip()
    return d
aln=rd(sys.argv[1]); cds=rd(sys.argv[2])
for n,ap in aln.items():
    s=cds[n].upper(); j=0; out=[]
    for a in ap:
        if a=="-": out.append("---")
        else: out.append(s[j:j+3] if j+3<=len(s) else "---"); j+=3
    print(">"+n); print("".join(out))
PY
echo "[$(date)] backtrans done ($(grep -c '>' c8_c9_codon.aln) seqs)" >> run.log

# 4) IQ-TREE: ModelFinder + ML + ultrafast bootstrap 1000 (unrooted; root on C8 at visualization)
iqtree -s c8_c9_codon.aln -m MFP -bb 1000 -nt AUTO -ntmax 16 -pre c8_c9_tree >> run.log 2>&1
echo "[$(date)] DONE" >> run.log
