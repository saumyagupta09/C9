#!/usr/bin/env bash
# asaccus_tblastn_genomewide.sh -- settle the Asaccus_caudivolvulus C9 "absent"
# call with a genome-wide tblastn (curated C9 PROTEINS vs translated genome).
#
# Why: 07_tblastn_rescue.py only rescues exons around an ALREADY-located locus.
# Asaccus was scored absent (no blastn locus), so the sensitive protein search
# was never run there. This does an unrestricted, whole-genome tblastn instead.
#
# Read-only on the shared dataset; writes ONLY inside $OUT.
set -euo pipefail

ROOT=/media/ashutosh/disk3/Open_seminar/C9_related/Squamata/Squamata_genomes
ANALYSIS=$ROOT/C9_copy_number_2026
REFDIR=$ANALYSIS/runs/v2_ashu_C9A/refs
GENOME=$ROOT/genome_dir/Asaccus_caudivolvulus/GCA_042257475.1_ASM4225747v1_genomic.fna
OUT=$ANALYSIS/runs/v2_ashu_C9A/asaccus_tblastn_check
mkdir -p "$OUT"

# 1) Build a protein query: translate the curated C9A (22) + C9B (10) CDS.
#    Simple frame-1 translation of each curated CDS (they are full-length ORFs).
QUERY=$OUT/ashu_C9_proteins.fa
python3 - "$REFDIR/ashu_C9A.fa" "$REFDIR/ashu_C9B.fa" > "$QUERY" <<'PY'
import sys
codon = {}
bases = "TCAG"
aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i,a in enumerate(aa):
    codon["".join(bases[(i>>s)&3] for s in (4,2,0))] = a
def translate(seq):
    seq = seq.upper().replace("\n","")
    p = "".join(codon.get(seq[i:i+3],"X") for i in range(0,len(seq)-2,3))
    return p.rstrip("*")
def reads(fp):
    h=None; s=[]
    for ln in open(fp):
        ln=ln.rstrip()
        if ln.startswith(">"):
            if h: yield h,"".join(s)
            h=ln[1:]; s=[]
        else: s.append(ln)
    if h: yield h,"".join(s)
for fp in sys.argv[1:]:
    for h,s in reads(fp):
        print(">"+h); print(translate(s))
PY
echo "[query] $(grep -c '^>' "$QUERY") C9 proteins -> $QUERY"

# 2) tblastn against the genome. Use the existing nucleotide blastdb if present
#    (tblastn translates it on the fly); otherwise build one inside $OUT.
DB=$GENOME
if ! ls "${GENOME}".n?? >/dev/null 2>&1; then
  echo "[db] no blastdb beside genome; building a local copy in $OUT"
  makeblastdb -in "$GENOME" -dbtype nucl -out "$OUT/Asaccus_db" >/dev/null
  DB=$OUT/Asaccus_db
fi

# 3) Run. Lenient e-value so a real-but-divergent exon can surface; we judge by
#    bitscore/length/contig clustering, not by passing a strict cutoff.
OUTTSV=$OUT/asaccus_C9_tblastn.outfmt6.tsv
tblastn -query "$QUERY" -db "$DB" -evalue 1 -seg no -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out "$OUTTSV"

echo "[done] hits: $(wc -l < "$OUTTSV")  -> $OUTTSV"
echo "=== top hits by bitscore ==="
sort -t$'\t' -k12,12gr "$OUTTSV" | head -20
echo "=== distinct contigs hit ==="
cut -f2 "$OUTTSV" | sort -u