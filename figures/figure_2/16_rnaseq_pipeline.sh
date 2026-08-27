#!/usr/bin/env bash
# 16_rnaseq_pipeline.sh -- one merged BAM per species of C9-copy RNA-seq expression.
#
# Per species (resumable; run in tmux): tries BioProjects largest-first until the
# primary copy (C9A) is expressed; if none express after MAX_PROJECTS, reports it.
# For each tried project:
#   * select <=4 liver/lung Illumina runs by size (rule below),
#   * download fastq from ENA,
#   * build a STAR index of OUR assembly, align all selected runs, merge -> 1 BAM,
#   * subset to the DAB2-FYB1 region (+index), drop the species exon BED beside it,
#   * quantify per copy (total / unique MAPQ255 / depth) and per-exon coverage,
#   * read-vs-assembly sequence check (bcftools variants over the C9 exons) and
#     splice-site check (STAR junctions vs our exon boundaries) -- REPORT ONLY:
#     never edits exon_seqs_final; discrepancies go to parsed/rnaseq_seq_confirm.tsv
#     and parsed/rnaseq_splice_check.tsv for later review,
#   * DELETE the STAR index + fastq + full BAM; keep only the region BAM (never deleted).
# Stops cleanly if free disk < MIN_FREE_GB. Tissue selection (<=4 runs, largest by size):
#   both tissues: liver>=lung -> 3 liver+1 lung ; lung>liver -> 2 lung+2 liver
#   only liver -> 4 liver ; only lung -> 4 lung ; <=4 total -> all.
set -uo pipefail
cd "$(dirname "$0")"; source ./00_config.sh

THREADS="${THREADS:-8}"                         # used for STAR genome index build
STAR_THREADS="${STAR_THREADS:-$(( THREADS>16 ? 16 : THREADS ))}"  # align: cap (high counts crash STAR)
SORT_RAM="${SORT_RAM:-40000000000}"             # STAR BAM-sort RAM cap (40G)
MAX_RUNS="${MAX_RUNS:-4}"
MIN_FREE_GB="${MIN_FREE_GB:-150}"
REGION_BUF="${REGION_BUF:-200000}"
EXPR_MIN="${EXPR_MIN:-10}"        # min C9A unique reads to call the primary copy "expressed"
MAX_PROJECTS="${MAX_PROJECTS:-2}" # if C9A not expressed, try up to this many BioProjects
OUTDIR="$RUN/rnaseq"; mkdir -p "$OUTDIR"
WORK="$RUN/rnaseq_work"; mkdir -p "$WORK"
EXP="$PARSED/rnaseq_expression_v2.tsv"
[ -s "$EXP" ] || echo -e "species\tcopy\tscaffold\tregion\ttotal_reads\tunique_q255\tmean_depth\texon_cov_frac\truns_used\ttissues" > "$EXP"
SEQ="$PARSED/rnaseq_seq_confirm.tsv"
[ -s "$SEQ" ] || echo -e "species\tcopy\tregion\tvariants_in_exons\tnote" > "$SEQ"
BAMRUNS="$PARSED/rnaseq_bam_runs.tsv"
[ -s "$BAMRUNS" ] || echo -e "species\tbam_file\tn_runs\trun_ids\ttissues\tbioproject\tprimary_expressed" > "$BAMRUNS"
NOEXPR="$PARSED/rnaseq_no_expression.tsv"
[ -s "$NOEXPR" ] || echo -e "species\tn_projects_tried\tnote" > "$NOEXPR"

free_gb() { df -BG --output=avail "$ROOT" | tail -1 | tr -dc '0-9'; }

# runs_all <species> -> "Run<TAB>size_MB<TAB>BioProject<TAB>tissue" for liver+lung
# (esearch -> temp file -> efetch; the direct pipe is unreliable in nested PS.)
runs_all() {
  local sp="$1" tis es org
  for tis in liver lung; do
    org="\"${sp//_/ }\"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform] AND ${tis}[All Fields]"
    es=$(mktemp); esearch -db sra -query "$org" </dev/null >"$es" 2>/dev/null
    efetch -format runinfo <"$es" 2>/dev/null \
      | awk -F, -v t="$tis" 'NR>1 && $1 ~ /^[SED]RR/ {print $1"\t"($8+0)"\t"$22"\t"t}'
    rm -f "$es"
  done
}

ena_fastq() {   # echo ftp urls (https) for a run
  curl -s -m 40 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" \
    | awk -F'\t' 'NR>1{print $NF}' | tr ';' '\n' | sed 's#^#https://#' | grep -i 'fastq.gz$'
}

# ---- main: iterate target species that have liver/lung data --------------
avail="${AVAIL:-$ANALYSIS/config/rnaseq_availability.tsv}"   # AVAIL= to reorder (e.g. miocene-first)
ONLY="${ONLY:-}"   # set ONLY=Species to process a single species (testing)
# ashu species + their genera: skip a miocene species whose genus is already
# represented in the ashu set (e.g. Anolis_sagrei <- Anolis_carolinensis).
ASHU_SP=$(grep -v '^[[:space:]]*$' "$ROOT/ashu_list" 2>/dev/null | sort -u)
ASHU_GENERA=$(printf '%s\n' "$ASHU_SP" | sed 's/_.*//' | sort -u)
while IFS=$'\t' read -r sp clade local nll nany; do
  [ "$sp" = species ] && continue
  [ -n "$ONLY" ] && [ "$sp" != "$ONLY" ] && continue
  [ "${nll:-0}" -gt 0 ] || continue
  genus="${sp%%_*}"
  if printf '%s\n' "$ASHU_GENERA" | grep -qx "$genus" && ! printf '%s\n' "$ASHU_SP" | grep -qx "$sp"; then
    echo "[skip: genus $genus already in ashu] $sp"; continue; fi
  if ls "$OUTDIR/$sp/"*.C9region.bam >/dev/null 2>&1; then echo "[skip done] $sp"; continue; fi
  fg=$(free_gb); if [ "${fg:-0}" -lt "$MIN_FREE_GB" ]; then echo "[STOP] free ${fg}GB < ${MIN_FREE_GB}GB -- clear space and re-run"; break; fi

  genome=$(awk -F'\t' -v s="$sp" '$1==s{print $6}' "$TARGETS")
  [ -s "$genome" ] || { echo "[no genome] $sp"; continue; }
  echo "=== $sp (free ${fg}GB) ==="
  sw="$WORK/$sp"; rm -rf "$sw"; mkdir -p "$sw/fastq" "$sw/idx" "$sw/aln"
  # all liver/lung runs grouped by BioProject, planned (project order + tissue rule)
  runs_all "$sp" > "$sw/allruns.tsv"
  python3 "$ANALYSIS/scripts/select_runs.py" < "$sw/allruns.tsv" > "$sw/plan.tsv"
  [ -s "$sw/plan.tsv" ] || { echo "[no runs] $sp"; rm -rf "$sw"; continue; }

  # DAB2-FYB1 region + primary copy (C9A) -- computed once
  read -r scaf dab fyb < <(awk -F'\t' -v s="$sp" '$1==s && $7!="NA" && $8!="NA"{print $3,$7,$8; exit}' "$PARSED/synteny_order.tsv")
  if [ -n "${dab:-}" ] && [ -n "${fyb:-}" ]; then
    lo=$(( (dab<fyb?dab:fyb) - REGION_BUF )); hi=$(( (dab>fyb?dab:fyb) + REGION_BUF ))
  else
    read -r scaf lo hi < <(awk -F'\t' -v s="$sp" '$1==s{if(mn==""||$5<mn)mn=$5; if($6>mx)mx=$6; sc=$4} END{print sc, mn-500000, mx+500000}' "$PARSED/master_copies.tsv")
  fi
  (( lo<1 )) && lo=1; region="${scaf}:${lo}-${hi}"
  prim=$(awk -F'\t' -v s="$sp" 'NR>1 && $1==s && $4=="C9A"{print "C9A"; exit}' "$PARSED/master_copies.tsv")
  [ -z "${prim:-}" ] && prim=$(awk -F'\t' -v s="$sp" 'NR>1 && $1==s{print $4; exit}' "$PARSED/master_copies.tsv")
  mkdir -p "$OUTDIR/$sp"; cp -f "$RUN/beds/${sp}.C9_exons.bed" "$OUTDIR/$sp/" 2>/dev/null || true

  # STAR index once per species (reused across projects), deleted at the end
  STAR --runMode genomeGenerate --genomeDir "$sw/idx" --genomeFastaFiles "$genome" \
       --runThreadN "$THREADS" --genomeSAindexNbases 14 --limitGenomeGenerateRAM 64000000000 >/dev/null 2>"$sw/star_index.log" || { echo "  STAR index FAIL"; rm -rf "$sw"; continue; }

  expressed=0; ntried=0
  for rank in $(cut -f1 "$sw/plan.tsv" | sort -un); do
    ntried=$((ntried+1)); [ "$ntried" -gt "$MAX_PROJECTS" ] && break
    proj=$(awk -F'\t' -v r="$rank" '$1==r{print $2; exit}' "$sw/plan.tsv")
    mapfile -t PRUNS < <(awk -F'\t' -v r="$rank" '$1==r{print $3"\t"$4}' "$sw/plan.tsv")
    echo "  project $ntried [$proj]: ${#PRUNS[@]} runs"
    rm -rf "$sw"/aln/*                                      # clear previous project's alignments (incl _STARtmp)
    used=(); utis=(); bams=()
    for pr in "${PRUNS[@]}"; do
      run=$(echo "$pr"|cut -f1); tis=$(echo "$pr"|cut -f2)
      # download with up to 3 tries; verify each fastq.gz is intact (network is flaky)
      ok=0
      for try in 1 2 3; do
        rm -f "$sw/fastq/${run}"*.fastq.gz
        mapfile -t urls < <(ena_fastq "$run")
        [ "${#urls[@]}" -gt 0 ] || break
        for u in "${urls[@]}"; do wget -q -c -t 3 -T 60 -P "$sw/fastq" "$u" || true; done
        ok=1; for f in "$sw/fastq/${run}"*.fastq.gz; do [ -f "$f" ] && gzip -t "$f" 2>/dev/null || ok=0; done
        [ "$ok" = 1 ] && break || echo "    download try $try incomplete/corrupt for $run"
      done
      r1=$(ls "$sw/fastq/${run}"*_1.fastq.gz 2>/dev/null | head -1)
      r2=$(ls "$sw/fastq/${run}"*_2.fastq.gz 2>/dev/null | head -1)
      [ -z "$r1" ] && r1=$(ls "$sw/fastq/${run}"*.fastq.gz 2>/dev/null | head -1)
      if [ "$ok" = 1 ] && [ -n "$r1" ]; then
        rf=("$r1"); [ -n "$r2" ] && rf=("$r1" "$r2")
        STAR --genomeDir "$sw/idx" --readFilesIn "${rf[@]}" --readFilesCommand zcat \
             --runThreadN "$STAR_THREADS" --limitBAMsortRAM "$SORT_RAM" --outBAMsortingThreadN 4 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "$sw/aln/${run}." >/dev/null 2>"$sw/aln/${run}.log" \
          && { bams+=("$sw/aln/${run}.Aligned.sortedByCoord.out.bam"); used+=("$run"); utis+=("$tis"); } \
          || echo "    STAR failed for $run (see $sw/aln/${run}.log)"
      else echo "    no usable fastq for $run (ENA download)"; fi
      rm -f "$sw/fastq/${run}"*.fastq.gz
    done
    [ "${#bams[@]}" -gt 0 ] || { echo "    no alignments (project $proj)"; continue; }
    full="$sw/full.bam"
    if [ "${#bams[@]}" -eq 1 ]; then cp "${bams[0]}" "$full"; else samtools merge -f -@ "$THREADS" "$full" "${bams[@]}"; fi
    samtools index "$full"
    used_join=$(IFS=_; echo "${used[*]}"); used_str=$(IFS=,; echo "${used[*]}"); utis_str=$(IFS=,; echo "${utis[*]}")
    done_bam="$OUTDIR/$sp/${sp}.${used_join}.C9region.bam"
    samtools view -b -@ "$THREADS" "$full" "$region" > "$done_bam"; samtools index "$done_bam"
    rm -f "$full" "$full.bai"
    # region splice junctions kept beside the BAM (small)
    cat "$sw"/aln/*SJ.out.tab 2>/dev/null | awk -F'\t' -v sc="$scaf" -v lo="$lo" -v hi="$hi" 'BEGIN{OFS="\t"}
      $1==sc && $2>=lo && $3<=hi {k=$1 FS $2 FS $3; if(!(k in s))s[k]=$0; u[k]+=$7; m[k]+=$8}
      END{for(k in s){split(s[k],a,FS); print a[1],a[2],a[3],a[4],a[5],a[6],u[k],m[k],a[9]}}' \
      > "$OUTDIR/$sp/${sp}.${used_join}.SJ.tab"

    # per-copy quantification (capture primary-copy unique reads)
    c9a_uniq=0
    while IFS=$'\t' read -r csp name cscaf cs ce; do
      reg2="${cscaf}:${cs}-${ce}"
      tot=$(samtools view -c -F 0x104 "$done_bam" "$reg2" 2>/dev/null)
      uq=$(samtools view -c -q 255 -F 0x104 "$done_bam" "$reg2" 2>/dev/null)
      dep=$(samtools depth -r "$reg2" "$done_bam" 2>/dev/null | awk '{s+=$3;n++} END{if(n)printf "%.1f",s/n; else print 0}')
      ecov=$(awk -F'\t' -v n="${sp}_${name}_" '$4 ~ n{print $1"\t"$2"\t"$3}' "$RUN/beds/${sp}.C9_exons.bed" 2>/dev/null \
             | while read -r c a b; do d=$(samtools depth -r "$c:$((a+1))-$b" "$done_bam" 2>/dev/null | awk '$3>0{n++} END{print n+0}'); echo "$d $((b-a))"; done \
             | awk '{cov+=$1; t+=$2} END{if(t)printf "%.2f",cov/t; else print 0}')
      echo -e "${sp}\t${name}\t${cscaf}\t${reg2}\t${tot:-0}\t${uq:-0}\t${dep:-0}\t${ecov:-0}\t${used_str}\t${utis_str}" >> "$EXP"
      [ "$name" = "$prim" ] && c9a_uniq=${uq:-0}
    done < <(awk -F'\t' -v s="$sp" 'NR>1 && $1==s{print $1"\t"$4"\t"$5"\t"$6"\t"$7}' "$PARSED/master_copies.tsv")

    pexp=$([ "${c9a_uniq:-0}" -ge "$EXPR_MIN" ] && echo yes || echo no)
    echo -e "${sp}\t$(basename "$done_bam")\t${#used[@]}\t${used_str}\t${utis_str}\t${proj}\t${pexp}" >> "$BAMRUNS"

    if [ "$pexp" = yes ]; then
      echo "    [expressed] $prim uniq=$c9a_uniq (project $proj)"
      python3 "$ANALYSIS/scripts/17_splice_check.py" "$RUN" "$sp" "$OUTDIR/$sp/${sp}.${used_join}.SJ.tab" || true
      while IFS=$'\t' read -r csp name cscaf cs ce; do
        nv=$(bcftools mpileup -f "$genome" -r "${cscaf}:${cs}-${ce}" "$done_bam" 2>/dev/null | bcftools call -mv 2>/dev/null | awk -F'\t' '!/^#/ && $6>20' | wc -l)
        echo -e "${sp}\t${name}\t${cscaf}:${cs}-${ce}\t${nv:-NA}\tQUAL>20" >> "$SEQ"
      done < <(awk -F'\t' -v s="$sp" 'NR>1 && $1==s{print $1"\t"$4"\t"$5"\t"$6"\t"$7}' "$PARSED/master_copies.tsv")
      expressed=1; break
    else
      echo "    [no expression] $prim uniq=$c9a_uniq < $EXPR_MIN (project $proj) -- trying next project"
    fi
  done
  [ "$expressed" -eq 1 ] || { echo "  [no_expression_found] $sp after $ntried project(s)"; echo -e "${sp}\t${ntried}\tC9A unique reads < ${EXPR_MIN} in all tried projects" >> "$NOEXPR"; }
  rm -rf "$sw"                                              # delete index + all transient
  echo "  [done] $sp"
done < "$avail"
echo "[16] finished. BAMs+BEDs in $OUTDIR/<species>/ ; expression -> $EXP"
