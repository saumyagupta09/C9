#!/bin/bash
P=/media/ashutosh/disk3/Open_seminar/C9_related/Squamata/Squamata_genomes/C9_copy_number_2026
SIF=$P/tools/hyphy.sif
SR=$P/runs/v2_ashu_C9A/Selection_analysis/RELAX
OUT=$P/runs/v2_ashu_C9A/Selection_analysis
mkdir -p $OUT/MEME $OUT/FEL
LOG=$OUT/meme_fel_progress.log; echo "MEME+FEL START $(date)" > "$LOG"
declare -A ALN=( [snakes_C9A]=snakes_C9A_relaxready.fasta [duplicated_snakes_C9A_C9B]=dupsnakes_C9A_C9B_relaxready.fasta [colubrid_C9A_C9B_C9C]=colubrid_ABC_relaxready.fasta )
declare -A TRE=( [snakes_C9A]=snakes_C9A.treefile [duplicated_snakes_C9A_C9B]=dupsnakes_C9A_C9B.treefile [colubrid_C9A_C9B_C9C]=colubrid_ABC.treefile )
for ds in snakes_C9A duplicated_snakes_C9A_C9B colubrid_C9A_C9B_C9C; do
  aln=$SR/$ds/${ALN[$ds]}; tree=$SR/$ds/${TRE[$ds]}
  apptainer exec --bind $P $SIF hyphy CPU=40 meme --alignment "$aln" --tree "$tree" --output $OUT/MEME/${ds}.MEME.json > $OUT/MEME/${ds}.meme.log 2>&1
  echo "MEME $ds $(date)" >> "$LOG"
  apptainer exec --bind $P $SIF hyphy CPU=40 fel --alignment "$aln" --tree "$tree" --output $OUT/FEL/${ds}.FEL.json > $OUT/FEL/${ds}.fel.log 2>&1
  echo "FEL  $ds $(date)" >> "$LOG"
done
echo "ALL_MEMEFEL_DONE $(date)" >> "$LOG"
