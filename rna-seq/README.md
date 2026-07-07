# rna-seq/ — C9-locus RNA-seq (expression)
Per-species IGV minilocus bundles: `<sp>.minilocus.bam`(+.bai) = RNA-seq reads over the C9 locus, with
`<sp>.minilocus.fa`(+.fai), `.exons.bed`, `.offset.txt` for direct IGV loading. 32 species (mostly liver,
some lung). Quantified per copy in `../tables/c9_orf_rnaseq_table.csv` (reads + locus coverage). Result:
C9A is the dominant transcript in snakes (C9B 1-3 orders of magnitude lower, C9C near-silent); in Anolis
the pattern is reversed (C9B dominant). Minor-copy counts are upper bounds (tandem-paralog read cross-mapping).
