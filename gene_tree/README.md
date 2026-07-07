# gene_tree/ — ML gene trees per clade (for RELAX / codeml)
IQ-TREE 1.6.12 (`-m MFP -bb 1000`) trees built from the per-clade MACSE codon alignments
(`../alignment/<clade>/all_seq_<clade>_macse.aln`). One subfolder per clade: `<clade>.treefile` (ML tree,
used as the topology for selection tests), `.iqtree` (model/report), `.contree` (bootstrap consensus, if
>=4 taxa), plus the `<clade>_codon.aln` input. Includes the combined **Iguania_Anguimorpha** tree; no
separate Iguania/Anguimorpha and no global trees (per analysis design).
