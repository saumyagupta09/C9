# Alt-fragment assembly context (`c9_alt_fragments_assembly_context.csv`)

Assembly-level audit of the 6 "C9_alt / fragment" copies (the `fragment` category in
`c9_all_copies_status_table.csv`), to test whether they are genuine partial C9 genes or unanchored
assembly pieces.

## Columns
| column | meaning |
|---|---|
| assembly_used | the genome assembly we analysed (from config/targets.tsv) |
| assembly_level | NCBI level: Contig / Scaffold / Chromosome |
| n_scaffolds | number of scaffolds in that assembly |
| chromosome_level_assembly_available | is any chromosome-level assembly available for the species on NCBI |
| fragment_scaffold | the scaffold/contig the fragment hit sits on |
| scaffold_type | placed chromosome vs unplaced scaffold/contig |
| interpretation | genuine copy vs assembly artifact |

## Conclusion
All 6 fragments are on **unplaced / non-chromosomal scaffolds**. Only **Naja** and **Notechis** have
chromosome-level assemblies, and in both the fragment is still on an **unplaced scaffold** => the
fragments are **assembly artifacts** (unanchored duplicate/haplotig or segmental-dup pieces of the real
C9A/C9B loci), not independent partial C9 genes. They stay excluded from copy-number totals and analyses.
Notechis C9_alt is specifically the split C9B locus (same orphan scaffold NW_020719678 as C9B exons 2-6).
See `CORRECTION_LOG.md` (2026-07-02).
