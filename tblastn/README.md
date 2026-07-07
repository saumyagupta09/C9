# tblastn/

Per-species (`<species>/`) protein→translated-genome searches. outfmt 6 columns
throughout: `qseqid sseqid qlen length qstart qend sstart send pident evalue bitscore`.

| file | what it is |
|---|---|
| `ashu_C9_prot_on_<sp>_tblastn.out` | **Chromosome scan** — curated C9 proteins (`ashu_C9_prot.fa`, the C9A+C9B reference set) tblastn'd against the **whole genome** of `<sp>` (`-evalue 1e-5 -seg no`). Locates every C9-like locus genome-wide, independent of our copy calls. |
| `<sp>_on_<sp>_tblastn.out` | **Self tblastn** — the species' own **finalized** peptides (pooled `*.pep.fa` from `exonwise_sequence_final/<sp>/`) tblastn'd back against its **own** genome. Self-hits are ~100 % identity at the called loci; confirms each finalized copy maps cleanly back to the assembly. |
| `<sp>_rescue.tblastn.out` | tblastn protein rescue of exons missed by blastn, within the locus window (only present for copies that needed rescue; empties — copies already complete — were removed). Each block header records `subject=<seqid> region_offset=<lo>` so the searched region is reproducible from the genome. |
| `<sp>_C9Aprot_vs_C9Bregion.tblastn.tsv` | C9A-protein vs C9B-region audit (elapid/viper subset) used to confirm C9B identity. |
| `asaccus_C9_tblastn.outfmt6.tsv` | Asaccus genome-wide C9 search (C9-loss evidence; see also `asaccus_C9_loss/`). |

Notes:
- *Asaccus* (C9 lost) and *Dopasia* (C9A pseudogene) have a chromosome scan but **no**
  self tblastn — they carry no finalized peptide.
- Query/reference FASTAs: `ashu_C9_prot.fa` (in `miniprot/QUERY_ashu_C9_prot.faa`); per-copy
  finalized peptides are under `sequences/<species>/`.
