# Gene-conversion check on FINALIZED C9 sequences (redo, 2026-06-25)

Input: the 82 genome-validated, synteny-named finalized copies
(`exonwise_sequence_final/`), pooled in `all_finalized_cds.fa`. 21 species carry
≥2 intact copies (29 paralog pairs incl. triplicates; remnant/pseudogene copies
included in the divergence scan, excluded from the global alignment if <1500 nt).

## Tests
1. **Per-species paralog divergence + tract scan** (`paralog_divergence.tsv`,
   `conversion_tracts.tsv`): mafft-align each species' copies; background pairwise
   identity + sliding 90-nt windows; a candidate tract = contiguous windows ≥99 %
   identity inside a pair whose CDS background is <96 %.
2. **Conversion-vs-conservation discriminator** (`conversion_vs_conservation.tsv`):
   one global alignment (80 copies); per window, compare a species' own-paralog
   identity to the best identity of that copy to ANY OTHER species' copy. A true
   conversion tract = own paralogs ≥98 % AND ≥6 % closer to each other than to any
   other species' copy (i.e. locally each other's closest relative despite global
   divergence). A conserved domain instead has some other species' copy equally close.
3. **Topology test** (`global_fasttree.nwk`): FastTree GTR on the global alignment;
   per multi-copy species, are its copies a sister clade (by-species = concerted
   evolution / very recent dup) or interspersed (by-copy/clade = orthologous)?

## Results
- **Background divergence:** most paralog pairs sit at **~81–85 % identity** (normal
  orthologous paralog divergence). A few colubrid within-triplicate A/B pairs
  (*Arizona, Hemorrhois, Pituophis, Pantherophis*) are ~94–95 %, and *Bothrops* A/B
  ~98 % — i.e. more recently duplicated, but still no conversion (below).
- **Tract scan** flagged localized ≥99 % windows in 9 pairs, BUT these **recur at the
  same relative CDS positions across unrelated species** (~510, ~900, ~1410 nt) =
  conserved C9 functional regions, not species-specific tracts.
- **Conversion-vs-conservation: 0 / 29 pairs** have a single conversion window. At
  every high-identity window some *other* species' copy is equally close → the windows
  are **conserved domains (purifying selection), not gene conversion.**
- **Topology: 21 / 21 multi-copy species have paralogs INTERSPERSED (by-copy), 0
  sister.** A species' C9A's nearest neighbour is almost always another species' copy
  (e.g. *Vipera berus* C9A ↔ *V. seoanei* C9A; *Pituophis* C9A ↔ *Pantherophis* C9A;
  *Anolis carolinensis* C9A ↔ *A. sagrei* C9B). Even the ~95 % colubrid A/B pairs are
  *less* similar to each other than to their cross-genus orthologs → the C9A/C9B split
  **predates** the colubrid/viperid/elapid radiations.

## Conclusion
**Gene conversion / concerted evolution is REJECTED on the finalized sequences.** The
two C9 copies are ancient, orthologously-inherited paralogs that diverge cleanly and
group by copy-type across species; high-identity windows are conserved functional
domains. The Sharma et al. 2022 "copies do not separate" (concerted-evolution) signal
was inflated by the tandem-copy extraction chimeras (see §1.10, §4e–4g); the
locus-anchored finalized set removes that artifact. Consistent with the earlier
corrected-set result (`gene_conversion_corrected/`, conversion rejected), now confirmed
with the more specific cross-species window test and the full finalized set.

Deliverables in this folder: all_finalized_cds.fa · intact_for_global.fa · global.aln
· global_fasttree.nwk · paralog_divergence.tsv · conversion_tracts.tsv ·
conversion_vs_conservation.tsv · tract_scan.py · this summary.
