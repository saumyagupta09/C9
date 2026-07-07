# RepeatMasker screen of the finalized C9 CDS (2026-06-26)

**Input:** the 82 genome-validated finalized CDS (`all_finalized_cds.fa`, pooled from
`exonwise_sequence_final/`). Total 142,429 bp, GC 43.5 %.

**Tool:** RepeatMasker 4.2.4 (rmblastn 2.17.1+, TRF, HMMER) via a self-contained
Apptainer image installed **inside the project folder**:
`C9_copy_number_2026/tools/repeatmasker.sif` (from `docker://dfam/tetools:latest`,
bundles Dfam 4.0 / 30,659 curated families). Runs were fully contained — image, cache
(`repeatmasker_final/rm_home/.RepeatMaskerCache`) and temp dirs all kept in-folder.

Two runs:
- **`*` (default, `-noint`)** — low-complexity + simple/tandem repeats only (library-
  independent). Library taxonomy irrelevant.
- **`squamata_*` (`-species squamata`)** — adds interspersed/TE screening. NOTE: curated
  Dfam 4.0 has **0 squamate-specific families** (only 218 ancestral amniote-level), so
  this detects conserved TEs only, not squamate-specific ones.

## Result — the CDS are essentially repeat-free
| metric | `-noint` | `-species squamata` |
|---|---|---|
| bases masked | 326 bp (**0.23 %**) | 543 bp (**0.38 %**) |
| SINE / LTR / DNA transposon / unclassified | 0 | 0 |
| LINE (CR1) | 0 | 4 hits, 217 bp (0.15 %) |
| simple repeats | 2 (74 bp) | 2 |
| low-complexity | 5 (252 bp) | 7 |

### Interpretation
- **No SINEs, LTRs, DNA transposons, or unclassified repeats** — the finalized CDS
  carry no transposable-element content / assembly contamination. They are bona fide
  coding sequence.
- **The 4 "X7B_LINE / CR1" hits are NOT real TE insertions.** They are short (31–80 bp),
  moderately diverged (10–24 %) segments at the **same conserved position (~nt 965–1020,
  exon 7/8)** in 4 unrelated snakes (*Chrysopelea* C9A, *Ptyas* C9A, *Crotalus_adamanteus*
  C9B, *Malpolon* C9B), all in copies with **clean intact ORFs**. A genuine TE insertion
  would be copy-specific and ORF-disrupting; an orthologous, conserved, ORF-preserving
  segment shared across species is instead coincidental similarity to the CR1 consensus
  (CR1 is the dominant reptile LINE) and/or an ancient exapted fragment — not a quality
  problem.
- **Trace low-complexity:** a recurrent A-rich stretch at ~nt 1144 in several colubrid
  minor copies (*Arizona* C9C, *Pituophis* C9C, *Pantherophis* C9C, *Masticophis* C9B,
  *Crotalus_adamanteus* C9B), plus 2 short simple repeats (*Bothrops* C9A (TGAC)n,
  *Crotalus_tigris* C9B (GAT)n).
- **The exon-6 Gly/Asp insertions** (*Acritoscincus, Anniella, Elgaria*) were **NOT
  flagged** as low-complexity/simple repeats — they fall below RepeatMasker's threshold,
  i.e. they are genuine (compositionally biased) coding sequence, not maskable repeats.

**Conclusion:** the finalized C9 CDS pass a clean repeat/TE QC — no transposon content,
only trace low-complexity, and the divergent exon-6 insertions are confirmed as real
coding sequence rather than repeat artefacts.

## Files
`all_finalized_cds.fa.{tbl,out,out.gff}` (default/-noint) ·
`squamata_all_finalized_cds.fa.{tbl,out,out.gff}` (-species squamata) · `rm_run.log`.
Server originals: `…/C9_copy_number_2026/repeatmasker_final/`.

## Re-running later (image persists in-folder, no re-download)
```
PROJECT=…/C9_copy_number_2026; SIF=$PROJECT/tools/repeatmasker.sif; OUT=$PROJECT/repeatmasker_final
apptainer exec --bind $PROJECT --home $OUT/rm_home --pwd $OUT $SIF \
  /opt/RepeatMasker/RepeatMasker -noint -xsmall -gff -pa 8 -dir $OUT <input.fa>
```