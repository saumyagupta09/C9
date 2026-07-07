# C9 ORF + RNA-seq status table (`c9_orf_rnaseq_table.csv`)

One row per **finalized C9 copy** (82 copies across 55 species + degenerate
copies). Built from the genome-validated, synteny-named copies in
`exonwise_sequence_final/<species>/<copy>.cds.fa` (+ `.exons.fa`).

## Columns
| column | meaning |
|---|---|
| `clade` | one of the 9 broad clades used for the alignments (Gekkota, Scincoidea, Lacertoidea, Anguimorpha, Iguania, Boidae/Pythonidae, Colubridae, Viperidae, Elapidae) |
| `species` | species |
| `copy` | C9A (DAB2-proximal) / C9B / C9C (synteny naming, Sharma et al. 2022) |
| `copy_id` | exact sequence/file name (incl. `_pseudogene` / `_remnant` suffixes) |
| `cds_len_nt`, `protein_aa`, `n_exons` | CDS length, translated length, coding-exon count |
| `ORF_status` | computed from the CDS: `intact` = ATG start + terminal stop + 0 internal stops; otherwise `DEGRADED(...)` with the reason (noATG / no-stop / N internal stops) |
| `has_RNAseq` | yes/no — does the species have a C9-region RNA-seq BAM on the server |
| `RNAseq_tissue` | tissue(s) of the RNA-seq runs (liver/lung) |
| `RNAseq_status_per_copy` | per-copy expression: ORF call from reads + unique read count + locus coverage fraction |
| `RNAseq_seq_exact_match` | sequence-level confirmation for the 16 representative species (1–2/clade); see below |

## ORF notes (the only non-`intact` copies, all expected/documented)
- `Anolis_sagrei_C9A` (10 exons, no-stop) — C-terminally truncated copy (parallel to *A. carolinensis* C9A).
- `Pantherophis_obsoletus_C9B` (10 exons, no-stop) — degrading triplicate copy.
- `Protobothrops_mucrosquamatus_C9B_remnant` (4 exons, noATG) — 5′-truncated remnant (exons 8–11).
- `Naja_naja_C9B` — **corrected 2026-07-02: functional gene** (11 exons, 1752 nt / 583 aa, 0 internal
  stops). SRA-read blast (Illumina WGS) showed the prior "25 internal stops" were a frameshift artifact
  from ~5 bp of assembly deletion (exon 7 −4 bp, exon 10 −1 bp). Old copy kept as
  `Naja_naja_C9B_pseudogene.*`; see `CORRECTION_LOG.md` (2026-07-02).
- `Notechis_scutatus_C9B` — **re-finalized 2026-06-30: a near-intact functional duplicate, NOT a
  truncated remnant.** Exons 2–11 (558 aa, 0 internal stops, terminal stop); only **exon 1 (signal
  peptide) missing to an assembly gap** (exons 2–6 sit on orphan scaffold `NW_020719678`, exons
  7–11 on the main scaffold). The earlier "5 coding exons / 308 aa" call only saw the main block.
  See `CORRECTION_LOG.md` (2026-06-30) and `Notechis_scutatus_C9B.PROVENANCE.txt`. Pending long-read
  (SRR30595366) / Illumina (ERR2714264) validation of exon 1.

## RNA-seq sequence exact-match (`RNAseq_seq_exact_match`) — method + the chimera caveat
For 1–2 representative species per clade we rebuilt a **spliced RNA-seq consensus** per copy
(per-exon `samtools consensus` over the genome exon intervals, concatenated, strand-corrected)
and compared it to the finalized CDS by blastn.

**Important (tandem-paralog read cross-mapping):** for near-identical tandem copies, a *minor*
(low-expression) copy's RNA-seq consensus is frequently a **read-swap chimera** — reads from the
dominant paralog multi-map into the minor locus and dominate its consensus. To avoid being
fooled, each consensus was **tiled (100 bp)** and every tile assigned to its best-matching
finalized copy. A copy is only called confirmed when **every tile maps to its own CDS**
(`clean self-only`) at full length.

Verdicts:
- `CONFIRMED_exact` — clean self-only, ≥98% identity over ≥95% of the CDS. The dominant/single
  copies. (Residual 0–1.4% = heterozygous SNPs in the wild RNA-seq individual vs the reference
  assembly, not sequence errors.)
- `CONFIRMED_partial` — clean self-only but partial coverage (low-expression but not chimeric).
- `CHIMERIC_consensus:NOT_rnaseq-confirmable` — the RNA-seq consensus is a cross-mapping chimera
  (e.g. *Thamnophis* C9B = 100% dominant-C9A reads; *Notechis* C9B = C9A + C9B-remnant). These
  minor copies **cannot** be confirmed from RNA-seq and instead rest on genome + miniprot + NCBI
  validation (documented in `CORRECTION_LOG.md` / `ashu_sequence_audit.md`).
- `—` = species has RNA-seq but was not part of the representative seq-match subset.
- `NA` = no RNA-seq.

The chimerism is only in the RNA-seq *confirmation probe* — the finalized CDS were built
locus-anchored from the genome (`miniprot_corrected`), never from this consensus.

## Source files (server: `…/C9_copy_number_2026/runs/v2_ashu_C9A/`)
- consensus FASTA: `rnaseq_exactmatch/rnaseq_consensus.fa`
- per-copy match + chimera verdict: `rnaseq_exactmatch/rnaseq_seq_match_final.tsv`
  (mirrored locally at `results/v2_ashu_C9A/rnaseq_seq_match_final.tsv`)
- expression / splice / coverage: `results/v2_ashu_C9A/rnaseq_*.tsv`