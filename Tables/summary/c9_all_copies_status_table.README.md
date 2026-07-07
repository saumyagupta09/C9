# C9 all-copies status table (`c9_all_copies_status_table.csv`)

One row per **C9 copy/locus** across **all 57 squamate species** — the complete inventory,
including the non-functional copies that are *excluded* from the functional MACSE
alignments. Complements `c9_orf_rnaseq_table.csv` (which lists only the 82 finalized
copies + RNA-seq evidence); this table adds the lost copy, pseudogenes, and assembly
fragments so every locus is accounted for.

## Columns
- **clade** — Gekkota, Scincoidea, Lacertoidea, Anguimorpha, Iguania, Boidae/Pythonidae, Colubridae, Viperidae, Elapidae
- **species**, **copy** — synteny-based copy name (C9A = DAB2-proximal, C9B = FYB1-proximal, C9C = colubrid third copy; `C9_alt` / `(fragment)` = partial)
- **n_exons**, **protein_aa**, **ORF_status** — as called from the finalized CDS (miniprot)
- **category** — one of:
  - `functional` — intact ORF (ATG + stop + 0 internal stops) **and** 11 exons → used in selection analyses
  - `degraded/pseudogene` — a real copy with a broken ORF (no stop, internal stops, or a multi-exon deletion / short remnant)
  - `fragment` — only a partial C9 segment recoverable (assembly gap / other-scaffold hit); not a complete copy and not necessarily a pseudogene
  - `lost` — no C9 locus detectable genome-wide
- **in_MACSE_alignment** — whether the copy is in the functional MACSE codon alignments (`Alignment/`)
- **notes** — copy-specific detail

## Summary (90 rows / 57 species)
| category | n |
|---|---|
| functional | 79 |
| degraded/pseudogene | 4 |
| fragment | 6 |
| lost | 1 |

*Update 2026-07-01: **Notechis_scutatus C9B** is now a **complete functional gene** (11 exons, 583 aa) —
exons 1 and 2 were recovered from Illumina WGS reads (SRA ERR2714264), resolving the assembly gap. It
was previously called a 5-exon remnant, then a near-intact copy; it is now a full functional copy and
enters the MACSE alignment. Notechis is a genuine two-copy species (C9A + C9B both functional).*

*Update 2026-07-02: **Naja_naja C9B** is now a **complete functional gene** (11 exons, 583 aa) — SRA-read
blast (Illumina WGS) shows the prior "25 internal stops" were a frameshift artifact from ~5 bp of
assembly deletion errors (exon 7 −4 bp, exon 10 −1 bp). Reclassified `degraded/pseudogene` → `functional`;
old erroneous copy kept as `Naja_naja_C9B_pseudogene.*`. Naja is a genuine two-copy species.*

## Notable non-functional calls (characterized in detail)
- **Asaccus_caudivolvulus** — C9 **lost** (no locus genome-wide; see `asaccus_C9_loss/`).
- **Dopasia_gracilis C9A** — heavily degraded **pseudogene** (frameshift, 39 internal stops).
- **Pantherophis_obsoletus C9B** — genuinely **truncated pseudogene**: intact ATG + 10/11 exons,
  but **exon 11 (C-terminal TSP1 + stop) is deleted** at the locus (confirmed by tblastn of an
  intact Pituophis C9B protein — coverage stops at aa ~533 — and zero blastn hit for an intact
  C9B exon 11, despite ~10 kb of downstream sequence present).
- **Chrysopelea_ornata C9B** — genuine C9B ortholog but only an **assembly fragment** (~120 aa,
  exons ~6–7) on a 1221-bp contig; no start/stop recoverable, no internal stops where it aligns —
  unassembled, **not** a pseudogene (fragmented genome, no RNA-seq to rescue).
- **Notechis_scutatus C9B** — **COMPLETE functional gene (finalized 2026-07-01).** 11 exons,
  1752 nt / 583 aa, ATG start, terminal stop, 0 internal stops. Exons 1 and 2 were **recovered from
  Illumina WGS reads (SRA ERR2714264)**, resolving the assembly gap (orphan scaffold NW_020719678
  holds exons 3–6; main scaffold NW_020716626.1:11,074,675–11,083,116 holds exons 7–11). Earlier
  called a 5-exon remnant (308 aa), then near-intact (558 aa). C9A (Copy-1) = the complete contiguous
  gene at NW_020716626.1:11,106,491–11,139,037. **Notechis is a genuine two-copy species (C9A + C9B).**
- **Naja_naja C9B** — **corrected 2026-07-02: functional gene, NOT a pseudogene.** SRA-read blast shows
  the 25 "internal stops" were a frameshift artifact from ~5 bp assembly deletion (exon 7 −4 bp, exon 10
  −1 bp); the read consensus is a complete intact ORF (11 exons, 583 aa). Old copy kept as `*_pseudogene.*`.
- `C9_alt (fragment)` rows (Bothrops, Naja, Notechis, Pantherophis×2) — partial C9 hits on
  separate scaffolds; assembly fragments, not complete copies.
