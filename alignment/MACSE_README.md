# C9 alignments — finalized set (rebuilt 2026-07-03)

Rebuilt from the **finalized** C9 CDS in `exonwise_sequence_final/` after the read-validation /
reorganization. **81 copies** total: **55 C9A + 22 C9B + 4 C9C**.

## Inclusion / exclusion
- **Included:** all functional copies, the read-corrected copies (*Naja naja* C9B, *Notechis scutatus*
  C9B — now intact 583-aa genes), and the kept **truncated-but-meaningful** copies that lost only exon 11
  and terminate at an in-frame read-through stop (*Anolis sagrei* C9A 537 aa, *Pantherophis obsoletus*
  C9B 533 aa, *Anolis carolinensis* C9A 538 aa).
- **Excluded (pseudogenes / lost / fragments):** *Dopasia gracilis* C9A, *Protobothrops mucrosquamatus*
  C9B, *Asaccus* (lost) — all now in `Pseudogenes/`; the 6 assembly fragments; the superseded
  `Notechis_scutatus_C9B_remnant`; and the `Naja_naja_C9B_pseudogene` record.

## Global alignments (MAFFT only) — 5, in `Alignment/`
`C9A_global_anolisC9A` (55), `C9A_global_anolisC9B` (55), `C9B_global_anolisC9A` (22),
`C9B_global_anolisC9B` (22), `C9C_global` (4). Each = `<name>.fa` + `<name>_mafft.aln`.
(Anole variants: the two *Anolis* copies represented as C9A or as C9B.)

## Per-clade alignments (MAFFT + MACSE) — 9, in `Alignment/<clade>/`
`all_seq_<clade>.fa` (input, every paralog in that clade) + `_mafft.aln` (nucleotide) +
`_macse.aln` (codon-aware NT) + `_macse_AA.fasta` (protein).

| clade | copies | species | MACSE aln len | notes |
|---|---|---|---|---|
| Gekkota | 6 | 6 | 1875 | exon-6 length variation |
| Scincoidea | 6 | 6 | 1827 | exon-6 length variation |
| Lacertoidea | 6 | 6 | 1779 | |
| Anguimorpha | 2 | **2** | 1941 | **<3 species — too small for a standalone clade selection test**; long anguimorph exon 6 (Anniella 420 bp / Elgaria 357 bp) |
| Iguania | 12 | 10 | 1818 | includes both anoles' C9A+C9B |
| Boidae_Pythonidae | 3 | 3 | 1752 | |
| Colubridae | 21 | 9 | 1755 | |
| Viperidae | 11 | 6 | 1752 | |
| Elapidae | 14 | 7 | 1761 | incl. corrected Naja/Notechis C9B |

## Codon QC (all 9 MACSE alignments)
Alignment length % 3 = 0 · **0 frameshift markers (`!`)** · **0 real internal stops** → valid codon
alignments, ready for codeml / HyPhy RELAX after stripping the terminal-stop column.

## Tools
MAFFT v7.525 (`--auto`), MACSE v2.07 (`alignSequences`). Server build:
`…/C9_copy_number_2026/runs/v2_ashu_C9A/Alignment_rebuild/`.
