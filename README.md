# C9   — squamate complement C9 copy-number, phylogeny & selection study

Consolidated, reproducible result set. All files are **copies** of originals under
`results/v2_ashu_C9A/`, `exonwise_sequence_final/`, `Alignment/`, `Selection_analysis/`, `gene_tree/`,
`Pseudogenes/`, and the server run dir (originals untouched). Naming: **C9A = DAB2-proximal, C9B =
FYB1-proximal, C9C = colubrid third copy** (Sharma et al. 2022).

## Folder index (18)
| folder | contents |
|---|---|
| `tables/` | master CSVs: `c9_all_copies_status_table.csv` (every copy, ORF status, category), `c9_orf_rnaseq_table.csv` (per-copy RNA-seq expression), `c9_nonfunctional_copies.csv`, `c9_alt_fragments_assembly_context.csv`, `ashu_species_assembly_chimera_table.csv` (+ READMEs). |
| `sequences/<species>/` | **finalized** C9 CDS / peptide / exon FASTAs per copy (81 functional copies = 55 C9A + 22 C9B + 4 C9C; all clean ORFs) + PROVENANCE notes. |
| `Pseudogenes/<species>_<copy>/` | the genuinely non-functional loci excluded from analyses: *Dopasia* C9A (pseudogene), *Protobothrops* C9B (5'-remnant), *Asaccus* (C9 lost) — sequences (where any) + DETAILS.md. |
| `paralog_history/` | ancestral-copy & duplication-order reconstruction (3 tests): `PARALOG_HISTORY_SUMMARY.md`, `ANOLIS_PARALOG_HISTORY.md`, `test1_synteny_distances.csv`, all-copy alignment + tree. **C9A ancestral → C9B → C9C**; anole duplication independent. |
| `alignment/` | MAFFT + MACSE codon alignments: 5 globals (MAFFT; the `*_macse.aln` globals are older pre-correction builds kept only for reference) + 10 per-clade `all_seq_<clade>_{mafft,macse}.aln` (incl. Iguania_Anguimorpha) + `MACSE_README.md`. |
| `gene_tree/<clade>/` | IQ-TREE ML gene trees (`.treefile`/`.iqtree`) per clade + combined Iguania_Anguimorpha, built from the MACSE codon alignments — the topologies used for RELAX/codeml. |
| `Selection_analysis/` | **RELAX + codeml + aBSREL + MEME + FEL** results: summary CSVs (`RELAX_results.csv`, `codeml_results.csv`, `aBSREL_results.csv`, `MEME_sites.csv`, `FEL_sites.csv`, `sites_by_domain.csv`), narrative `RESULTS.md`, raw JSONs + human-readable tables (`.sites.tsv`/`.branches.tsv`/`.summary.txt`), codeml runs, labeled trees, and `domains/` (InterProScan). See `RESULTS.md`. |
| `synteny/<species>/` | DAB2/FYB1 synteny tblastn per species + `synteny_order.tsv` (gene order & copy naming). |
| `coordinates/` | `<species>.C9_exons.bed` — genomic exon coordinates of each C9 copy per species. |
| `phylogeny/` | C9 species/gene tree (`C9_tree.treefile`/`.iqtree`) and taxon list. |
| `gene-conversion/` | gene-conversion tract analysis (REJECTED; uses the pre-correction set — superseded by `paralog_history/`). |
| `rna-seq/<species>/` | IGV minilocus BAM bundles over the C9 locus (32 species) used for the expression analysis in `tables/c9_orf_rnaseq_table.csv`. |
| `sra_blast/` | SRA read-validation SAMs (Naja/Notechis/Pantherophis/Protobothrops/Dopasia) + consensus + README — the read-level checks that corrected/confirmed degraded copies. |
| `blast/` | blastn detection outputs (curated C9A/C9B CDS vs each genome; self-blasts; cross-genome checks). |
| `tblastn/<species>/` | tblastn results per species (whole-genome C9-protein scans, self-tblastn, exon rescue) + README. |
| `miniprot/<species>/` | miniprot protein-to-genome C9 models per species (locus-anchored gene models used for finalization). |
| `asaccus_C9_loss/` | evidence for the *Asaccus* C9 loss (genome-wide tblastn negative + SRA pooled-read check). |
| `repeatmasker/` | RepeatMasker screen of the finalized CDS (confirms repeat-free) + README. |

## Key result documents
- **Selection**: `Selection_analysis/RESULTS.md` (+ CSVs). **Paralog history**: `paralog_history/PARALOG_HISTORY_SUMMARY.md`.
- **Copy status / expression**: `tables/`. **Thesis-style write-up**: `../theisis/C9_Results_thesis_draft.md`.
- **Methods & full log**: `../METHODS_AND_RESULTS.md`, `../CORRECTION_LOG.md`. **Scripts**: `../scripts/`.
