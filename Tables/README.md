# Tables/ — all C9 project tables (single consolidated location)

Every project table lives here, organised by category. Copies of tables that also belong to an analysis
folder (`../Selection_analysis/`, `../paralog_history/`, `../gene-conversion/`, `../figures/`) are kept there
too; this folder is the one-stop reference. Copy naming: **C9A = DAB2-proximal, C9B = FYB1-proximal,
C9C = colubrid third copy.**

## Subfolders
| folder | contents |
|---|---|
| **`summary/`** | copy number, ORF status, gene structure, assembly: `c9_all_copies_status_table` (master), `c9_nonfunctional_copies`, `c9_alt_fragments_assembly_context`, `exon_length_clade_table`, `master_copies`, `master_species`, `ashu_species_assembly_chimera_table`, `table_all_squamata_genomes` (NCBI genome metadata, 96 squamata), `loci_detail.json` (+ per-table READMEs). |
| **`synteny/`** | `synteny_order` (gene order DAB2–C9A–C9B–C9C–FYB1, flank positions) + `test1_synteny_distances` (copy→DAB2/FYB1 distances). |
| **`expression/`** | RNA-seq: `c9_orf_rnaseq_table` (per-copy reads + coverage) + `rnaseq_species_summary` (wide, figure input). |
| **`selection_analyses/`** | `Table_3_1_Selection_Analyses_Summary` (main-text), `RELAX_results`, `codeml_results`, `aBSREL_results`, `MEME_sites`, `FEL_sites`, `sites_by_domain`, `Pituophis_C9A_interproscan`. |
| **`gene_conversion/`** | `paralog_divergence`, `conversion_tracts`, `conversion_vs_conservation`. |
| **`QC/`** | pipeline QC / intermediate tables (not result tables): `contiguity_check_refined`, `rescue_report`, `termini_fix_report`, `exon_coords`, `miniprot_copies`, `clade_legend`. |
| **`appendix_main_tables/`** | tables destined for the thesis. `Table_3_1_Selection_Analyses_Summary` (main-text selection summary); `appendix_present_study_genomes` — the 40-species C9 dataset with genome provenance: **38 high-quality reference genomes analysed in the present study (35 unique + 3 shared with Sharma et al. 2022)**, plus **Vipera_berus** and **Chrysopelea_ornata**, whose assemblies did not meet the HQ criteria (contig N50 < 5 Mb) and were taken from **Sharma et al. 2022**. NCBI assembly metrics included; Chrysopelea_ornata metrics are from our assembly QC table. |

## Clade legend (group → clade)
1 Gekkota · 2 Scincoidea · 3 Lacertoidea · 4 Anguimorpha · 5 Iguania · 6 Boidae/Pythonidae · 7 Colubridae · 8 Viperidae · 9 Elapidae

*Replaces the former `tables/` and `all_tables/` folders (consolidated here on 2026-07-07).*
