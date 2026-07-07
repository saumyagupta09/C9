# figures/ — C9 RNA-seq expression figure (self-contained: script + inputs + outputs)

Everything needed to regenerate the C9 transcriptional-expression figure lives in this folder.

## Make the figure
```bash
Rscript figures/plot_rnaseq_expression.R      # run from the repo root
```
Requires R (≥4) with the **ape** package only (base graphics otherwise). Writes both
`c9_rnaseq_expression_figure.png` and `.pdf`.

### `plot_rnaseq_expression.R` — what it draws
Three aligned panels sharing the tree tip order:
1. **species tree**, branches coloured by copy number (single = light blue, 2 = green, 3 = pink);
2. **Read counts (uniquely mapped)** — per-copy log10-stacked bars (C9A steelblue3, C9B darkgoldenrod1,
   C9C mediumseagreen), copy key above; *Bothrops* shown as "no expression detected";
3. **Fractional locus** coverage — gradient cells (white→aquamarine→turquoise4→purple4, 0→1) with a
   colourbar, plus a **Tissue** column (brown3 = liver, plum = lung).

### Inputs (both in this folder — do not move)
| file | role |
|---|---|
| `rna_seq_sps.nwk` | Newick species tree of the 32 RNA-seq species (tip `Bassiana_duperreyi` is reconciled to `Acritoscincus_duperreyi` in-script) |
| `rnaseq_species_summary.csv` | per-species data: clade, tissue, n_copies, copies_present, `reads_C9A/B/C`, `cov_C9A/B/C` |

### Outputs
`c9_rnaseq_expression_figure.png`, `c9_rnaseq_expression_figure.pdf`.

## Where the numbers come from — calculation scripts
Standalone tools reproducing the per-copy quantification in `../scripts/16_rnaseq_pipeline.sh`.

| script | computes | method |
|---|---|---|
| `count_unique_reads_per_copy.sh <bam> <scaffold:start-end> [name]` | total + **uniquely-mapped** reads per copy | `samtools view -c -F 0x104`; unique adds `-q 255` (STAR unique-mapping MAPQ) — separates near-identical tandem paralogues |
| `calc_fractional_coverage.sh <bam> <exon_bed> [region]` | **fractional locus coverage** (0–1) | fraction of exonic bases with depth ≥ 1 = Σ(covered exon bases)/Σ(exon lengths) |

**Validation** (Naja naja): script gives C9A 128224 reads / cov 0.998, C9B 76 / 0.324 — matches the
master table (128227 / 1.00, 76 / 0.32). Official values were computed on the full-genome STAR BAMs.

## Other data files
| file | contents |
|---|---|
| `rnaseq_expression_tidy.csv` | long form (1 row/species×copy): species, clade, copy, reads, coverage, tissue, orf — 49 rows |

Source table: `../exonwise_sequence_final/c9_orf_rnaseq_table.csv`.
Copy naming: C9A = DAB2-proximal, C9B = FYB1-proximal, C9C = colubrid third copy.
