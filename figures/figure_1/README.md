# figure_1/ — C9 status & selection across squamate representatives (self-contained)

Figure 1: a 23-species representative tree with, per species, the **DAB2–C9A–C9B–C9C–FYB1 synteny**
(copy number + ORF status) and **selection** results.

## Make the figure
```bash
python3 figure_1/build_status_data.py       # -> c9_status_data.csv
Rscript figure_1/plot_figure1_scaled.R       # -> figure1_c9_status_scaled.png / .pdf   (main / final)
Rscript figure_1/plot_figure1.R              # -> figure1_c9_status.png / .pdf           (earlier purple variant)
```
Requires R + the **ape** package (base graphics otherwise).

## Encoding (final variant, `plot_figure1_scaled.R`)
- **Tree branches:** thin gray10; **seagreen4 (thick) = aBSREL-positive** lineage (Naja, Hemorrhois). aBSREL only.
- **Synteny boxes:** DAB2 & FYB1 fixed flanks (gray70); C9 copies **evenly distributed** between them
  (single-copy C9 centred). Fill = identity — C9A cadetblue4, C9B mediumseagreen, C9C steelblue3.
- **ORF status:** functional = solid; **red (firebrick2)** = pseudogene w/ multiple stops (Dopasia C9A);
  **faded + dotted** = remnant (Protobothrops C9B); **"locus deleted"** (red) = C9 gone (Asaccus).
- **Box outline:** **orange (darkorange2) = C9B relaxation** (RELAX K<1); others plain grey.
- **Events:** ◆ = C9B / C9C duplication; ✕ = loss / pseudogenisation.

## Files
| file | role |
|---|---|
| `rep_species.nwk` | INPUT — 23-species representative time tree |
| `c9_status_data.csv` | per-species copy status (c9a/b/c) + relax + absrel (built) |
| `build_status_data.py` | builds the table from the master copy table + RELAX/aBSREL results |
| `plot_figure1_scaled.R` | **final** figure script (fixed flanks, even distribution) |
| `plot_figure1.R` | earlier fixed-column purple variant |
| `figure1_c9_status_scaled.png/.pdf` | **final** figure |
| `figure1_c9_status.png/.pdf` | earlier variant |

## Data provenance
Copy/ORF status from `../Deliverables/tables/c9_all_copies_status_table.csv`; RELAX/aBSREL from
`../Selection_analysis/`. **Asaccus C9 loss** confirmed by whole-genome alignment to Correlophus
(`../Deliverables/asaccus_C9_loss/genomic_alignment/`): DAB2/FYB1 flanks align, C9 exons absent → clean deletion.
