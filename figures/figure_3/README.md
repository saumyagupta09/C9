# figure_3/ — C8 / C9 maximum-likelihood gene tree (Sharma-2022 Fig-5 style)

A nucleotide ML tree of **C8α, C8β and C9** to place the squamate C9 duplication in the context of its
sister gene family — reproducing (and extending, with C9C) Fig 5 of Sharma et al. 2022.

## Sequence set (114)
- **81 squamate C9** — **our curated intact copies** (55 C9A + 22 C9B + 4 C9C); pseudogenes/fragments excluded.
- **24 C8** (C8α + C8β) — reused **exactly from Sharma 2022** (`../data/…C8A_C8B.fa`).
- **9 non-squamate C9 outgroups** — reused from Sharma 2022: *Xenopus, Rana* (frogs), *Struthio, Dromaius,
  Taeniopygia* (birds), *Chrysemys, Trachemys* (turtles), *Alligator, Gavialis* (crocs).

## Pipeline (reproduce)
```bash
python3 figure_3/build_input.py     # -> c8_c9_cds.fna  (assemble the 114 CDS)
bash    figure_3/run_tree.sh        # (server) codon-aware align + IQ-TREE  [needs mafft + iqtree]
Rscript figure_3/plot_tree.R        # -> figure3_c8c9_tree.png / .pdf
```
Tree: **nucleotide, codon-aware alignment** (protein-guided: translate → MAFFT → back-translate) →
**IQ-TREE 1.6.12**, ModelFinder best model **GTR+F+R5**, **ultrafast bootstrap ×1000**, rooted on C8.

## Files
| file | role |
|---|---|
| `c8_c9_cds.fna` | input — 114 CDS | 
| `c8_c9_codon.aln` | codon-aware nucleotide alignment (tree input) |
| `c8_c9_tree.treefile` / `.contree` / `.iqtree` | ML tree / bootstrap consensus / model+log |
| `build_input.py`, `run_tree.sh`, `plot_tree.R` | scripts |
| `figure3_c8c9_tree.png` / `.pdf` | figure |

## Result
**C8 (α+β) is monophyletic** with strong support (clean outgroup), but the **C9 copies (C9A/C9B/C9C) are
not reciprocally monophyletic** — the same *weak internal C9 resolution* Sharma reported. Notably the copies
are **not homogenised** (a species' own C9A vs C9B are clearly distinct), so the weak resolution is not gene
conversion but genuinely poor deep signal among the young, sequence-similar copies.
