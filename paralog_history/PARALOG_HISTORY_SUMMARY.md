# C9 paralog history — ancestral copy & duplication order (2026-07-03)

Three independent tests to establish which C9 copy is ancestral and the duplication/triplication order.

## Test 1 — synteny position (single-copy C9 vs DAB2 / FYB1)
Source: `results/v2_ashu_C9A/synteny_order.tsv`. For every single-copy species (lizards + boids), the
single C9 was scored as DAB2-proximal or FYB1-proximal (distance of copy midpoint to each flank).
**Result: 23/23 single-copy species are DAB2-proximal (0 FYB1-proximal)** = the C9A position.
Snake reference confirms order **DAB2 — C9A — C9B — C9C — FYB1** (C9C most FYB1-proximal).

## Test 2 — ML-tree clustering (patristic distance)
Source: `results/v2_ashu_C9A/allC9_ABC_ml.treefile` (IQ-TREE, all A/B/C copies).
- Single-copy lizard/boid C9: mean patristic 0.729 to snake **C9A** vs 0.836 to C9B → **clusters with C9A**.
- All 4 C9C: mean 0.224 to **C9B** vs 0.345 to C9A → **C9C clusters with C9B**.
- Paralog clades distinct: between A–B 0.360 > within-A 0.308 / within-B 0.309.

## Test 3 — pairwise %identity (all-copy MAFFT alignment, 81 seqs)
- Single-copy C9 mean identity: **C9A 69.0%** > C9B 68.3% > C9C 67.7% (29/32 best-match C9A).
- C9C: **C9B 86–87%** ≫ C9A 81–82% (all 4 colubrid C9C).

## Conclusion
- **C9A = ancestral copy** (ancestral DAB2-proximal position; orthologous lineage of the single ancestral
  gene; highest identity to single-copy C9). Single-copy lizards/boids retain C9A.
- **C9B = derived tandem duplicate of C9A** (positioned toward FYB1).
- **C9C = triplication derived from C9B** (clusters with C9B, ~87% identity, most FYB1-proximal).
- **Order of origin = chromosomal order: DAB2 — C9A(ancestral) — C9B(from C9A) — C9C(from C9B) — FYB1.**

Caveat: single-copy vs C9A/C9B identity gap is small (~0.7%) due to deep divergence (~65-69%) and likely
C9A<->C9B gene conversion; but all three tests agree in direction, and the C9C<-C9B signal is strong.

Inputs/reproducibility: synteny_order.tsv, allC9_ABC_ml.treefile, allC9_81_mafft.aln (on server
`…/Alignment_rebuild/`). Driver: `scripts/paralog_history_tests.py`.
