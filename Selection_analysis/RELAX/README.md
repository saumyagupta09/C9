# Selection_analysis / RELAX — datasets, foregrounds, files (2026-07-03)

MACSE codon alignments + IQ-TREE gene trees + foreground-labeled trees for HyPhy RELAX and PAML codeml.
Each dataset folder has: `*_input.fa` (CDS), `*_macse.aln` (codon NT) + `*_macse_AA.fasta`,
`*_codon.aln` (! -> - sanitized, tree input), `*.treefile` (+ .iqtree/.contree/.log), and one
`*_FG<scheme>_relax.nwk` (HyPhy `{Test}`) + `*_FG<scheme>_codeml.nwk` (PAML `#1`) per foreground scheme.
Tools: MACSE v2.07, IQ-TREE 1.6.12 (`-m MFP -bb 1000`).

## Universal exclusions (NOT in any RELAX/codeml analysis)
Dopasia_gracilis, Elgaria_multicarinata, Asaccus_caudivolvulus, Anniella_stebbinsi (all copies) +
Protobothrops_mucrosquamatus_C9B (that copy only; Protobothrops C9A is retained).

## 1. snakes_C9A/ (25 seqs) — single-copy boid/python C9A vs advanced-snake C9A
Boidae/Pythonidae(3) + Colubridae(9) + Elapidae(7) + Viperidae(6, incl Protobothrops C9A) C9A.
Foreground schemes (each: that clade = Test, rest = Reference):
`snakes_C9A_FGboidae`, `_FGcolubridae`, `_FGelapidae`, `_FGviperidae`.

## 2. duplicated_snakes_C9A_C9B/ (42 seqs) — does C9B differ from C9A after duplication
Colubridae+Elapidae+Viperidae C9A(22) + C9B(20). Boid/python C9A, C9C, Protobothrops C9B excluded.
Schemes: `dupsnakes_FGC9B` (all C9B = Test, all C9A = Ref) and `dupsnakes_FGC9A` (all C9A = Test, all C9B = Ref).

## 3. colubrid_C9A_C9B_C9C/ (12 seqs) — does C9C differ after triplication
Triplicated colubrids only: Arizona, Hemorrhois, Pituophis, Pantherophis (each C9A+C9B+C9C).
Schemes: `colubrid_ABC_FGC9C` (Test=C9C), `_FGC9B` (Test=C9B), `_FGC9A` (Test=C9A); rest = background.

## 4. iguania_anolis/ (12 seqs) — Anolis duplicated copies
Anolis carolinensis + sagrei C9A/C9B (4) + 8 other-iguanian C9A (background).
- `iguania_anolis_4A_FGanoleC9B`: Test = anole C9B, Reference = C9A copies.
- `iguania_anolis_4B_FGanoleC9A`: Test = anole C9A, Reference = other-iguanian C9A + anole C9B.

## Codon QC
All 4 MACSE alignments: length %3 = 0, 0 frameshift markers, 0 internal stops.
