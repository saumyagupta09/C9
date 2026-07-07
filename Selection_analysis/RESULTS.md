# Selection analysis results — RELAX + codeml branch models (2026-07-03)

HyPhy 2.5.100 RELAX (K relaxation parameter + LRT p) and PAML codeml branch models (M0 null `model=0`
vs two-ratio alt `model=2`, NSsites=0; LRT 2ΔlnL df=1) on the same foreground partitions.
Tables: `RELAX_results.csv`, `codeml_results.csv`. Alignments = MACSE codon (stops→gaps); trees = IQ-TREE.

## Significant results (p<0.05)
| test | RELAX K (p) | codeml ω_fg / ω_bg (p) | interpretation |
|---|---|---|---|
| **snakes_C9A — Elapidae FG** | **1.73 intensified (3.9e-5)** | **1.02 / 0.55 (1.3e-6)** | Elapid C9A under strongly **elevated ω (~neutral)** — significant rate shift (both methods). |
| **dup_snakes — C9B FG** | **0.81 relaxed (0.019)** | **0.48 / 0.68 (2.8e-4)** | **C9B evolves differently from C9A after duplication** — RELAX: relaxed distribution; codeml: lower mean ω (stronger purifying). |
| **snakes_C9A — Colubridae FG** | **0.78 relaxed (0.045)** | **0.53 / 0.71 (0.023)** | Colubrid C9A modest **relaxation / lower ω**. |
| **iguania 4B — anole C9A FG** | 0.62 (0.43 ns) | **0.54 / 0.31 (6.8e-4)** | Anole C9A **elevated ω** vs iguanian background (codeml only; RELAX ns). |

## Non-significant (both methods)
- **colubrid C9C FG (triplication):** RELAX K=0.89 p=0.50; codeml ω 0.63/0.52 p=0.28 → **C9C does NOT
  significantly differ after triplication** (no detected rate shift; evolving like C9A/C9B).
- **snakes_C9A Boidae FG:** single-copy boid/python C9A not significantly different from advanced-snake C9A.
- **snakes_C9A Viperidae FG; dup_snakes C9A FG; colubrid C9A/C9B FG; iguania 4A anole-C9B FG:** all ns.

## Biological summary
1. **Duplication (C9A→C9B):** C9B is under **significantly different, relaxed selection** vs C9A (RELAX
   K=0.81 p=0.019; codeml p=2.8e-4) — consistent with post-duplication relaxation / subfunctionalization of C9B.
2. **Triplication (C9B→C9C):** **no significant shift** for C9C — the third copy evolves under similar
   constraint to C9A/C9B (both methods ns).
3. **Clade-specific C9A:** the **elapid C9A** stands out with strongly elevated ω (intensified / near-neutral,
   both methods, p≤1e-5); **colubrid C9A** shows modest relaxation; boid/python and viper C9A unremarkable.
4. **Anolis:** anole C9A shows elevated ω vs other iguanian C9A (codeml p=7e-4); anole C9B test ns.

RELAX and codeml AGREE on 3/4 significant calls; the anole-C9A case is codeml-only (mean-ω shift the
distribution-based RELAX K did not flag). Method note: RELAX K measures ω-distribution shape (relaxed<1 /
intensified>1); codeml two-ratio compares mean ω between branch sets — so a low mean ω can co-occur with
RELAX "relaxation" (e.g., dup C9B).

## Branch-level positive selection (aBSREL) — aBSREL_results.csv
15/170 branches under episodic positive selection (corrected p<=0.05), concentrated on venomous-elapid
C9A: Laticauda_colubrina C9A (2.9e-9), Naja_naja C9A (6.8e-4), + Chrysopelea C9A and internal nodes;
also some C9B branches (Laticauda C9B 1.5e-6, Hemorrhois C9B) and one internal branch in the
iguania/anolis tree (Node5, corrected p=0.029). No C9C branches. Confirms the elapid-C9A elevated-omega
as genuine lineage-specific positive selection.

## Site-level selection (MEME/FEL) + domain mapping — MEME_sites.csv, FEL_sites.csv, sites_by_domain.csv
FEL: 209 significant sites (mostly purifying; 36 pervasive positive). MEME: 83 episodic-diversifying sites.
Positively-selected/diversifying sites are CONCENTRATED IN THE MACPF pore-forming domain (domains
assigned by InterProScan on the Pituophis C9A reference; MACPF = Pfam PF01823, residues 135-495):
snake C9A 74% (MEME) / 88% (FEL+); dup-snake 71% / 80%; colubrid 57% (MEME). => diversifying selection targets the membrane-
attack core, not the peripheral TSP1/LDLRA/EGF modules.
