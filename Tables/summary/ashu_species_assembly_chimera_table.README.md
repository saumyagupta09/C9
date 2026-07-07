# Ashu (reference-study) species — assembly quality & chimera assessment

`ashu_species_assembly_chimera_table.csv` — the **22 Ashu reference species**
(`role = ashu` [17] or `both` [5] in `targets.tsv`; = the species curated in the
Sharma et al. study and present in `data/ashu_C9A.fa` / `ashu_C9B.fa`). These were taken
**as-is** from the prior study; the assembly-quality selection criteria were **not** applied
to them (those apply only to the 35 newly-selected `miocene` species).

## Columns
- **ncbi_assembly_level** — NCBI's label (Contig / Scaffold / Chromosome). **Misleading on its
  own** — see below.
- **contig_N50 / scaffold_N50 / n_scaffolds** — assembly contiguity (bp).
- **effectively_contig_level** — YES if **scaffold N50 < ~25 kb** (a C9 gene ≈ 20–25 kb cannot
  sit on a single scaffold → the gene is physically shattered, regardless of NCBI's label).
- **copy_number_in_species** — number of real C9 gene copies WE identified (functional +
  degraded; assembly-fragment "C9_alt" hits not counted).
- **copies_in_github_ashu_alignment** — number of copies the **curated study** included for
  this species (`data/ashu_C9A.fa` + `ashu_C9B.fa`).
- **chimera_likelihood** — Yes / Possible / No (see rule below).

## Why NCBI level is misleading (the key point)
- **Chrysopelea_ornata** is labeled **"Scaffold"** but scaffold N50 = **1,564 bp** with
  **1.12 M scaffolds** — it is *effectively a contig assembly* and the only genuinely shattered
  one. Its C9A had to be stitched from ~10 contigs (→ chimera risk); its C9B is unrecoverable.
- **Crotalus_tigris** is labeled **"Contig"** but has contig N50 **2.1 Mb** — perfectly fine;
  the C9 gene sits on one piece.
- Many snake assemblies (Pantherophis, Ptyas, Pseudonaja, Laticauda, …) have small **contig**
  N50 (~30–50 kb) but large **scaffold** N50 (Mb) — the C9 gene sits intact on one scaffold
  (contigs joined by intronic N-gaps). Verified: their copies are single-scaffold.

## chimera_likelihood rule
- **Yes** — scaffold N50 < gene length **and** the copy is reconstructed across multiple
  contigs **and** a paralog exists to mis-stitch → **only Chrysopelea_ornata**.
- **Possible** — the copy had exons on >1 scaffold in the self-blast QC and the species has a
  paralog (needs per-exon verification) → **Notechis_scutatus C9A** (3/11 exons off the
  dominant scaffold despite a 6 Mb scaffold N50).
- **No** — single contiguous-scaffold gene, **or** a single-copy species (no paralog to
  chimerize with, e.g. *Python_bivittatus*, whose 1 off-scaffold exon cannot be a paralog mix).

## Note: where we found MORE copies than the study
Our `copy_number_in_species` > `copies_in_github_ashu_alignment` for Thamnophis, Vipera_berus,
Naja_naja, Protobothrops (we recovered a C9B the study lacked) and Pantherophis, Pituophis
(we recovered a third copy, C9C). These extra copies are on well-assembled genomes (chimera = No).

## Outstanding verification
Even though our Chrysopelea C9A = the curated reference 100%, the **curated sequence itself**
could be a mosaic (the study faced the same shattered assembly). A sliding-window per-paralog
scan over diagnostic sites is needed to confirm/refute chimerism of Chrysopelea C9A (and to
resolve Notechis C9A).