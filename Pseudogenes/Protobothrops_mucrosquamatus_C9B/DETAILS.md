# Protobothrops_mucrosquamatus C9B — genuine 5'-truncated REMNANT (removed from analysis)

Clade: Viperidae. Genome: GCA_018340635.1.

- C9B is only the C-terminal **exons 8-11** (a 5'-truncated remnant; **exons 1-7 deleted** at the locus).
- **SRA-read confirmed genuine (not an assembly artifact):** Crotalus tigris C9B (a close viper C9B) blasted
  vs Protobothrops Illumina reads (DRR049668/9) yields abundant candidate-C9B reads at exons 8-11 (e.g.
  57/57 at exon 8; 24 biallelic sites at exon 10) but **ZERO candidate-C9B reads at exons 1-7** (all C9A).
  A bait that detects C9B at 8-11 but not 1-7 => exons 1-7 are truly absent.
- Because it lost 7 exons (not <=1), it cannot encode a meaningful protein -> pseudogene/remnant.
- Sequences here: Protobothrops_mucrosquamatus_C9B_remnant.{cds,pep,exons}.fa (exons 8-11).
- SAM: Deliverables/sra_blast/Crotalus_tigris_C9B_on_Protobothrops_sra.sam.
