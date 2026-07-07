# sra_blast/

BLAST-in-SAM results of **short-read SRA (Illumina WGS) reads vs C9 exon queries**, used to
validate/correct assembly-derived C9 sequences from raw reads (resolving assembly gaps and
frameshift/indel errors that masqueraded as pseudogenisation).

| file | reads | query | used for |
|---|---|---|---|
| `Naja_naja_C9B_sra.sam` | Naja Illumina WGS | Naja C9B exons | **Corrected Naja C9B**: read consensus is an intact 583-aa ORF → C9B is functional, not a 25-stop pseudogene (the stops were a ~5 bp assembly-deletion frameshift, exon 7 −4 bp / exon 10 −1 bp). |
| `Naja_naja_C9A_sra.sam` | Naja Illumina WGS | Naja C9A exons | C9A read validation / phasing companion to the C9B correction (Naja C9A + C9B locus reads). |
| `Protobothrops_C9A_sra.sam` | Protobothrops Illumina (DRR049668, DRR049669) | Protobothrops C9A exons | Protobothrops C9 read validation (C9A + C9B locus reads; C9B is flagged a 4-exon remnant, pending consensus check). |
| `Pituophis_C9B_on_Pantherophis_sra.sam` + `Pantherophis_exon11_read_consensus.fa` | Pantherophis Illumina (SRR10405238) | Pituophis C9B exons | Testing Pantherophis C9B exon 11: per-read classification: 15 reads=C9A, 2 reads=C9C, **0 distinct C9B** (query is 97% to C9A / 81% to C9C -> C9A capture bias). Supports genuine C9B exon-11 deletion. |
| `Protobothrops_C9A_lowstringency_sra.sam` | Protobothrops Illumina (DRR049668/9) | Protobothrops C9A exons (low-stringency) | Detects the divergent C9B 2nd copy at exon 10 (9 diagnostic sites) but NOT at exons 1-7 -> supports genuine 4-exon C9B remnant (exons 8-11), not an artifact. |
| `Crotalus_tigris_C9B_on_Protobothrops_sra.sam` | Protobothrops Illumina (DRR049668/9) | Crotalus tigris C9B exons (viper C9B bait) | **Confirms Protobothrops C9B is a 5'-truncated remnant**: candidate-C9B reads abundant at exons 8-11, ZERO at exons 1-7 (all C9A) -> exons 1-7 genuinely absent, not an assembly artifact. |
| `Dopasia_on_Dopasia_sra.sam` | Dopasia gracilis Illumina WGS (SRR1914353) | Dopasia C9A exons (12-exon, 39-stop annotation) | Re-checking the C9A pseudogene call with Dopasia's OWN reads (previously only Ophisaurus was used); testing the exon-4 boundary/frameshift. |
| `Dopasia_extended_SRA.sam` | Dopasia gracilis WGS (SRR1914353) | Dopasia C9A extended exons (+/-50bp, incl. introns) | Splice-site check: all acceptor-AG/donor-GT are read-confirmed at every exon (incl. 3,4,6); no frame-fixing indels -> canonical splice sites are REAL, frameshifts genuine. |
| `notechis_scutatus_sra.sam`, `notechis_scutatus_sra_again.sam` | Notechis Illumina WGS (ERR2714264) | Notechis C9A/C9B exons | **Recovered Notechis C9B exon 1** (and exon 2) from reads → completed C9B (assembly gap). |
| `notechis_exon1_extended_sra.sam` | Notechis Illumina WGS | Notechis C9B exon-1 extended region | Confirming the recovered C9B exon 1. |

Format: blastn `-outfmt SAM`. @SQ headers = the per-exon query references; alignment records = reads.
See `CORRECTION_LOG.md` (Notechis 2026-06-30/07-01, Naja 2026-07-02) and each copy's `*.PROVENANCE.txt`.
