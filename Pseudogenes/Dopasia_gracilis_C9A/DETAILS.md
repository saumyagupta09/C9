# Dopasia_gracilis C9A — genuine PSEUDOGENE (removed from analysis)

Clade: Anguimorpha (Anguidae, Anguinae). Genome: GCA_052054735.1, CM125122.1 (minus strand).

- Annotated as 12 exons (a v3 extraction split of the long anguimorph exon 6), total CDS 1750 nt, **%3=1
  (frameshift), 39 internal stop codons**.
- **All 12 exon boundaries are canonical GT...AG splice sites** (checked vs genome; both tblastn & miniprot
  agree) -> the exon structure is correct; the frameshifts are genuine internal coding disruptions, NOT
  mis-annotation.
- **Own-species reads (Illumina WGS SRR1914353)** were blasted vs the exons (Dopasia_on_Dopasia_sra.sam):
  the consensus reproduces the assembly bases -> the disruptions are NOT sequencing/assembly errors.
- exon-4 +1bp shift rescues exons 4-5 (39->26 stops) but frameshifts resume at exon 6 -> multiple
  independent disruptions.
- Conclusion: bona-fide pseudogene, pseudogenization likely ancestral to Anguinae (shared with Ophisaurus).
- Sequences here: Dopasia_gracilis_C9A.{cds,pep,exons}.fa (the 12-exon, 39-stop annotation).
- Supporting extractions/analyses: ../../Dopasia_extraction/ (splice sites, extended exons, read consensus).
