#!/bin/bash
# InterProScan on the reference C9A (Pituophis_catenifer_C9A, 583 aa) via EBI REST API.
# Reproduces the domain intervals used to map MEME/FEL sites. Requires internet.
SEQ=$(python3 -c "print(''.join(l.strip() for l in open('exonwise_sequence_final/Pituophis_catenifer/Pituophis_catenifer_C9A.pep.fa') if not l.startswith('>')).rstrip('*'))")
JOB=$(curl -s -X POST --data "email=rohitgupta.dev@gmail.com" --data "stype=p" --data-urlencode "sequence=$SEQ" \
      https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run)
echo "job: $JOB"
until [ "$(curl -s https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/$JOB)" = "FINISHED" ]; do sleep 15; done
curl -s https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/$JOB/tsv > Selection_analysis/domains/Pituophis_C9A_interproscan.tsv
echo "wrote Selection_analysis/domains/Pituophis_C9A_interproscan.tsv"
# Curated domain intervals extracted (reference codon positions):
# signal 1-20 | TSP1_N 38-91 | LDLRA 96-133 | MACPF 135-495 | EGF 496-536 | TSP1_C 537-583
