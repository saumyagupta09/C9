#Script to make a multifasta

for i in `ls -1`
do
echo "$i"
dos2unix "$i"
cat "$i" >> multifasta.fa
done
dos2unix multifasta
cat multifasta.fa #check the gaps
awk 'NF' multifasta.fa | sed -i 's/>/\n>/g' multifasta.fa   #remove all the empty lines #add an extra line adter each fasta 

#make everything in one line
sed -e '/>/s/$/#/' multifasta.fa | tr -d '\n' | sed 's/#/\n/g' | sed 's/>/\n\n>/g' | sed '1,2d' | sed -e '$a\' > multifasta12.fa

#running guidance 
perl /home/ceglab09/programs/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile multifasta12.fa --msaProgram PRANK --seqType codon --outDir prank_100 --genCode 1 --bootstraps 100 --proc_num 8



#to make a multifasta with exons having name of species
for i in `ls *`
do 
echo $i
n=`ls $i | cut -c7- | rev | cut -c10- | rev`
sed 's/exon/'$n'_exon/g' $i >> multifasta.fa
done

#to extract a portion from a file
awk '/Query_7     61/,/Query_7     121/' pipra_motacilla_passermontanus_on_zebrafinch_eval0.05_srr.out

awk '/Query_7     61/,/Query_7     121/' pipra_motacilla_passermontanus_on_zebrafinch_eval0.05_srr.out | cut -c26- |awk '/^\./{a++}END{print a}' #to count the number of changes in sra blast output result.Calculates the first character of each line.
#or
awk '/Query_7     61/,/Query_7     121/' pipra_motacilla_passermontanus_on_zebrafinch_eval0.05_srr.out | cut -c26 | grep "G" | wc -l



# to extract a consensus sequence from a particular region in bamfile using mpileup

bcftools mpileup -f <genome> <bamfile> -r chromosome:regionfrom-regionto | bcftools call -c | vcfutils.pl vcf2fq > output_filename.fa

seqtk seq -q20 -n N rawread_from_mpileup_exon1_zwint.fa | tr -d 'N/!'| sed -n '2p' > output_filename.fa

