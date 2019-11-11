#!/bin/sh
#requires activating "conda activate PromoterArchitecturePipeline"
#$1 is promoter gff3 location. $2 is the genome.fasta file location.
#identify input file name without path or extension
xbase=${1##*/}
xpref=${xbase%.*}

#convert gff file to bed file
gff2bed < $1 > ../../data/FIMO/${xpref}.bed
#create fasta file of promoters from genome fasta file and from the promoters_renamedChr.bed file
bedtools getfasta -fi $2 -bed ../../data/FIMO/${xpref}.bed -fo ../../data/FIMO/${xpref}.fasta -name
