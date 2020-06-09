#!/bin/sh
#requires activating "conda activate PromoterArchitecturePipeline"
#$1 is promoter gff3 location. $2 is the genome.fasta file location. $3 is the foldername

#identify input file name without path or extension
xbase=${1##*/}
xpref=${xbase%.*}

#make directory if doesn't already exist
mkdir -p ../../data/FIMO/$3
#convert gff file to bed file
gff2bed < $1 > ../../data/FIMO/$3/${xpref}.bed
#create fasta file of promoters from genome fasta file and from the promoters_renamedChr.bed file
bedtools getfasta -fi $2 -bed ../../data/FIMO/$3/${xpref}.bed -fo ../../data/FIMO/$3/${xpref}.fasta -name
