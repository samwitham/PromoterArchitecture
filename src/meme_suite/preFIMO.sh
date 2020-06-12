#!/bin/sh
#requires activating "conda activate PromoterArchitecturePipeline"
#$1 is promoter gff3 location. $2 is the genome.fasta file location. $3 is the foldername

#identify input file name without path or extension
xbase=${1##*/}
xpref=${xbase%.*}

#make directory if doesn't already exist
mkdir -p ../../data/output/$3
mkdir -p ../../data/output/$3/FIMO
#convert gff file to bed file
gff2bed < $1 > ../../data/output/$3/FIMO/${xpref}.bed
#remove 'gene:' from column 4
awk 'BEGIN {FS = OFS = "\t"} { gsub(/gene:/,"", $4); print}' ../../data/output/$3/FIMO/${xpref}.bed > ../../data/output/$3/FIMO/${xpref}.bed.tmp && mv ../../data/output/$3/FIMO/${xpref}.bed.tmp ../../data/output/$3/FIMO/${xpref}.bed
#create fasta file of promoters from genome fasta file and from the promoters_renamedChr.bed file
bedtools getfasta -fi $2 -bed ../../data/output/$3/FIMO/${xpref}.bed -fo ../../data/output/$3/FIMO/${xpref}.fasta -name