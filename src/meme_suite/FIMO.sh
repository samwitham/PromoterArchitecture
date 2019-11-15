#!/bin/sh
#requires activating "conda activate MemeSuite"
#$1 is the promoter fasta file location
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file 

#retrieve file name from file path
xbase=${1##*/}

xpref=${xbase%.*}
#create background file for FIMO
fasta-get-markov -m 3 -dna $1 ../../data/FIMO/${xpref}.bfile
#run FIMO
fimo --bfile ../../data/FIMO/${xpref}.bfile --o ../../data/FIMO/output/${xpref}_FIMO --thresh $2 --max-stored-scores $3 $4 $1
