#!/bin/sh
#requires activating "conda activate MemeSuite"
#$1 is the promoter fasta file location
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file 
#$5 is the folder name

#retrieve file name from file path
xbase=${1##*/}

xpref=${xbase%.*}
#make directory if doesn't already exist
mkdir -p ../../data/FIMO/$5
#create background file for FIMO
fasta-get-markov -dna $1 ../../data/FIMO/$5/${xpref}.bfile
#run FIMO
fimo --bfile ../../data/FIMO/$5/${xpref}.bfile --o ../../data/FIMO/$5/output/${xpref}_FIMO --thresh $2 --max-stored-scores $3 $4 $1
