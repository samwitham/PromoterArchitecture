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
mkdir -p "../../data/output/$5"
mkdir -p "../../data/output/$5/FIMO"
mkdir -p "../../data/output/$5/FIMO/output"
#create background file for FIMO
fasta-get-markov -dna "$1" "../../data/output/$5/FIMO/${xpref}.bfile"
#run FIMO
fimo --bfile "../../data/output/$5/FIMO/${xpref}.bfile" --o "../../data/output/$5/FIMO/output/${xpref}_plantpan" --thresh "$2" --max-stored-scores "$3" "$4" "$1"
