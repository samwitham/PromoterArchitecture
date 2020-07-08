#!/bin/sh
#requires activating "conda activate MemeSuite"
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the promoter fasta file location of constitutive, variable and random genes only
#$3 is meme motif file 
#$4 is the folder name

#retrieve file name from file path
xbaseconstitutive=${1##*/}
xprefconstitutive=${xbase%.*}


#run centrimo
#--local = compute enrichment of all regions
centrimo --neg $2 --o ../../data/output/$4/centrimo/output/${xbaseconstitutive}_DAPseq --local $1 $3