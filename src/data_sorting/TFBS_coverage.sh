#!/bin/bash -e
#% of nucleotides falling within motifs - number of base pairs covered by at least one motif in a given sequence
#$1 is promoter bed file
#$2 is the folder name

#make directory if doesn't already exist
mkdir -p ../../data/output/$2/TFBS_coverage/
#retrieve filename from file path
xbase=${1##*/}
xpref=${xbase%.*}
coverageBed -a $1 -b ../../data/output/$2/FIMO/${xpref}_motifs.bed > ../../data/output/$2/TFBS_coverage/${xpref}.bp_covered.txt