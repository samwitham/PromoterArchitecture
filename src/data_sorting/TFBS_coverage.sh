#!/bin/bash -e
#% of nucleotides falling within motifs - number of base pairs covered by at least one motif in a given sequence
#$1 is promoter bed file
xbase=${1##*/}
xpref=${xbase%.*}
coverageBed -a $1 -b ../../data/FIMO/${xpref}_motifs.bed > ../../data/promoter_analysis/${xpref}.bp_covered.txt

