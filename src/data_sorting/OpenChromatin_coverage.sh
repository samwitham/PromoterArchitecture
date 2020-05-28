#!/bin/bash -e
#% of nucleotides falling within open chromatin
#$1 is promoter bed file

#retrieve filename from file path
xbase=${1##*/}
xpref=${xbase%.*}
coverageBed -a $1 -b ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed > ../../data/promoter_analysis/${xpref}RootOpenChrom.bp_covered.txt
coverageBed -a $1 -b ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed > ../../data/promoter_analysis/${xpref}ShootOpenChrom.bp_covered.txt