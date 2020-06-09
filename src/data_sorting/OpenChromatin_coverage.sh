#!/bin/bash -e
#% of nucleotides falling within open chromatin
#$1 is promoter bed file
#$2 is the folder name

#retrieve filename from file path
xbase=${1##*/}
xpref=${xbase%.*}
coverageBed -a $1 -b ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed > ../../data/promoter_analysis/$2/${xpref}RootOpenChrom.bp_covered.txt
coverageBed -a $1 -b ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed > ../../data/promoter_analysis/$2/${xpref}ShootOpenChrom.bp_covered.txt
coverageBed -a $1 -b ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed > ../../data/promoter_analysis/$2/${xpref}ShootRootIntersectOpenChrom.bp_covered.txt