#!/bin/bash -e
#run gat enrichment
#the following input files are those created by src/data_sorting/gat_enrichment.py
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name as prefix)

#constitutive gene enrichment
gat-run.py --ignore-segment-tracks --segments=$4 `#eg. TATA box annotations` \
    --annotations=$2 `#constitutive promoter annotations` \
    --workspace=$1 `#constitutive and variable promoters` \
    --num-samples=1000 --log=$5/gat_$6_constitutive.log > $5/gat_$6_constitutive.out
    
#variable gene enrichment
gat-run.py --ignore-segment-tracks --segments=$4 `#eg. TATA box annotations` \
    --annotations=$3 `#variable promoter annotations` \
    --workspace=$1 `#constitutive and variable promoters` \
    --num-samples=1000 --log=$5/gat_$6_variable.log > $5/gat_$6_variable.out