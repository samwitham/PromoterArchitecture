#!/bin/sh
#generate a background file for use with MEME_suite. Markov model of order 3 is used.
#Follow this page to add interactive mode: http://linuxcommand.org/lc3_wss0120.php
#look into document root as base location for code - you set this once?
#Run ./background.sh At_allpromoters.fasta
if [ "$1" -eq 0 ]
  then
    echo "No arguments supplied"
fi
fasta-get-markov -m 3 -dna "$1" /ei/workarea/group-eg/project_PromoterArchitecturePipeline/data_output/meme_suite/background.bfile