#!/bin/sh
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the promoter fasta file location of constitutive, variable and random genes only
#$3 is jaspar motif file
#$4 is the promoterpref folder name
#$5 is the CiiiDER output folder name

#make CiiiDER folder
mkdir -p ../../data/output/$4
mkdir -p ../../data/output/$4/CiiiDER
mkdir -p ../../data/output/$4/CiiiDER/$5





java -jar ../../software/CiiiDER/CiiiDER.jar -n ../../software/CiiiDER/config.ini