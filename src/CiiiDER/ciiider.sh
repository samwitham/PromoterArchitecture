#!/bin/sh
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the background gene fasta file
#$3 is jaspar motif file
#$4 is the promoterpref folder name
#$5 is the CiiiDER output folder name

#make CiiiDER folder
mkdir -p ../../data/output/"$4"
mkdir -p ../../data/output/"$4"/CiiiDER
mkdir -p ../../data/output/"$4"/CiiiDER/"$5"

echo "[General parameters]
STARTPOINT = 1
ENDPOINT = 2
OUTPUTFOLDER = ../../data/output/$4/CiiiDER/$5
PROJECTOUTPUTFILE = project.cdr

[Scan parameters]
GENELISTFILENAME = $1
MATRIXFILE = $3
DEFICIT = 0.15
GENESCANRESULTS = binding_site_data.bsl

[Enrichment parameters]
BGGENELISTFILENAME = $2
ENRICHMENTCOVERAGEPVALUE = 0.05
ENRICHMENTSITEPVALUE = 1.0
ENRICHMENTOUTPUTFILE = enrichmentoutput.txt" > ../../data/output/"$4"/CiiiDER/"$5"/config.ini




java -jar ../../software/CiiiDER/CiiiDER.jar -n ../../data/output/"$4"/CiiiDER/"$5"/config.ini