#!/bin/sh
#requires activating "conda activate MemeSuite"
#$1 is the promoter bed file location of constitutive genes only
#$2 is the promoter bed file location of constitutive and variable
#$3 is the folder name
#$4 is the genome fasta file
#$5 is the constitutive promoter fasta file output location
#$6 is the constitutive, variable and random promoters fasta file output location

#make directories for centrimo if they don't already exist
mkdir -p ../../data/output/$3
mkdir -p ../../data/output/$3/centrimo
mkdir -p ../../data/output/$3/centrimo/output
#make fasta files for just constitutive promoters, and also for constitutive variable and control promoters

bedtools getfasta -fi $4 -bed $1 -fo $5 -name
bedtools getfasta -fi $4 -bed $2 -fo $6 -name