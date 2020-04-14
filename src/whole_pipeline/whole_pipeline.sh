#!/bin/sh
source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate # need to change this to hardlink
conda activate PromoterArchitecturePipeline
##extract promoters from a genome file
python ../data_sorting/extract_promoter.py

##extract promoters

## optional - merge meme files together into one file ready for FIMO
#$1 directory containing lots of sub directories containing .txt files containing motifs
#$2 output file name
#../meme_suite/./preFIMO.sh ../../data/FIMO/motif_data/dap_data_v4/motifs dap_combined.meme


##run preFIMO.sh script. $1 is promoter gff3 location. $2 is the genome.fasta file location.
promoter_gff_file=../../data/genomes/promoters_renamedChr.gff3
genome_fasta=../../data/genomes/TAIR10_chr_all.fas

../meme_suite/./preFIMO.sh $promoter_gff_file $genome_fasta

## create FIMO background file:
#activate correct conda env
conda activate MemeSuite2
#identify the output filename created by preFIMO.sh
promoterbase=${promoter_gff_file##*/}
promoterpref=${promoterbase%.*}
#run FIMO.sh
#$1 is promoter fasta file. 
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file 

../meme_suite/./FIMO.sh ../../data/FIMO/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/dap_combined.meme

## filter the FIMO output
conda activate PromoterArchitecturePipeline

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bedfile file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
python ../data_sorting/./FIMO_filter.py ../../data/FIMO/output/${promoterpref}_FIMO/fimo.tsv ../../data/FIMO/${promoterpref}.bed ../../data/FIMO/${promoterpref}_motifs.bed 0.05

## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
../data_sorting/./TFBS_coverage.sh ../../data/FIMO/${promoterpref}.bed


