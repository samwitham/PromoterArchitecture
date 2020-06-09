#!/bin/sh
#source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate # need to change this to hardlink
eval "$(conda shell.bash hook)"
conda activate PromoterArchitecturePipeline
#directory_path=ei/workarea/group-eg/project_PromoterArchitecturePipeline
directory_path=/home/witham/Documents/pipeline_new/PromoterArchitecture
file_names=non-overlapping_includingbidirectional_all_genes

##extract promoters from a genome file
#python ../data_sorting/extract_promoter.py $directory_path $file_names --remove_bidirectional --prevent_overlapping_genes --fiveUTR
python ../data_sorting/extract_promoter.py $directory_path $file_names --prevent_overlapping_genes --fiveUTR
##extract promoters

## optional - merge meme files together into one file ready for FIMO
#$1 directory containing lots of sub directories containing .txt files containing motifs
#$2 output file name
#../meme_suite/./preFIMO.sh ../../data/FIMO/motif_data/dap_data_v4/motifs dap_combined.meme


##run preFIMO.sh script. $1 is promoter gff3 location. $2 is the genome.fasta file location.
#remember to change location depending on if promoter file used or UTR file
#promoter_gff_file=../../data/genomes/$file_names/promoters_renamedChr.gff3
promoter_gff_file=../../data/genomes/$file_names/promoters_5UTR_renamedChr.gff3
genome_fasta=../../data/genomes/TAIR10_chr_all.fas

../meme_suite/./preFIMO.sh $promoter_gff_file $genome_fasta $file_names

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
#$5 is the folder name

../meme_suite/./FIMO.sh ../../data/FIMO/$file_names/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/dap_combined.meme $file_names

## filter the FIMO output
conda activate PromoterArchitecturePipeline

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bed file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
python ../data_sorting/./FIMO_filter.py ../../data/FIMO/$file_names/output/${promoterpref}_FIMO/fimo.tsv ../../data/FIMO/$file_names/${promoterpref}.bed ../../data/FIMO/$file_names/${promoterpref}_motifs.bed 0.05

## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
../data_sorting/./TFBS_coverage.sh ../../data/FIMO/$file_names/${promoterpref}.bed

#%coverage of open chromatin
../data_sorting/./OpenChromatin_coverage.sh ../../data/FIMO/$3/$file_names.bed $file_names
