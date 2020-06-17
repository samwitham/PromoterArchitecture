#!/bin/sh
#source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate # need to change this to hardlink
eval "$(conda shell.bash hook)"
conda activate PromoterArchitecturePipeline
#directory_path=ei/workarea/group-eg/project_PromoterArchitecturePipeline
directory_path=/home/witham/Documents/pipeline_new/PromoterArchitecture #remove this, better to use relative path
file_names=non-overlapping_includingbidirectional_all_genes_newannotation

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
promoter_gff_file=../../data/output/$file_names/promoters_5UTR.gff3
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

../meme_suite/./FIMO.sh ../../data/output/$file_names/FIMO/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/dap_combined.meme $file_names

## filter the FIMO output
conda activate PromoterArchitecturePipeline

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bed file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_FIMO/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed 0.05

#If using DAP seq cistrome motifs without AGI names, map motifs
#arg1 is the input location of the motif bed file
#arg2 is the input location of the geneIDtable file derived from the .json (see .ipynb notebook)
#arg3 is the output location of the motif bed file
python ../meme_suite/./map_motif_ids.py ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed



## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
#$2 is the folder name
../data_sorting/./TFBS_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}.bed $file_names



##calculate promoter GC content
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of the promoter fasta
#arg3 is the output location of the promoters GC_content tsv file
python ../data_sorting/./promoter_GC_content.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.fasta ../../data/output/$file_names/GC_content/${promoterpref}_GC_content.tsv

#%coverage of open chromatin
#$1 is promoter bed file
#$2 is the folder name
../data_sorting/./OpenChromatin_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}.bed $file_names


#create sliding windows
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed 100 50



#Create Czechowski et al 2005 ranked cv dataset gene categories filtering out any genes not in the extracted promoters from this pipeline. The create subsets of N constitutive, variable or control genes
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the location of Czechowski et al 2005 ranked cv dataset reanalysed by Will Nash
#arg4 is the size, N, of the gene subsets
#arg5 is the gene category output file containing the selected gene subsets of size N
python ../data_sorting/./choose_genes_cv.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/genes/AtGE_dev_gcRMA__all_probes__CV.tsv 100 ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt


#TFBS coverage sliding window
#arg1 is the promoter extraction output folder name
#arg2 is the input location of motifs bed file
#arg3 is the output location of rolling window % coverage bed file
#arg4 is the input location of rolling window bed file
python ../rolling_window/./TFBScoverage_rw.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/rolling_window/TFBS_coverage_rw/${promoterpref}_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed


#GC content sliding window
#arg1 is the promoter extraction output folder name
#arg2 is the output location of the rolling window GC content tsv
#arg3 is the input location of the rolling window bed file
#arg4 is the input location of the genome fasta file
#arg5 is the output location of the rolling window fasta file
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/rolling_window/GC_content_rw/${promoterpref}_GCcontent_rw.tsv ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed $genome_fasta ../../data/output/$file_names/rolling_window/${promoterpref}_windows.fasta 

