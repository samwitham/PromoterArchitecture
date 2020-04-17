#!/bin/zsh
source /home/witham/opt/anaconda3/bin/activate
conda activate PromoterArchitecturePipeline
##extract promoters from a genome file. Deactivated for now until argparse implemented
#python ../data_sorting/extract_promoter.py

##extract promoters

## optional - merge meme files together into one file ready for FIMO
#$1 directory containing lots of sub directories containing .txt files containing motifs
#$2 output file name
#../meme_suite/./preFIMO.sh ../../data/FIMO/motif_data/dap_data_v4/motifs dap_combined.meme


##run preFIMO.sh script. $1 is promoter gff3 location. $2 is the genome.fasta file location.
#promoter_bed_file=../../data/FIMO/responsivepromoters.bed
#genome_fasta=../../data/genomes/TAIR10_chr_all.fas
promoter_fasta=../../../../meme_analysis/NLP7/promoters/NLP7promoters.fasta
#retrieve file name from file path
promoterbase=${promoter_fasta##*/}
promoterpref=${promoterbase%.*}

#create bed file from fasta
python ../meme_suite/fasta2bed.py ../../../../meme_analysis/NLP7/promoters/${promoterpref}.fasta ../../../../meme_analysis/NLP7/promoters/${promoterpref}.bed




## create FIMO background file:
#activate correct conda env
conda activate MemeSuite2
#identify the output filename created by preFIMO.sh

#run FIMO.sh
#$1 is promoter fasta file. 
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file 
#pvalue 1
#../meme_suite/./FIMO.sh ../../data/FIMO/${promoterpref}.fasta 1 5000000 ../../data/FIMO/motif_data/dap_combined.meme
#plantpan scan
../meme_suite/./FIMO.sh ../../../../meme_analysis/NLP7/promoters/${promoterpref}.fasta 1 5000000 ../../../../meme_analysis/NLP7/meme_m1.txt

## filter the FIMO output
conda activate PromoterArchitecturePipeline

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bedfile file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
#qvalue 1
#python ../data_sorting/./FIMO_filter.py ../../data/FIMO/output/${promoterpref}_FIMO/fimo.tsv ../../data/FIMO/${promoterpref}.bed ../../data/FIMO/${promoterpref}_motifs.bed 0.05
python ../data_sorting/./FIMO_filter.py ../../data/FIMO/output/${promoterpref}_FIMO/fimo.tsv ../../data/FIMO/${promoterpref}.bed ../../data/FIMO/${promoterpref}_motifs.bed 1

## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
#../data_sorting/./TFBS_coverage.sh ../../data/FIMO/${promoterpref}.bed
#../data_sorting/./TFBS_coverage.sh ../../data/FIMO/responsivepromoters.bed



##map gene IDs - need to use generic names before activating this section, have to edit map_motif_ids.py for this
python ../meme_suite/map_motif_ids.py ../../data/FIMO/${promoterpref}_motifs.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/FIMO/${promoterpref}_motifs_mapped.bed

