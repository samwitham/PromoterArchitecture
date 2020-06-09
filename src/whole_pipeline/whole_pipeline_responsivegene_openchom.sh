#!/bin/sh
#source /ei/projects/fa26d297-7b5b-4a43-be1f-bd87b73ae0d2/data/witham/PromoterArchitecture/software/miniconda3/bin/activate
#source /home/witham/opt/anaconda3/bin/activate
#eval "$(conda shell.bash hook)"
#conda activate PromoterArchitecturePipeline
##extract promoters from a genome file. Deactivated for now until argparse implemented
#python ../data_sorting/extract_promoter.py

##extract promoters

## optional - merge meme files together into one file ready for FIMO
#$1 directory containing lots of sub directories containing .txt files containing motifs
#$2 output file name
#../meme_suite/./preFIMO.sh ../../data/FIMO/motif_data/dap_data_v4/motifs dap_combined.meme


##run preFIMO.sh script. $1 is promoter gff3 location. $2 is the genome.fasta file location.
#promoter_bed_file=../../data/FIMO/responsivepromoters.bed
#promoter_bed_file=../../data/FIMO/promoters_5UTR.bed
promoter_bed_file=../../data/promoter_analysis/non-overlapping_includingbidirectional_variable_constitutive.bed
genome_fasta=../../data/genomes/TAIR10_chr_all.fas

#retrieve file name from file path
promoterbase=${promoter_bed_file##*/}
promoterpref=${promoterbase%.*}


#deactivated preFIMO.sh as don't need bedtool creation part
#../meme_suite/./preFIMO.sh $promoter_gff_file $genome_fasta

#create fasta file of promoters from genome fasta file and from the responsivepromoters.bed file
#bedtools getfasta -fi $genome_fasta -bed $promoter_bed_file -fo ../../data/FIMO/${promoterpref}.fasta -name


## create FIMO background file:
#activate correct conda env
#conda activate MemeSuite2
#identify the output filename created by preFIMO.sh

#run FIMO.sh
#$1 is promoter fasta file. 
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file 

#../meme_suite/./FIMO.sh ../../data/FIMO/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/dap_combined.meme

## filter the FIMO output
#conda activate PromoterArchitecturePipeline

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bedfile file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
#python ../data_sorting/./FIMO_filter.py ../../data/FIMO/output/${promoterpref}_FIMO/fimo.tsv ../../data/FIMO/${promoterpref}.bed ../../data/FIMO/${promoterpref}_motifs.bed 0.05
#python ../data_sorting/./FIMO_filter.py ../../data/FIMO/output/responsivepromoters_FIMO/fimo.tsv ../../data/FIMO/responsivepromoters.bed ../../data/FIMO/responsivepromoters_motifs.bed 0.05
### need to add the map_motif_ids.py option to this for if using DAP-seq cistrome motifs. Alternatively can use motif file containing gene IDs already.
## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
#../data_sorting/./TFBS_coverage.sh ../../data/FIMO/${promoterpref}.bed
#../data_sorting/./TFBS_coverage.sh ../../data/FIMO/responsivepromoters.bed

#%coverage of open chromatin
../data_sorting/./OpenChromatin_coverage.sh ../../data/promoter_analysis/${promoterpref}.bed


##map gene IDs - need to use generic names before activating this section, have to edit map_motif_ids.py for this
#python ../meme_suite/map_motif_ids.py ../../data/FIMO/responsivepromoters_motifs.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/FIMO/responsivepromoters_motifs_mapped.bed

