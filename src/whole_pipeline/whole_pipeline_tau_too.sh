#!/bin/sh
#source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate # need to change this to hardlink
eval "$(conda shell.bash hook)"
conda activate PromoterArchitecturePipeline
#directory_path=ei/workarea/group-eg/project_PromoterArchitecturePipeline
#directory_path=/home/witham/Documents/pipeline_new/PromoterArchitecture #remove this, better to use relative path
file_names=non-overlapping_includingbidirectional_all_genes_newannotation
#arabidopsis_selected_tissues
#non-overlapping_includingbidirectional_all_genes_newannotation

##extract promoters from a genome file
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the max length of promoters to be extracted upstream of the Araport TSS
#arg3 is an option to exclude potentially bidirectional promoters with upstream TSS >2000bp from TSS in opposite direction
#arg4 is an option to reduce the size of promoters if needed until they are not overlapping other genes.
#arg5 is an option to extend promoters to the first start codon.

#python ../data_sorting/extract_promoter.py $directory_path $file_names --remove_bidirectional --prevent_overlapping_genes --fiveUTR
python ../data_sorting/extract_promoter.py $file_names 1000 --prevent_overlapping_genes --fiveUTR
##extract promoters

## optional - merge meme files together into one file ready for FIMO
#$1 directory containing lots of sub directories containing .txt files containing motifs
#$2 output file name
#../meme_suite/./preFIMO.sh ../../data/FIMO/motif_data/dap_data_v4/motifs dap_combined.meme

#locations of promoter gff3 file and genome_fasta file
#remember to change location depending on if promoter file used or UTR file
#promoter_gff_file=../../data/genomes/$file_names/promoters_renamedChr.gff3
promoter_gff_file=../../data/output/$file_names/promoters_5UTR.gff3
genome_fasta=../../data/genomes/TAIR10_chr_all.fas

#identify the base filename
promoterbase=${promoter_gff_file##*/}
promoterpref=${promoterbase%.*}

#shortened promoter length
promoter_length=400

#window size and offset
window_size=100
window_offset=50


##run preFIMO.sh script. $1 is promoter gff3 location. $2 is the genome.fasta file location.
../meme_suite/./preFIMO.sh $promoter_gff_file $genome_fasta $file_names

## create FIMO background file:
#activate correct conda env
conda activate MemeSuite3

#run FIMO.sh
#$1 is promoter fasta file. 
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file
#$5 is the folder name
#scan for motifs using DAP-seq dataset
../meme_suite/./FIMO.sh ../../data/output/$file_names/FIMO/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/dap_combined.meme $file_names
#scan for motifs using plantpan dataset
../meme_suite/./FIMO_plantpan.sh ../../data/output/$file_names/FIMO/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/At_plantpan_BEEML.meme $file_names

## filter the FIMO output
conda activate PromoterArchitecturePipeline

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bed file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_DAPseq/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed 0.05

#filter plantpan motifs
python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_plantpan/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_plantpan.bed 0.05

#If using DAP seq cistrome motifs without AGI names, map motifs
#arg1 is the input location of the motif bed file
#arg2 is the input location of the geneIDtable file derived from the .json (see .ipynb notebook)
#arg3 is the output location of the motif bed file
python ../meme_suite/./map_motif_ids.py  ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed 

#map plantpan motifs
#python ../meme_suite/./map_motif_ids.py ../../data/output/$file_names/FIMO/${promoterpref}_motifs_plantpan.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs__plantpan_mapped.bed



#filter microarray to include only conditions of interest for TAU tissue specificity and CV values - not filtering any more, just generating TAU table
#arg1 is the input location of expression data table (schmid et al 2005 microarray AtGE_dev_gcRMA.txt.newline)
#arg2 is the output location of filtered microarray
#arg3 is the output location of TAU table
#arg4 is the promoter extraction output folder name

python ../data_sorting/./filter_microarray_conditions.py ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline.filtered ../../data/output/$file_names/genes/tissue_specific/${promoterpref}_schmid_allfilteredgenes_TAU.txt $file_names


#Prepare raw Czechowski et al 2005 microarray data filtering out genes that aren't expressed in >=80% of tissues/conditions
#arg1 is the input file
#arg2 is the mas5 affymetrix presence/absence data
#arg3 is the current Arabidopsis thaliana annotation data
#arg4 is the Araport housekeeping genes from Data S4 from Cheng et al. 2016
#arg5 is the output directory
#arg6 is the promoter extraction output folder name
#arg7 is the filter_threshold (percentage of conditions genes have to be expressed in otherwise they are filtered)

python ../data_sorting/./expressionVar_AtGE_dev_gcRMA_editedbySam.py ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline ../../data/genes/AtGeneExpress_CV_2020/E-TABM-17_affy_mas5_calls_20200529.txt ../../data/genomes/Arabidopsis_thaliana/annotation/Arabidopsis_thaliana.TAIR10.47.gff3 ../../data/genes/AtGeneExpress_CV_2020/Araport11_housekeeping_genes_fromDataS4_Chengetal2016.txt ../../data/output/$file_names/genes/filtered80 $file_names 80

#remove the incorrect header line of the file - no longer needed fixed in the expressionVar_AtGE_dev_gcRMA.py and expressionVar_AtGE_dev_gcRMA_editedbySam.py scripts
#tail -n +2 ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv > ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV_noheader.tsv

#Prepare raw Czechowski et al 2005 microarray data this time not filtering out genes that aren't expressed in >80% of tissues/conditions (output file now currently used)
#arg1 is the input file
#arg2 is the mas5 affymetrix presence/absence data
#arg3 is the current Arabidopsis thaliana annotation data
#arg4 is the Araport housekeeping genes from Data S4 from Cheng et al. 2016
#arg5 is the output directory
#arg6 is the promoter extraction output folder name
#arg7 is the filter_threshold (percentage of conditions genes have to be expressed in otherwise they are filtered)

python ../data_sorting/./expressionVar_AtGE_dev_gcRMA_editedbySam.py ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline ../../data/genes/AtGeneExpress_CV_2020/E-TABM-17_affy_mas5_calls_20200529.txt ../../data/genomes/Arabidopsis_thaliana/annotation/Arabidopsis_thaliana.TAIR10.47.gff3 ../../data/genes/AtGeneExpress_CV_2020/Araport11_housekeeping_genes_fromDataS4_Chengetal2016.txt ../../data/output/$file_names/genes/tissue_specific $file_names 0

#remove the incorrect header line of the file - no longer needed fixed in the expressionVar_AtGE_dev_gcRMA.py and expressionVar_AtGE_dev_gcRMA_editedbySam.py scripts
#tail -n +2 ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv > ../../data/output/$file_names/genes/tissue_specific/AtGE_dev_gcRMA__all_probes__CV_noheader.tsv



#Create Czechowski et al 2005 ranked cv dataset gene categories filtering out any genes not in the extracted promoters from this pipeline. The create subsets of N constitutive, variable or control genes
#Also filters out promoters which have 100% overlapping promoters with other genes (where only a 5UTR is present that's not overlapping)
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the location of the Czechowski et al 2005 ranked cv dataset reanalysed by Will Nash
#arg4 is the input location of the Mergner et al 2020 ranked cv dataset
#arg5 is the size, N, of the gene subsets
#arg6 is the czechowski gene category output file containing the selected gene subsets of size N
#arg7 is the mergner gene category output file containing the selected gene subsets of size N
#arg8 is the input location of the promoter mapped motifs bed file
#arg9 is the output location of the promoter bed file filtered so that each promoter contains at least one TFBS
#arg10 is the output location of all filtered microarray genes
#arg11 is the output location of all filtered RNAseq genes
#arg12 is the input location of promoters gff3 file
python ../data_sorting/./choose_genes_cv.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV_noheader.tsv ../../data/genes/RNA_CVs.csv 100 ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/${file_names}/genes/${promoterpref}_mergner_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_filtered_contain_motifs.bed ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt ../../data/output/$file_names/genes/${promoterpref}_mergner_allfilteredgenes.txt ../../data/output/$file_names/promoters.gff3

#run again with 300 genes per promoter category for gene ontology enrichment analysis
python ../data_sorting/./choose_genes_cv.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV_noheader.tsv ../../data/genes/RNA_CVs.csv 300 ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random_300.txt ../../data/output/${file_names}/genes/${promoterpref}_mergner_constitutive_variable_random_300.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_filtered_contain_motifs.bed ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt ../../data/output/$file_names/genes/${promoterpref}_mergner_allfilteredgenes.txt ../../data/output/$file_names/promoters.gff3

#Create Scmid et al 2005 ranked tau dataset gene categories filtering out any genes not in the extracted promoters from this pipeline. The create subsets of N constitutive, tissue_specific or control genes
#Also filters out promoters which have 100% overlapping promoters with other genes (where only a 5UTR is present that's not overlapping)
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the location of the Czechowski et al 2005 ranked tau dataset
#arg4 is the size, N, of the gene subsets
#arg5 is the schmid gene category output file containing the selected gene subsets of size N
#arg6 is the input location of the promoter mapped motifs bed file
#arg7 is the output location of the promoter bed file filtered so that each promoter contains at least one TFBS
#arg8 is the output location of all filtered microarray genes
#arg9 is the input location of promoters gff3 file
#arg10 is the location of the gene categories ranked by coefficient of variation
#arg11 is the optional input location of coefficient of variation all genes after filtering to include only genes present in 80% of condition. This is used to choose control genes that are also present in this file and then replace the CV control genes with the tau control genes

python ../data_sorting/./choose_genes_tau.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt 100 ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_filtered_contain_motifs.bed ../../data/output/$file_names/genes/${promoterpref}_schmid_allfilteredgenes.txt ../../data/output/$file_names/promoters.gff3 ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt

#run again with 300 genes per promoter category for gene ontology enrichment analysis
python ../data_sorting/./choose_genes_tau.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt 300 ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random_300.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_filtered_contain_motifs.bed ../../data/output/$file_names/genes/${promoterpref}_schmid_allfilteredgenes.txt ../../data/output/$file_names/promoters.gff3 ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random_300.txt ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt


## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
#$2 is the folder name
#$3 is the motifs file location
#$4 is the output file location 
../data_sorting/./TFBS_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}.bed $file_names ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/TFBS_coverage/${promoterpref}.bp_covered.txt

##calculate promoter GC content
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of the promoter fasta
#arg3 is the output location of the promoters GC_content tsv file
python ../data_sorting/./promoter_GC_content.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.fasta ../../data/output/$file_names/GC_content/${promoterpref}_GC_content.tsv

#%coverage of open chromatin
#also generates a motifs file containing only TFBSs which overlap with open chromatin
#$1 is promoter bed file
#$2 is the folder name
#$3 is the TFBS motifs bed file
#$4 is the TFBS motifs_mapped bed file
../data_sorting/./OpenChromatin_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}.bed $file_names ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed 





#create sliding windows
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ${window_size} ${window_offset} ../../data/output/$file_names/overlapping_promoters_5UTR.bed
#convert promoters.gff3 into bed file
#convert gff file to bed file
gff2bed < ../../data/output/$file_names/promoters.gff3 > ../../data/output/$file_names/FIMO/promoters.bed
#remove 'gene:' from column 4
awk 'BEGIN {FS = OFS = "\t"} { gsub(/gene:/,"", $4); print}' ../../data/output/$file_names/FIMO/promoters.bed > ../../data/output/$file_names/FIMO/promoters.bed.tmp && mv ../../data/output/$file_names/FIMO/promoters.bed.tmp ../../data/output/$file_names/FIMO/promoters.bed
#remove lines not beginning with numbers (ie. remove mt and pt chromosomes)
grep '^[0-9]' ../../data/output/$file_names/FIMO/promoters.bed > ../../data/output/$file_names/FIMO/promoters.bed.tmp && mv ../../data/output/$file_names/FIMO/promoters.bed.tmp ../../data/output/$file_names/FIMO/promoters.bed


#create eukaryotic promoter database promoter (promoters going up until EPD annotated TSS)
#Also create 500 bp eukaryotic promoter database promoter + 100 bp 5'UTR bed file
#arg1 is the promoter extraction output folder name
#arg2 is the input location of promoter-5UTR bedfile
#arg3 is the input location of Eukaryotic Promoter Database TSS bed file
#arg4 is the output location of the EPD_promoters bedfile
#arg5 is the output location of the EPD_promoters_5UTR bedfile
#arg6 is the output location of flagged promoters which are not in the Eukaryotic promoter database
#arg7 is the output location of flagged EPD_promoters which have TSSs in coding regions
#arg8 is the output location of flagged EPD_promoters which are overlapping other genes so they are only a shortened 5\'UTR
#arg9 is the length of promoter (bp) to use in the promoter-5UTR bedfile output
#arg10 is the length of 5'UTR (bp) to use in the promoter_5UTR bedfile output
#arg11 is the output location of flagged specified length EPD_promoters_5UTR which overlap other genes

python ../data_sorting/./create_EPD_promoters.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/EPD_promoters.bed ../../data/output/$file_names/FIMO/EPD_promoters_5UTR.bed ../../data/output/$file_names/FIMO/flagged_proms_not_in_EPD.bed ../../data/output/$file_names/FIMO/flagged_EPD_TSS_in_CDS.bed ../../data/output/$file_names/FIMO/flagged_EPD_overlappingprom_so_only5UTRs.bed 500 100 ../../data/output/$file_names/FIMO/flagged_EPD_overlapping_specified_length.bed

#create sliding windows starting from Araport TSS instead (use only promoters not 5'UTR)
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/promoters.bed ../../data/output/$file_names/rolling_window/promoters_windows.bed ${window_size} ${window_offset} ../../data/output/$file_names/overlapping_promoters.bed

#create sliding windows starting from EPD TSS instead (use only promoters not 5'UTR)
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/EPD_promoters.bed ../../data/output/$file_names/rolling_window/EPD_promoters_windows.bed ${window_size} ${window_offset} ../../data/output/$file_names/EPD_overlapping_promoters.bed




#create Araport11 5'UTR bed file and artificially swap the strand from + to - and vice versa. 
#This is so the window number correctly starts from the TSS and goes downstream towards the ATG)
#arg1 is the promoter extraction output folder name
#arg2 is the input location of promoter bedfile
#arg3 is the input location of promoter-5UTR bedfile
#arg4 is the output location of the 5UTR bedfile
#arg5 is the output location of the 5UTR bedfile with artificially switched strands for sliding window analysis
#arg6 is the output location of genes which have no 5UTR

python ../data_sorting/./create_5UTRs.py $file_names ../../data/output/$file_names/FIMO/promoters.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/Araport11_5UTR.bed ../../data/output/$file_names/FIMO/Araport11_5UTR_swapped_strands.bed ../../data/output/$file_names/FIMO/genes_no_5UTR.bed

#create Eukaryotic promoter database (EPD) 5'UTR bed file and artificially swap the strand from + to - and vice versa. 
#This is so the window number correctly starts from the TSS and goes downstream towards the ATG)
#arg1 is the promoter extraction output folder name
#arg2 is the input location of promoter bedfile
#arg3 is the input location of promoter-5UTR bedfile
#arg4 is the output location of the 5UTR bedfile
#arg5 is the output location of the 5UTR bedfile with artificially switched strands for sliding window analysis
#arg6 is the output location of genes which have no 5UTR

python ../data_sorting/./create_5UTRs.py $file_names ../../data/output/$file_names/FIMO/EPD_promoters.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/EPD_5UTR.bed ../../data/output/$file_names/FIMO/EPD_5UTR_swapped_strands.bed ../../data/output/$file_names/FIMO/EPD_genes_no_5UTR.bed

#create sliding windows starting from the Araport TSS for the 5'UTRs with numbering going downstream towards the ATG start codon
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/Araport11_5UTR_swapped_strands.bed ../../data/output/$file_names/rolling_window/5UTR_windows_swapped_strands.bed ${window_size} ${window_offset} ../../data/output/$file_names/overlapping_5UTRs.bed

#create sliding windows starting from the EPD TSS for the 5'UTRs with numbering going downstream towards the ATG start codon
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/EPD_5UTR_swapped_strands.bed ../../data/output/$file_names/rolling_window/EPD_5UTR_windows_swapped_strands.bed ${window_size} ${window_offset} ../../data/output/$file_names/EPD_overlapping_5UTRs.bed


#merge the promoters and Araport11 5UTR rolling windows, with window numbers going head to head around the Araport TSS (so window 1 of promoters becomes -1, while 5UTR window numbers stay the same)
#arg1 is the promoter extraction output folder name
#arg2 is the input location of 5UTR sliding windows
#arg3 is the input location of the promoter 5UTR bed file to get strand information
#arg4 is the input location of the promoter sliding windows file
#arg5 is the output location of the promoter + 5UTR sliding window numbers going outwards from the TSS
#arg6 is the size of the rolling window in bp
#arg7 is the size of the window offset in bp
python ../rolling_window/./TSS_outward_rw.py $file_names ../../data/output/$file_names/rolling_window/5UTR_windows_swapped_strands.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed ${window_size} ${window_offset}

#merge the promoters and EPD 5UTR rolling windows, with window numbers going head to head around the Araport TSS (so window 1 of promoters becomes -1, while 5UTR window numbers stay the same)
#arg1 is the promoter extraction output folder name
#arg2 is the input location of 5UTR sliding windows
#arg3 is the input location of the promoter 5UTR bed file to get strand information
#arg4 is the input location of the promoter sliding windows file
#arg5 is the output location of the promoter + 5UTR sliding window numbers going outwards from the TSS
#arg6 is the size of the rolling window in bp
#arg7 is the size of the window offset in bp
python ../rolling_window/./TSS_outward_rw.py $file_names ../../data/output/$file_names/rolling_window/EPD_5UTR_windows_swapped_strands.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/rolling_window/EPD_promoters_windows.bed ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed ${window_size} ${window_offset}






#TFBS coverage sliding window promoters/5'UTR
#arg1 is the promoter extraction output folder name
#arg2 is the name of the output directory
#arg3 is the input location of rolling window bed file
#arg4 is the input location of bed file for which % coverage is being calculated for
#arg4 is the output location of the rolling window % coverage bed file
python ../rolling_window/./coverage_rw.py $file_names TFBS_coverage_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/rolling_window/TFBS_coverage_rw/${promoterpref}_bpcovered_rw.bed

#TFBS coverage sliding window promoters only (starting at Araport11 TSS)
#arg1 is the promoter extraction output folder name
#arg2 is the name of the output directory
#arg3 is the input location of rolling window bed file
#arg4 is the input location of bed file for which % coverage is being calculated for
#arg4 is the output location of the rolling window % coverage bed file
python ../rolling_window/./coverage_rw.py $file_names TFBS_coverage_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/rolling_window/TFBS_coverage_rw_promoterno5UTR/promoters_bpcovered_rw.bed

#TFBS coverage sliding window promoters/5UTR TSS outward windows (Araport 11 TSS)
#arg1 is the promoter extraction output folder name
#arg2 is the name of the output directory
#arg3 is the input location of rolling window bed file
#arg4 is the input location of bed file for which % coverage is being calculated for
#arg4 is the output location of the rolling window % coverage bed file
python ../rolling_window/./coverage_rw.py $file_names TFBS_coverage_rw_Araport11_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/rolling_window/TFBS_coverage_rw_Araport11_TSS_outward_promoter5UTR/Araport11_TSS_outward_promoter5UTR_bpcovered_rw.bed


#TFBS coverage sliding window promoters/5UTR TSS outward windows (EPD TSS)
#arg1 is the promoter extraction output folder name
#arg2 is the name of the output directory
#arg3 is the input location of rolling window bed file
#arg4 is the input location of bed file for which % coverage is being calculated for
#arg4 is the output location of the rolling window % coverage bed file
python ../rolling_window/./coverage_rw.py $file_names TFBS_coverage_rw_EPD_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/rolling_window/TFBS_coverage_rw_EPD_TSS_outward_promoter5UTR/EPD_TSS_outward_promoter5UTR_bpcovered_rw.bed



#Root open chromatin coverage sliding window promoters_5'UTR
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed

#Shoot open chromatin coverage sliding window promoters_5'UTR
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed

#Root/Shoot intersect coverage sliding window promoters_5'UTR
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters only
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_promoterno5UTR/promoters_root_bpcovered_rw.bed

#Shoot open chromatin coverage sliding window promoters only
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_promoterno5UTR/promoters_shoot_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters only
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_promoterno5UTR/promoters_rootshootintersect_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters/5UTR TSS outward windows (Araport11)
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_Araport11_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_Araport11_TSS_outward_promoter5UTR/Araport11_TSS_outward_promoter5UTR_root_bpcovered_rw.bed

#Shoot open chromatin coverage sliding window promoters/5UTR TSS outward windows (Araport11)
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_Araport11_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_Araport11_TSS_outward_promoter5UTR/Araport11_TSS_outward_promoter5UTR_shoot_bpcovered_rw.bed

#Root/Shoot intersect coverage sliding window promoters/5UTR TSS outward windows (Araport11)
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_Araport11_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_Araport11_TSS_outward_promoter5UTR/Araport11_TSS_outward_promoter5UTR_rootshootintersect_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters/5UTR TSS outward windows (EPD)
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_EPD_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_EPD_TSS_outward_promoter5UTR/EPD_TSS_outward_promoter5UTR_root_bpcovered_rw.bed

#Shoot open chromatin coverage sliding window promoters/5UTR TSS outward windows (EPD)
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_EPD_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_EPD_TSS_outward_promoter5UTR/EPD_TSS_outward_promoter5UTR_shoot_bpcovered_rw.bed

#Root/Shoot intersect coverage sliding window promoters/5UTR TSS outward windows (EPD)
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_EPD_TSS_outward_promoter5UTR ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_EPD_TSS_outward_promoter5UTR/EPD_TSS_outward_promoter5UTR_rootshootintersect_bpcovered_rw.bed




#GC content sliding window promoters & 5'UTR
#arg1 is the promoter extraction output folder name
#arg2 is the output location of the rolling window GC content tsv
#arg3 is the input location of the rolling window bed file
#arg4 is the input location of the genome fasta file
#arg5 is the output location of the rolling window fasta file
#arg6 is the output folder name
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/rolling_window/GC_content_rw/${promoterpref}_GCcontent_rw.tsv ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed $genome_fasta ../../data/output/$file_names/rolling_window/${promoterpref}_windows.fasta GC_content_rw

#GC content sliding window promoters only
#arg1 is the promoter extraction output folder name
#arg2 is the output location of the rolling window GC content tsv
#arg3 is the input location of the rolling window bed file
#arg4 is the input location of the genome fasta file
#arg5 is the output location of the rolling window fasta file
#arg6 is the output folder name
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/rolling_window/GC_content_rw_promoterno5UTR/promoters_GCcontent_rw.tsv ../../data/output/$file_names/rolling_window/promoters_windows.bed $genome_fasta ../../data/output/$file_names/rolling_window/promoters_windows.fasta GC_content_rw_promoterno5UTR

#GC content sliding window promoters/5UTR TSS outward windows (Araport11)
#arg1 is the promoter extraction output folder name
#arg2 is the output location of the rolling window GC content tsv
#arg3 is the input location of the rolling window bed file
#arg4 is the input location of the genome fasta file
#arg5 is the output location of the rolling window fasta file
#arg6 is the output folder name
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/rolling_window/GC_content_rw_Araport11_TSS_outward_promoter5UTR/Araport11_TSS_outward_promoter5UTR_GCcontent_rw.tsv ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed $genome_fasta ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.fasta  GC_content_rw_Araport11_TSS_outward_promoter5UTR

#GC content sliding window promoters/5UTR TSS outward windows (EPD)
#arg1 is the promoter extraction output folder name
#arg2 is the output location of the rolling window GC content tsv
#arg3 is the input location of the rolling window bed file
#arg4 is the input location of the genome fasta file
#arg5 is the output location of the rolling window fasta file
#arg6 is the output folder name
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/rolling_window/GC_content_rw_EPD_TSS_outward_promoter5UTR/EPD_TSS_outward_promoter5UTR_GCcontent_rw.tsv ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed $genome_fasta ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.fasta  GC_content_rw_EPD_TSS_outward_promoter5UTR



#TF_diversity sliding window promoters & 5'UTR
#arg1 is the promoter extraction output folder name
#arg2 is the input location of the rolling window bed file
#arg3 is the input location of the promoters mapped motif bed
#arg4 is the output location of windows_motifs intersect bed
#arg5 is the output location of the window TF_diversity_bed
#arg6 is the output folder name
python ../rolling_window/./TF_diversity_rw.py $file_names ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/rolling_window/${promoterpref}_windows_motifs.bed ../../data/output/$file_names/rolling_window/TF_diversity_rw/${promoterpref}_TF_diversity.bed TF_diversity_rw

#TF_diversity sliding window promoters only
#arg1 is the promoter extraction output folder name
#arg2 is the input location of the rolling window bed file
#arg3 is the input location of the promoters mapped motif bed
#arg4 is the output location of windows_motifs intersect bed
#arg5 is the output location of the window TF_diversity_bed
#arg6 is the output folder name
python ../rolling_window/./TF_diversity_rw.py $file_names ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/rolling_window/promoters_windows_motifs.bed ../../data/output/$file_names/rolling_window/TF_diversity_rw_promoterno5UTR/promoters_TF_diversity.bed TF_diversity_rw_promoterno5UTR

#TF_diversity sliding window promoters/5UTR TSS outward windows (Araport11)
#arg1 is the promoter extraction output folder name
#arg2 is the input location of the rolling window bed file
#arg3 is the input location of the promoters mapped motif bed
#arg4 is the output location of windows_motifs intersect bed
#arg5 is the output location of the window TF_diversity_bed
#arg6 is the output folder name
python ../rolling_window/./TF_diversity_rw.py $file_names ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/rolling_window/Araport11_TSS_outward_promoter5UTR_windows_motifs.bed ../../data/output/$file_names/rolling_window/TF_diversity_rw_Araport11_TSS_outward_promoter5UTR/Araport11_TSS_outward_promoter5UTR_TF_diversity.bed TF_diversity_rw_Araport11_TSS_outward_promoter5UTR

#TF_diversity sliding window promoters/5UTR TSS outward windows (EPD)
#arg1 is the promoter extraction output folder name
#arg2 is the input location of the rolling window bed file
#arg3 is the input location of the promoters mapped motif bed
#arg4 is the output location of windows_motifs intersect bed
#arg5 is the output location of the window TF_diversity_bed
#arg6 is the output folder name
python ../rolling_window/./TF_diversity_rw.py $file_names ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/rolling_window/EPD_TSS_outward_promoter5UTR_windows_motifs.bed ../../data/output/$file_names/rolling_window/TF_diversity_rw_EPD_TSS_outward_promoter5UTR/EPD_TSS_outward_promoter5UTR_TF_diversity.bed TF_diversity_rw_EPD_TSS_outward_promoter5UTR



#prepare files for gat analysis TATA enrichment
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the Input location of TATAbox_location bed file (from Eukaryotic promoter database)
#use ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed (stringent TATAs))
#Used the following search parameters for download:
## FindM Genome Assembly : A. thaliana (Feb 2011 TAIR10/araTha1)
##Series : EPDnew, the Arabidopsis Curated Promoter Database
##Sample : TSS from EPDnew rel 004
##Repeat masking: off
##5' border: -50     3' border: 0 
##Search mode: forward
##Selection mode : all matches)

#or ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_renamed.bed (Downloaded TATAbox_location.bed from EPD to data/EPD_promoter_analysis/EPDnew_promoters

#This is less stringent than before, selects TATA boxes +-100bp from most common EPD TSS
#Used the following search parameters for download:
## FindM Genome Assembly : A. thaliana (Feb 2011 TAIR10/araTha1)
##Series : EPDnew, the Arabidopsis Curated Promoter Database
##Sample : TSS from EPDnew rel 004
##Repeat masking: off
##5' border: -100     3' border: 100
##Search mode: forward
##Selection mode : all matches)
python ../data_sorting/./TATA_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Czechowski ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed

#run gat (Genomic association tester) enrichment for TATA boxes using Czechowski gene categories
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#use ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed (stringent TATAs) or ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_renamed.bed (Downloaded TATAbox_location.bed from EPD to data/EPD_promoter_analysis/EPDnew_promoters

#This is less stringent than before, selects TATA boxes +-100bp from most common EPD TSS
#Used the following search parameters for download:
## FindM Genome Assembly : A. thaliana (Feb 2011 TAIR10/araTha1)
##Series : EPDnew, the Arabidopsis Curated Promoter Database
##Sample : TSS from EPDnew rel 004
##Repeat masking: off
##5' border: -100     3' border: 100
##Search mode: forward
##Selection mode : all matches)
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#7 is the optional --ignore-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TATA/gat_analysis/Czechowski_${promoterpref}_workspace.bed ../../data/output/$file_names/TATA/gat_analysis/Czechowski_${promoterpref}_constitutive.bed ../../data/output/$file_names/TATA/gat_analysis/Czechowski_${promoterpref}_variable.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed ../../data/output/$file_names/TATA/gat_analysis ${promoterpref}_Czechowski_TATA --ignore-segment-tracks



##PLOTS
#Whole promoter plots:

#Plot GC content
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters GC_content tsv file
#arg4 is the output folder name ending in a forward slash
python ../plotting/./GC_content_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/GC_content/${promoterpref}_GC_content.tsv 


#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}.bp_covered.txt

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
python ../plotting/./TF_diversity_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed

#TATA enrichment plot Czechowski gene categories
#arg1 is the promoter extraction output folder name
#arg2 is the Location of constitutive promoter gat analysis output
#arg3 is the Location of variable promoter gat analysis output
#arg4 is the Output prefix to add to plot file name
python ../plotting/./TATA_enrichment_plots.py $file_names ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Czechowski_TATA_constitutive.out ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Czechowski_TATA_variable.out Czechowski_${promoterpref}

#Open chromatin coverage shoot-root intersect using Potter et al 2018 ATAC-seq negative controls for root/shoot open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of open chromatin bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash - open chromatin tissue
#arg5 is the optional output folder name ending in a forward slash - promoter set name 
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt RootShootIntersect/ ${promoterpref}/
#Open chromatin coverage root
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}RootOpenChrom.bp_covered.txt Root/ ${promoterpref}/
#Open chromatin coverage shoot
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootOpenChrom.bp_covered.txt Shoot/ ${promoterpref}/



#rerun analyses at shorter promoter length

#optionally run post-FIMO analysis at specific promoter length

#This shortens to promoter to 400bp upstream of the ATG
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is promoterpref - prefix name for the promoter files
#arg3 is the location of the extracted promoters/5UTRs
#arg4 is the length to shorten promoters to

python ../data_sorting/./shorten_promoters.py $file_names ${promoterpref} ../../data/output/$file_names/FIMO/${promoterpref}.bed $promoter_length


#TFBS coverage of shortened promoters
#arg1 is the promoter extraction output folder name
#arg2 is the name of the output directory
#arg3 is the input location of shortened promoter bed
#arg4 is the input location of bed file for which % coverage is being calculated for
#arg4 is the output location of the % coverage bed file
python ../rolling_window/./coverage_rw.py $file_names TFBS_coverage ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bp.bp_covered.txt


#GC content of shortened promoters
#arg1 is the promoter extraction output folder name
#arg2 is the output location of the GC content tsv
#arg3 is the input location of the shortened promoter bed file
#arg4 is the input location of the genome fasta file
#arg5 is the output location of the shortened promoter fasta file
#arg6 is the output folder name
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/GC_content/${promoterpref}_${promoter_length}bp_GC_content.tsv ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed $genome_fasta ../../data/output/$file_names/${promoterpref}_${promoter_length}bp.fasta 

#Create shortened promoter motifs_mapped file
#arg1 is the input location of the shortened promoter file
#arg2 is the input location of the promoters mapped motif bed
#arg3 is the output location of the shortened mapped_motifs bed
python ../data_sorting/create_motif_mapped_from_intersect.py ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed

#%coverage of open chromatin
#$1 is promoter bed file
#$2 is the folder name
#$3 is the motifs bed file location
#$4 is the TFBS motifs_mapped bed file
../data_sorting/./OpenChromatin_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed $file_names ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed 




##PLOT SLIDING WINDOWS

#GC content sliding window
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of GC content tsv file
#arg4 is the input location of eukaryotic promoter database transcription start site bed file
#arg5 is the input location of promoter bed file
#arg6 is the input location of promoter no 5UTR bed file
#arg7 is the output folder name prefix to use
#arg8 is the input location of root chromatin bed file
#arg9 is the input location of shoot chromatin bed file
#arg10 is the input location of rootshootintersect chromatin bed file
#arg11 is the optional replacement colour palette for plots',default = None
python ../plotting/rolling_window/./GC_content_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/rolling_window/GC_content_rw/${promoterpref}_GCcontent_rw.tsv ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 GC_content_rw ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed deep

#TF diversity sliding window plot
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of the window bed file
#arg4 is the input location of the TF diversity bed file
#arg5 is the input location of eukaryotic promoter database transcription start site bed file
#arg6 is the input location of promoter bed file
#arg7 is the input location of promoter no 5UTR bed file
#arg8 is the output folder name prefix to use
#arg9 is the input location of root chromatin bed file
#arg10 is the input location of shoot chromatin bed file
#arg11 is the input location of rootshootintersect chromatin bed file
#arg12 is the optional replacement colour palette for plots',default = None
#arg13 is the optional author name to add to output file names',default = 'Czechowski
#arg14 is the optional variable 2 name eg. tissue_specific',default = 'variable'

python ../plotting/rolling_window/./TF_diversity_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/output/$file_names/rolling_window/TF_diversity_rw/${promoterpref}_TF_diversity.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 TF_diversity_rw ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed deep

#TFBS coverage sliding window plot

#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of the TFBS base pairs covered file
#arg4 is the input location of eukaryotic promoter database transcription start site bed file
#arg5 is the input location of promoter bed file
#arg6 is the input location of promoter no 5UTR bed file
#arg7 is the output folder name prefix to use
#arg8 is the input location of root chromatin bed file
#arg9 is the input location of shoot chromatin bed file
#arg10 is the input location of rootshootintersect chromatin bed file
#arg11 is the optional replacement colour palette for plots',default = None
#arg12 is the optional author name to add to output file names',default = 'Czechowski
#arg13 is the optional variable 2 name eg. tissue_specific',default = 'variable'

python ../plotting/rolling_window/./TFBScoverage_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/rolling_window/TFBS_coverage_rw/${promoterpref}_bpcovered_rw.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 TFBS_coverage_rw ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed deep

#open chromatin coverage sliding window plot

#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of the GC content TSV file
#arg4 is the input location of eukaryotic promoter database transcription start site bed file
#arg5 is the input location of promoter bed file
#arg6 is the input location of promoter no 5UTR bed file
#arg7 is the output folder name prefix to use
#arg8 is the input location of root chromatin bed file
#arg9 is the input location of shoot chromatin bed file
#arg10 is the input location of rootshootintersect chromatin bed file
#arg11 is the optional replacement colour palette for plots',default = None
#arg12 is the optional author name to add to output file names',default = 'Czechowski
#arg13 is the optional variable 2 name eg. tissue_specific',default = 'variable'
python ../plotting/rolling_window/./openchromatin_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/rolling_window/GC_content_rw/${promoterpref}_GCcontent_rw.tsv ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 OpenChromatin_rw ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed deep

#PLOTTING

#Shortened promoter plots:

#Plot GC content
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters GC_content tsv file
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./GC_content_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/GC_content/${promoterpref}_${promoter_length}bp_GC_content.tsv ${promoterpref}_${promoter_length}bp/


#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bp.bp_covered.txt ${promoterpref}_${promoter_length}bp/


#Create shortened promoter motifs_mapped file in correct format for the TF_diversity_plots.py script
#arg1 is the input location of the shortened promoter file
#arg2 is the input location of the promoters mapped motif bed
#arg3 is the output location of the shortened mapped_motifs bed
#python ../data_sorting/create_motif_mapped_from_intersect.py ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TF_diversity_plots_shortenedprom.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp/

#Open chromatin coverage shoot-root intersect using Potter et al 2018 ATAC-seq negative controls for root/shoot open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of open chromatin bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash - open chromatin tissue
#arg5 is the optional output folder name ending in a forward slash - promoter set name 
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt RootShootIntersect/ ${promoterpref}_${promoter_length}bp/
#Open chromatin coverage root
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpRootOpenChrom.bp_covered.txt Root/ ${promoterpref}_${promoter_length}bp/
#Open chromatin coverage shoot
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootOpenChrom.bp_covered.txt Shoot/ ${promoterpref}_${promoter_length}bp/


#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, whole promoters open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Czechowski ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_openchromrootshootintersect_TFfamily.bed

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, whole promoters
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Czechowski ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_TFfamily.bed

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, 400bp promoters, openchromatin
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed Czechowski_400bp ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect_TFfamily.bed

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, 400bp promoters
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed Czechowski_400bp ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_TFfamily.bed



#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - whole promoters open chromatin only
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_openchromrootshootintersect_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoterpref}_Czechowski_TFBS_openchrom --with-segment-tracks

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - whole promoters 
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoterpref}_Czechowski_TFBS --with-segment-tracks

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - 400bp promoters open chromatin only
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks flag or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoter_length}bp_${promoterpref}_Czechowski_TFBS_openchrom --with-segment-tracks

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - 400bp promoters
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks flag or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoter_length}bp_${promoterpref}_Czechowski_TFBS --with-segment-tracks

#run precentrimo.sh - 400bp promoters CV
#$1 is the promoter bed file location of constitutive genes only
#$2 is the promoter bed file location of constitutive and variable genes (background genes)
#$3 is the folder name
#$4 is the genome fasta file location
#$5 is the constitutive promoter fasta file output location
#$6 is the constitutive, variable and random promoters fasta file output location
../meme_suite/./precentrimo.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable_gat.bed $file_names $genome_fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta

#run precentrimo.sh - 400bp promoters CV
#$1 is the promoter bed file location of variable genes only
#$2 is the promoter bed file location of constitutive and variable genes (background genes)
#$3 is the folder name
#$4 is the genome fasta file location
#$5 is the constitutive promoter fasta file output location
#$6 is the constitutive, variable and random promoters fasta file output location
../meme_suite/./precentrimo.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable_gat.bed $file_names $genome_fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_variable.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta


#activate memesuite conda environment (install memesuite 5.0.2 using conda install 'meme=5.0.2' 'icu=58.2'. Centrimo is broken on memesuite versions higher than 5.0.2)
#MemeSuite4 env has 
conda activate MemeSuite4
#run centrimo to look for enriched motifs in constitutive promoters (using whole promoters)
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the promoter fasta file location of the background genes
#$3 is meme motif file 
#$4 is the folder name
../meme_suite/./centrimo.sh ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta ../../data/FIMO/motif_data/dap_combined.meme $file_names 

../meme_suite/./centrimo.sh ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_variable.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta ../../data/FIMO/motif_data/dap_combined.meme $file_names 


#run ciiider software - constitutive genes
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the background gene fasta file
#$3 is jaspar motif file
#$4 is the folder name
#$5 is the CiiiDER output folder name (Czechowski_400bp_promoters_5UTR_constitutive)
../CiiiDER/ciiider.sh ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta ../../data/FIMO/motif_data/DAP_seq_motifs.jaspar $file_names Czechowski_${promoter_length}bp_${promoterpref}_constitutive

#run ciiider software - variable genes
#$1 is the promoter fasta file location of variable genes only
#$2 is the background gene fasta file
#$3 is jaspar motif file
#$4 is the folder name
#$5 is the CiiiDER output folder name (Czechowski_400bp_promoters_5UTR_constitutive)
../CiiiDER/ciiider.sh ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_variable.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta ../../data/FIMO/motif_data/DAP_seq_motifs.jaspar $file_names Czechowski_${promoter_length}bp_${promoterpref}_variable

#change conda env
conda activate PromoterArchitecturePipeline

#sort CiiiDER output adding TF_family and AGI - constitutive genes
#arg1 is the input location of the geneIDtable
#arg2 is the input location of the motifs_csv output from CiiiDER output
#arg3 is the output location of tsv file containing significantly enriched motifs with AGI and TF family name
python ../data_sorting/./CiiiDER_sort_output.py ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_constitutive/enrichmentoutput_MostSigDeficit.txt ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_constitutive/enrichmentoutput_MostSigDeficit_mapped.tsv

#sort CiiiDER output adding TF_family and AGI - variable genes
#arg1 is the input location of the geneIDtable
#arg2 is the input location of the motifs_csv output from CiiiDER output
#arg3 is the output location of tsv file containing significantly enriched motifs with AGI and TF family name
python ../data_sorting/./CiiiDER_sort_output.py ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_variable/enrichmentoutput_MostSigDeficit.txt ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_variable/enrichmentoutput_MostSigDeficit_mapped.tsv

#create list of enriched TF families 
#arg1 is the input location of enriched TFs with mapped AGI
#arg2 is the output location of enriched TF family coutns text file
#constitutive genes
python ../data_sorting/./enriched_TFfamilies.py ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_constitutive/enrichmentoutput_MostSigDeficit_mapped.tsv ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_constitutive/enriched_TFBS_family_counts.txt
#variable genes
python ../data_sorting/./enriched_TFfamilies.py ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_variable/enrichmentoutput_MostSigDeficit_mapped.tsv ../../data/output/$file_names/CiiiDER/Czechowski_${promoter_length}bp_${promoterpref}_variable/enriched_TFBS_family_counts.txt

#flag genes from the czechowski constitutive/variable/control gene_set which are transcription factors
#arg1 is the promoter extraction output folder name
#arg2 is the gene_categories input file
#arg3 is the input location of the Arabidopsis transcription factor list
#arg4 is the output location of the flagged TF genes
python ../data_sorting/./flag_TF_genes.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/genes/Ath_TF_list.txt ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random_variable_TFs_only.txt

## gene ontology enrichment analysis

#arg1 is the promoter extraction output folder name
#arg2 is the directory location of go term enrichment files
#arg3 is the location of background gene set
#arg4 is the location of NCBI gene list
#For this I downloaded all Arabidopsis protein coding genes from NCBI using the following search (instructions shown here https://github.com/tanghaibao/goatools/blob/master/notebooks/backround_genes_ncbi.ipynb):

#"3702"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]

#I downloaded the tabular file format of all the genes (27562) on 19/08/20

#arg5 is the location of genes of interest
#analysis of top 300 constitutive and top 300 variable genes
python ../data_sorting/./go_term_enrichment.py $file_names ../../data/output/$file_names/genes/gene_ontology ../${promoterpref}_czechowski_allfilteredgenes.txt ../../../../genes/gene_result.txt ../${promoterpref}_czechowski_constitutive_variable_random_300.txt





######Rerun TFBS coverage and diversity analyses and plotting just using TFBSs falling within open chromatin######





## run coverageBed to find TFBS % nucleotide coverage of a promoter
#$1 is promoter bed file
#$2 is the folder name
#$3 is the motifs file location
#$4 is the output file location 
../data_sorting/./TFBS_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}.bed $file_names ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsintersect.bed ../../data/output/$file_names/TFBS_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt


#TFBS coverage of shortened promoters
#arg1 is the promoter extraction output folder name
#arg2 is the name of the output directory
#arg3 is the input location of shortened promoter bed
#arg4 is the input location of bed file for which % coverage is being calculated for
#arg5 is the output location of the % coverage bed file
python ../rolling_window/./coverage_rw.py $file_names TFBS_coverage ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsintersect.bed ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt


###Whole promoter plots###
#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt whole_prom_openchrom/

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TF_diversity_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed whole_prom_openchrom/


###Shortened promoter plots###

#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt ${promoterpref}_${promoter_length}bp_openchrom/


#Create shortened promoter motifs_mapped file in correct format for the TF_diversity_plots.py script
#arg1 is the input location of the shortened promoter file
#arg2 is the input location of the promoters mapped motif bed
#arg3 is the output location of the shortened mapped_motifs bed
#python ../data_sorting/create_motif_mapped_from_intersect.py ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom_motifs_mapped.bed

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TF_diversity_plots_shortenedprom.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_${promoter_length}bp_openchrom/

#need to add KEGG analysis to this too

####
####
####
#### TAU ANALYSIS


#Ran filter_microarray_conditions.ipynb to filter conditions of interest, outputting AtGE_dev_gcRMA.txt.newline.filtered
#filter microarray to include only conditions of interest for TAU tissue specificity
#arg1 is the input location of expression data table (schmid et al 2005 microarray AtGE_dev_gcRMA.txt.newline)
#arg2 is the output location of filtered microarray
#arg3 is the output location of TAU table

#python ../data_sorting/./filter_microarray_conditions.py ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline.filtered ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt




#prepare files for gat analysis TATA enrichment
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Schmid gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the Input location of TATAbox_location bed file (from Eukaryotic promoter database)
#arg7 is the optional replacement name for first variable (eg. non-specific)
#arg8 is the optional replacement name for the second variable (eg. tissue_specific)
#use ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed (stringent TATAs))
#Used the following search parameters for download:
## FindM Genome Assembly : A. thaliana (Feb 2011 TAIR10/araTha1)
##Series : EPDnew, the Arabidopsis Curated Promoter Database
##Sample : TSS from EPDnew rel 004
##Repeat masking: off
##5' border: -50     3' border: 0 
##Search mode: forward
##Selection mode : all matches)

#or ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_renamed.bed (Downloaded TATAbox_location.bed from EPD to data/EPD_promoter_analysis/EPDnew_promoters

#This is less stringent than before, selects TATA boxes +-100bp from most common EPD TSS
#Used the following search parameters for download:
## FindM Genome Assembly : A. thaliana (Feb 2011 TAIR10/araTha1)
##Series : EPDnew, the Arabidopsis Curated Promoter Database
##Sample : TSS from EPDnew rel 004
##Repeat masking: off
##5' border: -100     3' border: 100
##Search mode: forward
##Selection mode : all matches)
python ../data_sorting/./TATA_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Schmid ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed non-specific tissue_specific


#run gat (Genomic association tester) enrichment for TATA boxes using Schmid gene categories
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#use ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed (stringent TATAs) or ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_renamed.bed (Downloaded TATAbox_location.bed from EPD to data/EPD_promoter_analysis/EPDnew_promoters

#This is less stringent than before, selects TATA boxes +-100bp from most common EPD TSS
#Used the following search parameters for download:
## FindM Genome Assembly : A. thaliana (Feb 2011 TAIR10/araTha1)
##Series : EPDnew, the Arabidopsis Curated Promoter Database
##Sample : TSS from EPDnew rel 004
##Repeat masking: off
##5' border: -100     3' border: 100
##Search mode: forward
##Selection mode : all matches)
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#7 is the optional --ignore-segment-tracks flag
../data_sorting/./gat_enrichment_tau.sh ../../data/output/$file_names/TATA/gat_analysis/Schmid_${promoterpref}_workspace.bed ../../data/output/$file_names/TATA/gat_analysis/Schmid_${promoterpref}_non-specific.bed ../../data/output/$file_names/TATA/gat_analysis/Schmid_${promoterpref}_tissue_specific.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_-50to0_renamed.bed ../../data/output/$file_names/TATA/gat_analysis ${promoterpref}_Schmid_TATA --ignore-segment-tracks


##PLOTS TAU
#Whole promoter plots:

#Plot GC content
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters GC_content tsv file
#arg4 is the output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional variable2 name (default is 'variable')
#arg7 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg8 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./GC_content_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/GC_content/${promoterpref}_GC_content.tsv tau/ non-specific tissue_specific Schmid tab10_r


#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional variable2 name (default is 'variable')
#arg7 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg8 is the optional seaborn colour palette for plots, default is None (sns.color_palette())
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}.bp_covered.txt tau/ non-specific tissue_specific Schmid tab10_r

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional variable2 name (default is 'variable')
#arg7 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg8 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./TF_diversity_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed tau/ non-specific tissue_specific Schmid tab10_r

#TATA enrichment plot Czechowski gene categories
#arg1 is the promoter extraction output folder name
#arg2 is the Location of constitutive promoter gat analysis output
#arg3 is the Location of variable promoter gat analysis output
#arg4 is the Output prefix to add to plot file name
#arg5 is the output folder name ending in a forward slash
#arg6 is the optional variable1 name (default is 'constitutive')
#arg7 is the optional variable2 name (default is 'variable')
#arg8 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./TATA_enrichment_plots.py $file_names ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Schmid_TATA_non-specific.out ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Schmid_TATA_tissue_specific.out Schmid_${promoterpref} tau/ non-specific tissue_specific tab10_r

#Open chromatin coverage shoot-root intersect using Potter et al 2018 ATAC-seq negative controls for root/shoot open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of open chromatin bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash - open chromatin tissue
#arg5 is the optional output folder name ending in a forward slash - promoter set name
#arg6 is the optional variable1 name (default is 'constitutive')
#arg7 is the optional variable2 name (default is 'variable')
#arg8 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg9 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt RootShootIntersect/ ${promoterpref}_tau/ non-specific tissue_specific Schmid tab10_r
#Open chromatin coverage root
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}RootOpenChrom.bp_covered.txt Root/ ${promoterpref}_tau/ non-specific tissue_specific Schmid tab10_r
#Open chromatin coverage shoot
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootOpenChrom.bp_covered.txt Shoot/ ${promoterpref}_tau/ non-specific tissue_specific Schmid tab10_r



##PLOT SLIDING WINDOWS TAU

#GC content sliding window
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of GC content tsv file
#arg4 is the input location of eukaryotic promoter database transcription start site bed file
#arg5 is the input location of promoter bed file
#arg6 is the input location of promoter no 5UTR bed file
#arg7 is the output folder name prefix to use
#arg8 is the input location of root chromatin bed file
#arg9 is the input location of shoot chromatin bed file
#arg10 is the input location of rootshootintersect chromatin bed file
#arg11 is the optional replacement colour palette for plots',default = None
#arg13 is the optional author name to add to output file names',default = 'Czechowski
#arg14 is the optional variable1 name (default is 'constitutive')
#arg15 is the optional variable 2 name eg. tissue_specific',default = 'variable'
python ../plotting/rolling_window/./GC_content_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/rolling_window/GC_content_rw/${promoterpref}_GCcontent_rw.tsv ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 GC_content_rw_tau ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed tab10_r Schmid non-specific tissue_specific

#TF diversity sliding window plot
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of the window bed file
#arg4 is the input location of the TF diversity bed file
#arg5 is the input location of eukaryotic promoter database transcription start site bed file
#arg6 is the input location of promoter bed file
#arg7 is the input location of promoter no 5UTR bed file
#arg8 is the output folder name prefix to use
#arg9 is the input location of root chromatin bed file
#arg10 is the input location of shoot chromatin bed file
#arg11 is the input location of rootshootintersect chromatin bed file
#arg12 is the optional replacement colour palette for plots',default = None
#arg13 is the optional author name to add to output file names',default = 'Czechowski
#arg14 is the optional variable1 name (default is 'constitutive')
#arg15 is the optional variable 2 name eg. tissue_specific',default = 'variable'

python ../plotting/rolling_window/./TF_diversity_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/output/$file_names/rolling_window/TF_diversity_rw/${promoterpref}_TF_diversity.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 TF_diversity_rw_tau ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed tab10_r Schmid non-specific tissue_specific

#TFBS coverage sliding window plot

#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of the TFBS base pairs covered file
#arg4 is the input location of eukaryotic promoter database transcription start site bed file
#arg5 is the input location of promoter bed file
#arg6 is the input location of promoter no 5UTR bed file
#arg7 is the output folder name prefix to use
#arg8 is the input location of root chromatin bed file
#arg9 is the input location of shoot chromatin bed file
#arg10 is the input location of rootshootintersect chromatin bed file
#arg11 is the optional replacement colour palette for plots',default = None
#arg12 is the optional author name to add to output file names',default = 'Czechowski
#arg13 is the optional variable1 name (default is 'constitutive')
#arg14 is the optional variable 2 name eg. tissue_specific',default = 'variable'

python ../plotting/rolling_window/./TFBScoverage_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/rolling_window/TFBS_coverage_rw/${promoterpref}_bpcovered_rw.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 TFBS_coverage_rw_tau ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed tab10_r Schmid non-specific tissue_specific

#open chromatin coverage sliding window plot

#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of the GC content TSV file
#arg4 is the input location of eukaryotic promoter database transcription start site bed file
#arg5 is the input location of promoter bed file
#arg6 is the input location of promoter no 5UTR bed file
#arg7 is the output folder name prefix to use
#arg8 is the input location of root chromatin bed file
#arg9 is the input location of shoot chromatin bed file
#arg10 is the input location of rootshootintersect chromatin bed file
#arg11 is the optional replacement colour palette for plots',default = None
#arg12 is the optional author name to add to output file names',default = 'Czechowski
#arg13 is the optional variable1 name (default is 'constitutive')
#arg14 is the optional variable 2 name eg. tissue_specific',default = 'variable'
python ../plotting/rolling_window/./openchromatin_rw_plots_single.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/rolling_window/GC_content_rw/${promoterpref}_GCcontent_rw.tsv ../../data/EPD_promoter_analysis/EPDnew_promoters/At_EPDnew.bed ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/promoters.gff3 OpenChromatin_rw_tau ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed tab10_r Schmid non-specific tissue_specific


#rerun analyses at shorter promoter length


#PLOTTING

#Shortened promoter plots:

#Plot GC content
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters GC_content tsv file
#arg4 is the output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional variable2 name (default is 'variable')
#arg7 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg8 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())

python ../plotting/./GC_content_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/GC_content/${promoterpref}_${promoter_length}bp_GC_content.tsv ${promoterpref}_${promoter_length}bp_tau/ non-specific tissue_specific Schmid tab10_r



#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional variable2 name (default is 'variable')
#arg7 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg8 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bp.bp_covered.txt ${promoterpref}_${promoter_length}bp_tau/ non-specific tissue_specific Schmid tab10_r

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional variable2 name (default is 'variable')
#arg7 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg8 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./TF_diversity_plots_shortenedprom.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp_tau/ non-specific tissue_specific Schmid tab10_r



#Open chromatin coverage shoot-root intersect using Potter et al 2018 ATAC-seq negative controls for root/shoot open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of open chromatin bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash - open chromatin tissue
#arg5 is the optional output folder name ending in a forward slash - promoter set name
#arg6 is the optional variable1 name (default is 'constitutive')
#arg7 is the optional variable2 name (default is 'variable')
#arg8 is the optional replacement name for author in reference to the geneset (default Czechowski)
#arg9 is the optional seaborn colour palette for plots, default is "" (sns.color_palette())
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt RootShootIntersect/ ${promoterpref}_${promoter_length}bp_tau/ non-specific tissue_specific Schmid tab10_r
#Open chromatin coverage root
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpRootOpenChrom.bp_covered.txt Root/ ${promoterpref}_${promoter_length}bp_tau/ non-specific tissue_specific Schmid tab10_r
#Open chromatin coverage shoot
python ../plotting/./OpenChromatin_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootOpenChrom.bp_covered.txt Shoot/ ${promoterpref}_${promoter_length}bp_tau/ non-specific tissue_specific Schmid tab10_r



## gene ontology enrichment analysis

#arg1 is the promoter extraction output folder name
#arg2 is the directory location of go term enrichment files
#arg3 is the location of background gene set
#arg4 is the location of NCBI gene list
#For this I downloaded all Arabidopsis protein coding genes from NCBI using the following search (instructions shown here https://github.com/tanghaibao/goatools/blob/master/notebooks/backround_genes_ncbi.ipynb):

#"3702"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]

#I downloaded the tabular file format of all the genes (27562) on 19/08/20

#arg5 is the location of genes of interest
#arg6 is the optional variable1 name (default is 'constitutive')
#arg7 is the optional replacement name for 2nd variable eg. tissue_specific',default = 'variable'
#arg8 is the optional replacement name for author in reference to the geneset',default = 'Czechowski'

#analysis of top 300 non-specific and top 300 tissue_specific genes
python ../data_sorting/./go_term_enrichment.py $file_names ../../data/output/$file_names/genes/gene_ontology ../${promoterpref}_schmid_allfilteredgenes.txt ../../../../genes/gene_result.txt ../${promoterpref}_schmid_non-specific_tissue_specific_random_300.txt non-specific tissue_specific Schmid


###TAU analysis only in open chromatin###

###Whole promoter plots###
#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional replacement name for 2nd variable eg. tissue_specific
#arg7 is the optional replacement name for author in reference to the geneset
#arg8 is the optional replacement colour palette for plots
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt whole_prom_openchrom_tau/ non-specific tissue_specific Schmid tab10_r

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the optional output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional replacement name for 2nd variable eg. tissue_specific
#arg7 is the optional replacement name for author in reference to the geneset
#arg8 is the optional replacement colour palette for plots
python ../plotting/./TF_diversity_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed whole_prom_openchrom_tau/ non-specific tissue_specific Schmid tab10_r


###Shortened promoter plots###

#Plot TFBS coverage
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters bp_covered txt file
#arg4 is the optional output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional replacement name for 2nd variable eg. tissue_specific
#arg7 is the optional replacement name for author in reference to the geneset
#arg8 is the optional replacement colour palette for plots
python ../plotting/./TFBS_coverage_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt ${promoterpref}_${promoter_length}bp_openchrom_tau/ non-specific tissue_specific Schmid tab10_r




#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the optional output folder name ending in a forward slash
#arg5 is the optional variable1 name (default is 'constitutive')
#arg6 is the optional replacement name for 2nd variable eg. tissue_specific
#arg7 is the optional replacement name for author in reference to the geneset
#arg8 is the optional replacement colour palette for plotss

python ../plotting/./TF_diversity_plots_shortenedprom.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_${promoter_length}bp_openchrom_tau/ non-specific tissue_specific Schmid tab10_r

#Tau TFBS family enrichment

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, whole promoters open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
#arg8 is the optional variable1 name (default is 'constitutive')
#arg9 is the optional variable2 name eg. tissue_specific (default is 'variable')
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Schmid ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_openchromrootshootintersect_TFfamily.bed non-specific tissue_specific

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, whole promoters
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
#arg8 is the optional variable1 name (default is 'constitutive')
#arg9 is the optional variable2 name eg. tissue_specific (default is 'variable')
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Schmid ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_TFfamily.bed non-specific tissue_specific

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, 400bp promoters, openchromatin
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
#arg8 is the optional variable1 name (default is 'constitutive')
#arg9 is the optional variable2 name eg. tissue_specific (default is 'variable')
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed Schmid_400bp ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect_TFfamily.bed non-specific tissue_specific

#prepare files for gat analysis TFs enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, 400bp promoters
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
#arg8 is the optional variable1 name (default is 'constitutive')
#arg9 is the optional variable2 name eg. tissue_specific (default is 'variable')
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed Schmid_400bp ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_TFfamily.bed non-specific tissue_specific

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - whole promoters open chromatin only
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks or --with-segment-tracks flag
../data_sorting/./gat_enrichment_tau.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoterpref}_non-specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoterpref}_tissue_specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_openchromrootshootintersect_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoterpref}_Schmid_TFBS_openchrom --with-segment-tracks

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - whole promoters 
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks or --with-segment-tracks flag
../data_sorting/./gat_enrichment_tau.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoterpref}_non-specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoterpref}_tissue_specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoterpref}_Schmid_TFBS --with-segment-tracks

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - 400bp promoters open chromatin only
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks flag or --with-segment-tracks flag
../data_sorting/./gat_enrichment_tau.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_non-specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoter_length}bp_${promoterpref}_Schmid_TFBS_openchrom --with-segment-tracks

#run gat (Genomic association tester) enrichment for TFs using Czechowski gene categories - 400bp promoters
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks flag or --with-segment-tracks flag
../data_sorting/./gat_enrichment_tau.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_non-specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoter_length}bp_${promoterpref}_Schmid_TFBS --with-segment-tracks

#run precentrimo.sh - 400bp promoters Tau - non-specific
#$1 is the promoter bed file location of non-specific genes only
#$2 is the promoter bed file location of non-specific and tissue_specific genes (background genes)
#$3 is the folder name
#$4 is the genome fasta file location
#$5 is the constitutive promoter fasta file output location
#$6 is the constitutive, variable and random promoters fasta file output location 

../meme_suite/./precentrimo.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_non-specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/Schmid_${promoter_length}bp_${promoterpref}_non-specifictissue_specific_gat.bed $file_names $genome_fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific.fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific_tissue_specific.fasta


#run precentrimo.sh - 400bp promoters Tau - tissue_specific
#$1 is the promoter bed file location of non-specific genes only
#$2 is the promoter bed file location of non-specific and tissue_specific genes (background genes)
#$3 is the folder name
#$4 is the genome fasta file location
#$5 is the constitutive promoter fasta file output location
#$6 is the constitutive, variable and random promoters fasta file output location 

../meme_suite/./precentrimo.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific_gat.bed ../../data/output/$file_names/TFBS_enrichment/Schmid_${promoter_length}bp_${promoterpref}_non-specifictissue_specific_gat.bed $file_names $genome_fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific.fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific_tissue_specific.fasta

#activate memesuite conda environment (install memesuite 5.0.2 using conda install 'meme=5.0.2' 'icu=58.2'. Centrimo is broken on memesuite versions higher than 5.0.2)
#MemeSuite4 env has 
conda activate MemeSuite4
#run centrimo to look for enriched motifs in constitutive promoters (using whole promoters)
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the promoter fasta file location of the background genes
#$3 is meme motif file 
#$4 is the folder name
../meme_suite/./centrimo.sh ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific.fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific_tissue_specific.fasta ../../data/FIMO/motif_data/dap_combined.meme $file_names 

#run ciiider software - non-specific genes
#$1 is the promoter fasta file location of non-specific genes only
#$2 is the background gene fasta file
#$3 is jaspar motif file
#$4 is the folder name
#$5 is the CiiiDER output folder name (Czechowski_400bp_promoters_5UTR_constitutive)
../CiiiDER/ciiider.sh ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific.fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific_tissue_specific.fasta ../../data/FIMO/motif_data/DAP_seq_motifs.jaspar $file_names Schmid_${promoter_length}bp_${promoterpref}_non-specific

#run ciiider software - tissue_specific genes
#$1 is the promoter fasta file location of tissue_specific genes only
#$2 is the background gene fasta file
#$3 is jaspar motif file
#$4 is the folder name
#$5 is the CiiiDER output folder name (Czechowski_400bp_promoters_5UTR_constitutive)
../CiiiDER/ciiider.sh ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific.fasta ../../data/output/$file_names/centrimo/Schmid_${promoter_length}bp_${promoterpref}_non-specific_tissue_specific.fasta ../../data/FIMO/motif_data/DAP_seq_motifs.jaspar $file_names Schmid_${promoter_length}bp_${promoterpref}_tissue_specific

#run ciiider software - non-specific genes open chromatin

#run ciiider software - tissue_specific genes open chromatin
conda activate PromoterArchitecturePipeline

#sort CiiiDER output adding TF_family and AGI - non-specific genes
#arg1 is the input location of the geneIDtable
#arg2 is the input location of the motifs_csv output from CiiiDER output
#arg3 is the output location of tsv file containing significantly enriched motifs with AGI and TF family name
python ../data_sorting/./CiiiDER_sort_output.py ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_non-specific/enrichmentoutput_MostSigDeficit.txt ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_non-specific/enrichmentoutput_MostSigDeficit_mapped.tsv


#sort CiiiDER output adding TF_family and AGI - tissue_specific genes
#arg1 is the input location of the geneIDtable
#arg2 is the input location of the motifs_csv output from CiiiDER output
#arg3 is the output location of tsv file containing significantly enriched motifs with AGI and TF family name
python ../data_sorting/./CiiiDER_sort_output.py ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific/enrichmentoutput_MostSigDeficit.txt ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific/enrichmentoutput_MostSigDeficit_mapped.tsv


#create list of enriched TF families 
#arg1 is the input location of enriched TFs with mapped AGI
#arg2 is the output location of enriched TF family coutns text file
#non-specific genes
python ../data_sorting/./enriched_TFfamilies.py ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_non-specific/enrichmentoutput_MostSigDeficit_mapped.tsv ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_non-specific/enriched_TFBS_family_counts.txt
#tissue_specific genes
python ../data_sorting/./enriched_TFfamilies.py ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific/enrichmentoutput_MostSigDeficit_mapped.tsv ../../data/output/$file_names/CiiiDER/Schmid_${promoter_length}bp_${promoterpref}_tissue_specific/enriched_TFBS_family_counts.txt



#flag genes from the czechowski constitutive/variable/control gene_set which are transcription factors
#arg1 is the promoter extraction output folder name
#arg2 is the gene_categories input file
#arg3 is the input location of the Arabidopsis transcription factor list
#arg4 is the input location of Czechowski et al 2005 microarray all genes data
#arg5 is the output location of the flagged TF genes
python ../data_sorting/./flag_TF_genes_tau.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/genes/Ath_TF_list.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_allfilteredgenes.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random_variable_TFs_only.txt



#generate heatmap of gene expression across each tissue/condition for genes in each gene category

#arg1 is the promoter extraction output folder name
#arg2 is the input location of the coefficient of variation gene categories text file
#arg3 is the input location of the tau tissue/condition specificity gene categories text file
#arg4 is the input location of the log2_transformed microarray expression data in triplicates for each condition/tissue
python ../plotting/./gene_group_expression.py $file_names ../../data/output/non-overlapping_includingbidirectional_all_genes_newannotation/genes/promoters_5UTR_czechowski_constitutive_variable_random.txt ../../data/output/non-overlapping_includingbidirectional_all_genes_newannotation/genes/promoters_5UTR_schmid_non-specific_tissue_specific_random.txt ../../data/genes/AtGeneExpress_CV_2020/AtGE_dev_gcRMA.txt.newline



#Map CV of cognate TFs binding promoters
#arg1 is the promoter extraction output folder name
#arg2 is the input location of the gene categories text file
#arg3 is the input location of the all genes ranking file
#arg4 is the input location of the mapped motifs bed file
#arg5 is the output folder name ending in a forward slash
#arg6 is the optional replacement name for 2nd variable eg. non-specific
#arg7 is the optional replacement name for 2nd variable eg. tissue_specific
#arg8 is the optional replacement name for the author in reference to the genes. Default = Czechowski
#arg9 is the optional replacement colour palette for plot
#whole promoters
python ../plotting/./map_TF2CV.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ${promoterpref}_TF_CV/

#whole promoters using TFBS falling within open chromatin
python ../plotting/./map_TF2CV.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_TF_CV_openchrom/


#400bp promoters
python ../plotting/./map_TF2CV.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp_TF_CV/

#400bp promoters open chromatin
python ../plotting/./map_TF2CV.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_${promoter_length}bp_TF_CV_openchrom/


#Map tau of cognate TFs binding promoters
#arg1 is the promoter extraction output folder name
#arg2 is the input location of the gene categories text file
#arg3 is the input location of the all genes ranking file
#arg4 is the input location of the mapped motifs bed file
#arg5 is the output folder name ending in a forward slash
#arg6 is the optional replacement name for 2nd variable eg. non-specific
#arg7 is the optional replacement name for 2nd variable eg. tissue_specific
#arg8 is the optional replacement name for the author in reference to the genes. Default = Czechowski
#arg9 is the optional replacement colour palette for plot
#whole promoters
python ../plotting/./map_TF2Tau.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ${promoterpref}_TF_Tau/ non-specific tissue_specific Schmid tab10_r

#whole promoters using TFBS falling within open chromatin
python ../plotting/./map_TF2Tau.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_TF_Tau_openchrom/ non-specific tissue_specific Schmid tab10_r


#400bp promoters
python ../plotting/./map_TF2Tau.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp_TF_Tau/ non-specific tissue_specific Schmid tab10_r

#400bp promoters open chromatin
python ../plotting/./map_TF2Tau.py $file_names ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_${promoter_length}bp_TF_Tau_openchrom/ non-specific tissue_specific Schmid tab10_r


#choose top 100 constitutive and top 100 variable TFs based on CVs
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the number of genes in each category to subset
#arg3 is the list of Arabidopsis transcription factors
#arg4 is the input location of all genes ranked by CV produced by choose_genes_cv.py
#arg5 is the output location of the TF gene categories
#arg6 is the output location of all ranked TFs by CV value
python ../data_sorting/./choose_TFs_cv.py $file_names 100 ../../data/genes/Ath_TF_list.txt ../../data/output/${file_names}/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/${file_names}/genes/tfs_czechowski_constitutive_variable_random_100.txt ../../data/output/${file_names}/genes/tfs_czechowski_all_cv.txt 


#choose top 100 non-specific and top 100 tissue_specific TFs based on tau
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the number of genes in each category to subset
#arg3 is the list of Arabidopsis transcription factors
#arg4 is the Input location of all filtered microarray genes ranked by tau (output from choose_genes_tau.py)
#arg5 is the output location of TF gene categories
#arg6 is the input location of Czechowski coefficient of variation TF gene categories
#arg7 is the optional input location of coefficient of variation all genes after filtering to include only genes present in 80% of conditions
python ../data_sorting/./choose_TFs_tau.py $file_names 100 ../../data/genes/Ath_TF_list.txt ../../data/output/${file_names}/genes/${promoterpref}_schmid_allfilteredgenes.txt ../../data/output/${file_names}/genes/tfs_schmid_non-specific_tissue_specific_random_100.txt ../../data/output/$file_names/genes/tfs_czechowski_constitutive_variable_random_100.txt ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt ../../data/output/${file_names}/genes/tfs_schmid_all_tau.txt 


#Do descriptive plots - distribution of CV and tau in the TFs
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the table of interest showing all TFs with their cv value
#arg3 is the table of interest showing all TFs with their tau value
#arg4 is the cv tf gene categories
#arg5 is the tau tf gene categories
python ../plotting/./distribution_tf_ranking.py $file_names ../../data/output/${file_names}/genes/tfs_czechowski_all_cv.txt ../../data/output/${file_names}/genes/tfs_schmid_all_tau.txt ../../data/output/${file_names}/genes/tfs_czechowski_constitutive_variable_random_100.txt ../../data/output/${file_names}/genes/tfs_schmid_non-specific_tissue_specific_random_100.txt


#plot and analyse the TF-targets - are constitutive TFs more likely to bind genes with lower CVs?#

#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the Input location of transcription factor gene categories
#arg3 is the Input location of all ranked genes
#arg4 is the Input location of mapped TF-promoter target file
#arg5 is the dependent variable name [tau or expression_cv]
#arg6 is the output folder name ending with a forward slash
#arg7 is the optional replacement name for 1st variable eg. non-specific',default = 'constitutive'
#arg8 is the optional replacement name for 2nd variable eg. tissue&condition_specific',default = 'variable'
#arg9 is the optional replacement name for author in reference to the geneset',default = 'Czechowski'
#arg10 is the optional replacement colour palette for plots',default = None
#CV ranked genes
python ../plotting/./TF_target_genetype.py $file_names ../../data/output/${file_names}/genes/tfs_czechowski_constitutive_variable_random_100.txt ../../data/output/${file_names}/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/${file_names}/TFBS_TF_class/${promoterpref}_TF_CV/promoter_TF_CV.txt expression_CV TF_target/


#TAU ranked genes
python ../plotting/./TF_target_genetype.py $file_names ../../data/output/${file_names}/genes/tfs_schmid_non-specific_tissue_specific_random_100.txt ../../data/output/${file_names}/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/${file_names}/TFBS_TF_class/${promoterpref}_TF_Tau/promoter_TF_TAU.txt TAU TF_target/ non-specific tissue_specific Schmid tab10_r 


#make combined CV tau GC content plots
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of coefficient of variation gene categories text file
#arg3 is the input location of tau tissue specific gene categories text file
#arg4 is the input location of promoters GC_content tsv file
#arg5 is the optional output folder name ending in a forward slash
#arg7 is the optional replacement colour palette for cv
#arg8 is the optional replacement colour palette for tau

#make combined CV tau GC content plots whole promoter
python ../plotting/./GC_content_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/GC_content/${promoterpref}_GC_content.tsv ${promoterpref}_combined_cv_tau/ deep tab10_r
# shortened promoter combined GC content plot
python ../plotting/./GC_content_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/GC_content/${promoterpref}_${promoter_length}bp_GC_content.tsv ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r

#make combined CV tau % bp open chromatin plots

#Combined cv and tau open chromatin plots whole promoters
#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the input location of coefficient of variation gene categories text file
#arg3 is the input location of tau tissue specific gene categories text file
#arg4 is the input location of bp covered of chromatin in promoters text file
#arg5 is the optional output folder name ending in a forward slash
#arg6 is the optional output folder name ending in a forward slash (name this after promoter set name
#arg7 is the optional replacement colour palette for cv
#arg8 is the optional replacement colour palette for tau

#Open chromatin coverage root-shoot intersect
python ../plotting/./OpenChromatin_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt RootShootIntersect/ ${promoterpref}_combined_cv_tau/ deep tab10_r
#Open chromatin coverage root
python ../plotting/./OpenChromatin_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}RootOpenChrom.bp_covered.txt Root/ ${promoterpref}_combined_cv_tau/ deep tab10_r 
#Open chromatin coverage shoot
python ../plotting/./OpenChromatin_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootOpenChrom.bp_covered.txt Shoot/ ${promoterpref}_combined_cv_tau/ deep tab10_r 


#Combined cv and tau open chromatin plots shortened promoters

#Open chromatin coverage root-shoot intersect
python ../plotting/./OpenChromatin_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt RootShootIntersect/ ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r
#Open chromatin coverage root
python ../plotting/./OpenChromatin_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpRootOpenChrom.bp_covered.txt Root/ ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r 
#Open chromatin coverage shoot
python ../plotting/./OpenChromatin_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootOpenChrom.bp_covered.txt Shoot/ ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r 


#ATA enrichment combined cv and tau plot

#arg1 is the promoter extraction output folder name
#arg2 is the location of the constitutive promoter gat analysis output
#arg3 is the location of the variable promoter gat analysis output
#arg4 is the location of the non-specific promoter gat analysis output
#arg5 is the location of the tissue_specific promoter gat analysis output
#arg6 is the Output prefix to add to plot file name
#arg7 is the optional output_folder_name
#arg8 is the colour palette for the cv gene categories
#arg8 is the colour palette for the tau gene categories

python ../plotting/./TATA_enrichment_plots_combined.py $file_names ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Czechowski_TATA_constitutive.out ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Czechowski_TATA_variable.out ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Schmid_TATA_non-specific.out ../../data/output/$file_names/TATA/gat_analysis/gat_${promoterpref}_Schmid_TATA_tissue_specific.out TATA_enrichment ${promoterpref}_combined_cv_tau/ deep tab10_r


#TF diversity combined plots cv and tau

#arg1 is the Name of folder and filenames for the promoters extracted
#arg2 is the input location of coefficient of variation gene categories text file
#arg3 is the Input location of tau tissue specific gene categories text file
#arg4 is the Input location of promoters mapped motif bed
#arg5 is the Optional output folder name ending in a forward slash
#arg6 is the Optional replacement colour palette for cv
#arg7 is the Optional replacement colour palette for tau

#whole promoter
python ../plotting/./TF_diversity_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ${promoterpref}_combined_cv_tau/ deep tab10_r

#whole promoter open chromatin
python ../plotting/./TF_diversity_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_combined_cv_tau_openchrom/ deep tab10_r

#shortened promoter
python ../plotting/./TF_diversity_plots_shortenedprom_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r

#shortened promoter open chromatin
python ../plotting/./TF_diversity_plots_shortenedprom_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_${promoter_length}bp_combined_cv_tau_openchrom/ deep tab10_r

#TFBS % coverage combined plots cv and tau
#arg1 is the promoter extraction output folder name
#arg2 is the Input location of coefficient of variation gene categories text file
#arg3 is the Input location of tau tissue specific gene categories text file
#arg4 is the Input location of promoters bp_covered txt file
#arg5 is the Optional output folder name ending in a forward slash
#arg6 is the Optional replacement colour palette for cv categories
#arg7 is the Optional replacement colour palette for tau categories

#whole promoter combined cv & tau plot
python ../plotting/./TFBS_coverage_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}.bp_covered.txt ${promoterpref}_combined_cv_tau/ deep tab10_r

#whole promoter open chromatin
python ../plotting/./TFBS_coverage_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}ShootRootIntersectOpenChrom.bp_covered.txt ${promoterpref}_combined_cv_tau_openchrom/ deep tab10_r

#shortened promoter
python ../plotting/./TFBS_coverage_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bp.bp_covered.txt ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r 

#shortened promoter open chromatin
python ../plotting/./TFBS_coverage_plots_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/TFBS_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.bp_covered.txt ${promoterpref}_${promoter_length}bp_combined_cv_tau_openchrom/ deep tab10_r 


#map_TF2CVandTAU combined plots

#arg1 is the Name of folder and filenames for the promoters extracted
#arg2 is the Input location of coefficient of variation gene categories text file
#arg3 is the Input location of tau tissue specific gene categories text file
#arg4 is the Input location of all genes ranking cv categories"
#arg5 is the Input location of all genes ranking tau categories
#arg6 is the Input location of mapped motifs bed file
#arg7 is the Optional replacement colour palette for cv categories
#arg8 is the Optional replacement colour palette for tau categories

#whole promoters
python ../plotting/./map_TF2CV_TF2Tau_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ${promoterpref}_combined_cv_tau/ deep tab10_r

#whole promoters open chromatin
python ../plotting/./map_TF2CV_TF2Tau_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}ShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_combined_cv_tau_openchrom/ deep tab10_r

#shortened promoters
python ../plotting/./map_TF2CV_TF2Tau_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp_combined_cv_tau/ deep tab10_r

#shortened promoters open chromatin
python ../plotting/./map_TF2CV_TF2Tau_combined.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/genes/${promoterpref}_schmid_non-specific_tissue_specific_random.txt ../../data/output/$file_names/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/$file_names/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/$file_names/chromatin_coverage/${promoterpref}_${promoter_length}bpShootRootIntersectOpenChrom.motifsmappedintersect.bed ${promoterpref}_${promoter_length}bp_combined_cv_tau_openchrom/ deep tab10_r


#plot and analyse the TF-targets - are constitutive TFs more likely to bind genes with lower CVs?#

#arg1 is the name of folder and filenames for the promoters extracted
#arg2 is the Input location of transcription factor gene categories
#arg3 is the Input location of all ranked genes
#arg4 is the Input location of mapped TF-promoter target file
#arg5 is the dependent variable name [tau or expression_cv]
#arg6 is the output folder name ending with a forward slash
#arg7 is the optional replacement name for 1st variable eg. non-specific',default = 'constitutive'
#arg8 is the optional replacement name for 2nd variable eg. tissue&condition_specific',default = 'variable'
#arg9 is the optional replacement name for author in reference to the geneset',default = 'Czechowski'
#arg10 is the optional replacement colour palette for plots',default = None
#CV ranked genes
python ../plotting/./TF_target_genetype.py $file_names ../../data/output/${file_names}/genes/tfs_czechowski_constitutive_variable_random_100.txt ../../data/output/${file_names}/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/${file_names}/TFBS_TF_class/${promoterpref}_TF_CV/promoter_TF_CV.txt expression_CV TF_target/


#TAU ranked genes
python ../plotting/./TF_target_genetype.py $file_names ../../data/output/${file_names}/genes/tfs_schmid_non-specific_tissue_specific_random_100.txt ../../data/output/${file_names}/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/${file_names}/TFBS_TF_class/${promoterpref}_TF_Tau/promoter_TF_TAU.txt TAU TF_target/ non-specific tissue_specific Schmid tab10_r 

#plot and analyse the TF-targets combined cv and tau plot
#arg1 is the Name of folder and filenames for the promoters extracted
#arg2 is the Input location of transcription factor gene categories cv
#arg3 is the Input location of transcription factor gene categories tau
#arg4 is the Input location of all ranked genes cv categories
#arg5 is the Input location of all ranked genes tau categories
#arg6 is the Input location of cv mapped TF-promoter target file
#arg7 is the Input location of tau mapped TF-promoter target file
#arg8 is the Output folder name ending in forward slash
#arg9 is the Optional replacement colour palette for cv categories
#arg10 is the Optional replacement colour palette for tau categories
python ../plotting/./TF_target_genetype_combined.py $file_names ../../data/output/${file_names}/genes/tfs_czechowski_constitutive_variable_random_100.txt ../../data/output/${file_names}/genes/tfs_schmid_non-specific_tissue_specific_random_100.txt ../../data/output/${file_names}/genes/filtered80/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/output/${file_names}/genes/tissue_specific/promoters_5UTR_schmid_allfilteredgenes_TAU.txt ../../data/output/${file_names}/TFBS_TF_class/${promoterpref}_TF_CV/promoter_TF_CV.txt ../../data/output/${file_names}/TFBS_TF_class/${promoterpref}_TF_Tau/promoter_TF_TAU.txt TF_target/ deep tab10_r