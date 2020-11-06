#!/bin/sh
#source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate # need to change this to hardlink
eval "$(conda shell.bash hook)"
conda activate PromoterArchitecturePipeline
#directory_path=ei/workarea/group-eg/project_PromoterArchitecturePipeline
directory_path=/home/witham/Documents/pipeline_new/PromoterArchitecture #remove this, better to use relative path
file_names=non-overlapping_includingbidirectional_all_genes_newannotation_3KB

##extract promoters from a genome file
#arg1 is the location of the base directory
#arg2 is the name of folder and filenames for the promoters extracted
#arg3 is the max length of promoters to be extracted upstream of the Araport TSS
#arg4 is an option to exclude potentially bidirectional promoters with upstream TSS >2000bp from TSS in opposite direction
#arg5 is an option to reduce the size of promoters if needed until they are not overlapping other genes.
#arg6 is an option to extend promoters to the first start codon.

#python ../data_sorting/extract_promoter.py $directory_path $file_names --remove_bidirectional --prevent_overlapping_genes --fiveUTR
python ../data_sorting/extract_promoter.py $directory_path $file_names 3000 --prevent_overlapping_genes --fiveUTR
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
conda activate MemeSuite3
#identify the output filename created by preFIMO.sh
promoterbase=${promoter_gff_file##*/}
promoterpref=${promoterbase%.*}
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
python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_DAPseq/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_q1.bed 1

python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_DAPseq/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_q0_05.bed 0.05

#filter plantpan motifs
python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_plantpan/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_plantpan.bed 0.05

#If using DAP seq cistrome motifs without AGI names, map motifs
#arg1 is the input location of the motif bed file
#arg2 is the input location of the geneIDtable file derived from the .json (see .ipynb notebook)
#arg3 is the output location of the motif bed file
python ../meme_suite/./map_motif_ids.py  ../../data/output/$file_names/FIMO/${promoterpref}_motifs_q1.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped_q1.bed 

python ../meme_suite/./map_motif_ids.py  ../../data/output/$file_names/FIMO/${promoterpref}_motifs_q0_05.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped_q0_05.bed 

#map plantpan motifs
#python ../meme_suite/./map_motif_ids.py ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs__plantpan_mapped.bed




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
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed 100 50 ../../data/output/$file_names/overlapping_promoters_5UTR.bed
#convert promoters.gff3 into bed file
#convert gff file to bed file
gff2bed < ../../data/output/$file_names/promoters.gff3 > ../../data/output/$file_names/FIMO/promoters.bed
#remove 'gene:' from column 4
awk 'BEGIN {FS = OFS = "\t"} { gsub(/gene:/,"", $4); print}' ../../data/output/$file_names/FIMO/promoters.bed > ../../data/output/$file_names/FIMO/promoters.bed.tmp && mv ../../data/output/$file_names/FIMO/promoters.bed.tmp ../../data/output/$file_names/FIMO/promoters.bed
#remove lines not beginning with numbers (ie. remove mt and pt chromosomes)
grep '^[0-9]' ../../data/output/$file_names/FIMO/promoters.bed > ../../data/output/$file_names/FIMO/promoters.bed.tmp && mv ../../data/output/$file_names/FIMO/promoters.bed.tmp ../../data/output/$file_names/FIMO/promoters.bed

#create sliding windows starting from Araport TSS instead (use only promoters not 5'UTR)
#arg1 is the promoter extraction output folder name
#arg2 is the promoter bed file
#arg3 is the output location of rolling window bed file
#arg4 is the size of the rolling window in bp
#arg5 is the size of the window offset in bp
#arg6 is the output location of the overlapping promoters bed file
python ../rolling_window/./rolling_window.py $file_names ../../data/output/$file_names/FIMO/promoters.bed ../../data/output/$file_names/rolling_window/promoters_windows.bed 100 50 ../../data/output/$file_names/overlapping_promoters.bed

#Create Czechowski et al 2005 ranked cv dataset gene categories filtering out any genes not in the extracted promoters from this pipeline. The create subsets of N constitutive, variable or control genes
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
python ../data_sorting/./choose_genes_cv.py $file_names ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/genes/AtGE_dev_gcRMA__all_probes__CV.tsv ../../data/genes/RNA_CVs.csv 100 ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/${file_names}/genes/${promoterpref}_mergner_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_filtered_contain_motifs.bed ../../data/output/$file_names/genes/${promoterpref}_czechowski_allfilteredgenes.txt ../../data/output/$file_names/genes/${promoterpref}_mergner_allfilteredgenes.txt

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

#Root open chromatin coverage sliding window promoters_5'UTR
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_root_bpcovered_rw.bed

#Shoot open chromatin coverage sliding window promoters_5'UTR
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_shoot_bpcovered_rw.bed

#Root/Shoot intersect coverage sliding window promoters_5'UTR
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw ../../data/output/$file_names/rolling_window/${promoterpref}_windows.bed ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw/${promoterpref}_rootshootintersect_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters only
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/ATAC-seq/potter2018/Roots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_promoterno5UTR/promoters_root_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters only
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/ATAC-seq/potter2018/Shoots_NaOH_peaks_all.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_promoterno5UTR/promoters_shoot_bpcovered_rw.bed

#Root open chromatin coverage sliding window promoters only
python ../rolling_window/./coverage_rw.py $file_names OpenChromatin_rw_promoterno5UTR ../../data/output/$file_names/rolling_window/promoters_windows.bed ../../data/ATAC-seq/potter2018/intersectRootsShoots_PeaksInBoth.bed ../../data/output/$file_names/rolling_window/OpenChromatin_rw_promoterno5UTR/promoters_rootshootintersect_bpcovered_rw.bed



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

#prepare files for gat analysis TATA enrichment
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the Input location of TATAbox_location bed file (from Eukaryotic promoter database)
python ../data_sorting/./TATA_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Czechowski ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_renamed.bed

#run gat (Genomic association tester) enrichment for TATA boxes using Czechowski gene categories
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#7 is the optional --ignore-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TATA/gat_analysis/Czechowski_${promoterpref}_workspace.bed ../../data/output/$file_names/TATA/gat_analysis/Czechowski_${promoterpref}_constitutive.bed ../../data/output/$file_names/TATA/gat_analysis/Czechowski_${promoterpref}_variable.bed ../../data/EPD_promoter_analysis/EPDnew_promoters/TATAbox_location_renamed.bed ../../data/output/$file_names/TATA/gat_analysis ${promoterpref}_Czechowski_TATA --ignore-segment-tracks



##PLOTS
#Whole promoter plots:

#Plot GC content
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters GC_content tsv file
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
promoter_length=400
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
python ../rolling_window/./GC_content_rw.py $file_names ../../data/output/$file_names/GC_content/${promoterpref}_${promoter_length}bp_GC_content.tsv ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed $genome_fasta ../../data/output/$file_names/${promoterpref}_${promoter_length}bp.fasta

#%coverage of open chromatin
#$1 is promoter bed file
#$2 is the folder name
../data_sorting/./OpenChromatin_coverage.sh ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed $file_names

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
python ../data_sorting/create_motif_mapped_from_intersect.py ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed

#Plot TF and TF family diversity, and do PCA and Kmeans clustering
#arg1 is the promoter extraction output folder name
#arg2 is the input location of Czechowski gene categories text file
#arg3 is the input location of promoters mapped motif bed
#arg4 is the optional output folder name ending in a forward slash
python ../plotting/./TF_diversity_plots.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ${promoterpref}_${promoter_length}bp/

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


#prepare files for gat analysis TATA enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, whole promoters open chromatin
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Czechowski ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped_openchromrootshootintersect.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_openchromrootshootintersect_TFfamily.bed

#prepare files for gat analysis TATA enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, whole promoters
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}.bed Czechowski ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_TFfamily.bed

#prepare files for gat analysis TATA enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, 400bp promoters, openchromatin
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed Czechowski_400bp ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect_TFfamily.bed

#prepare files for gat analysis TATA enrichment - TFBSs only in open chromatin (root-shoot intersect) - just looking at constitutive genes, 400bp promoters
#arg1 is the promoter extraction output folder name
#arg2 is the promoter prefix name
#arg3 is the Input location of the Czechowski gene categories text file
#arg4 is the Input location of promoters bed file
#arg5 is the Gene category prefix (eg. Czechowski)
#arg6 is the input location of the mapped_motifs bed file
#arg7 is the output location of the mapped motifs bed file with TF family in column 4
python ../data_sorting/./prepare_TFBS_enrichment.py $file_names ${promoterpref} ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp.bed Czechowski_400bp ../../data/output/$file_names/FIMO/${promoterpref}_${promoter_length}bp_motifs_mapped.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_TFfamily.bed



#run gat (Genomic association tester) enrichment for TATA boxes using Czechowski gene categories - whole promoters open chromatin only
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_openchromrootshootintersect_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoterpref}_Czechowski_TFBS_openchrom --with-segment-tracks

#run gat (Genomic association tester) enrichment for TATA boxes using Czechowski gene categories - whole promoters 
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_motifs_mapped_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoterpref}_Czechowski_TFBS --with-segment-tracks

#run gat (Genomic association tester) enrichment for TATA boxes using Czechowski gene categories - 400bp promoters open chromatin only
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks flag or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_openchromrootshootintersect_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoter_length}bp_${promoterpref}_Czechowski_TFBS_openchrom --with-segment-tracks

#run gat (Genomic association tester) enrichment for TATA boxes using Czechowski gene categories - 400bp promoters
#$1 is the workspace file containing all promoters of interest
#$2 is the constitutive promoters file for annotations
#$3 is the variable promoters file for annotations
#$4 is the segments bed file
#$5 is the output folder location
#$6 is the prefix to add the the output files (recommend using segment name + gene categories author as prefix)
#$7 is the optional --ignore-segment-tracks flag or --with-segment-tracks flag
../data_sorting/./gat_enrichment.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_workspace.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_constitutive_gat.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_variable_gat.bed ../../data/output/$file_names/TFBS_enrichment/${promoterpref}_${promoter_length}bp_motifs_mapped_TFfamily.bed ../../data/output/$file_names/TFBS_enrichment/gat_analysis ${promoter_length}bp_${promoterpref}_Czechowski_TFBS --with-segment-tracks

#run precentrimo.sh
#$1 is the promoter bed file location of constitutive genes only
#$2 is the promoter bed file location of constitutive and variable genes (background genes)
#$3 is the folder name
#$4 is the genome fasta file location
#$5 is the constitutive promoter fasta file output location
#$6 is the constitutive, variable and random promoters fasta file output location
../meme_suite/./precentrimo.sh ../../data/output/$file_names/TFBS_enrichment/gat_analysis/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.bed ../../data/output/$file_names/TFBS_enrichment/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.bed $file_names $genome_fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta

#activate memesuite conda environment (install memesuite 5.0.2 using conda install 'meme=5.0.2' 'icu=58.2'. Centrimo is broken on memesuite versions higher than 5.0.2)
#MemeSuite4 env has 
conda activate MemeSuite4
#run centrimo to look for enriched motifs in constitutive promoters (using whole promoters)
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the promoter fasta file location of constitutive, variable and random genes only
#$3 is meme motif file 
#$4 is the folder name
../meme_suite/./centrimo.sh ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_allgenetypes.fasta ../../data/FIMO/motif_data/dap_combined.meme $file_names

#run ciiider software
#$1 is the promoter fasta file location of constitutive genes only
#$2 is the promoter fasta file location of constitutive, variable and random genes only
#$3 is jaspar motif file
#$4 is the folder name
#$5 is the CiiiDER output folder name (Czechowski_400bp_promoters_5UTR_constitutive)
../CiiiDER/ciiider.sh ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutive.fasta ../../data/output/$file_names/centrimo/Czechowski_${promoter_length}bp_${promoterpref}_constitutivevariable.fasta ../../data/FIMO/motif_data/DAP_seq_motifs.jaspar $file_names Czechowski_400bp_promoters_5UTR_constitutive


#flag genes from the czechowski constitutive/variable/control gene_set which are transcription factors
#arg1 is the promoter extraction output folder name
#arg2 is the gene_categories input file
#arg3 is the input location of the Arabidopsis transcription factor list
#arg4 is the output location of the flagged TF genes
python ../data_sorting/./flag_TF_genes.py $file_names ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random.txt ../../data/genes/Ath_TF_list.txt ../../data/output/$file_names/genes/${promoterpref}_czechowski_constitutive_variable_random_variable_TFs_only.txt