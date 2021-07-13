#!/bin/sh
eval "$(conda shell.bash hook)"
conda activate PromoterArchitecturePipeline
file_names=300bpsynprom
#promoterpref=TRAMP_random_300bp
promoterpref=TRAMP_random_300bp_manual

#convert fasta file to bed file
#arg1 is "Input location of promoter.fasta file",
#arg2 is Output location of promoter bed file
python ../meme_suite/./fasta2bed.py ../../data/output/$file_names/FIMO/${promoterpref}.fasta ../../data/output/$file_names/FIMO/${promoterpref}.bed


conda activate MemeSuite3

#run FIMO.sh
#$1 is promoter fasta file. 
#$2 is pvalue threshold. 
#$3 is max stored sequences.
#$4 is meme motif file
#$5 is the folder name
#scan for motifs using DAP-seq dataset
../meme_suite/./FIMO.sh ../../data/output/$file_names/FIMO/${promoterpref}.fasta 0.0001 5000000 ../../data/FIMO/motif_data/dap_combined.meme $file_names

conda activate PromoterArchitecturePipeline

#convert fasta file to bed file
#arg1 is "Input location of promoter.fasta file",
#arg2 is Output location of promoter bed file
python ../meme_suite/./fasta2bed.py ../../data/output/$file_names/FIMO/${promoterpref}.fasta ../../data/output/$file_names/FIMO/${promoterpref}.bed

#$1 is the fimo.tsv file location
#arg1 is Input location of FIMO.tsv file
#arg2 is Input location of promoter bed file
#arg3 is Output location of motifs bed file
#arg4 is q_value threshold for filtering
#arg5 is option to prevent shortening of sequence name up to the first colon character
python ../data_sorting/./FIMO_filter.py ../../data/output/$file_names/FIMO/output/${promoterpref}_DAPseq/fimo.tsv ../../data/output/$file_names/FIMO/${promoterpref}.bed ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed 1 --prevent_shorten_sequence_name

python ../meme_suite/./map_motif_ids.py  ../../data/output/$file_names/FIMO/${promoterpref}_motifs.bed ../../data/FIMO/motif_data/motif_map_IDs.txt ../../data/output/$file_names/FIMO/${promoterpref}_motifs_mapped.bed 