#!/bin/sh
#requires activating "conda activate pacbio"
conda activate pacbio
# bwa-mem2 index ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/References/genes.fa
# #Align one read to reference promoter
# bwa-mem2 ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/References/genes.fa ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/bc1017_ARF9_SW346_SW442.trim.fasta > ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/bc1017_ARF9_SW346_SW442.trim.aligned.sam
# CRISPRessoPooled -r1 ../../data/CRISPR_library/pacbio/CCS_HiFi/Data_Package_Batch_15_03_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Batch_15_03_2022-Demultiplex_Barcodes_CCS_Analysis/pb_demux_ccs_r64036e_20220308_172304-2_B01/outputs/demultiplex.bc1017_BAK8B_OA--bc1017_BAK8B_OA.hifi_reads.fastq.gz -f ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/References/genes.tsv --name library1

#convert trimmed demultiplexed fasta files into fastq files with dummy quality scores
# seqtk seq -F '#' ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/bc1017_ARF9_SW346_SW442.trim.fasta > ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/fastq/bc1017_ARF9_SW346_SW442.trim.fastq

#on multiple fasta files:
for filename in ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/*.fasta
do
    seqtk seq -F '#' "../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/${filename}.fasta" > "../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/fastq/${filename}.fastq"
    
done

CRISPRessoPooled -r1 ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/fastq/bc1017_ARF9_SW346_SW442.trim.fastq -f ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/References/genes.tsv --name library1

#use suppress report as  receive an error image size too large
CRISPResso --fastq_r1 ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/fastq/bc1017_ARF9_SW346_SW442.trim.fastq --amplicon_seq TTTGTTTGGCAACGATCAACATATATCATGACTCAAATAAATTTTGATGTTTCTAATTAGTAAAATGTGGAGTTGGAAAAAAGTTCTTGCAAAAGGTGAAAAGAACAAAAATCAGAGATGGAAACACTTCCCACATGCTTATGTTAGTCCCACCATTTAGTGGTAGAAATGTAATTTATCCAAAGAAGGAGGTCAAAGCTGCTTCGTCTCTAATCCAAGTCTAGCCAGTATTGTTGCCTGCATAAGAGCAAACGGTAAAAGTGACGGACTCTTAACGATGACGGGTGTTACTGTTTAATGTCTCTTGATCGACTTTTGTTACCGTAACGACAACGGCACGCGCCCTCACGCTCCGTACCAGCCTTAAAAGCCAACACATTTATGATTTATCCAGTAACTAATACTAATTCCTTCATCCACTAACGCAATTGCATTTAATTGTTTAGTACGTATAGAAAAAAAGTAAAAAGCTATTTATAATTATAACTTACAAAGTCTCTTTCAATTAGTACAGTGTAAAACAACCAATTATTTTCTATATTATACAGTATAAAGGTGTTTATCATAGACCCATTATTTATTTGTTTTTGACTAAAATATCAAAGTGTATTAAAACATAAACCAAATAAATTTGAAATTAAAGTTTTTTTTTTTTTTGTAAAATTAAAACTGATTTTTGGTTTAAATAGAGAAGCAAAGCTACGCTTTTTTTTTTGCTAGCTCTTTTTTTATTTAATCAAACAAAAATTGAGGAAAACAATTATTGTTGTGCTAAGCATTTCTATCATATATTTTTTTATAAACAATACATTTAATTTCAAAATAACTTTCTACTTTCTTCAAATCAAGCTATTTCCTAGAAATTTGGTACTGTCTGGATAACAACATGTCTTAGAAAGTATAGTAAGGAATTTACAATTGGTTATACAATTAATGTATAGCGTATTATTATATACGAAAAGGTGTCCAATTTTAGACATACATATAAAATCAATTTATTAAAAATAAAAGTATTGAAAAAAGAGAGAGTGAAACGGAGATAAAAGCAAACGGAGGTGTTAACTTAAAGACGGGAGAGAAGCGAGTGAGAGATTCTTTCTTGTAACGGCGTTGTTTGACGATTATAGTGACGAAGAGAGACAAAAAAAAAGAGTTTAGTCTTCTTTGTTTGGGAGACGAAGACGATAGAAGCTGGTTTTACGGGAGAGAAATTAGAGAGGTTACAGTTACAGAGCAGGAAGGATTGCGTTGCGTTACGATTCAGATCTCTCTCGATCGTTACTTACATGCATTAGAGAACAAGAAAACGTGGAGACGGAATCCAAGGGACCTTCTTCTTCAATTTAATTCCTCACTTTATTACGCACCCTTCTTTTATTTTCGACTGTTTCAGCTTCCGGTTCT --suppress_report

#enter fastq directory
cd ../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/fastq

#batch run ARF9

CRISPRessoBatch --batch_settings ../../../../References/ARF9.batch --amplicon_seq TTTGTTTGGCAACGATCAACATATATCATGACTCAAATAAATTTTGATGTTTCTAATTAGTAAAATGTGGAGTTGGAAAAAAGTTCTTGCAAAAGGTGAAAAGAACAAAAATCAGAGATGGAAACACTTCCCACATGCTTATGTTAGTCCCACCATTTAGTGGTAGAAATGTAATTTATCCAAAGAAGGAGGTCAAAGCTGCTTCGTCTCTAATCCAAGTCTAGCCAGTATTGTTGCCTGCATAAGAGCAAACGGTAAAAGTGACGGACTCTTAACGATGACGGGTGTTACTGTTTAATGTCTCTTGATCGACTTTTGTTACCGTAACGACAACGGCACGCGCCCTCACGCTCCGTACCAGCCTTAAAAGCCAACACATTTATGATTTATCCAGTAACTAATACTAATTCCTTCATCCACTAACGCAATTGCATTTAATTGTTTAGTACGTATAGAAAAAAAGTAAAAAGCTATTTATAATTATAACTTACAAAGTCTCTTTCAATTAGTACAGTGTAAAACAACCAATTATTTTCTATATTATACAGTATAAAGGTGTTTATCATAGACCCATTATTTATTTGTTTTTGACTAAAATATCAAAGTGTATTAAAACATAAACCAAATAAATTTGAAATTAAAGTTTTTTTTTTTTTTGTAAAATTAAAACTGATTTTTGGTTTAAATAGAGAAGCAAAGCTACGCTTTTTTTTTTGCTAGCTCTTTTTTTATTTAATCAAACAAAAATTGAGGAAAACAATTATTGTTGTGCTAAGCATTTCTATCATATATTTTTTTATAAACAATACATTTAATTTCAAAATAACTTTCTACTTTCTTCAAATCAAGCTATTTCCTAGAAATTTGGTACTGTCTGGATAACAACATGTCTTAGAAAGTATAGTAAGGAATTTACAATTGGTTATACAATTAATGTATAGCGTATTATTATATACGAAAAGGTGTCCAATTTTAGACATACATATAAAATCAATTTATTAAAAATAAAAGTATTGAAAAAAGAGAGAGTGAAACGGAGATAAAAGCAAACGGAGGTGTTAACTTAAAGACGGGAGAGAAGCGAGTGAGAGATTCTTTCTTGTAACGGCGTTGTTTGACGATTATAGTGACGAAGAGAGACAAAAAAAAAGAGTTTAGTCTTCTTTGTTTGGGAGACGAAGACGATAGAAGCTGGTTTTACGGGAGAGAAATTAGAGAGGTTACAGTTACAGAGCAGGAAGGATTGCGTTGCGTTACGATTCAGATCTCTCTCGATCGTTACTTACATGCATTAGAGAACAAGAAAACGTGGAGACGGAATCCAAGGGACCTTCTTCTTCAATTTAATTCCTCACTTTATTACGCACCCTTCTTTTATTTTCGACTGTTTCAGCTTCCGGTTCT -p 7 --suppress_report -bo ../../../../Variant_call/ARF9 --skip_failed