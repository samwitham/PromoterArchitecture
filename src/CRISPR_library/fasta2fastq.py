# %%
"""
Convert FASTA to FASTQ file with a static
Usage:
$ ./fasta_to_fastq NAME.fasta NAME.fastq
"""

import os
import glob
from Bio import SeqIO
from pathlib import Path


# %%
input_directory = "../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene"
output_directory = "../../data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/Sam_Witham_EI_SW_ENQ-5142_A_01_Additional_Barcode_Analysis/Cutadapt/Trim/Bygene/fastq"


# %%
def fasta_2_fastq(fasta):
    """ Function to read and convert fasta file to fastq file using a fixed quality score of 40."""
    filename = os.path.basename(fasta)
    removed_extension = os.path.splitext(filename)[0]
    # find parent directory to the one the fasta files are in
    path = Path(
        fasta
    ).parent
    fastq_path = f'{path}/fastq'
    # make fastq
    with open(f"{path}/{removed_extension}.fasta", "r", encoding="utf-8") as fasta, open(f'{fastq_path}/{removed_extension}.fastq', "w", encoding="utf-8") as fastq:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [40] * len(record)
            SeqIO.write(record, fastq, "fastq")


# %%
# find all fasta files in the folder
fasta_filenames = glob.glob(
    f"{input_directory}/*.fasta", recursive=False
)


# run the fasta_2_fastq function across all fasta file in to_be_sorted folder
list(map(fasta_2_fastq, fasta_filenames))

# %%
