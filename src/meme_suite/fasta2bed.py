import argparse
import os

import pandas as pd
from pybedtools import BedTool
from pyfaidx import Fasta


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="Fasta2Bed")
    parser.add_argument(
        "promoter_fasta",
        type=str,
        help="Input location of promoter.fasta file",
    )
    parser.add_argument(
        "promoter_bed", type=str, help="Output location of promoter bed file"
    )

    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def fasta_chromsizes(genome, output_file):
    """extracts chromosome sizes in a form compatible with BedTools.flank (as a dictionary)"""
    # output = open(output_file, 'w')
    the_genome = Fasta(
        genome
    )  # using pyfaidx Fasta, parse the genome file as a fasta

    chromsizes = {}  # make dictionary called chromsizes
    for key in the_genome.keys():

        # chromsizes[key] = f'({0}, {len(the_genome[key])})' #add the chromosome name and length to dictionary
        chromsizes[
            f"{key}"
        ] = f"({len(the_genome[key])})"  # add the chromosome name and length to dictionary
        # chromsizes[key] = f'({len(the_genome[key])})' #add the chromosome name and length to dictionary
    # output_file.write(chrom)
    chromsizes_string = ""

    for k, v in chromsizes.items():
        chromsizes_string = chromsizes_string + f"{k}\t{v}\n"

    with open(output_file, "w") as output:
        output.write(chromsizes_string.replace("(", "").replace(")", ""))


def chromsize2bed(chromsize, bed_file):
    chrom_df = pd.read_table(chromsize, sep="\t", header=None)
    cols = ["chr", "stop"]
    chrom_df.columns = cols
    chrom_df["start"] = 1
    chrom_df = chrom_df[["chr", "start", "stop"]]
    # add extra columns so compatible with FIMO_filter.py
    chrom_df["gene"] = chrom_df.chr
    chrom_df["dot"] = "."
    chrom_df["strand"] = "+"
    chrom_df["source"] = "manual"
    chrom_df["type"] = "promoter"
    chrom_df["dot2"] = "."
    chrom_df["details"] = "none"
    sorted_proms = chrom_df.sort_values(["chr", "start"])
    BedTool.from_dataframe(sorted_proms).saveas(bed_file)


def main(args):
    # parse arguments
    args = parse_args(args)
    # create temporary file
    temp = "../../data/FIMO/Fasta2Bed.tmp"

    fasta_chromsizes(args.promoter_fasta, temp)
    chromsize2bed(temp, args.promoter_bed)
    os.remove(temp)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
