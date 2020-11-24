import argparse
import os

from Bio import SeqIO
from Bio.SeqUtils import GC
from pybedtools import BedTool


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="GC_content_rw")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "GC_content_tsv",
        type=str,
        help="Output location of rolling window GC content tsv",
    )
    parser.add_argument(
        "window_bed",
        type=str,
        help="Input location of rolling window bed file",
    )
    parser.add_argument(
        "genome_fasta", type=str, help="Input location of genome fasta file"
    )
    parser.add_argument(
        "window_fasta",
        type=str,
        help="Output location of rolling window fasta file",
    )
    parser.add_argument(
        "outputfolder",
        type=str,
        help="Output folder name",
        default="",
        nargs="?",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def make_window_fasta(window_bed, genome_fasta, window_fasta):
    """Make a fasta file containing all of the promoter windows"""
    # make bedtools object of window_bed
    # windows = BedTool(window_bed)
    # create fasta file of promoters from genome fasta file and from the promoters_renamedChr.bed file
    BedTool().sequence(
        bed=window_bed, fi=genome_fasta, fo=window_fasta, name=True
    )


def GC_content(window_fasta, output_file):
    """function to calculate GC content % of sequences in a fasta file and output a tsv file"""

    fasta_sequences = SeqIO.parse(open(window_fasta), "fasta")
    with open(output_file, "w") as fout:

        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

            fout.write(f"{name}\t{GC(sequence)}\n")


def main(args):
    # parse arguments
    args = parse_args(args)
    # make directory for the output files to be exported to
    dirName = f"../../data/output/{args.file_names}/rolling_window/{args.outputfolder}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    make_window_fasta(args.window_bed, args.genome_fasta, args.window_fasta)
    GC_content(args.window_fasta, args.GC_content_tsv)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
