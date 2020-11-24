import argparse
import os

from pybedtools import BedTool


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="coverage_rw")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "output_directory", type=str, help="Name of output directory"
    )
    parser.add_argument(
        "window_bed",
        type=str,
        help="Input location of rolling window bed file",
    )
    parser.add_argument(
        "input_coverage_bed",
        type=str,
        help="Input location of bed file for which % coverage is being calculated for",
    )
    parser.add_argument(
        "output_coverage_bed",
        type=str,
        help="Output location of rolling window % coverage bed file",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def coverage_bed(window_bed, input_coverage_bed, output_coverage_bed):
    """function to calculate TFBS %coverage on sliding windows"""
    windows = BedTool(window_bed)
    motifs = BedTool(input_coverage_bed)
    # calculate TFBS coverage and save the output file
    windows.coverage(motifs, output=output_coverage_bed)


def main(args):
    # parse arguments
    args = parse_args(args)
    # make directory for the output files to be exported to
    # dirName = f'{args.directory_path}/data/output/{args.file_names}'
    dirName = f"../../data/output/{args.file_names}/rolling_window/{args.output_directory}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    coverage_bed(
        args.window_bed, args.input_coverage_bed, args.output_coverage_bed
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
