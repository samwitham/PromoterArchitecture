import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="TATA_enrichment_plots")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "gat_TATA_constitutive_output",
        type=str,
        help="Location of constitutive promoter gat analysis output",
    )
    parser.add_argument(
        "gat_TATA_variable_output",
        type=str,
        help="Location of variable promoter gat analysis output",
    )
    parser.add_argument(
        "output_prefix",
        type=str,
        help="Output prefix to add to plot file name",
    )
    parser.add_argument(
        "output_folder_name",
        type=str,
        help="Optional output folder name ending in a forward slash",
        default="",
        nargs="?",
    )
    parser.add_argument(
        "variable1_name",
        type=str,
        help="Optional replacement name for 2nd variable eg. non-specific",
        default="constitutive",
        nargs="?",
    )
    parser.add_argument(
        "variable2_name",
        type=str,
        help="Optional replacement name for 2nd variable eg. tissue_specific",
        default="variable",
        nargs="?",
    )
    parser.add_argument(
        "palette",
        type=str,
        help="Optional replacement colour palette for plots",
        default=None,
        nargs="?",
    )

    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def create_plot(
    file_names,
    output_folder_name,
    output_prefix,
    gat_TATA_constitutive_output,
    gat_TATA_variable_output,
    palette,
    variable1_name,
    variable2_name,
):
    """import and process the raw outputs after running gat (Genomic association tester). Then create barplot of constitutive and variable gene TATA enrichment"""
    # import gat output files as dfs
    constitutive = pd.read_table(
        gat_TATA_constitutive_output, sep="\t", header=0
    )
    variable = pd.read_table(gat_TATA_variable_output, sep="\t", header=0)
    # merge dfs
    merged = pd.concat([constitutive, variable], ignore_index=True)
    # set style to ticks
    sns.set(style="ticks", color_codes=True)

    # set colour palette
    colours = sns.color_palette(palette)

    # bar chart, 95% confidence intervals
    plot = sns.barplot(
        x="annotation",
        y="l2fold",
        data=merged,
        order=[variable1_name, variable2_name],
        palette=colours,
    )
    plot.axhline(0, color="black")
    plt.xlabel("Gene type")
    plt.ylabel("Log2-fold enrichment over background").get_figure().savefig(
        f"../../data/output/{file_names}/TATA/{output_folder_name}plots/{output_prefix}_log2fold.pdf",
        format="pdf",
    )


def main(args):
    # parse arguments
    args = parse_args(args)

    # make directory for the plots to be exported to
    dirName = (
        f"../../data/output/{args.file_names}/TATA/{args.output_folder_name}"
    )
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/TATA/{args.output_folder_name}/plots"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # Create barplot
    create_plot(
        args.file_names,
        args.output_folder_name,
        args.output_prefix,
        args.gat_TATA_constitutive_output,
        args.gat_TATA_variable_output,
        args.palette,
        args.variable1_name,
        args.variable2_name,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
