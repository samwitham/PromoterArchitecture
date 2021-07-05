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
        "gat_TATA_nonspecific_output",
        type=str,
        help="Location of non-specific promoter gat analysis output",
    )
    parser.add_argument(
        "gat_TATA_tissue_specific_output",
        type=str,
        help="Location of tissue_specific promoter gat analysis output",
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
        "palette_cv",
        type=str,
        help="Optional replacement colour palette for cv",
        default=None,
        nargs="?",
    )
    parser.add_argument(
        "palette_tau",
        type=str,
        help="Optional replacement colour palette for tau",
        default=None,
        nargs="?",
    )

    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def create_plot(
    gat_TATA_constitutive_output,
    gat_TATA_variable_output,
    gat_TATA_nonspecific_output,
    gat_TATA_tissuespecific_output,
    output_prefix,
    palette_cv,
    palette_tau,
    output_folder_name,
    file_names,
):
    """import and process the raw outputs after running gat (Genomic association tester). Then create barplot of constitutive and variable gene TATA enrichment"""
    # import gat output files as dfs
    def merge_gatfiles(gat_output1, gat_output2, palette):
        df1 = pd.read_table(gat_output1, sep="\t", header=0)
        df2 = pd.read_table(gat_output2, sep="\t", header=0)

        # merge dfs
        merged = pd.concat([df1, df2], ignore_index=True)
        # set style to ticks
        sns.set(style="ticks", color_codes=True)

        # set colour palette
        colours = sns.color_palette(palette)

        return merged, colours

    merged_cv, colours_cv = merge_gatfiles(
        gat_TATA_constitutive_output, gat_TATA_variable_output, palette_cv
    )
    merged_tau, colours_tau = merge_gatfiles(
        gat_TATA_nonspecific_output,
        gat_TATA_tissuespecific_output,
        palette_tau,
    )
    # make subplots
    f, (ax1, ax2) = plt.subplots(
        1,
        2,
    )
    f.set_size_inches(11, 5)

    # bar chart, 95% confidence intervals
    sns.barplot(
        x="annotation",
        y="l2fold",
        data=merged_cv,
        order=["constitutive", "variable"],
        palette=colours_cv,
        ax=ax1,
    )
    sns.barplot(
        x="annotation",
        y="l2fold",
        data=merged_tau,
        order=["non-specific", "tissue_specific"],
        palette=colours_tau,
        ax=ax2,
    )
    # add horizontal black line
    ax1.axhline(0, color="black")
    ax2.axhline(0, color="black")
    # add graph titles
    ax1.set_title("A", x=0.02, fontsize=16)
    ax2.set_title("B", x=0.02, fontsize=16)
    # set y axis range
    # first get y limits for both plots
    ax1ymin, ax1ymax = ax1.get_ylim()
    ax2ymin, ax2ymax = ax2.get_ylim()
    minimum = min(ax1ymin, ax2ymin)
    maximum = max(ax1ymax, ax2ymax)
    ax1.set_ylim(minimum, maximum)
    ax2.set_ylim(minimum, maximum)
    # tight layout
    # plt.tight_layout()
    # label axes
    ax1.set_xlabel("Gene type")
    ax1.set_ylabel("Log2-fold enrichment over background")
    ax2.set_xlabel("Gene type")
    ax2.set_ylabel("Log2-fold enrichment over background")
    f.savefig(
        f"../../data/output/{file_names}/TATA/{output_folder_name}plots/{output_prefix}_log2fold.pdf",
        format="pdf",
        # bbox_inches="tight",
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
        args.gat_TATA_constitutive_output,
        args.gat_TATA_variable_output,
        args.gat_TATA_nonspecific_output,
        args.gat_TATA_tissue_specific_output,
        args.output_prefix,
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        args.file_names,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
