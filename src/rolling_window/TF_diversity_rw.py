import argparse
import os

import numpy as np
import pandas as pd
import skbio
from pybedtools import BedTool


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="TF_diversity_rw")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "window_bed",
        type=str,
        help="Input location of rolling window bed file",
    )
    parser.add_argument(
        "mapped_motif_bed",
        type=str,
        help="Input location of promoters mapped motif bed",
    )
    parser.add_argument(
        "window_motifs_bed",
        type=str,
        help="Output location of windows_motifs intersect bed",
    )
    parser.add_argument(
        "TF_diversity_bed",
        type=str,
        help="Output location of the window TF_diversity_bed",
    )
    parser.add_argument("outputfolder", type=str, help="Output folder name")
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def merge_windows_motifs(window_bed, mapped_motif_bed, output_file):
    """perform bedtools intersect on the two dfs"""
    windows = BedTool(window_bed)
    motifs = mapped_motif_bed
    # -wao =Write the original A and B entries plus the number of base pairs of overlap between the two features.
    # However, A features w/o overlap are also reported with a NULL B feature and overlap = 0
    intersect = windows.intersect(motifs, wao=True)
    # Write to output_file
    with open(output_file, "w") as output:
        # Each line in the file contains bed entry a and bed entry b that it overlaps plus the number of bp in the overlap so 19 columns
        output.write(str(intersect))


def calculate_shannon_diversity(window_motifs_bed):
    """read in mapped motifs_bed, window_bed, merge them and then calculate Shannon diversity-like calculations for each window"""
    window_motifs_df = pd.read_table(window_motifs_bed, sep="\t", header=None)
    cols = [
        "chr",
        "start",
        "stop",
        "name",
        "motif_chr",
        "motif_start",
        "motif_stop",
        "name_rep",
        "score",
        "strand",
        "promoter_AGI",
        "p-value",
        "q-value",
        "matched_sequence",
        "TF_name",
        "TF_family",
        "TF_AGI",
        "bp_overlap",
    ]
    window_motifs_df.columns = cols
    # filter columns
    window_motifs_df = window_motifs_df[
        [
            "chr",
            "start",
            "stop",
            "name",
            "motif_chr",
            "motif_start",
            "motif_stop",
            "TF_name",
            "TF_family",
            "TF_AGI",
            "bp_overlap",
        ]
    ]
    # add AGI column
    window_motifs_df["promoter_AGI"] = window_motifs_df.name.str.split(
        "_", expand=True
    )[0]
    # make window number column
    window_motifs_df["window_number"] = window_motifs_df.name.str.split(
        "_", expand=True
    )[1]
    # add window length column
    window_motifs_df = window_motifs_df.assign(
        window_length=window_motifs_df.stop - window_motifs_df.start
    )
    # add motif length column
    window_motifs_df["motif_length"] = 0
    window_motifs_df.motif_length = (
        window_motifs_df.motif_stop - window_motifs_df.motif_start
    )
    # if bp_covered is greater than 0 but less than half of the motif length, remove line
    df = window_motifs_df[
        ~(
            (window_motifs_df["bp_overlap"] > 0)
            & (
                window_motifs_df["bp_overlap"]
                < 0.5 * window_motifs_df.motif_length
            )
        )
    ].copy()
    # replace columns with only dots in them with NaN
    no_motifs = df[df.bp_overlap == 0].copy()
    no_motifs["motif_chr"] = np.NaN
    no_motifs["motif_start"] = np.NaN
    no_motifs["motif_stop"] = np.NaN
    no_motifs["TF_name"] = np.NaN
    no_motifs["TF_family"] = np.NaN
    no_motifs["TF_AGI"] = np.NaN
    # make df without the lines with no motifs
    df_reduced = df[~(df.bp_overlap == 0)].copy()
    # merge the two dfs
    df = pd.concat([df_reduced, no_motifs])
    # sort on name
    df.sort_values("name", inplace=True, ignore_index=True)
    # count no. of each TF binding in each window
    groupby_promoter_counts = (
        df.groupby("name")["TF_AGI"]
        .value_counts(dropna=True)
        .unstack(fill_value=0)
    )
    # count no. of TF families binding in each promoter
    groupby_promoter_counts_family = (
        df.groupby("name")["TF_family"]
        .value_counts(dropna=True)
        .unstack(fill_value=0)
    )
    # Individual TF shannon diversity using arbitrary log2 base
    shannon_div_df = groupby_promoter_counts.apply(
        pd.Series(lambda x: skbio.diversity.alpha.shannon(x, base=2)), axis=1
    )
    # shannon diversity for TF family
    shannon_div_TF_family_df = groupby_promoter_counts_family.apply(
        pd.Series(lambda x: skbio.diversity.alpha.shannon(x, base=2)), axis=1
    )
    # convert rownames into column
    cols = ["name", "shannon"]
    shannon_div_df.index.name = "name"
    shannon_div_df.reset_index(inplace=True)
    shannon_div_TF_family_df.index.name = "name"
    shannon_div_TF_family_df.reset_index(inplace=True)
    # rename column
    shannon_div_df.rename(
        columns={"<lambda>": "Shannon_diversity_TF"}, inplace=True
    )
    shannon_div_TF_family_df.rename(
        columns={"<lambda>": "Shannon_diversity_TF_family"}, inplace=True
    )
    # merge individual TF and TF family diversity dfs
    diversity_df = pd.merge(
        shannon_div_df, shannon_div_TF_family_df, on="name"
    )

    # calculate unique TF counts
    # groupby promoter, and include only unique TFs within each promoter group. Preserve column names.
    unique_TF_count = df.groupby(by="name", as_index=False).agg(
        {"TF_AGI": pd.Series.nunique}
    )
    # rename column
    unique_TF_count.rename(columns={"TF_AGI": "unique_TF_count"}, inplace=True)

    # calculate total TF counts
    total_TF_count = df.groupby(by="name", as_index=False).agg(
        {"TF_AGI": pd.Series.count}
    )
    # rename column
    total_TF_count.rename(columns={"TF_AGI": "total_TF_count"}, inplace=True)

    # calculate total TF family counts
    total_TF_family_count = df.groupby(by="name", as_index=False).agg(
        {"TF_family": pd.Series.nunique}
    )
    # rename column
    total_TF_family_count.rename(
        columns={"TF_family": "TF_family_count"}, inplace=True
    )

    # merge diversity df with unique_TF_count
    diversity_df = pd.merge(diversity_df, unique_TF_count, on="name")
    # then merge with total_TF_count
    diversity_df = pd.merge(diversity_df, total_TF_count, on="name")
    # then merge with TF_family_count
    diversity_df = pd.merge(diversity_df, total_TF_family_count, on="name")

    # if name isn't in diversity_df, add it and change all other columns to 0
    # select all unique names that aren't in diversity_df (will select all values but good to be sure)
    missing_names = no_motifs[
        ~no_motifs.name.isin(diversity_df.name)
    ].name.unique()
    # turn into df
    missing_names_diversity = pd.DataFrame(missing_names)
    # make columns
    missing_names_diversity.columns = ["name"]
    missing_names_diversity["Shannon_diversity_TF"] = -0
    missing_names_diversity["Shannon_diversity_TF_family"] = -0
    missing_names_diversity["unique_TF_count"] = 0
    missing_names_diversity["total_TF_count"] = 0
    missing_names_diversity["TF_family_count"] = 0
    # concatenate missing_names_diversity with missing_names_diversity
    diversity_df = pd.concat([diversity_df, missing_names_diversity])
    return diversity_df


def main(args):
    # parse arguments
    args = parse_args(args)

    # make directory for the output files to be exported to
    # dirName = f'{args.directory_path}/data/output/{args.file_names}'
    dirName = f"../../data/output/{args.file_names}/rolling_window/"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    dirName = f"../../data/output/{args.file_names}/rolling_window/{args.outputfolder}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # merge windows with motifs and create output bed file
    merge_windows_motifs(
        args.window_bed, args.mapped_motif_bed, args.window_motifs_bed
    )
    # make shannon df
    shannon_df = calculate_shannon_diversity(args.window_motifs_bed)
    # write out df
    shannon_df.to_csv(args.TF_diversity_bed, index=False, sep="\t", header=1)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
