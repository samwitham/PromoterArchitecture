import argparse
import os

import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(description="choose_TFs_cv")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "no_of_genes",
        type=int,
        help="Number of genes in each category to subset",
    )
    parser.add_argument(
        "TF_list", type=str, help="List of Arabidopsis transcription factors"
    )
    parser.add_argument(
        "Czechowski_rankedcv",
        type=str,
        help="Input location of all genes ranked by CV",
    )
    parser.add_argument(
        "Czechowski_gene_categories_tfs",
        type=str,
        help="Output location of the TF gene categories",
    )
    parser.add_argument(
        "Czechowski_all_tfs",
        type=str,
        help="Output location of all ranked TFs by CV value",
    )

    return parser.parse_args(args)


def filter_genes_czechowski(TF_list, select_genes_file, Czechowski_all_tfs):
    """filter out genes from the microarray data which aren't in the promoter_bed"""
    # select_genes = pd.read_table(select_genes_file, sep='\t', header=0)
    select_genes = pd.read_table(select_genes_file, sep="\t", header=None)
    cols = [
        "rank",
        "probe_id",
        "AGI",
        "expression_mean",
        "expression_SD",
        "expression_CV",
        "proportion_of_values_present_in_mas5",
        "presence_in_araport11",
        "constitutive_in_araport11",
    ]
    select_genes.columns = cols
    # make AGI uppercase
    select_genes.AGI = select_genes.AGI.str.upper()

    # read in TF list
    TF_list_df = pd.read_table(TF_list, sep="\t", header=1)

    merged = pd.merge(
        TF_list_df, select_genes, left_on="Gene_ID", right_on="AGI", how="left"
    )
    # remove NaNs in expression_CV column
    filtered_df = merged[merged.expression_CV.notnull()].copy()

    # sort by TF ID
    filtered_df.sort_values(["TF_ID"], inplace=True, ignore_index=True)

    # remove duplicates based on Gene_ID column
    filtered_df.drop_duplicates(subset=["Gene_ID"], inplace=True)

    # sort by CV value
    filtered_df.sort_values("expression_CV", inplace=True, ignore_index=True)
    # save df
    filtered_df.to_csv(
        Czechowski_all_tfs, sep="\t", columns=filtered_df.columns, index=False
    )

    return filtered_df


def subSet_onCV(in_df, out_dir, no_of_genes):
    """
    Extract the constitutive, variable, and control subsets based on CV values
    """
    # filtering based on presence in the Araport11 annotation, define the first
    # n rows as the constitutive set and add label
    constitutive = in_df[in_df.presence_in_araport11 == 1][0:no_of_genes]
    constitutive["state"] = "constitutive"

    # define the last n rows as the variable set and add label
    variable = in_df[in_df.presence_in_araport11 == 1][-no_of_genes:]
    variable["state"] = "variable"

    # extract the rest of the rows as the control search space
    mid_range = in_df[in_df.presence_in_araport11 == 1][
        (no_of_genes + 1) : -(no_of_genes + 1)
    ]

    # create 10 labelled bins
    mid_range["bins"] = pd.Series(
        pd.qcut(mid_range["expression_CV"], q=10, precision=2)
    )

    # extract 10 random rows from these bins and label as the control set
    sample = no_of_genes / 10
    sample_integar = int(
        str(sample).replace(".0", "")
    )  # convert sample to an integar
    samples_from_bins = mid_range.groupby("bins").apply(
        pd.DataFrame.sample, sample_integar, random_state=2
    )
    samples_from_bins["state"] = "control"

    # concatenate and write as output
    output_set = pd.concat(
        [
            constitutive[["AGI", "state"]],
            variable[["AGI", "state"]],
            samples_from_bins[["AGI", "state"]],
        ],
        ignore_index=True,
    )
    output_set.to_csv(out_dir, sep="\t", index=False, header=False)

    # function from expressionVar_subsets_plot.py
    # __author__ = "Will Nash"
    # __copyright__ = "Copyright 2020, The Earlham Institute"
    # __credits__ = ["Will Nash", "Wilfried Haerty"]
    # __license__ = "GPL"
    # __version__ = "1.0"
    # __maintainer__ = "Will Nash"
    # __email__ = "will.nash@earlham.ac.uk"
    # __status__ = "Testing"
    # __modified_by__ "Sam Witham"


def main(args):
    # parse arguments
    args = parse_args(args)
    # make directory for the output files to be exported to
    # dirName = f'{args.directory_path}/data/output/{args.file_names}'
    dirName = f"../../data/output/{args.file_names}/genes"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    filtered_czechowski = filter_genes_czechowski(
        args.TF_list, args.Czechowski_rankedcv, args.Czechowski_all_tfs
    )
    # filtered_mergner = filter_genes_mergner(args.promoters_filtered_contain_motifs,args.Mergner_rankedcv)
    # czechowksi subset
    subSet_onCV(
        filtered_czechowski,
        args.Czechowski_gene_categories_tfs,
        args.no_of_genes,
    )
    # mergner subset
    # subSet_onCV(filtered_mergner,args.Mergner_gene_categories,args.no_of_genes)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
