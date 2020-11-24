import argparse

import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(description="choose_TFs_tau")
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
        "Schmid_allgenes",
        type=str,
        help="Input location of all filtered microarray genes ranked by tau",
    )
    parser.add_argument(
        "Schmid_gene_categories_tfs",
        type=str,
        help="Output location of TF gene categories",
    )
    parser.add_argument(
        "CV_TF_gene_categories",
        type=str,
        help="Input location of Czechowski coefficient of variation gene categories",
    )
    parser.add_argument(
        "CV_all_genes",
        type=str,
        help="Optional input location of coefficient of variation all genes after filtering to include only genes present in 80% of conditions",
        default=None,
        nargs="?",
    )
    parser.add_argument(
        "Schmid_all_tfs",
        type=str,
        help="Output location of all ranked TFs by tau value",
    )
    return parser.parse_args(args)


def filter_genes_schmid(TF_list, select_genes_file, Schmid_all_tfs):
    """filter out genes from the microarray data which aren't in the promoter_bed"""
    select_genes = pd.read_table(select_genes_file, sep="\t", header=0)

    # make AGI capitalised
    # select_genes.AGI = select_genes.AGI.str.upper()

    TF_list_df = pd.read_table(TF_list, sep="\t", header=1)
    merged = pd.merge(
        TF_list_df, select_genes, left_on="Gene_ID", right_on="AGI", how="left"
    )
    # remove NaNs in tau column
    filtered_df = merged[merged.tau.notnull()].copy()

    # sort by TF ID
    filtered_df.sort_values(["TF_ID"], inplace=True, ignore_index=True)

    # remove duplicates based on Gene_ID column
    filtered_df.drop_duplicates(subset=["Gene_ID"], inplace=True)

    # sort by tau value
    filtered_df.sort_values("tau", inplace=True, ignore_index=True)
    # save df
    filtered_df.to_csv(
        Schmid_all_tfs, sep="\t", columns=filtered_df.columns, index=False
    )

    return filtered_df


def subSet_ontau(
    in_df, out_dir, no_of_genes, CV_gene_categories, CV_all_genes
):
    """
    Extract the non-specific, tissue_specific, and control subsets based on tau values
    """
    # filtering based on presence in the Araport11 annotation, define the first
    # n rows as the non-specific set and add label
    nonspecific = in_df[0:no_of_genes].copy()
    nonspecific["state"] = "non-specific"

    # define the last n rows as the tissue_specific set and add label
    tissue_specific = in_df[-no_of_genes:].copy()
    tissue_specific["state"] = "tissue_specific"

    # read in czechowski gene categories so the midrange samples exclude the constitutive and variable genes from there
    CV_categories = pd.read_table(CV_gene_categories, sep="\t", header=None)
    # Name columns
    cols = ["AGI", "gene_type"]
    CV_categories.columns = cols
    # exclude control genes from CV gene categories
    CV_categories_nocontrol = CV_categories[
        ~(CV_categories.gene_type == "control")
    ]

    # extract the rest of the rows as the control search space
    mid_range_full = in_df[(no_of_genes + 1) : -(no_of_genes + 1)].copy()

    # exclude the constitutive and variable genes
    mid_range = mid_range_full[
        ~(mid_range_full.AGI.isin(CV_categories_nocontrol.AGI))
    ]

    # exclude genes which aren't found within the CV all gene set
    # read in CV all genes
    if CV_all_genes is None:
        mid_range_reduced = mid_range.copy()
    else:
        CV_all_genes_df = pd.read_table(CV_all_genes, header=0, sep="\t")
        mid_range_reduced = mid_range[
            mid_range.AGI.isin(CV_all_genes_df.AGI)
        ].copy()

    # create 10 labelled bins
    mid_range_reduced["bins"] = pd.Series(
        pd.qcut(mid_range_reduced["tau"], q=10, precision=2)
    ).copy()

    # extract 10 random rows from these bins and label as the control set
    sample = no_of_genes / 10
    sample_integar = int(
        str(sample).replace(".0", "")
    )  # convert sample to an integar
    samples_from_bins = (
        mid_range_reduced.groupby("bins")
        .apply(pd.DataFrame.sample, sample_integar, random_state=2)
        .copy()
    )
    samples_from_bins["state"] = "control"

    # concatenate and write as output
    output_set = pd.concat(
        [
            nonspecific[["AGI", "state"]],
            tissue_specific[["AGI", "state"]],
            samples_from_bins[["AGI", "state"]],
        ],
        ignore_index=True,
    )
    output_set.to_csv(out_dir, sep="\t", index=False, header=False)

    # replace control genes in CV categories with same control genes as tau
    if CV_all_genes is None:
        pass
    else:
        CV_categories[CV_categories.gene_type == "control"] = output_set[
            output_set.state == "control"
        ]
        # save CV gene categories
        CV_categories.to_csv(
            CV_gene_categories, sep="\t", header=None, index=False
        )

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
    # filter out genes from the microarray data which aren't in the promoter_bed
    filtered_schmid = filter_genes_schmid(
        args.TF_list, args.Schmid_allgenes, args.Schmid_all_tfs
    )
    # Extract the non-specific, tissue_specific, and control subsets based on tau values
    subSet_ontau(
        filtered_schmid,
        args.Schmid_gene_categories_tfs,
        args.no_of_genes,
        args.CV_TF_gene_categories,
        args.CV_all_genes,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
