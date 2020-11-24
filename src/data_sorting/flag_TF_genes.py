import argparse

import pandas as pd


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="flag_TF_genes")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "gene_categories",
        type=str,
        help="Input location of Czechowski gene categories text file",
    )
    parser.add_argument(
        "TF_list",
        type=str,
        help="Input location of Arabidopsis transcription factor list",
    )
    parser.add_argument(
        "all_genes",
        type=str,
        help="Input location of Czechowski et al 2005 microarray all genes data",
    )
    parser.add_argument(
        "output", type=str, help="Output location of flagged TF genes"
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def flagTFs(genes_of_interest, TF_list, all_genes, output):
    """Function to export a file with all genes which are TFs from a list of genes of interest"""
    genes_of_interest_df = pd.read_table(
        genes_of_interest,
        sep="\t",
        header=None,
        names=["Gene_ID", "gene_type"],
    )
    TFs = pd.read_table(TF_list, sep="\t", header=1)
    # merge the dfs so that only TFs are included
    TFs_of_interest = pd.merge(genes_of_interest_df, TFs, on="Gene_ID")

    # Get the CV and expression data
    all_genes_df = pd.read_table(all_genes, sep="\t", header=0)
    # merge with TFs_of_interest
    TFs_of_interest = pd.merge(
        TFs_of_interest,
        all_genes_df,
        how="left",
        left_on="Gene_ID",
        right_on="AGI",
    )

    # select only columns of interest

    TFs_of_interest = TFs_of_interest[
        [
            "Gene_ID",
            "gene_type",
            "Family",
            "expression_mean",
            "expression_CV",
            "constitutive_in_araport11",
        ]
    ]
    # drop duplicates
    TFs_of_interest = TFs_of_interest.drop_duplicates("Gene_ID")
    # make output file containing only genes from genes of interest which are transcription factors
    TFs_of_interest.to_csv(output, sep="\t", index=False)


def main(args):
    # parse arguments
    args = parse_args(args)
    flagTFs(args.gene_categories, args.TF_list, args.all_genes, args.output)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
