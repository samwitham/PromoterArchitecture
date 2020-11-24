import argparse

import pandas as pd


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="enriched_TFfamilies")
    parser.add_argument(
        "input_tsv",
        type=str,
        help="Input location of enriched TFs with mapped AGI",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output location of enriched TF family coutns text file",
    )

    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def family_count(input_tsv, output):
    """count the number of enriched promoters in each family and save to file"""
    input_tsv_df = pd.read_table(input_tsv, sep="\t", header=0)
    # input_tsv_df.groupby(by='TF_family').agg('count')
    # include only positively enriched TFs
    input_tsv_df_pos = input_tsv_df[input_tsv_df.site_representation == "Up"]
    # count number of enriched genes in each family
    familycounts = input_tsv_df_pos["TF_family"].value_counts().reset_index()
    cols = ["TF_family", "number_of_enriched_TFs"]
    familycounts.columns = cols
    # save as file
    familycounts.to_csv(output, index=False)

    return input_tsv_df_pos["TF_family"].value_counts()


def main(args):
    # parse arguments
    args = parse_args(args)

    family_count(args.input_tsv, args.output)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
