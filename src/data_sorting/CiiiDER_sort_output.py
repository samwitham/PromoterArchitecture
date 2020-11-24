import argparse

import pandas as pd


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="CiiiDER_sort_output")
    parser.add_argument(
        "geneIDtable", type=str, help="location of geneIDtable"
    )
    parser.add_argument(
        "motifs_csv",
        type=str,
        help="location of motifs_csv output from CiiiDER output",
    )
    parser.add_argument(
        "output_tsv",
        type=str,
        help="Output location of tsv file containing significantly enriched motifs with AGI and TF family name",
    )

    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def rename_motif(motifs_csv):
    """function to add a TF name column in the motifs.bed file (for DAP-seq cistrome motifs only)"""
    # read in motifs_bed_file
    motifs = pd.read_csv(motifs_csv, header=0)
    cols = [
        "TF_ID",
        "TF_name",
        "deficit_no",
        "search_sites",
        "mean_sites_per_search_gene",
        "no_background_sites",
        "mean_sites_per_background_gene",
        "site_representation",
        "Mann_whitney_U",
        "site_p_value",
        "total_no_search_genes",
        "no_TF_search_genes",
        "total_background_genes",
        "no_TF_background_genes",
        "gene_representation",
        "gene_p_value",
        "average_log2_proportion_bound",
        "log2_enrichment",
        "significance_score",
    ]
    motifs.columns = cols
    motifs["TF"] = motifs.TF_name

    # define lambda function
    def capitalise(x):
        return x.upper()

    # capitalise = lambda x: x.upper()
    motifs.TF = motifs.TF.apply(capitalise)
    # replace characters upto and including the '.' in the TF column
    motifs.TF = motifs.TF.replace(r"^[^\.]*\.", "", regex=True)
    # replace characters after the '_'
    motifs.TF = motifs.TF.replace(r"_.*$", "", regex=True)
    return motifs


def TF_family(motifs_df):
    """function to add a TF family column in the motifs.bed file (for DAP-seq cistrome motifs only)"""
    motifs_df["TF_family"] = motifs_df.TF_name

    def capitalise(x):
        return x.upper()

    # capitalise = lambda x: x.upper()
    motifs_df.TF_family = motifs_df.TF_family.apply(capitalise)
    # replace characters inluding and after the '_' in the TF_family column
    motifs_df.TF_family = motifs_df.TF_family.replace("_.*$", "", regex=True)
    return motifs_df


def map_ID2(motifs, geneIDtable, output_tsv):
    """function to rename the TF column values in a motifs.bed file to the Arabidopsis gene ID nomenclature using geneIDtable. (for DAP-seq cistrome motifs only). Outputs a bed file."""
    # read in motifID table
    geneIDtable = pd.read_table(geneIDtable, sep="\t", header=None)
    # name columns
    cols = ["TF_ID", "AGI"]
    geneIDtable.columns = cols
    # remove '_m' from end of name_rep value in motifs
    motifs.TF_ID = motifs.TF_ID.str.replace("_m1", "")
    merged = pd.merge(motifs, geneIDtable, on="TF_ID")
    # print(merged.shape)TF_ID
    # make bed file
    # sorted_motifs = merged.sort_values(['chr','start'])
    # bed = BedTool.from_dataframe(sorted_motifs).saveas(output_bed)
    # save output file
    merged.to_csv(output_tsv, sep="\t", header=True, index=None)
    return merged


def main(args):
    # parse arguments
    args = parse_args(args)

    # rename motifs
    motifs_renamed = rename_motif(args.motifs_csv)

    # add motif family name
    motifs_renamed_family = TF_family(motifs_renamed)

    # add AGI and save output file
    map_ID2(motifs_renamed_family, args.geneIDtable, args.output_tsv)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
