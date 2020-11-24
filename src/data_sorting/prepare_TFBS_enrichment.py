import argparse
import io
import os

import pandas as pd
from pybedtools import BedTool


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="prepare_TFBS_enrichment")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument("promoterpref", type=str, help="Promoter prefix name")
    parser.add_argument(
        "gene_categories",
        type=str,
        help="Input location of gene categories text file",
    )
    parser.add_argument(
        "promoter_bed_file",
        type=str,
        help="Input location of promoters bed file",
    )
    parser.add_argument(
        "output_genecat_prefix",
        type=str,
        help="Gene category prefix (eg. Czechowski)",
    )
    parser.add_argument(
        "mapped_motifs",
        type=str,
        help="input location of mapped_motifs bed file",
    )  #
    parser.add_argument(
        "mapped_motifs_out",
        type=str,
        help="output location of mapped_motifs bed file with TF family in column 4",
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
        help="Optional replacement name for 2nd variable eg. tissue&condition_specific",
        default="variable",
        nargs="?",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def sort_files(
    promoter_bed_file, mapped_motifs_file, gene_categories, mapped_motifs_out
):
    """sort the promoters, motifs and genes types data ready for gat TATA enrichment"""
    promoters = pd.read_table(promoter_bed_file, sep="\t", header=None)
    col = [
        "chr",
        "start",
        "stop",
        "promoter_AGI",
        "dot1",
        "strand",
        "source",
        "type",
        "dot2",
        "attributes",
    ]
    promoters.columns = col
    # read in gene categories
    gene_cats = pd.read_table(gene_categories, sep="\t", header=None)
    cols = ["promoter_AGI", "gene_type"]
    gene_cats.columns = cols
    # merge with promoters
    promoters = pd.merge(gene_cats, promoters, on="promoter_AGI", how="left")
    # read in mapped_motifs file
    motifs = pd.read_table(mapped_motifs_file, sep="\t", header=None)

    # select the correct columns depending on the number of columns in the file
    if len(motifs.columns) == 24:
        cols = [
            "chr",
            "start",
            "stop",
            "promoter_AGI",
            "dot1",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
            "motif_chr",
            "motif_start",
            "motif_stop",
            "name_rep",
            "score",
            "motif_strand",
            "promoter_AGI2",
            "p-value",
            "q-value",
            "matched_sequence",
            "TF_name",
            "TF_family",
            "TF_AGI",
            "bp_overlap",
        ]
        motifs.columns = cols
        # filter columns
        motifs = motifs[
            [
                "motif_chr",
                "motif_start",
                "motif_stop",
                "name_rep",
                "score",
                "motif_strand",
                "promoter_AGI2",
                "p-value",
                "q-value",
                "matched_sequence",
                "TF_name",
                "TF_family",
                "TF_AGI",
            ]
        ]
        # rename columns
        cols = [
            "chr",
            "start",
            "stop",
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
        ]
        motifs.columns = cols

    elif len(motifs.columns) == 13:
        cols = [
            "chr",
            "start",
            "stop",
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
        ]
        motifs.columns = cols

    elif len(motifs.columns) == 17:
        cols = [
            "chr",
            "start",
            "stop",
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
            "chr_openchrom",
            "start_openchrom",
            "stop_openchrom",
            "bp_overlap",
        ]
        motifs.columns = cols

    # swap columns around so TF family is in column 4
    motifs = motifs[
        [
            "chr",
            "start",
            "stop",
            "TF_family",
            "score",
            "strand",
            "promoter_AGI",
            "p-value",
            "q-value",
            "matched_sequence",
            "TF_name",
            "name_rep",
            "TF_AGI",
        ]
    ]
    # write out the file
    BedTool.from_dataframe(motifs).saveas(mapped_motifs_out)

    return promoters


def prepare_gat(
    df,
    file_names,
    output_genecat_prefix,
    promoterpref,
    variable1_name,
    variable2_name,
):
    """prepare files for running gat analysis - outputs a workspace file containing all promoters, a variable promoter file and a constitutive promoter file"""
    # make buffer to save promoters
    buffer = io.StringIO()
    df.to_csv(buffer, sep="\t", header=None, index=False)
    buffer.seek(0)
    # select only constitutive and variable genes
    # df = df[(df.gene_type == 'constitutive') | (df.gene_type == 'variable')]
    # with gene_type in fourth column for gat enrichment analysis
    df_reordered_gene_type = df[
        [
            "chr",
            "start",
            "stop",
            "gene_type",
            "strand",
            "source",
            "attributes",
            "promoter_AGI",
        ]
    ]
    df_reordered_AGI = df[
        [
            "chr",
            "start",
            "stop",
            "promoter_AGI",
            "strand",
            "source",
            "attributes",
            "gene_type",
        ]
    ]
    # with promoter_AGI in fourth column for other enrichment analysis such as CiiiDER
    # sort by chromosome and start
    sorted_motifs_gene_type = df_reordered_gene_type.sort_values(
        ["chr", "start"]
    )
    sorted_motifs_AGI = df_reordered_AGI.sort_values(["chr", "start"])
    # save bed file with all gene types
    BedTool.from_dataframe(sorted_motifs_gene_type).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/{output_genecat_prefix}_{promoterpref}_allgenetypes_gat.bed"
    )
    BedTool.from_dataframe(sorted_motifs_gene_type).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/{output_genecat_prefix}_{promoterpref}_allgenetypes_AGI.bed"
    )

    # save bed file with only constitutive and variable genetypes
    sorted_motifs_gene_type2 = sorted_motifs_gene_type[
        (sorted_motifs_gene_type.gene_type == variable1_name)
        | (sorted_motifs_gene_type.gene_type == variable2_name)
    ]
    sorted_motifs_AGI2 = sorted_motifs_AGI[
        (sorted_motifs_AGI.gene_type == variable1_name)
        | (sorted_motifs_AGI.gene_type == variable2_name)
    ]
    BedTool.from_dataframe(sorted_motifs_gene_type2).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/{output_genecat_prefix}_{promoterpref}_{variable1_name}{variable2_name}_gat.bed"
    )
    BedTool.from_dataframe(sorted_motifs_AGI2).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/{output_genecat_prefix}_{promoterpref}_{variable1_name}{variable2_name}_AGI.bed"
    )
    # make a new gat workspace file with all promoters (first 3 columns)
    BedTool.from_dataframe(
        sorted_motifs_gene_type[["chr", "start", "stop"]]
    ).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/gat_analysis/{output_genecat_prefix}_{promoterpref}_workspace.bed"
    )
    # select only variable promoters
    variable_promoters_extended = sorted_motifs_gene_type[
        sorted_motifs_gene_type["gene_type"] == variable2_name
    ]
    sorted_variable = variable_promoters_extended.sort_values(["chr", "start"])
    BedTool.from_dataframe(sorted_variable).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/gat_analysis/{output_genecat_prefix}_{promoterpref}_{variable2_name}_gat.bed"
    )
    # same again for other enrichment analyses
    variable_promoters_extended_agi = sorted_motifs_AGI[
        sorted_motifs_AGI["gene_type"] == variable2_name
    ]
    sorted_variable_AGI = variable_promoters_extended_agi.sort_values(
        ["chr", "start"]
    )
    BedTool.from_dataframe(sorted_variable_AGI).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/other_analysis/{output_genecat_prefix}_{promoterpref}_{variable2_name}_AGI.bed"
    )
    # make a constitutive only file
    constitutive_promoters = sorted_motifs_gene_type[
        sorted_motifs_gene_type["gene_type"] == variable1_name
    ]
    sorted_constitutive = constitutive_promoters.sort_values(["chr", "start"])
    BedTool.from_dataframe(sorted_constitutive).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/gat_analysis/{output_genecat_prefix}_{promoterpref}_{variable1_name}_gat.bed"
    )
    # same again for other enrichment analyses with AGI in fourth column
    constitutive_promoters_AGI = sorted_motifs_AGI[
        sorted_motifs_AGI["gene_type"] == variable1_name
    ]
    sorted_constitutive_AGI = constitutive_promoters_AGI.sort_values(
        ["chr", "start"]
    )
    BedTool.from_dataframe(sorted_constitutive_AGI).saveas(
        f"../../data/output/{file_names}/TFBS_enrichment/other_analysis/{output_genecat_prefix}_{promoterpref}_{variable1_name}_AGI.bed"
    )


def main(args):
    # parse arguments
    args = parse_args(args)
    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/TFBS_enrichment"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

        # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/TFBS_enrichment/plots"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

        # make directory for the plots to be exported to
    dirName = (
        f"../../data/output/{args.file_names}/TFBS_enrichment/gat_analysis"
    )
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

        # make directory for the plots to be exported to
    dirName = (
        f"../../data/output/{args.file_names}/TFBS_enrichment/other_analysis"
    )
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")
    promoters = sort_files(
        args.promoter_bed_file,
        args.mapped_motifs,
        args.gene_categories,
        args.mapped_motifs_out,
    )
    # create gat files
    prepare_gat(
        promoters,
        args.file_names,
        args.output_genecat_prefix,
        args.promoterpref,
        args.variable1_name,
        args.variable2_name,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
