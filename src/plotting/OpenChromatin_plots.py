import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import seaborn as sns
from statannot import add_stat_annotation


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="OpenChromatin_plots")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "Czechowski_gene_categories",
        type=str,
        help="Input location of Czechowski gene categories text file",
    )
    parser.add_argument(
        "Chromatin_bp_covered",
        type=str,
        help="Input location of bp covered of chromatin in promoters text file",
    )
    parser.add_argument(
        "output_folder_name",
        type=str,
        help="Optional output folder name ending in a forward slash",
        default="",
        nargs="?",
    )
    parser.add_argument(
        "output_folder_name_promoter",
        type=str,
        help="Optional output folder name ending in a forward slash (name this after promoter set name)",
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
        "author_name",
        type=str,
        help="Optional replacement name for author in reference to the geneset",
        default="Czechowski",
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


def rep_sample(df, col, n, random_state):
    """function to return a df with equal sample sizes
    taken from here: https://stackoverflow.com/questions/39457762/python-pandas-conditionally-select-a-uniform-sample-from-a-dataframe"""
    # identify number of categories
    nu = df[col].nunique()
    # find number of rows
    # m = len(df)
    # integar divide total sample size by number of categories
    mpb = n // nu
    # multiply this by the number of categories and subtract from the number of samples to find the remainder
    mku = n - mpb * nu
    # make an array fileld with zeros corresponding to each category
    fills = np.zeros(nu)

    # make values in the array 1s up until the remainder
    fills[:mku] = 1

    # calculate sample sizes for each category
    sample_sizes = (np.ones(nu) * mpb + fills).astype(int)

    # group the df by categories
    gb = df.groupby(col)

    # define sample size function
    def sample(sub_df, i):
        return sub_df.sample(sample_sizes[i], random_state=random_state)

    # sample = lambda sub_df, i: sub_df.sample(
    #     sample_sizes[i], random_state=random_state
    # )
    # run sample size function on each category
    subs = [sample(sub_df, i) for i, (_, sub_df) in enumerate(gb)]
    # return concatenated sub dfs
    return pd.concat(subs)


def percent_coverage(bp_covered):
    """function to calculate the % coverage from the output file of bedtools coverage"""

    coverage_df = pd.read_table(bp_covered, sep="\t", header=None)
    col = [
        "chr",
        "start",
        "stop",
        "AGI",
        "dot",
        "strand",
        "source",
        "type",
        "dot2",
        "details",
        "no._of_overlaps",
        "no._of_bases_covered",
        "promoter_length",
        "fraction_bases_covered",
    ]
    coverage_df.columns = col
    # add % bases covered column
    coverage_df["percentage_bases_covered"] = (
        coverage_df.fraction_bases_covered * 100
    )

    # remove unnecessary columns
    coverage_df_reduced_columns = coverage_df[
        [
            "chr",
            "start",
            "stop",
            "AGI",
            "strand",
            "no._of_overlaps",
            "no._of_bases_covered",
            "promoter_length",
            "fraction_bases_covered",
            "percentage_bases_covered",
        ]
    ]
    return coverage_df_reduced_columns


def merge_genecategories(df, gene_categories):
    """merge df with gene categories"""
    # read in gene categories
    gene_cats = pd.read_csv(gene_categories, sep="\t", header=None)
    gene_cats.columns = ["AGI", "gene_type"]
    # merge to limit to genes of interest
    df_categories = pd.merge(gene_cats, df, how="left", on="AGI")
    return df_categories


def dunn_posthoc_test(df, dependent_variable, between):
    """dunn_posthoc tests with bonferroni multiple correction"""
    return sp.posthoc_dunn(
        df,
        val_col=dependent_variable,
        group_col=between,
        p_adjust="bonferroni",
    )


def make_plot(
    file_names,
    output_folder_name_promoter,
    output_folder_name,
    df,
    x_variable,
    y_variable,
    x_label,
    y_label,
    output_prefix,
    plot_kind,
    palette,
    variable1_name,
    variable2_name,
    dependent_variable,
):
    """function to make and save plot"""
    # allow colour codes in seaborn
    sns.set(color_codes=True)
    sns.set_style("whitegrid")
    # plot
    x = x_variable
    y = y_variable
    order = [variable1_name, variable2_name, "control"]
    # set colour palette
    colours = sns.color_palette(palette)

    # make copy of df
    merged2_unique = df.copy()
    # make sample sizes equal for comparison
    # identify sample size of the minimum category
    minimum_sample_size = merged2_unique.gene_type.value_counts().min()
    # print this
    print(f"sample size in each category = {minimum_sample_size}")
    # save sample size as file
    with open(
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name_promoter}{output_folder_name}plots/number_of_genes_in_each_category_{dependent_variable}.txt",
        "w",
    ) as file:
        file.write(
            "number_of_genes_in_each_category=" + str(minimum_sample_size)
        )

    # multiply this by the number of categories
    total_sample_size = minimum_sample_size * len(
        merged2_unique.gene_type.unique()
    )
    # select equal sample sizes of each category with a random state of 1 so it's reproducible
    equal_samplesizes = rep_sample(
        merged2_unique, "gene_type", total_sample_size, random_state=1
    )
    # now filter out genes which were not selected using the minimum sample size
    to_remove = merged2_unique[~merged2_unique.AGI.isin(equal_samplesizes.AGI)]
    df = df[~df.AGI.isin(to_remove.AGI)]

    # if violin plot don't extend past datapoints
    if plot_kind == "violin":
        sns.catplot(
            x=x,
            y=y,
            data=df,
            kind=plot_kind,
            order=order,
            palette=colours,
            cut=0,
        )
    else:
        sns.catplot(
            x=x, y=y, data=df, kind=plot_kind, order=order, palette=colours
        )
    # plot points
    ax = sns.swarmplot(x=x, y=y, data=df, color=".25", order=order)
    # add significance if necessary - dunn's posthocs with multiple Bonferroni correction
    stat = dunn_posthoc_test(df, y_variable, x_variable)
    # label box pairs
    box_pairs = [
        (variable1_name, variable2_name),
        (variable1_name, "control"),
        (variable2_name, "control"),
    ]
    # make empty list of p_values
    p_values = []
    # populate the list of p_values accoridng to the box_pairs
    for pair in box_pairs:
        print(pair)
        # select p value for each pair
        p = stat.loc[pair[0], pair[1]]
        p_values.append(p)

    # add stats annotation to the plot
    add_stat_annotation(
        ax,
        data=df,
        x=x,
        y=y,
        order=order,
        box_pairs=box_pairs,
        text_format="star",
        loc="outside",
        verbose=2,
        perform_stat_test=False,
        pvalues=p_values,
        test_short_name="Dunn",
    )

    # change axes labels
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    # tight layout
    plt.tight_layout()
    # save figure
    ax.get_figure().savefig(
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name_promoter}{output_folder_name}plots/{output_prefix}_{plot_kind}.pdf",
        format="pdf",
    )


def all_prom_distribution(
    df,
    x_variable,
    x_label,
    output_prefix,
    file_names,
    output_folder_name_promoter,
    output_folder_name,
    dependent_variable,
):
    """function to return distribution plot of all promoters and the chromatin % bp covered"""

    dist_plot = df[x_variable]
    # create figure with no transparency
    dist_plot_fig = sns.distplot(dist_plot).get_figure()
    plt.xlabel(x_label)

    # save to file
    dist_plot_fig.savefig(
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name_promoter}{output_folder_name}plots/{output_prefix}_distribution.pdf",
        format="pdf",
    )


def main(args):
    # parse arguments
    args = parse_args(args)

    dependent_variable = "chromatin_coverage"

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{dependent_variable}/{args.output_folder_name_promoter}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{dependent_variable}/{args.output_folder_name_promoter}{args.output_folder_name}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{dependent_variable}/{args.output_folder_name_promoter}{args.output_folder_name}plots"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # add % coverage column
    chromatin_coverage = percent_coverage(args.Chromatin_bp_covered)

    # add gene categories to the df
    cats = merge_genecategories(
        chromatin_coverage, args.Czechowski_gene_categories
    )

    # all promoter distribution plot
    all_prom_distribution(
        cats,
        "percentage_bases_covered",
        "% bp covered",
        f"{dependent_variable}_allproms",
        args.file_names,
        args.output_folder_name_promoter,
        args.output_folder_name,
        dependent_variable,
    )

    # Czechowski_gene_categories box plot
    make_plot(
        args.file_names,
        args.output_folder_name_promoter,
        args.output_folder_name,
        cats,
        "gene_type",
        "percentage_bases_covered",
        "Gene type",
        "% bp covered",
        f"{args.author_name}_{dependent_variable}",
        "box",
        args.palette,
        args.variable1_name,
        args.variable2_name,
        dependent_variable,
    )
    make_plot(
        args.file_names,
        args.output_folder_name_promoter,
        args.output_folder_name,
        cats,
        "gene_type",
        "percentage_bases_covered",
        "Gene type",
        "% bp covered",
        f"{args.author_name}_{dependent_variable}",
        "violin",
        args.palette,
        args.variable1_name,
        args.variable2_name,
        dependent_variable,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
