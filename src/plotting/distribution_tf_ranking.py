import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args(args):
    # arg1 is the name of folder and filenames for the promoters extracted
    # arg2 is the table of interest showing all TFs with their cv value
    # arg3 is the table of interest showing all TFs with their tau value
    # arg4 is the cv tf gene categories
    # arg5 is the tau tf gene categories

    parser = argparse.ArgumentParser(description="GC_content_plots")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "cv_all_tfs",
        type=str,
        help="Input file containing all transcription factors with their CV values",
    )
    parser.add_argument(
        "tau_all_tfs",
        type=str,
        help="Input file containing all transcription factors with their TAU values",
    )
    parser.add_argument(
        "cv_tf_gene_categories",
        type=str,
        help="Input file containing CV transcription factor categories",
    )
    parser.add_argument(
        "tau_tf_gene_categories",
        type=str,
        help="Input file containing TAU transcription factor categories",
    )

    return parser.parse_args(args)


def read_files(
    cv_all_tfs, tau_all_tfs, cv_tf_gene_categories, tau_tf_gene_categories
):
    """Read in the files as dataframes"""
    cv_all_tfs_df = pd.read_table(cv_all_tfs, sep="\t", header=0)
    tau_all_tfs_df = pd.read_table(tau_all_tfs, sep="\t", header=0)
    cv_tf_gene_categories_df = pd.read_table(
        cv_tf_gene_categories, sep="\t", header=None
    )
    tau_tf_gene_categories_df = pd.read_table(
        tau_tf_gene_categories, sep="\t", header=None
    )
    cols = ["AGI", "gene_type"]
    cv_tf_gene_categories_df.columns = cols
    tau_tf_gene_categories_df.columns = cols
    return (
        cv_all_tfs_df,
        tau_all_tfs_df,
        cv_tf_gene_categories_df,
        tau_tf_gene_categories_df,
    )


def all_prom_distribution(
    file_names,
    output_prefix,
    palette_colours,
    df,
    x_variable,
    x_label,
    df2=pd.DataFrame(),
    df1_label="",
    df2_label="",
    labels=False,
    min_x_constitutive=False,
    max_x_constitutive=False,
    min_x_variable=False,
    max_x_variable=False,
    save=False,
):
    """function to return distribution plot of all promoters of variable of interest.
    df1_label and df2 labels are the names of the respective gene type subset in the df"""
    # if only 1 dataframe provided then create just 1 plot
    # set colour palette
    sns.set_palette(palette_colours)
    if df2.empty:
        dist_plot = df[x_variable]
        # create figure with no transparency
        dist_plot_fig = sns.distplot(dist_plot)
        plt.xlabel(x_label)
        plt.ylabel("Kernel density")
    # else if 2 dataframes provided plot them on the same axes
    else:
        dist_plot1 = df[x_variable]
        dist_plot2 = df2[x_variable]
        dist_plot_fig = sns.distplot(
            dist_plot1, hist=False, rug=True, label=df1_label
        )
        sns.distplot(dist_plot2, hist=False, rug=True, label=df2_label)
        # create legend
        plt.legend()
        plt.ylabel("Kernel density")
    if labels is True:
        # get axes
        ax = plt.axes()
        # constitutive annotation
        ax.annotate(
            "top 100 constitutive range",
            xy=(max_x_constitutive, 0.2),
            xycoords="data",
            ha="left",
            xytext=(50, 100),
            textcoords="offset points",
            arrowprops=dict(
                arrowstyle="->", connectionstyle="arc3,rad=0.4", color="black"
            ),
        )
        ax.annotate(
            "",
            xy=(max_x_constitutive, 0.2),
            xytext=(min_x_constitutive, 0.2),
            xycoords="data",
            textcoords="data",
            arrowprops={
                "arrowstyle": "|-|,widthA=0.2,widthB=0.2",
                "color": "blue",
            },
        )
        # Variable annotation
        ax.annotate(
            "top 100 variable range",
            xy=(max_x_variable, 0.2),
            xycoords="data",
            ha="right",
            xytext=(0, -20),
            textcoords="offset points",
        )
        ax.annotate(
            "",
            xy=(max_x_variable, 0.2),
            xytext=(min_x_variable, 0.2),
            xycoords="data",
            textcoords="data",
            arrowprops={
                "arrowstyle": "|-|,widthA=0.2,widthB=0.2",
                "color": "orange",
            },
        )

    # save to file
    if save is True:
        dist_plot_fig.get_figure().savefig(
            f"../../data/output/{file_names}/genes/TFs/plots/{output_prefix}_distribution.pdf",
            format="pdf",
        )
    plt.clf()
    return dist_plot_fig


def main(args):
    args = parse_args(args)

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/genes/TFs/"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")
        # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/genes/TFs/plots"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    (
        cv_all_tfs_df,
        tau_all_tfs_df,
        cv_tf_gene_categories_df,
        tau_tf_gene_categories_df,
    ) = read_files(
        args.cv_all_tfs,
        args.tau_all_tfs,
        args.cv_tf_gene_categories,
        args.tau_tf_gene_categories,
    )
    # all czechowski CV distribution plot
    all_prom_distribution(
        args.file_names,
        "cv_dist_plot_czechowski",
        "deep",
        cv_all_tfs_df,
        "expression_CV",
        "Expression CV",
        save=True,
    )
    # all TAU distribution plot
    all_prom_distribution(
        args.file_names,
        "tau_dist_plot_czechowski",
        "tab10_r",
        tau_all_tfs_df,
        "tau",
        "Tau",
        save=True,
    )

    # plot TAU for top variable and top constitutive genes from CV ranking'
    merged = pd.merge(
        cv_tf_gene_categories_df,
        tau_all_tfs_df,
        left_on="AGI",
        right_on="Gene_ID",
        how="left",
    )

    # plot CV from CV ranking for top tissue_specific and top constitutive genes from CV ranking'

    # merge dfs
    merged2 = pd.merge(
        tau_tf_gene_categories_df,
        cv_all_tfs_df,
        left_on="AGI",
        right_on="Gene_ID",
        how="left",
    )
    # ONLY ? GENES FROM TAU SET ARE IN THE CV SET
    # constitutive genes from RNA-seq in microarray
    # variable genes from RNA-seq in microarray
    # control genes from RNA-seq in microarray
    all_prom_distribution(
        args.file_names,
        "cv_taucategories_dist_czechowski",
        "tab10_r",
        merged2[merged2.gene_type == "non-specific"],
        "expression_CV",
        "Expression CV",
        merged2[merged2.gene_type == "tissue_specific"],
        df1_label="non-specific",
        df2_label="tissue_specific",
        save=True,
    )

    # plot the tau of the CV ranking variable and constitutive gene sets '

    all_prom_distribution(
        args.file_names,
        "tau_cvcategories_dist_czechowski",
        "deep",
        merged[merged.gene_type == "constitutive"],
        "tau",
        "Tau",
        merged[merged.gene_type == "variable"],
        df1_label="constitutive",
        df2_label="variable",
        save=True,
    )
    # Indeed constitutive genes set from the microarray have lower TAU than the variable gene set

    # czechovski CV distribution
    # using CVs and gene categories from only the microarray, plot CVs of the constitutive and variable gene sets
    merged_czechowski_cv = pd.merge(
        cv_tf_gene_categories_df, cv_all_tfs_df, on="AGI", how="left"
    )
    all_prom_distribution(
        args.file_names,
        "cv_categories_czechowski",
        "deep",
        merged_czechowski_cv[merged_czechowski_cv.gene_type == "constitutive"],
        "expression_CV",
        "Expression CV",
        merged_czechowski_cv[merged_czechowski_cv.gene_type == "variable"],
        df1_label="constitutive",
        df2_label="variable",
        save=True,
    )

    # TAU distribution
    # using tau and gene categories from only the tau ranking data, plot tau of the constitutive and tissue_specific gene sets
    merged_czechowski_tau = pd.merge(
        tau_tf_gene_categories_df, tau_all_tfs_df, on="AGI", how="left"
    )
    all_prom_distribution(
        args.file_names,
        "tau_categories_czechowski",
        "tab10_r",
        merged_czechowski_tau[
            merged_czechowski_tau.gene_type == "non-specific"
        ],
        "tau",
        "Tau",
        merged_czechowski_tau[
            merged_czechowski_tau.gene_type == "tissue_specific"
        ],
        df1_label="non-specific",
        df2_label="tissue_specific",
        save=True,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
