import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import seaborn as sns
from pingouin import kruskal
from statannot import add_stat_annotation


def parse_args(args):
    parser = argparse.ArgumentParser(description="GC_content_plots")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "cv_gene_categories",
        type=str,
        help="Input location of coefficient of variation gene categories text file",
    )
    parser.add_argument(
        "tau_gene_categories",
        type=str,
        help="Input location of tau tissue specific gene categories text file",
    )

    parser.add_argument(
        "GC_content_tsv",
        type=str,
        help="Input location of promoters GC_content tsv file",
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
        help="Optional replacement colour palette for cv categories",
        default=None,
        nargs="?",
    )
    parser.add_argument(
        "palette_tau",
        type=str,
        help="Optional replacement colour palette for tau categories",
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


def read_GC_file(GC_content_tsv):
    """read in GC file and make extra columns"""
    # read in GC content tsv
    GC_content = pd.read_table(GC_content_tsv, sep="\t", header=None)
    cols2 = ["name", "percentage_GC_content"]
    GC_content.columns = cols2

    # Make AGI column
    GC_content["AGI"] = GC_content.name.str.split(":", expand=True)[0]
    # make window number column
    # GC_content = GC_content.assign(window_number=GC_content.name.str.extract(r'_(.*?)\:'))
    # make chr column
    GC_content = GC_content.assign(
        chr=GC_content.name.str.split(":", n=3, expand=True)[2]
    )
    # make start column
    GC_content = GC_content.assign(
        start=GC_content.name.str.split(":", n=3, expand=True)[3].str.split(
            "-", expand=True
        )[0]
    )
    # make stop column
    GC_content = GC_content.assign(
        stop=GC_content.name.str.split(":", n=3, expand=True)[3].str.split(
            "-", expand=True
        )[1]
    )
    return GC_content


def mergeGC_genecategories(GC_content_df, gene_categories):
    """merged GC content df with gene categories"""
    # read in gene categories
    gene_cats = pd.read_csv(gene_categories, sep="\t", header=None)
    gene_cats.columns = ["AGI", "gene_type"]
    # merge to limit to genes of interest
    GC_content_categories = pd.merge(
        gene_cats, GC_content_df, how="left", on="AGI"
    )
    return GC_content_categories


def dunn_posthoc_test(df, dependent_variable, between):
    """dunn_posthoc tests with bonferroni multiple correction"""
    return sp.posthoc_dunn(
        df,
        val_col=dependent_variable,
        group_col=between,
        p_adjust="bonferroni",
    )


def make_plot(
    df_cv,
    df_tau,
    x_variable,
    y_variable,
    x_label,
    y_label,
    output_prefix,
    plot_kind,
    palette_cv,
    palette_tau,
    output_folder_name,
    dependent_variable,
    file_names,
):

    """function to make and save plot"""
    # allow colour codes in seaborn
    sns.set(color_codes=True)
    sns.set_style("ticks")
    # plot
    x = x_variable
    y = y_variable

    def equalise_samples_sizes(
        df, variable1_name, variable2_name, palette, categorisation_name
    ):
        # descriptive statistics
        def describe_stats(df, dependent_variable, between):
            """return descriptve statistics"""
            return df.groupby([between])[dependent_variable].describe()

        order = [variable1_name, variable2_name, "control"]

        # set colour palette
        colours = sns.color_palette(palette)

        # make copy of df
        merged2_unique = df.copy()

        # make sample sizes equal for comparison
        # identify sample size of the minimum category
        minimum_sample_size = merged2_unique.gene_type.value_counts().min()

        # print this
        print(f"sample size in each category cv = {minimum_sample_size}")

        # save sample size as file
        with open(
            f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/number_of_genes_in_each_category_{categorisation_name}_{dependent_variable}.txt",
            "w",
        ) as file:
            file.write(
                "number_of_genes_in_each_category=" + str(minimum_sample_size)
            )
        # descriptive stats

        # multiply this by the number of categories
        total_sample_size = minimum_sample_size * len(
            merged2_unique.gene_type.unique()
        )

        # select equal sample sizes of each category with a random state of 1 so it's reproducible
        equal_samplesizes = rep_sample(
            merged2_unique, "gene_type", total_sample_size, random_state=1
        )

        # now filter out genes which were not selected using the minimum sample size
        to_remove = merged2_unique[
            ~merged2_unique.AGI.isin(equal_samplesizes.AGI)
        ]
        df = df[~df.AGI.isin(to_remove.AGI)]

        describe = describe_stats(df, y, x)
        # save sample size as file
        with open(
            f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{dependent_variable}_descriptivestats_{categorisation_name}.txt",
            "w",
        ) as file:
            file.write(str(describe))

        return df, order, colours

    df_cv, order_cv, colours_cv = equalise_samples_sizes(
        df_cv, "constitutive", "variable", palette_cv, "cv"
    )
    df_tau, order_tau, colours_tau = equalise_samples_sizes(
        df_tau, "non-specific", "tissue_specific", palette_tau, "tau"
    )

    # if violin plot, don't extend past datapoints
    if plot_kind == "violin":
        # make subplots
        f, (ax1, ax2) = plt.subplots(
            1,
            2,
        )

        f.set_size_inches(10, 5)

        sns.violinplot(
            x=x,
            y=y,
            data=df_cv,
            order=order_cv,
            palette=colours_cv,
            cut=0,
            ax=ax1,
        )

        sns.violinplot(
            x=x,
            y=y,
            data=df_tau,
            order=order_tau,
            palette=colours_tau,
            cut=0,
            ax=ax2,
        )
    else:
        # make subplots
        f, (ax1, ax2) = plt.subplots(
            1,
            2,
        )
        f.set_size_inches(10, 5)

        sns.boxplot(
            x=x,
            y=y,
            data=df_cv,
            order=order_cv,
            palette=colours_cv,
            ax=ax1,
        )
        sns.boxplot(
            x=x,
            y=y,
            data=df_tau,
            order=order_tau,
            palette=colours_tau,
            ax=ax2,
        )
    # plot points
    sns.swarmplot(x=x, y=y, data=df_cv, color=".25", order=order_cv, ax=ax1)
    sns.swarmplot(x=x, y=y, data=df_tau, color=".25", order=order_tau, ax=ax2)

    # add significance if necessary - dunn's posthocs with multiple Bonferroni correction
    def add_stats(df, variable1, variable2, ax, order):

        stat = dunn_posthoc_test(df, y_variable, x_variable)
        # label box pairs
        box_pairs = [
            (variable1, variable2),
            (variable1, "control"),
            (variable2, "control"),
        ]

        # make empty list of p_values
        p_values = []
        # populate the list of p_values according to the box_pairs
        for pair in box_pairs:
            print(pair)
            # select p value for each pair
            p = stat.loc[pair[0], pair[1]]
            p_values.append(p)

        # add stats annotation to the plot
        add_stat_annotation(
            ax,
            # plot=plot_type,
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

    # run kruskall wallis test - if significant add significance annotations to plot
    def kruskal_test(df, dependent_variable, between):
        """Do Kruskal-Wallis analysis"""
        # Kruskal-Wallis one way analysis of variance
        return kruskal(data=df, dv=dependent_variable, between=between)

    # cv
    def statistics(df, variable1, variable2, ax, order, ranking):
        kruskal = kruskal_test(df, y, x)
        # add stats to plot if significant
        # make into df
        kruskal_df = pd.DataFrame(kruskal)
        # make column numeric
        kruskal_df = kruskal_df.astype({"p-unc": "float64"})
        # save kruskal table
        with open(
            f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{dependent_variable}_kruskal_{ranking}.txt",
            "w",
        ) as file:
            file.write(str(kruskal_df))

        if kruskal_df["p-unc"].iloc[0] < 0.05:
            add_stats(df, variable1, variable2, ax, order)

    statistics(df_cv, "constitutive", "variable", ax1, order_cv, "cv")
    statistics(
        df_tau, "non-specific", "tissue_specific", ax2, order_tau, "tau"
    )
    # change axes labels
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax2.set_xlabel(x_label)
    ax2.set_ylabel(y_label)
    # add graph titles
    ax1.set_title("A", x=0.02, fontsize=16)
    ax2.set_title("B", x=0.02, fontsize=16)

    # set yaxis range
    # set y axis range
    # first get y limits for both plots
    ax1ymin, ax1ymax = ax1.get_ylim()
    ax2ymin, ax2ymax = ax2.get_ylim()
    minimum = min(ax1ymin, ax2ymin)
    maximum = max(ax1ymax, ax2ymax)
    ax1.set_ylim(minimum, maximum)
    ax2.set_ylim(minimum, maximum)

    # tight layout
    plt.tight_layout()
    # save figure
    f.savefig(
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{output_prefix}_{plot_kind}.pdf",
        format="pdf",
        bbox_inches="tight",
    )


def main(args):
    # parse arguments
    args = parse_args(args)
    dependent_variable = "GC_content"

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{dependent_variable}/{args.output_folder_name}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{dependent_variable}/{args.output_folder_name}plots"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # read GC file
    GC_content_df = read_GC_file(args.GC_content_tsv)

    # merge with CV promoter categories
    GC_content_cv_gene_categories = mergeGC_genecategories(
        GC_content_df, args.cv_gene_categories
    )
    # merge with tau promoter categories
    GC_content_tau_gene_categories = mergeGC_genecategories(
        GC_content_df, args.tau_gene_categories
    )

    # Czechowski_gene_categories violin plot
    make_plot(
        GC_content_cv_gene_categories,
        GC_content_tau_gene_categories,
        "gene_type",
        "percentage_GC_content",
        "Gene type",
        "Mean % GC content",
        f"{dependent_variable}",
        "violin",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )

    # Czechowski_gene_categories box plot
    make_plot(
        GC_content_cv_gene_categories,
        GC_content_tau_gene_categories,
        "gene_type",
        "percentage_GC_content",
        "Gene type",
        "Mean % GC content",
        f"{dependent_variable}",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
