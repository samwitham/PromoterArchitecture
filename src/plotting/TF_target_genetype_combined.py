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
    """define arguments"""
    parser = argparse.ArgumentParser(description="TF_target_genetype")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "cv_tf_gene_categories",
        type=str,
        help="Input location of cv transcription factor gene categories",
    )
    parser.add_argument(
        "tau_tf_gene_categories",
        type=str,
        help="Input location of tau transcription factor gene categories",
    )
    parser.add_argument(
        "cv_all_genes_ranked",
        type=str,
        help="Input location of all ranked genes cv categories",
    )
    parser.add_argument(
        "tau_all_genes_ranked",
        type=str,
        help="Input location of all ranked genes tau categories",
    )
    parser.add_argument(
        "cv_promoter_TF",
        type=str,
        help="Input location of cv mapped TF-promoter target file",
    )
    parser.add_argument(
        "tau_promoter_TF",
        type=str,
        help="Input location of tau mapped TF-promoter target file",
    )
    parser.add_argument(
        "output_folder_name",
        type=str,
        help="Output folder name ending in forward slash",
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


def read_files(tf_gene_categories, promoter_TF):
    """read in the files as dfs"""
    # all_gene_categories_dfs = pd.read_table(all_gene_categories, sep='\t',header=0)
    promoter_TF_df = pd.read_table(promoter_TF, sep="\t", header=0)
    tf_gene_categories_df = pd.read_table(
        tf_gene_categories, sep="\t", header=None
    )
    cols = ["TF_AGI", "gene_type"]
    tf_gene_categories_df.columns = cols
    # merge the dfs
    merged = pd.merge(
        tf_gene_categories_df, promoter_TF_df, on="TF_AGI", how="left"
    )
    # remove NaN
    if "TF_AGI_allgenes" in merged.columns:
        merged_filtered = merged[merged.TF_AGI_allgenes.notna()]
    elif "AGI code" in merged.columns:
        merged_filtered = merged[merged["AGI code"].notna()]

    return tf_gene_categories_df, promoter_TF_df, merged_filtered


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


def dunn_posthoc_test(df, dependent_variable, between):
    """dunn_posthoc tests with bonferroni multiple correction"""
    return sp.posthoc_dunn(
        df,
        val_col=dependent_variable,
        group_col=between,
        p_adjust="bonferroni",
    )


def calculate_mean(df, dependent_variable):
    """calculate the mean coefficient of variation or tau of the promoters the TF binds to"""
    # group by TF and calculate mean for each TF
    means = df.groupby("TF_AGI")[dependent_variable].mean()
    # turn into a dataframe
    means_df = pd.DataFrame(means)
    # turn the index into a new column
    means_df.reset_index(level=0, inplace=True)
    # name columns
    cols = ["TF_AGI", f"mean_{dependent_variable}"]
    means_df.columns = cols

    # group by promoter and calculate SD (standard deviation) for each promoter
    # sd = df.groupby('promoter_AGI')['expression_CV'].std()
    # turn into a dataframe
    # sd_df = pd.DataFrame(sd)
    # turn the index into a new column
    # sd_df.reset_index(level=0, inplace=True)
    # name columns
    # cols = ['promoter_AGI', 'sd']
    # sd_df.columns = cols

    # merge the dfs
    # merged = pd.merge(means_df,sd_df)
    return means_df


def merge_genetype(df, gene_categories):
    """merge df with gene_categories file adding the genetype of the TFs (if in top 100 constitutive or top 100 variable TFs)"""
    gene_cats = pd.read_table(gene_categories, sep="\t", header=None)
    cols = ["TF_AGI", "gene_type"]
    gene_cats.columns = cols
    merged = pd.merge(gene_cats, df, on="TF_AGI", how="left")
    # drop NaN
    merged_filtered = merged.dropna()
    # reset index
    merged_filtered_index = merged_filtered.reset_index(drop=True)

    return merged_filtered_index


def make_plot(
    df_cv,
    df_tau,
    x_variable,
    y_variable,
    y_variable2,
    x_label,
    y_label,
    y_label2,
    output_prefix,
    plot_kind,
    palette_cv,
    palette_tau,
    output_folder_name,
    dependent_variable,
    file_names,
    mean=False,
):

    """function to make and save plot"""
    # allow colour codes in seaborn
    sns.set(color_codes=True)
    sns.set_style("ticks")
    # plot

    def equalise_samples_sizes(
        df,
        variable1_name,
        variable2_name,
        palette,
        categorisation_name,
        y_variable,
    ):
        # descriptive statistics
        def describe_stats(df, dependent_variable, between):
            """return descriptve statistics"""
            return df.groupby([between])[dependent_variable].describe()

        order = [variable1_name, variable2_name, "control"]

        # set colour palette
        colours = sns.color_palette(palette)

        # make copy of df
        # make copy of df
        if mean is False:
            merged2 = df.copy()

            merged2_unique = merged2.drop_duplicates(
                ["promoter_AGI"], keep="last"
            )
        elif mean is True:
            merged2_unique = df.copy()

        # make sample sizes equal for comparison
        # identify sample size of the minimum category
        minimum_sample_size = merged2_unique.gene_type.value_counts().min()

        # print this
        print(f"sample size in each category cv = {minimum_sample_size}")

        # save sample size as file
        with open(
            f"../../data/output/{file_names}/{output_folder_name}plots/number_of_genes_in_each_category_{categorisation_name}_{y_variable}.txt",
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
        to_remove = merged2_unique[
            ~merged2_unique.TF_AGI.isin(equal_samplesizes.TF_AGI)
        ]
        df = df[~df.TF_AGI.isin(to_remove.TF_AGI)]

        # descriptive stats
        describe = describe_stats(df, y_variable, x_variable)
        # save sample size as file
        with open(
            f"../../data/output/{file_names}/{output_folder_name}plots/{dependent_variable}_descriptivestats_{categorisation_name}.txt",
            "w",
        ) as file:
            file.write(str(describe))

        return df, order, colours

    df_cv, order_cv, colours_cv = equalise_samples_sizes(
        df_cv, "constitutive", "variable", palette_cv, "cv", y_variable
    )
    df_tau, order_tau, colours_tau = equalise_samples_sizes(
        df_tau,
        "non-specific",
        "tissue_specific",
        palette_tau,
        "tau",
        y_variable2,
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
            x=x_variable,
            y=y_variable,
            data=df_cv,
            order=order_cv,
            palette=colours_cv,
            cut=0,
            ax=ax1,
        )

        sns.violinplot(
            x=x_variable,
            y=y_variable2,
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
            x=x_variable,
            y=y_variable,
            data=df_cv,
            order=order_cv,
            palette=colours_cv,
            ax=ax1,
        )
        sns.boxplot(
            x=x_variable,
            y=y_variable2,
            data=df_tau,
            order=order_tau,
            palette=colours_tau,
            ax=ax2,
        )
    # plot points
    sns.swarmplot(
        x=x_variable,
        y=y_variable,
        data=df_cv,
        color=".25",
        order=order_cv,
        ax=ax1,
    )
    sns.swarmplot(
        x=x_variable,
        y=y_variable2,
        data=df_tau,
        color=".25",
        order=order_tau,
        ax=ax2,
    )

    # add significance if necessary - dunn's posthocs with multiple Bonferroni correction
    def add_stats(df, variable1, variable2, ax, order, y_variable):

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
            x=x_variable,
            y=y_variable,
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
    def kruskal_test(df, y_variable, between):
        """Do Kruskal-Wallis analysis"""
        # Kruskal-Wallis one way analysis of variance
        return kruskal(data=df, dv=y_variable, between=between)

    # cv
    def statistics(df, variable1, variable2, ax, order, y_variable, ranking):
        kruskal = kruskal_test(df, y_variable, x_variable)
        # add stats to plot if significant
        # make into df
        kruskal_df = pd.DataFrame(kruskal)
        # make column numeric
        kruskal_df = kruskal_df.astype({"p-unc": "float64"})

        # save kruskal table
        with open(
            f"../../data/output/{file_names}/{output_folder_name}plots/{dependent_variable}_kruskal_{ranking}.txt",
            "w",
        ) as file:
            file.write(str(kruskal_df))

        if kruskal_df["p-unc"].iloc[0] < 0.05:
            add_stats(df, variable1, variable2, ax, order, y_variable)

    statistics(
        df_cv, "constitutive", "variable", ax1, order_cv, y_variable, "cv"
    )
    statistics(
        df_tau,
        "non-specific",
        "tissue_specific",
        ax2,
        order_tau,
        y_variable2,
        "tau",
    )
    # change axes labels
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax2.set_xlabel(x_label)
    ax2.set_ylabel(y_label2)
    # add graph titles
    ax1.set_title("A", x=0.02, fontsize=16)
    ax2.set_title("B", x=0.02, fontsize=16)
    # tight layout
    plt.tight_layout()
    # save figure
    f.savefig(
        f"../../data/output/{file_names}/{output_folder_name}plots/{output_prefix}_{plot_kind}.pdf",
        format="pdf",
        bbox_inches="tight",
    )


def main(args):
    # parse arguments
    args = parse_args(args)
    dependent_variable = "TF_target_class"
    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{args.output_folder_name}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = (
        f"../../data/output/{args.file_names}/{args.output_folder_name}plots"
    )
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    cv_tf_gene_categories_dfs, cv_promoter_TF_df, cv_merged = read_files(
        args.cv_tf_gene_categories, args.cv_promoter_TF
    )
    (
        tau_tf_gene_categories_dfs,
        cv_promoter_TF_df,
        tau_merged,
    ) = read_files(args.tau_tf_gene_categories, args.tau_promoter_TF)

    # plot the CV for each promoter gene_type - whole promoter individual TF CVs
    # make_plot(merged,'gene_type', 'expression_CV','Gene type', 'Promoter target expression CV', 'Czechowski_CV', 'box',palette)

    # calculate mean CV/TAU of targets per TF
    cv_means = calculate_mean(cv_merged, "expression_CV")
    tau_means = calculate_mean(tau_merged, "TAU")

    # add TF gene_type categories
    cv_means_genetype = merge_genetype(cv_means, args.cv_tf_gene_categories)
    tau_means_genetype = merge_genetype(tau_means, args.tau_tf_gene_categories)

    # plot the mean CV for each promoter gene_type - whole promoter mean TF CVs

    make_plot(
        cv_means_genetype,
        tau_means_genetype,
        "gene_type",
        "mean_expression_CV",
        "mean_TAU",
        "Gene type",
        "Promoter target mean expression CV",
        "Promoter target mean tau",
        "Promoter_target_mean_category",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
        mean=True,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
