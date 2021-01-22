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
    parser = argparse.ArgumentParser(description="map_TF2CV")
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
        "cv_all_genes_ranking",
        type=str,
        help="Input location of all genes ranking cv categories",
    )
    parser.add_argument(
        "tau_all_genes_ranking",
        type=str,
        help="Input location of all genes ranking tau categories",
    )
    parser.add_argument(
        "mapped_motif_bed",
        type=str,
        help="Input location of mapped motifs bed file",
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


def map_cv(
    ranked_cvs_file,
    mapped_motifs_file,
    file_names,
    dependent_variable,
    output_folder_name,
):
    """function to map the CV value to the TFs which bind to each promoter"""
    # read in files
    cvs = pd.read_table(ranked_cvs_file, sep="\t", header=None)
    cols = [
        "rank",
        "probe_id",
        "TF_AGI",
        "expression_mean",
        "expression_SD",
        "expression_CV",
        "proportion_of_values_present_in_mas5",
        "presence_in_araport11",
        "constitutive_in_araport11",
    ]
    cvs.columns = cols
    # filter out any genes that aren't present in araport11
    cvs = cvs[cvs.presence_in_araport11 == 1]
    # read in mapped motifs
    mapped_motifs = pd.read_table(mapped_motifs_file, sep="\t", header=None)
    # if whole promoter, mapped motif will have 13 columns
    # if shortened promoter, mapped motif file will have 24 (as bp overlap is needed in TF_diversity_plots_shortenedprom.py to remove TFBSs where the middle isn't in the promoter)
    # if 24 columns, only select the subset of 13 columns
    # if 13 columns, keep them all
    # This is to make the dfs have identical column names
    if len(mapped_motifs.columns) == 24:
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
        mapped_motifs.columns = cols
        # filter columns
        mapped_motifs = mapped_motifs[
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
        mapped_motifs.columns = cols

    elif len(mapped_motifs.columns) == 13:
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
        mapped_motifs.columns = cols

    elif len(mapped_motifs.columns) == 17:
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
        mapped_motifs.columns = cols

    # merge CV df with mapped_motifs, adding the CVs to the respective AGIs
    merged = pd.merge(mapped_motifs, cvs, how="left", on="TF_AGI")
    merged_promotercv = pd.merge(
        mapped_motifs,
        cvs,
        how="left",
        left_on="promoter_AGI",
        right_on="TF_AGI",
        suffixes=["", "_allgenes"],
    )
    # Groupby promoter and then keep only unique TFs in each promoter
    # unique_CV_means = merged.groupby(['promoter_AGI', 'TF_AGI'])['expression_mean'].agg(lambda x: x.unique())
    unique_TFs = merged.drop_duplicates(
        ["promoter_AGI", "TF_AGI"]
    ).reset_index(drop=True)
    unique_TFs_promotercv = merged_promotercv.drop_duplicates(
        ["promoter_AGI", "TF_AGI"]
    ).reset_index(drop=True)
    # remove NaN
    unique_TFs_promotercv_filtered = unique_TFs_promotercv[
        unique_TFs_promotercv.TF_AGI_allgenes.notna()
    ]

    # save df as a file
    with open(
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}promoter_TF_CV.txt",
        "w",
    ) as f:
        unique_TFs_promotercv_filtered.to_csv(
            f, sep="\t", header=True, index=False
        )

    return merged, unique_TFs


def map_tau(
    ranked_cvs_file,
    mapped_motifs_file,
    file_names,
    dependent_variable,
    output_folder_name,
):
    """function to map the CV value to the TFs which bind to each promoter"""
    # read in files
    tau = pd.read_table(ranked_cvs_file, sep="\t", header=0)
    # make AGI uppercase
    tau["AGI code"] = tau["AGI code"].str.upper()
    # filter out any genes that aren't present in araport11
    # read in mapped motifs
    mapped_motifs = pd.read_table(mapped_motifs_file, sep="\t", header=None)
    # if whole promoter, mapped motif will have 13 columns
    # if shortened promoter, mapped motif file will have 24 (as bp overlap is needed in TF_diversity_plots_shortenedprom.py to remove TFBSs where the middle isn't in the promoter)
    # if 24 columns, only select the subset of 13 columns
    # if 13 columns, keep them all
    # This is to make the dfs have identical column names
    if len(mapped_motifs.columns) == 24:
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
        mapped_motifs.columns = cols
        # filter columns
        mapped_motifs = mapped_motifs[
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
        mapped_motifs.columns = cols

    elif len(mapped_motifs.columns) == 13:
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
        mapped_motifs.columns = cols

    elif len(mapped_motifs.columns) == 17:
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
        mapped_motifs.columns = cols
    # merge CV df with mapped_motifs, adding the CVs to the respective TF AGIs
    merged = pd.merge(
        mapped_motifs, tau, how="left", left_on="TF_AGI", right_on="AGI code"
    )
    merged_promotercv = pd.merge(
        mapped_motifs,
        tau,
        how="left",
        left_on="promoter_AGI",
        right_on="AGI code",
        suffixes=["", "_allgenes"],
    )
    # Groupby promoter and then keep only unique TFs in each promoter
    # unique_CV_means = merged.groupby(['promoter_AGI', 'TF_AGI'])['expression_mean'].agg(lambda x: x.unique())
    unique_TFs = merged.drop_duplicates(
        ["promoter_AGI", "TF_AGI"]
    ).reset_index(drop=True)
    unique_TFs_promotercv = merged_promotercv.drop_duplicates(
        ["promoter_AGI", "TF_AGI"]
    ).reset_index(drop=True)

    # remove NaN
    unique_TFs_promotercv_filtered = unique_TFs_promotercv[
        unique_TFs_promotercv.TF_AGI.notna()
    ]

    # save df as a file
    with open(
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}promoter_TF_TAU.txt",
        "w",
    ) as f:
        unique_TFs_promotercv_filtered.to_csv(
            f, sep="\t", header=True, index=False
        )

    return merged, unique_TFs


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


def merge_genetype(df, gene_categories):
    """merge df with gene_categories file adding the genetype of the promoters (if in top 100 constitutive or top 100 variable promoters)"""
    gene_cats = pd.read_table(gene_categories, sep="\t", header=None)
    cols = ["promoter_AGI", "gene_type"]
    gene_cats.columns = cols
    merged = pd.merge(gene_cats, df, on="promoter_AGI", how="left")
    # drop NaN
    merged_filtered = merged.dropna()
    # reset index
    merged_filtered_index = merged_filtered.reset_index(drop=True)

    return merged_filtered_index


def calculate_mean_SD_CV(df, ranking, mean_col_name):
    """calculate the mean coefficient of variation of the tFs binding to a promoter"""
    # group by promoter and calculate mean for each promoter
    means = df.groupby("promoter_AGI")[ranking].mean()
    # turn into a dataframe
    means_df = pd.DataFrame(means)
    # turn the index into a new column
    means_df.reset_index(level=0, inplace=True)
    # name columns
    cols = ["promoter_AGI", mean_col_name]
    means_df.columns = cols

    # group by promoter and calculate SD (standard deviation) for each promoter
    sd = df.groupby("promoter_AGI")[ranking].std()
    # turn into a dataframe
    sd_df = pd.DataFrame(sd)
    # turn the index into a new column
    sd_df.reset_index(level=0, inplace=True)
    # name columns
    cols = ["promoter_AGI", "sd"]
    sd_df.columns = cols

    # merge the dfs
    merged = pd.merge(means_df, sd_df)
    return merged


# def all_prom_distribution(
#     df,
#     x_variable,
#     x_label,
#     output_prefix,
#     file_names,
#     dependent_variable,
#     output_folder_name,
# ):
#     """function to return distribution plot of all promoters GC content"""

#     dist_plot = df[x_variable]
#     # create figure with no transparency
#     dist_plot_fig = sns.distplot(dist_plot).get_figure()
#     plt.xlabel(x_label)

#     # save to file
#     dist_plot_fig.savefig(
#         f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{output_prefix}_distribution.pdf",
#         format="pdf",
#     )
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
            f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/number_of_genes_in_each_category_{categorisation_name}_{y_variable}.txt",
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
            ~merged2_unique.promoter_AGI.isin(equal_samplesizes.promoter_AGI)
        ]
        df = df[~df.promoter_AGI.isin(to_remove.promoter_AGI)]

        # descriptive stats
        describe = describe_stats(df, y_variable, x_variable)
        # save sample size as file
        with open(
            f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{dependent_variable}_descriptivestats_{categorisation_name}_{y_variable}.txt",
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
            f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{dependent_variable}_kruskal_{ranking}_{y_variable}.txt",
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
        f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{output_prefix}_{plot_kind}.pdf",
        format="pdf",
        bbox_inches="tight",
    )


def main(args):
    # parse arguments
    args = parse_args(args)

    dependent_variable = "TFBS_TF_class"

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/{dependent_variable}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

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

    # map coefficient of variation (CV) values to each TF in the mapped_motifs file
    cv_merged, cv_unique_TF = map_cv(
        args.cv_all_genes_ranking,
        args.mapped_motif_bed,
        args.file_names,
        dependent_variable,
        args.output_folder_name,
    )
    # map tau tissue specificty values to each TF in the mapped_motifs file
    tau_merged, tau_unique_TF = map_tau(
        args.tau_all_genes_ranking,
        args.mapped_motif_bed,
        args.file_names,
        dependent_variable,
        args.output_folder_name,
    )

    # add gene_types for the promoters (eg. constitutive, variable or control)
    cv_genetypes = merge_genetype(cv_unique_TF, args.cv_gene_categories)
    tau_genetypes = merge_genetype(tau_unique_TF, args.tau_gene_categories)
    # calculate CV means per promoter df, ranking, mean_col_name
    cv_means_sd = calculate_mean_SD_CV(
        cv_genetypes, "expression_CV", "mean_cv"
    )
    tau_means_sd = calculate_mean_SD_CV(tau_genetypes, "TAU", "mean_tau")
    # merge dfs
    cv_means_sd_genetype = merge_genetype(cv_means_sd, args.cv_gene_categories)
    tau_means_sd_genetype = merge_genetype(
        tau_means_sd, args.tau_gene_categories
    )
    # remove promoters with no mean_cv or mean tau
    cv_means_sd_genetype = cv_means_sd_genetype[
        cv_means_sd_genetype.mean_cv.notnull()
    ]
    tau_means_sd_genetype = tau_means_sd_genetype[
        tau_means_sd_genetype.mean_tau.notnull()
    ]

    # check how many of each promoter type have mean_cv values
    # constitutive
    print("number of constitutive genes")
    print(
        len(
            cv_means_sd_genetype[
                cv_means_sd_genetype.gene_type == "constitutive"
            ]
        )
    )
    print("number of variable genes")
    print(
        len(cv_means_sd_genetype[cv_means_sd_genetype.gene_type == "variable"])
    )
    print("number of" + "control" + "genes")
    print(
        len(cv_means_sd_genetype[cv_means_sd_genetype.gene_type == "control"])
    )
    print("number of non-specific genes")
    print(
        len(
            tau_means_sd_genetype[
                tau_means_sd_genetype.gene_type == "non-specific"
            ]
        )
    )
    print("number of tissue_specific genes")
    print(
        len(
            tau_means_sd_genetype[
                tau_means_sd_genetype.gene_type == "tissue_specific"
            ]
        )
    )
    print("number of" + "control" + "genes")
    print(
        len(
            tau_means_sd_genetype[tau_means_sd_genetype.gene_type == "control"]
        )
    )

    # plot all promoter distribution of TF CV values
    # all_prom_distribution(
    #     Czechowski_merged,
    #     "expression_CV",
    #     "expression CV",
    #     "Czechowski_expressionCV",
    #     args.file_names,
    #     dependent_variable,
    #     args.output_folder_name,
    # )

    # plot the CV or tau for each promoter gene_type - whole promoter individual TF CVs/taus

    make_plot(
        cv_genetypes,
        tau_genetypes,
        "gene_type",
        "expression_CV",
        "TAU",
        "Gene type",
        "Cognate TF expression CV",
        "Cognate TF Tau tissue specificity",
        "TF_expression_categories",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )

    # plot the mean CV for each promoter gene_type - whole promoter mean TF CVs
    make_plot(
        cv_means_sd_genetype,
        tau_means_sd_genetype,
        "gene_type",
        "mean_cv",
        "mean_tau",
        "Gene type",
        "Mean cognate TF expression CV",
        "mean cognate TF Tau tissue specificity",
        "TF_expression_categories_mean",
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
