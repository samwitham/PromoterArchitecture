import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import seaborn as sns
import skbio
from pingouin import kruskal
from statannot import add_stat_annotation


# from sklearn.cluster import KMeans
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler
# import scipy.cluster.hierarchy as shc
def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="TF_diversity_plots")
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
        "mapped_motif_bed",
        type=str,
        help="Input location of promoters mapped motif bed",
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


def calculate_shannon_diversity(mapped_motif_bed):
    """read in mapped motifs_bed, calculate Shannon diversity"""

    # if whole promoter, mapped motif will have 13 columns
    # if shortened promoter, mapped motif file will have 24 (as bp overlap is needed in TF_diversity_plots_shortenedprom.py to remove TFBSs where the middle isn't in the promoter)
    # if 24 columns, only select the subset of 13 columns
    # if 13 columns, keep them all
    # This is to make the dfs have identical column names
    mapped_motifs = pd.read_table(mapped_motif_bed, sep="\t", header=None)
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

    df = mapped_motifs.copy()

    # count no. of each TF binding in each promoter
    groupby_promoter_counts = (
        df.groupby("promoter_AGI")["TF_AGI"]
        .value_counts()
        .unstack(fill_value=0)
    )
    # count no. of TF families binding in each promoter
    groupby_promoter_counts_family = (
        df.groupby("promoter_AGI")["TF_family"]
        .value_counts()
        .unstack(fill_value=0)
    )
    # Individual TF shannon diversity using arbitrary log2 base
    shannon_div_df = groupby_promoter_counts.apply(
        pd.Series(lambda x: skbio.diversity.alpha.shannon(x, base=2)), axis=1
    )
    # shannon diversity for TF family
    shannon_div_TF_family_df = groupby_promoter_counts_family.apply(
        pd.Series(lambda x: skbio.diversity.alpha.shannon(x, base=2)), axis=1
    )
    # convert rownames into column
    cols = ["promoter_AGI", "shannon"]
    shannon_div_df.index.name = "promoter_AGI"
    shannon_div_df.reset_index(inplace=True)
    shannon_div_TF_family_df.index.name = "promoter_AGI"
    shannon_div_TF_family_df.reset_index(inplace=True)
    # rename column
    shannon_div_df.rename(
        columns={"<lambda>": "Shannon_diversity_TF"}, inplace=True
    )
    shannon_div_TF_family_df.rename(
        columns={"<lambda>": "Shannon_diversity_TF_family"}, inplace=True
    )
    # merge individual TF and TF family diversity dfs
    diversity_df = pd.merge(
        shannon_div_df, shannon_div_TF_family_df, on="promoter_AGI"
    )

    # calculate unique TF counts
    # groupby promoter, and include only unique TFs within each promoter group. Preserve column names.
    unique_TF_count = df.groupby(by="promoter_AGI", as_index=False).agg(
        {"TF_AGI": pd.Series.nunique}
    )
    # rename column
    unique_TF_count.rename(columns={"TF_AGI": "unique_TF_count"}, inplace=True)

    # calculate total TF counts
    total_TF_count = df.groupby(by="promoter_AGI", as_index=False).agg(
        {"TF_AGI": pd.Series.count}
    )
    # rename column
    total_TF_count.rename(columns={"TF_AGI": "total_TF_count"}, inplace=True)

    # calculate total TF family counts
    total_TF_family_count = df.groupby(by="promoter_AGI", as_index=False).agg(
        {"TF_family": pd.Series.nunique}
    )
    # rename column
    total_TF_family_count.rename(
        columns={"TF_family": "TF_family_count"}, inplace=True
    )

    # merge diversity df with unique_TF_count
    diversity_df = pd.merge(diversity_df, unique_TF_count, on="promoter_AGI")
    # then merge with total_TF_count
    diversity_df = pd.merge(diversity_df, total_TF_count, on="promoter_AGI")
    # then merge with TF_family_count
    diversity_df = pd.merge(
        diversity_df, total_TF_family_count, on="promoter_AGI"
    )

    return diversity_df


def merge_shannon_genetype(shannon_df, gene_categories):
    """merge shannon diversity df with gene_categories file"""

    gene_cats = pd.read_table(gene_categories, sep="\t", header=None)
    cols = ["gene", "gene_type"]
    gene_cats.columns = cols
    merged = pd.merge(
        shannon_df, gene_cats, left_on="promoter_AGI", right_on="gene"
    )
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
    sns.set_style("white")
    # plot
    x = x_variable
    y = y_variable

    def equalise_samples_sizes(
        df, variable1_name, variable2_name, palette, categorisation_name
    ):
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
    def statistics(df, variable1, variable2, ax, order):
        kruskal = kruskal_test(df, y, x)
        # add stats to plot if significant
        # make into df
        kruskal_df = pd.DataFrame(kruskal)
        # make column numeric
        kruskal_df = kruskal_df.astype({"p-unc": "float64"})

        if kruskal_df["p-unc"].iloc[0] < 0.05:
            add_stats(df, variable1, variable2, ax, order)

    statistics(df_cv, "constitutive", "variable", ax1, order_cv)
    statistics(df_tau, "non-specific", "tissue_specific", ax2, order_tau)
    # change axes labels
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax2.set_xlabel(x_label)
    ax2.set_ylabel(y_label)
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
    plt.clf()


# def make_plot(
#     df,
#     x_variable,
#     y_variable,
#     x_label,
#     y_label,
#     output_prefix,
#     plot_kind,
#     palette,
#     file_names,
#     dependent_variable,
#     output_folder_name,
#     variable1_name,
#     variable2_name,
# ):
#     """function to make and save plot"""
#     # allow colour codes in seaborn
#     sns.set(color_codes=True)
#     sns.set_style("whitegrid")
#     # plot
#     x = x_variable
#     y = y_variable
#     order = [variable1_name, variable2_name, "control"]

#     # set colour palette
#     colours = sns.color_palette(palette)

#     # make copy of df
#     merged2_unique = df.copy()
#     # make sample sizes equal for comparison
#     # identify sample size of the minimum category
#     minimum_sample_size = merged2_unique.gene_type.value_counts().min()
#     # print this
#     print(f"sample size in each category = {minimum_sample_size}")
#     # save sample size as file
#     with open(
#         f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/number_of_genes_in_each_category_{y_variable}.txt",
#         "w",
#     ) as file:
#         file.write(
#             "number_of_genes_in_each_category=" + str(minimum_sample_size)
#         )

#         # multiply this by the number of categories
#     total_sample_size = minimum_sample_size * len(
#         merged2_unique.gene_type.unique()
#     )
#     # select equal sample sizes of each category with a random state of 1 so it's reproducible
#     equal_samplesizes = rep_sample(
#         merged2_unique, "gene_type", total_sample_size, random_state=1
#     )
#     # now filter out genes which were not selected using the minimum sample size
#     to_remove = merged2_unique[
#         ~merged2_unique.promoter_AGI.isin(equal_samplesizes.promoter_AGI)
#     ]
#     df = df[~df.promoter_AGI.isin(to_remove.promoter_AGI)]

#     # if violin plot don't extend past datapoints
#     if plot_kind == "violin":
#         sns.catplot(
#             x=x,
#             y=y,
#             data=df,
#             kind=plot_kind,
#             order=order,
#             palette=colours,
#             cut=0,
#         )
#     else:
#         sns.catplot(
#             x=x, y=y, data=df, kind=plot_kind, order=order, palette=colours
#         )
#     # plot points
#     ax = sns.swarmplot(x=x, y=y, data=df, color=".25", order=order)
#     # add significance if necessary - dunn's posthocs with multiple Bonferroni correction
#     stat = dunn_posthoc_test(df, y_variable, x_variable)
#     # label box pairs
#     box_pairs = [
#         (variable1_name, variable2_name),
#         (variable1_name, "control"),
#         (variable2_name, "control"),
#     ]
#     # make empty list of p_values
#     p_values = []
#     # populate the list of p_values accoridng to the box_pairs
#     for pair in box_pairs:
#         print(pair)
#         # select p value for each pair
#         p = stat.loc[pair[0], pair[1]]
#         p_values.append(p)

#     # add stats annotation to the plot
#     add_stat_annotation(
#         ax,
#         data=df,
#         x=x,
#         y=y,
#         order=order,
#         box_pairs=box_pairs,
#         text_format="star",
#         loc="outside",
#         verbose=2,
#         perform_stat_test=False,
#         pvalues=p_values,
#         test_short_name="Dunn",
#     )

#     # change axes labels
#     plt.ylabel(y_label)
#     plt.xlabel(x_label)
#     # tight layout
#     plt.tight_layout()
#     # save figure
#     ax.get_figure().savefig(
#         f"../../data/output/{file_names}/{dependent_variable}/{output_folder_name}plots/{output_prefix}_{plot_kind}.pdf",
#         format="pdf",
#     )
#     plt.clf()


# def run_PCA(mapped_motif_bed, shannon_Czechowski_gene_categories):
#     """perform a PCA"""
#     mapped_motifs = pd.read_table(mapped_motif_bed, sep="\t", header=None)
#     if len(mapped_motifs.columns) == 24:
#         cols = [
#             "chr",
#             "start",
#             "stop",
#             "promoter_AGI",
#             "dot1",
#             "strand",
#             "source",
#             "type",
#             "dot2",
#             "attributes",
#             "motif_chr",
#             "motif_start",
#             "motif_stop",
#             "name_rep",
#             "score",
#             "motif_strand",
#             "promoter_AGI2",
#             "p-value",
#             "q-value",
#             "matched_sequence",
#             "TF_name",
#             "TF_family",
#             "TF_AGI",
#             "bp_overlap",
#         ]
#         mapped_motifs.columns = cols
#         # filter columns
#         mapped_motifs = mapped_motifs[
#             [
#                 "motif_chr",
#                 "motif_start",
#                 "motif_stop",
#                 "name_rep",
#                 "score",
#                 "motif_strand",
#                 "promoter_AGI2",
#                 "p-value",
#                 "q-value",
#                 "matched_sequence",
#                 "TF_name",
#                 "TF_family",
#                 "TF_AGI",
#             ]
#         ]
#         # rename columns
#         cols = [
#             "chr",
#             "start",
#             "stop",
#             "name_rep",
#             "score",
#             "strand",
#             "promoter_AGI",
#             "p-value",
#             "q-value",
#             "matched_sequence",
#             "TF_name",
#             "TF_family",
#             "TF_AGI",
#         ]
#         mapped_motifs.columns = cols

#     elif len(mapped_motifs.columns) == 13:
#         cols = [
#             "chr",
#             "start",
#             "stop",
#             "name_rep",
#             "score",
#             "strand",
#             "promoter_AGI",
#             "p-value",
#             "q-value",
#             "matched_sequence",
#             "TF_name",
#             "TF_family",
#             "TF_AGI",
#         ]
#         mapped_motifs.columns = cols

#     elif len(mapped_motifs.columns) == 17:
#         cols = [
#             "chr",
#             "start",
#             "stop",
#             "name_rep",
#             "score",
#             "strand",
#             "promoter_AGI",
#             "p-value",
#             "q-value",
#             "matched_sequence",
#             "TF_name",
#             "TF_family",
#             "TF_AGI",
#             "chr_openchrom",
#             "start_openchrom",
#             "stop_openchrom",
#             "bp_overlap",
#         ]
#         mapped_motifs.columns = cols

#     df = mapped_motifs.copy()

#     # df.columns = cols
#     # count no. of TF families binding in each promoter
#     groupby_promoter_counts_family = (
#         df.groupby("promoter_AGI")["TF_family"]
#         .value_counts()
#         .unstack(fill_value=0)
#     )
#     # add gene type column
#     groupby_promoter_counts_family = pd.merge(
#         groupby_promoter_counts_family,
#         shannon_Czechowski_gene_categories[["promoter_AGI", "gene_type"]],
#         on="promoter_AGI",
#     )

#     # standardise the data - have to scale features before applying a PCA. Standardise to mean = 0, variance = 1
#     TF_Families = groupby_promoter_counts_family.columns.tolist()
#     # remove promoter_AGI and gene_type from column list
#     if "promoter_AGI" in TF_Families:
#         TF_Families.remove("promoter_AGI")
#     if "gene_type" in TF_Families:
#         TF_Families.remove("gene_type")
#     # separate out the families
#     x = groupby_promoter_counts_family.loc[:, TF_Families].values
#     # separate out the gene_type
#     # y = groupby_promoter_counts_family.loc[:, ["gene_type"]].values
#     # standardise the families
#     x = StandardScaler().fit_transform(x)
#     # run PCA, letting the algorithm decide on number of components such that 95% of the variation is maintained
#     # make instance of the model
#     pca = PCA(0.95)
#     # fit PCA to the data
#     principalComponents = pca.fit_transform(x)
#     # make into a dataframe
#     principalDf = pd.DataFrame(data=principalComponents)
#     # readd AGI and gene_type columns
#     finalDF_variable_promoters = pd.concat(
#         [
#             principalDf,
#             groupby_promoter_counts_family[["promoter_AGI", "gene_type"]],
#         ],
#         axis=1,
#     )
#     # calculate PCA variance
#     pca_variance = pca.explained_variance_ratio_

#     return finalDF_variable_promoters, pca_variance


# def hierarchical_clustering(PCA_df, file_names, output_folder_name):
#     """Run hierarchical clustering"""
#     # hierarchical clustering of PCA including 100 random genes

#     # separate out the families
#     x = PCA_df.drop(["gene_type"], axis=1)
#     x = x.set_index("promoter_AGI")
#     # separate out the gene_type
#     # y = PCA_df.loc[:, ["gene_type"]].values

#     # plot dendrograms to work out how many clusters to use
#     plt.figure(figsize=(10, 7))
#     shc.dendrogram(shc.linkage(x, method="ward"), leaf_rotation=45)
#     plt.gca()
#     # x_labels = ax.get_xmajorticklabels()
#     plt.savefig(
#         f"../../data/output/{file_names}/TF_diversity/{output_folder_name}plots/hierarchical_clustering_TF_family_counts.pdf"
#     )
#     # linkage matrix
#     z = shc.linkage(x, method="ward")

#     return x, z


# def elbow_method(z):
#     """run elbow method on hierachical clusters to decide how many clusters there are"""
#     # decide how many clusters there are
#     # try elbow method
#     last = z[-10:, 2]
#     last_rev = last[::-1]
#     idxs = np.arange(1, len(last) + 1)
#     plt.plot(idxs, last_rev)

#     acceleration = np.diff(last, 2)  # 2nd derivative of the distances
#     acceleration_rev = acceleration[::-1]
#     plt.plot(idxs[:-2] + 1, acceleration_rev)
#     plt.show()
#     k = (
#         acceleration_rev.argmax() + 2
#     )  # if idx 0 is the max of this we want 2 clusters
#     print("clusters:", k)
#     return k


# def kmeans_clustering(k, PCA_df, x_from_hierarchical_clustering):
#     """run kmeans clustering"""
#     PCA_kmeans = KMeans(n_clusters=k, random_state=0)
#     y_PCA_kmeans = PCA_kmeans.fit_predict(x_from_hierarchical_clustering)
#     PCA_df["Kmeans_PCA_cluster"] = y_PCA_kmeans

#     return PCA_df


# def plot_kmeans_clusters(
#     k,
#     PCA_df,
#     pca_variance,
#     file_names,
#     output_folder_name,
#     variable1_name,
#     variable2_name,
# ):
#     """make two subplots of the first 2 PCA components, the top subplot coloured by KMeans cluster, the bottom coloured by gene_type"""
#     # set seaborn graph background
#     sns.set(color_codes=True, font_scale=1)
#     sns.set_style("white")
#     # Create a figure instance, and the two subplots
#     fig = plt.figure(figsize=(6, 7))
#     ax1 = fig.add_subplot(211)
#     ax2 = fig.add_subplot(212)
#     # add custom palette size as sns doesnt like having numeric values for hue palette=sns.color_palette("Set1", 6)

#     sns.scatterplot(
#         x=0,
#         y=1,
#         hue="Kmeans_PCA_cluster",
#         data=PCA_df,
#         s=100,
#         palette=sns.color_palette("Set1", k),
#         ax=ax1,
#     )
#     sns.scatterplot(
#         x=0,
#         y=1,
#         hue="gene_type",
#         data=PCA_df,
#         s=100,
#         ax=ax2,
#         hue_order=[variable1_name, variable2_name, "control"],
#     )
#     # add graph titles
#     ax1.set(ylabel="", title="A")
#     ax2.set(xlabel="", ylabel="", title="B")
#     fig.tight_layout()

#     # Add axes labels
#     fig.text(
#         0.5,
#         0.01,
#         f"PC2 {(pca_variance[1]*100).round(1)}% of variance",
#         ha="center",
#         va="center",
#     )
#     fig.text(
#         0.0,
#         0.5,
#         f"PC1 {(pca_variance[0]*100).round(1)}% of variance",
#         ha="center",
#         va="center",
#         rotation="vertical",
#     )

#     fig.savefig(
#         f"../../data/output/{file_names}/TF_diversity/{output_folder_name}plots/PCA_Kmeans_TF_family_counts.pdf"
#     )
#     plt.clf()


def main(args):
    # parse arguments
    args = parse_args(args)

    dependent_variable = "TF_diversity"

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

    # make shannon df
    shannon_df = calculate_shannon_diversity(args.mapped_motif_bed)
    # merge shannon diversity df with Czechowski gene_categories file
    shannon_cv_gene_categories = merge_shannon_genetype(
        shannon_df, args.cv_gene_categories
    )
    shannon_tau_gene_categories = merge_shannon_genetype(
        shannon_df, args.tau_gene_categories
    )
    # # all promoter distribution plot - Shannon_diversity_TF
    # all_prom_distribution(
    #     shannon_df,
    #     "Shannon_diversity_TF",
    #     "TF Shannon diversity",
    #     "TF_diversity_shannon_allproms",
    #     args.file_names,
    #     dependent_variable,
    #     args.output_folder_name,
    # )
    # # all promoter distribution plot - Shannon_diversity_TF_family
    # all_prom_distribution(
    #     shannon_df,
    #     "Shannon_diversity_TF_family",
    #     "TF family Shannon diversity",
    #     "TFfamily_diversity_shannon_allproms",
    #     args.file_names,
    #     dependent_variable,
    #     args.output_folder_name,
    # )
    # # all promoter distribution plot - unique_TF_count
    # all_prom_distribution(
    #     shannon_df,
    #     "unique_TF_count",
    #     "unique TF count",
    #     "unique_TF_count_allproms",
    #     args.file_names,
    #     dependent_variable,
    #     args.output_folder_name,
    # )
    # # all promoter distribution plot - total_TF_count
    # all_prom_distribution(
    #     shannon_df,
    #     "total_TF_count",
    #     "total TF count",
    #     "total_TF_count_allproms",
    #     args.file_names,
    #     dependent_variable,
    #     args.output_folder_name,
    # )
    # # all promoter distribution plot - TF_family_count
    # all_prom_distribution(
    #     shannon_df,
    #     "TF_family_count",
    #     "TF family count",
    #     "TF_family_count_allproms",
    #     args.file_names,
    #     dependent_variable,
    #     args.output_folder_name,
    # )
    # Czechowski_gene_categories violin and boxplot

    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "Shannon_diversity_TF",
        "Gene type",
        "TF Shannon diversity",
        f"{dependent_variable}",
        "violin",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )
    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "Shannon_diversity_TF",
        "Gene type",
        "TF Shannon diversity",
        f"{dependent_variable}",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )
    # Czechowski_gene_categories violin and boxplot
    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "Shannon_diversity_TF_family",
        "Gene type",
        "TF family Shannon diversity",
        "TF_family_diversity",
        "violin",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )
    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "Shannon_diversity_TF_family",
        "Gene type",
        "TF family Shannon diversity",
        "TF_family_diversity",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )

    # Czechowski_gene_categories violin and boxplot
    # make_plot(shannon_Czechowski_gene_categories,'gene_type','unique_TF_count','Gene type','unique TF count', f'Czechowski_unique_TF_count', 'violin')
    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "unique_TF_count",
        "Gene type",
        "unique TF count",
        "unique_TF_count",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )
    # Czechowski_gene_categories violin and boxplot
    # make_plot(shannon_Czechowski_gene_categories,'gene_type','total_TF_count','Gene type','unique TF count', f'Czechowski_unique_TF_count', 'violin')
    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "total_TF_count",
        "Gene type",
        "total TF count",
        "total_TF_count",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )
    # Czechowski_gene_categories violin and boxplot
    # make_plot(shannon_Czechowski_gene_categories,'gene_type','TF_family_count','Gene type','TF family count', f'Czechowski_TF_family_count', 'violin')
    make_plot(
        shannon_cv_gene_categories,
        shannon_tau_gene_categories,
        "gene_type",
        "TF_family_count",
        "Gene type",
        "TF family count",
        "TF_family_count",
        "box",
        args.palette_cv,
        args.palette_tau,
        args.output_folder_name,
        dependent_variable,
        args.file_names,
    )

    # # Run PCA of TF family count
    # PCA_df_cv, pca_variance_cv = run_PCA(
    #     args.mapped_motif_bed, shannon_cv_gene_categories
    # )
    # PCA_df_tau, pca_variance_tau = run_PCA(
    #     args.mapped_motif_bed, shannon_tau_gene_categories
    # )

    # # Run hierarchical clustering
    # x_cv, z_cv = hierarchical_clustering(
    #     PCA_df_cv, args.file_names, args.output_folder_name
    # )
    # # decide how many clusters there are
    # # try elbow method
    # k = elbow_method(z)
    # # run kmeans clustering
    # kmeans_clustering(k, PCA_df, x)
    # # make two subplots of the first 2 PCA components, the top subplot coloured by KMeans cluster, the bottom coloured by gene_type
    # plot_kmeans_clusters(
    #     k,
    #     PCA_df,
    #     pca_variance,
    #     args.file_names,
    #     args.output_folder_name,
    #     args.variable1_name,
    #     args.variable2_name,
    # )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
