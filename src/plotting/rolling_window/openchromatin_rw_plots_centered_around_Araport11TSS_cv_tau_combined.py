import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args(args):
    parser = argparse.ArgumentParser(description="OpenChromatin_plots_rw")
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
        "GC_content_tsv", type=str, help="Input location of GC content tsv"
    )
    parser.add_argument(
        "EPD_TSS_bed",
        type=str,
        help="Input location of eukaryotic promoter database transcription start site bed file",
    )
    parser.add_argument(
        "promoter_bed", type=str, help="Input location of promoter bed file"
    )
    parser.add_argument(
        "promoter_5UTR_bed",
        type=str,
        help="Input location of promoter_5UTR bed file",
    )
    parser.add_argument(
        "foldername_prefix", type=str, help="Output folder name prefix to use"
    )
    parser.add_argument(
        "root_chrom_bp_covered",
        type=str,
        help="Input location of root chromatin bed file",
    )
    parser.add_argument(
        "shoot_chrom_bp_covered",
        type=str,
        help="Input location of shoot chromatin bed file",
    )
    parser.add_argument(
        "rootshootintersect_chrom_bp_covered",
        type=str,
        help="Input location of rootshootintersect chromatin bed file",
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


def process_input_files(GC_content_tsv, Czechowski_gene_categories):
    """process and merge the input files into a df"""
    promoters = pd.read_csv(Czechowski_gene_categories, sep="\t", header=None)
    cols = ["AGI", "gene_type"]
    promoters.columns = cols
    # read in GC content table just to get the window locations (GC content not wanted, just open chromatin)
    GC_content = pd.read_table(GC_content_tsv, sep="\t", header=None)
    cols2 = ["name", "percentage_GC_content"]
    GC_content.columns = cols2
    # Make AGI column
    GC_content["AGI"] = GC_content.name.str.split("_", expand=True)[0]
    # make window number column
    GC_content = GC_content.assign(
        window_number=GC_content.name.str.extract(r"_(.*?)\:")
    )
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
    # make df columns integars
    GC_content = GC_content.astype(
        {"stop": "int", "start": "int", "chr": "int"}
    )
    # add window length column
    GC_content = GC_content.assign(
        window_length=GC_content.stop - GC_content.start
    )

    # merge to limit to genes of interest
    GC_content = pd.merge(promoters, GC_content, how="left", on="AGI")

    # remove windows with fewer than 100 promoters extending to that location
    GC_content = GC_content[
        GC_content["window_number"].map(
            GC_content["window_number"].value_counts()
        )
        > 99
    ]

    # # allow colour codes in seaborn
    # sns.set(color_codes=True)
    # sns.set_style("ticks")
    # sns.set_palette(args.palette)

    return GC_content


def add_coverage(df, coverage_bed, suffix):
    """add % bp covered data from a bed file to the df. Prefix is a name added to any new columns"""
    # read in bed file
    coverage_df = pd.read_table(coverage_bed, sep="\t", header=None)
    cols = [
        "chr",
        "start",
        "stop",
        "name",
        "no._of_overlaps",
        "no._of_bases_covered",
        "window_length",
        "fraction_bases_covered",
    ]
    coverage_df.columns = cols
    # add % bases covered column
    coverage_df["percentage_bases_covered"] = (
        coverage_df["fraction_bases_covered"] * 100
    )
    # filter columns
    coverage_df = coverage_df[
        ["chr", "start", "stop", "name", "percentage_bases_covered"]
    ]
    # make df columns integars
    df = df.astype({"stop": "int", "start": "int", "chr": "int"})
    coverage_df = coverage_df.astype(
        {"stop": "int", "start": "int", "chr": "int"}
    )
    # merge the dfs
    merged = pd.merge(
        df,
        coverage_df,
        how="left",
        on=["chr", "start", "stop"],
        suffixes=("", f"_{suffix}"),
    )
    # remove NaN
    # merged = merged[merged['name'].notnull()]
    return merged


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


def windows_coords(
    file_names,
    foldername_prefix,
    output_prefix,
    variable_of_interest,
    variable_of_interest_name,
    variable_of_interest_df_cv,
    variable_of_interest_df_tau,
    promoter_bed,
    promoter_5UTR_bed,
    window_offset,
    EPD_TSS_bed,
    palette_cv,
    palette_tau,
    includeEPDTSS=False,
    includechromatin=False,
    chromatin_tissue_variable="percentage_bases_covered_rootshootintersect_chrom",
    chromatin_tissue_variable_name="% open chromatin root and shoot intersect",
    x_range=False,
    estimator="median",
    ci=95,
    n_boot=10000,
    genetype_cv=False,
    genetype_tau=False,
    genetype2_cv=False,
    genetype2_tau=False,
    genetype3_cv=False,
    genetype3_tau=False,
):
    """function to add the centre of each window corresponding to each window no. and return a lineplot. Also add promoter length distribution, Araport TSS distribution,
    EPD TSS distribution (add the most common TSS as documented on eukaryotic promoter database Arabidopsis last modified on EPD 06/06/2018)"""
    # function to equalise samples sizes
    def equalise_samples_sizes(
        df, categorisation_name, palette, genetype, genetype2, genetype3=False
    ):
        # make a subselection of categories so all sample sizes are equal
        # first select only the relevant genetypes
        # set colour palette
        colours = sns.color_palette(palette)

        if genetype3 is False:
            df = df[df.gene_type.isin([genetype, genetype2])]
        else:
            df = df[df.gene_type.isin([genetype, genetype2, genetype3])]

        # make each promoter unique
        df_unique = df.drop_duplicates("AGI")
        # identify sample size of the minimum category
        minimum_sample_size = df_unique.gene_type.value_counts().min()
        # print this
        print(f"sample size in each category = {minimum_sample_size}")
        # save sample size as file
        with open(
            f"../../data/output/{file_names}/rolling_window/{foldername_prefix}/plots/number_of_genes_in_each_category_{categorisation_name}.txt",
            "w",
        ) as file:
            file.write(
                "number_of_genes_in_each_category=" + str(minimum_sample_size)
            )

        # multiply this by the number of categories
        total_sample_size = minimum_sample_size * len(
            df_unique.gene_type.unique()
        )
        # select equal sample sizes of each category with a random state of 1 so it's reproducible
        equal_samplesizes = rep_sample(
            df_unique, "gene_type", total_sample_size, random_state=1
        )
        # now filter out genes which were not selected using the minimum sample size
        to_remove = df_unique[~df_unique.AGI.isin(equal_samplesizes.AGI)]
        df = df[~df.AGI.isin(to_remove.AGI)]
        return df, colours

    # read in bed file
    promoter_df = pd.read_table(promoter_bed, sep="\t", header=None)
    col = [
        "chr",
        "start",
        "stop",
        "AGI",
        "dot1",
        "strand",
        "source",
        "type",
        "dot2",
        "attributes",
    ]
    promoter_df.columns = col

    # merge promoter_bed with variable_of_interest_df on AGI
    promoter_5UTR_df = pd.read_table(promoter_5UTR_bed, sep="\t", header=None)
    promoter_5UTR_df.columns = col

    # add promoter length column
    promoter_df["promoter_length"] = promoter_df.stop - promoter_df.start

    # temporarily merge promoter_df with promoter_5UTR_bed
    temp_merged = pd.merge(
        promoter_df,
        promoter_5UTR_df,
        how="left",
        on="AGI",
        suffixes=("", "_promsUTR"),
    )
    # add 5UTR length column
    temp_merged["five_UTR_length"] = (
        temp_merged.stop_promsUTR - temp_merged.start_promsUTR
    ) - temp_merged.promoter_length
    # filter columns
    temp_merged = temp_merged[
        [
            "chr",
            "start",
            "stop",
            "AGI",
            "dot1",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
            "promoter_length",
            "five_UTR_length",
        ]
    ]
    # rename temp_merged back to promoter_df
    promoter_df = temp_merged.copy()

    def sort_df(
        variable_of_interest_df,
        promoter_df,
        genetypes_used,
    ):
        merged = pd.merge(
            variable_of_interest_df,
            promoter_df,
            on="AGI",
            how="left",
            suffixes=("", "_wholeprom"),
        )
        # remove NaN
        merged = merged[merged[variable_of_interest].notnull()]
        # make columns integars
        merged = merged.astype(
            {
                "stop_wholeprom": "int",
                "start_wholeprom": "int",
                "start": "int",
                "stop": "int",
            }
        )
        # split merged into 2 dfs by strand
        pos = merged[merged.strand == "+"].copy()
        neg = merged[merged.strand == "-"].copy()

        # add variable of interest position column where position is the middle of the window using the Araport TSS (end of promoter bed file) as a reference
        # this will lead to positive positions being in the 5'UTR and negative in the promoter region
        pos["position"] = (pos.stop_wholeprom) - (
            pos.start + 0.5 * (pos.stop - pos.start)
        )
        neg["position"] = (
            neg.start + 0.5 * (neg.stop - neg.start)
        ) - neg.start_wholeprom
        merged2 = pd.merge(pos, neg, how="outer")

        merged2 = merged2.astype({"position": "int64"})

        # make window number an integar
        variable_of_interest_df = variable_of_interest_df.astype(
            {"window_number": "float"}
        )

        # calculate promoter and 5UTR window length based on window cutoff
        # number_of_windows = len(variable_of_interest_df.window_number.unique())
        # promoter window number plus 1 because window 0 is excluded
        # promoter_window_number = (
        #     len(
        #         variable_of_interest_df[
        #             variable_of_interest_df.window_number < 0
        #         ].window_number.unique()
        #     )
        #     + 1
        # )
        # 5UTR window number plus 1
        # five_UTR_window_number = (
        #     len(
        #         variable_of_interest_df[
        #             variable_of_interest_df.window_number > 0
        #         ].window_number.unique()
        #     )
        #     + 1
        # )

        # max_promoter_length
        # window_length = variable_of_interest_df.window_length.max()
        # max_promoter_length = promoter_window_number * (
        #     window_length - window_offset
        # )
        # max_5UTR_length = five_UTR_window_number * (
        #     window_length - window_offset
        # )

        # make integars
        merged2 = merged2.astype(
            {
                f"{variable_of_interest}": "float64",
                f"{chromatin_tissue_variable}": "float64",
            }
        )

        # remove NaN (promoters in promoters.bed but not in promoters_5UTR.bed)
        merged2 = merged2[merged2.promoter_length.notnull()]
        # remove NaN (promoters in promoters_5UTR but not in promoters.bed - ie. only 5'UTRs)
        merged2 = merged2[merged2.chr.notnull()]

        # merged2.start_no_UTR = merged2.start_no_UTR - 1
        # # add Araport TSS location column
        # # merged2['TSS'] = int()
        # merged2.loc[merged2.strand == "+", "TSS"] = merged2.loc[
        #     merged2.strand == "+", "stop"
        # ]
        # merged2.loc[merged2.strand == "-", "TSS"] = (
        #     merged2.loc[merged2.strand == "-", "start"] - 1
        # )
        # # transform TSS location in the same way as the position column
        # merged2.loc[merged2.strand == "-", "TSS_transformed_Araport11"] = (
        #     merged2.loc[merged2.strand == "-", "TSS"]
        #     - merged2.loc[merged2.strand == "-", "start_wholeprom"]
        # )
        # merged2.loc[merged2.strand == "+", "TSS_transformed_Araport11"] = (
        #     merged2.loc[merged2.strand == "+", "stop_wholeprom"]
        #     - merged2.loc[merged2.strand == "+", "TSS"]
        # )

        # # make integars
        # merged2 = merged2.astype(
        #     {
        #         "start_no_UTR": "float64",
        #         "stop_no_UTR": "float64",
        #         "TSS": "float64",
        #         "TSS_transformed_Araport11": "float64",
        #         f"{variable_of_interest}": "float64",
        #         f"{chromatin_tissue_variable}": "float64",
        #     }
        # )
        # return merged2[['AGI','strand','start','stop','start_wholeprom','stop_wholeprom','start_no_UTR','stop_no_UTR','TSS','TSS_transformed','position','chr_no_UTR','window_number']]

        # set number of subplots so can easily change all output possibilities, where subplotA is the top
        # subplots = 2

        # make subplots
        #     if includeEPDTSS == True:
        #         subplots = subplots + 1
        #         f, axes = plt.subplots(subplots, figsize=(10,10))
        #         OpenChromplot = axes[subplots-subplots]
        #         Araport11TSSplot = axes[subplots-(subplots-1)]
        #         EPDTSSplot = axes[subplots-(subplots-2)]
        #         #promlengthsplot = axes[subplots-(subplots-3)]
        #         variableofinterestplot = axes[subplots-(subplots-3)]
        #     else:
        #         f, axes = plt.subplots(subplots, figsize=(10,8))
        #         OpenChromplot = axes[subplots-subplots]
        #         Araport11TSSplot = axes[subplots-(subplots-1)]
        #         #promlengthsplot = axes[subplots-(subplots-2)]
        #         variableofinterestplot = axes[subplots-(subplots-2)]

        # check the plot axes variables are there. If they are not, assign None to them
        # try:
        #     OpenChromplot
        # except NameError:
        #     OpenChromplot = None
        # try:
        #     Araport11TSSplot
        # except NameError:
        #     Araport11TSSplot = None
        # try:
        #     EPDTSSplot
        # except NameError:
        #     EPDTSSplot = None
        #     try:
        #         promlengthsplot
        #     except NameError:
        #         promlengthsplot = None
        # try:
        #     variableofinterestplot
        # except NameError:
        #     variableofinterestplot = None

        # If EPD TSS plot is present, filter promoters which aren't in EPD to remove NaNs
        # if EPDTSSplot is not None:
        #     # remove NaN (promoters in promoters_5UTR but not in promoters.gff - ie. only 5'UTRs)
        #     merged2 = merged2[merged2.TSS_transformed_EPD.notnull()]

        if genetypes_used is not False:
            # filter so only genetype subset present
            merged2 = merged2[merged2.gene_type.notnull()]
            # remove windows with fewer than 50 promoters extending to that location if looking at specific genetypes
            merged2 = merged2[
                merged2["window_number"].map(
                    merged2["window_number"].value_counts()
                )
                > 49
            ]
            # calculate promoter and 5UTR window length based on window cutoff
            # number_of_windows = len(merged2.window_number.unique())
            # promoter window number plus 1 because window 0 is excluded
            # promoter_window_number = (
            #     len(merged2[merged2.window_number < 0].window_number.unique())
            #     + 1
            # )
            # # 5UTR window number plus 1
            # five_UTR_window_number = (
            #     len(merged2[merged2.window_number > 0].window_number.unique())
            #     + 1
            # )
            # redefine max_promoter_length
            # window_length = merged2.window_length.max()
            # max_promoter_length = promoter_window_number * (
            #     window_length - window_offset
            # )
            # max_5UTR_length = five_UTR_window_number * (
            #     window_length - window_offset
            # )

        # make all values of interest negative as upstream from ATG
        # merged_positive = merged2.copy()
        merged2[["promoter_length", "position"]] = -merged2[
            ["promoter_length", "position"]
        ]
        return merged2

    df_cv = sort_df(
        variable_of_interest_df_cv,
        promoter_df,
        genetypes_used=genetype_cv,
    )
    df_tau = sort_df(
        variable_of_interest_df_tau,
        promoter_df,
        genetypes_used=genetype_tau,
    )
    # allow colour codes in seaborn
    sns.set(color_codes=True)
    sns.set_style("white")
    # set size
    plt.figure(figsize=(10, 5))

    # change estimator
    if estimator == "mean":
        new_estimator = estimator

    if estimator == "median":
        new_estimator = np.median

    if genetype_cv is False:

        # length_of_longest_promoter = merged_positive.length.max()
        # if openchromplot variable present, add that plot

        # if variableofinterestplot variable present, add that plot
        # variable of interest lineplot
        # make subplot
        sns.set_palette(sns.color_palette(palette_cv))
        f1 = plt.subplot(1, 2, 1)

        # f.set_size_inches(10, 5)
        sns.lineplot(
            y=df_cv[chromatin_tissue_variable],
            x=df_cv.position,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
        )
        # make subplot
        sns.set_palette(sns.color_palette(palette_tau))
        f2 = plt.subplot(1, 2, 2)

        # f.set_size_inches(10, 5)
        sns.lineplot(
            y=df_tau[chromatin_tissue_variable],
            x=df_tau.position,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
        )

    elif genetype2_cv is False:
        # filter so only genetype subset present
        df_cv = df_cv[df_cv.gene_type.notnull()]
        df_tau = df_tau[df_tau.gene_type.notnull()]

        # if openchromplot variable present, add that plot

        # if variableofinterestplot variable present, add that plot
        # variable of interest lineplot
        # make subplot
        sns.set_palette(sns.color_palette(palette_cv))
        f1 = plt.subplot(
            1,
            2,
            1,
        )

        # f.set_size_inches(10, 5)
        sns.lineplot(
            y=df_cv[df_cv.gene_type == genetype_cv][chromatin_tissue_variable],
            x=df_cv[df_cv.gene_type == genetype_cv].position,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
        )
        # make subplot
        sns.set_palette(sns.color_palette(palette_tau))
        f2 = plt.subplot(
            1,
            2,
            2,
        )

        # f.set_size_inches(10, 5)
        sns.lineplot(
            y=df_tau[df_tau.gene_type == genetype_tau][
                chromatin_tissue_variable
            ],
            x=df_tau[df_tau.gene_type == genetype_tau].position,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
        )

        # set y axis as maximum mean window % bp covered of all genetype subset
        # variableofinterestplot.set_ylim([0,merged2.groupby('window_number')[variable_of_interest].median().max()+20])
        # set x axis range if specified
    #         if x_range==False:
    #             pass
    #         else:
    #             length_of_longest_promoter = x_range

    #         #for all subplots:
    #         for n in axes:
    #             #remove grids
    #             n.grid(False)
    #             n.set_xlim([-length_of_longest_promoter,0])
    #         f.tight_layout()

    elif genetype3_cv is False:
        # filter so only genetype subset present
        df_cv = df_cv[df_cv.gene_type.notnull()]
        df_tau = df_tau[df_tau.gene_type.notnull()]

        # equalise samples sizes
        df_cv, colours_cv = equalise_samples_sizes(
            df_cv, "cv", palette_cv, "constitutive", "variable"
        )
        df_tau, colours_tau = equalise_samples_sizes(
            df_tau, "tau", palette_tau, "non-specific", "tissue_specific"
        )

        # if variableofinterestplot variable present, add that plot
        # lineplot variable of interest
        # make subplot
        sns.set_palette(sns.color_palette(palette_cv))
        f1 = plt.subplot(
            1,
            2,
            1,
        )
        sns.lineplot(
            y=df_cv[df_cv.gene_type == genetype_cv][chromatin_tissue_variable],
            x=df_cv[df_cv.gene_type == genetype_cv].position,
            label=genetype_cv,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_cv,
        )
        sns.lineplot(
            y=df_cv[df_cv.gene_type == genetype2_cv][
                chromatin_tissue_variable
            ],
            x=df_cv[df_cv.gene_type == genetype2_cv].position,
            label=genetype2_cv,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_cv,
        )

        # make subplot
        sns.set_palette(sns.color_palette(palette_tau))
        f2 = plt.subplot(
            1,
            2,
            2,
        )
        sns.set_palette(palette_tau)
        sns.lineplot(
            y=df_tau[df_tau.gene_type == genetype_tau][
                chromatin_tissue_variable
            ],
            x=df_tau[df_tau.gene_type == genetype_tau].position,
            label=genetype_tau,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_tau,
            # ax=ax2,
        )
        sns.lineplot(
            y=df_tau[df_tau.gene_type == genetype2_tau][
                chromatin_tissue_variable
            ],
            x=df_tau[df_tau.gene_type == genetype2_tau].position,
            label=genetype2_tau,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_tau,
            # ax=ax2,
        )
        # Create the legend
        f1.legend()
        f2.legend()

    #         if x_range==False:
    #             pass
    #         else:
    #             length_of_longest_promoter = x_range

    #         #for all subplots:
    #         for n in axes:
    #             #remove grids
    #             n.grid(False)
    #             n.set_xlim([-length_of_longest_promoter,0])
    #         f.tight_layout()
    else:
        # filter so only genetype subset present
        df_cv = df_cv[df_cv.gene_type.notnull()]
        df_tau = df_tau[df_tau.gene_type.notnull()]

        # equalise samples sizes
        df_cv, colours_cv = equalise_samples_sizes(
            df_cv, "cv", palette_cv, "constitutive", "variable", "control"
        )
        df_tau, colours_tau = equalise_samples_sizes(
            df_tau,
            "tau",
            palette_tau,
            "non-specific",
            "tissue_specific",
            "control",
        )

        # if variableofinterestplot variable present, add that plot
        # lineplot
        # make subplot
        sns.set_palette(sns.color_palette(palette_cv))
        f1 = plt.subplot(
            1,
            2,
            1,
        )
        sns.lineplot(
            y=df_cv[df_cv.gene_type == genetype_cv][chromatin_tissue_variable],
            x=df_cv[df_cv.gene_type == genetype_cv].position,
            label=genetype_cv,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_cv,
        )
        sns.lineplot(
            y=df_cv[df_cv.gene_type == genetype2_cv][
                chromatin_tissue_variable
            ],
            x=df_cv[df_cv.gene_type == genetype2_cv].position,
            label=genetype2_cv,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_cv,
        )
        sns.lineplot(
            y=df_cv[df_cv.gene_type == genetype3_cv][
                chromatin_tissue_variable
            ],
            x=df_cv[df_cv.gene_type == genetype3_cv].position,
            label=genetype3_cv,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_cv,
        )
        # make subplot
        sns.set_palette(sns.color_palette(palette_tau))
        f2 = plt.subplot(
            1,
            2,
            2,
        )

        sns.lineplot(
            y=df_tau[df_tau.gene_type == genetype_tau][
                chromatin_tissue_variable
            ],
            x=df_tau[df_tau.gene_type == genetype_tau].position,
            label=genetype_tau,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_tau,
        )
        sns.lineplot(
            y=df_tau[df_tau.gene_type == genetype2_tau][
                chromatin_tissue_variable
            ],
            x=df_tau[df_tau.gene_type == genetype2_tau].position,
            label=genetype2_tau,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_tau,
        )
        sns.lineplot(
            y=df_tau[df_tau.gene_type == genetype3_tau][
                chromatin_tissue_variable
            ],
            x=df_tau[df_tau.gene_type == genetype3_tau].position,
            label=genetype3_tau,
            estimator=new_estimator,
            ci=ci,
            n_boot=n_boot,
            palette=colours_tau,
        )
        # Create the legend
        f1.legend()
        f2.legend()

        # set x axis length
    #         if x_range==False:
    #             pass
    #         else:
    #             length_of_longest_promoter = x_range

    #         #for all subplots:
    #         for n in axes:
    #             #remove grids
    #             n.grid(False)
    #             n.set_xlim([-length_of_longest_promoter,0])
    #             leg = n.legend()
    #         f.tight_layout()
    # set x axis range if specified
    # if x_range is False:
    #     pass
    # else:
    #     length_of_longest_promoter = x_range
    # set titles & axes names
    f1.set_ylabel(f"{estimator} {variable_of_interest_name}")
    f1.set_xlabel("position relative to Araport 11 TSS")
    f2.set_ylabel(f"{estimator} {variable_of_interest_name}")
    f2.set_xlabel("position relative to Araport 11 TSS")

    # add graph titles
    f1.set_title("A", x=0.02, fontsize=16)
    f2.set_title("B", x=0.02, fontsize=16)
    # set x_axis range
    # if genetype_cv is not False:
    #     f1.set_xlim([1000, 700])
    #     f2.set_xlim([1000, 700])

    # remove grids
    # plt.grid(False)
    # plt.set_xlim([(-length_of_longest_promoter-50),0])
    # set a tight layout
    plt.tight_layout()

    # save figure
    plt.savefig(
        f"../../data/output/{file_names}/rolling_window/{foldername_prefix}/plots/{output_prefix}_{chromatin_tissue_variable}_{estimator}_sliding_window.pdf",
        format="pdf",
        bbox_inches="tight",
    )
    # remove plot
    plt.clf()


# def plot_length(
#     df, output_prefix, genetype=False, genetype2=False, genetype3=False
# ):
#     """function to plot length distribution of promoter"""
#     # make lengths positive by squaring and then square rooting
#     df.length = (df.length ** 2) ** (1 / 2)

#     if genetype is False:
#         dist_plot = df["length"]
#         # create figure with no transparency
#         sns.distplot(dist_plot, axlabel="length (bp)").get_figure()

#     elif genetype2 is False:
#         sns.distplot(
#             df[df.gene_type == genetype].length,
#             label=genetype,
#             axlabel="length (bp)",
#         )

#     elif genetype3 is False:
#         sns.distplot(
#             df[df.gene_type == genetype].length, hist=None, label=genetype
#         )
#         sns.distplot(
#             df[df.gene_type == genetype2].length,
#             hist=None,
#             label=genetype2,
#             axlabel="length (bp)",
#         ).get_figure()
#         plt.legend()
#     else:
#         sns.distplot(
#             df[df.gene_type == genetype].length, hist=None, label=genetype
#         )
#         sns.distplot(
#             df[df.gene_type == genetype2].length, hist=None, label=genetype2
#         )
#         sns.distplot(
#             df[df.gene_type == genetype3].length,
#             hist=None,
#             label=genetype3,
#             axlabel="length (bp)",
#         )
#         # plt.axlabel='length (bp)'
#         plt.legend()
#     # tight layout
#     plt.tight_layout()
#     # save figure
#     plt.savefig(
#         f"../../data/output/{args.file_names}/rolling_window/{args.foldername_prefix}/plots/{output_prefix}_promoter_lengths.pdf",
#         format="pdf",
#     )
#     # remove plot
#     plt.clf()


def add_genetype(df, gene_categories):
    """function to add gene type to the df, and remove random genes"""

    select_genes = pd.read_table(gene_categories, sep="\t", header=None)
    cols = ["AGI", "gene_type"]
    select_genes.columns = cols

    merged = pd.merge(select_genes, df, on="AGI", how="left")
    # no_random = merged_renamed[merged_renamed.gene_type != 'random']
    #  no_random.reset_index(drop=True, inplace=True)

    return merged


def main(args):
    # parse arguments
    args = parse_args(args)

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/rolling_window/"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/rolling_window/{args.foldername_prefix}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # make directory for the plots to be exported to
    dirName = f"../../data/output/{args.file_names}/rolling_window/{args.foldername_prefix}/plots"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # process input files
    GC_content_cv = process_input_files(
        args.GC_content_tsv, args.cv_gene_categories
    )
    # add root chromatin coverage data
    openchrom_cv = add_coverage(
        GC_content_cv, args.root_chrom_bp_covered, "root_chrom"
    )
    # add shoot chromatin coverage data
    openchrom_cv = add_coverage(
        openchrom_cv, args.shoot_chrom_bp_covered, "shoot_chrom"
    )
    # add rootshootintersect chromatin coverage data
    openchrom_cv = add_coverage(
        openchrom_cv,
        args.rootshootintersect_chrom_bp_covered,
        "rootshootintersect_chrom",
    )

    GC_content_tau = process_input_files(
        args.GC_content_tsv,
        args.tau_gene_categories,
    )

    # add root chromatin coverage data to the df
    openchrom_tau = add_coverage(
        GC_content_tau, args.root_chrom_bp_covered, "root_chrom"
    )
    # add shoot chromatin coverage data to the df
    openchrom_tau = add_coverage(
        openchrom_tau, args.shoot_chrom_bp_covered, "shoot_chrom"
    )
    # add rootshootintersect chromatin coverage data to the df
    openchrom_tau = add_coverage(
        openchrom_tau,
        args.rootshootintersect_chrom_bp_covered,
        "rootshootintersect_chrom",
    )

    # plot all promoters in genome
    windows_coords(
        args.file_names,
        args.foldername_prefix,
        "Araport11_allproms",
        "percentage_GC_content",
        "% open chromatin",
        openchrom_cv,
        openchrom_tau,
        args.promoter_bed,
        args.promoter_5UTR_bed,
        50,
        args.EPD_TSS_bed,
        args.palette_cv,
        args.palette_tau,
        estimator="median",
    )

    # add gene type
    openchrom_prom_types_cv = add_genetype(
        openchrom_cv, args.cv_gene_categories
    )
    openchrom_prom_types_tau = add_genetype(
        openchrom_tau, args.tau_gene_categories
    )
    # rename genetype column
    if "gene_type_x" in openchrom_prom_types_cv.columns:
        openchrom_prom_types_cv.rename(
            columns={"gene_type_x": "gene_type"}, inplace=True
        )
    if "gene_type_x" in openchrom_prom_types_tau.columns:
        openchrom_prom_types_tau.rename(
            columns={"gene_type_x": "gene_type"}, inplace=True
        )
    # plot with genetypes
    windows_coords(
        args.file_names,
        args.foldername_prefix,
        "Araport11_genetypenocontrol",
        "percentage_GC_content",
        "% open chromatin",
        openchrom_prom_types_cv,
        openchrom_prom_types_tau,
        args.promoter_bed,
        args.promoter_5UTR_bed,
        50,
        args.EPD_TSS_bed,
        args.palette_cv,
        args.palette_tau,
        includeEPDTSS=False,
        x_range=1350,
        estimator="median",
        ci=95,
        n_boot=10000,
        genetype_cv="constitutive",
        genetype_tau="non-specific",
        genetype2_cv="variable",
        genetype2_tau="tissue_specific",
    )

    # plot including control

    windows_coords(
        args.file_names,
        args.foldername_prefix,
        "Araport11_genetype",
        "percentage_GC_content",
        "% open chromatin",
        openchrom_prom_types_cv,
        openchrom_prom_types_tau,
        args.promoter_bed,
        args.promoter_5UTR_bed,
        50,
        args.EPD_TSS_bed,
        args.palette_cv,
        args.palette_tau,
        includeEPDTSS=False,
        x_range=1350,
        estimator="median",
        ci=95,
        n_boot=10000,
        genetype_cv="constitutive",
        genetype_tau="non-specific",
        genetype2_cv="variable",
        genetype2_tau="tissue_specific",
        genetype3_cv="control",
        genetype3_tau="control",
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
