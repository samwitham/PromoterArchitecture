# %%
# choose guides falling in my promoters of interest for the mutation library
# guide pairs falling around 100bp from each other are chosen to increase chances of 100bp deletion
# Installed DNA cauldron from Edinburgh genome foundry
# sudo pip install dnacauldron


# import dnacauldron as dc
import numpy as np

# %%
import pandas as pd
# from scipy.spatial.distance import pdist, squareform

# %%
ARF9_guides_file = "../../data/CRISPR_library/sgRNAs-ARF9.csv"
ARF18_guides_file = "../../data/CRISPR_library/sgRNAs-ARF18.csv"
DREB26_guides_file = "../../data/CRISPR_library/sgRNAs-DREB26.csv"
NLP7_guides_file = "../../data/CRISPR_library/sgRNAs-NLP7.csv"

# %%
# ARF9_guides_file = '../../data/CRISPR_library/500bp_region/sgRNAs-ARF9.csv'
# ARF18_guides_file = '../../data/CRISPR_library/500bp_region/sgRNAs-ARF18.csv'
# DREB26_guides_file = '../../data/CRISPR_library/500bp_region/sgRNAs-DREB26.csv'
# NLP7_guides_file = '../../data/CRISPR_library/500bp_region/sgRNAs-NLP7.csv'


# %%
def select_guides(
    input_file,
    region_start,
    region_end,
    lower_bound_distance,
    upper_bound_distance,
):
    """function to put select guide pairs with a specified lowerbound and upper bound range for the distance.
    Selects guides within a specified chromosome region too"""
    guides = pd.read_csv(input_file, header=2)
    # rename headers
    cols = [
        "Benchling_position",
        "distance_from_ATG_bp",
        "chromosome_position",
        "strand",
        "sequence",
        "PAM",
        "specificity_score_2013",
        "specificity_score_2016",
        "efficiency_score",
    ]
    guides.columns = cols
    # add old guide no. column
    guides["oldname"] = guides.index + 1
    # prepend string to name
    guides["oldname"] = "guide" + guides["oldname"].astype(str)
    # select guides only between specified region
    guides = guides[
        (guides.chromosome_position >= region_start)
        & (guides.chromosome_position <= region_end)
    ].copy()

    # remove if efficiency less than 40%
    guides = guides[guides.efficiency_score >= 40]
    # rename guides to match index and add to a new column
    # reset index first
    guides = guides.reset_index(drop=True)
    # rename the new index column
    guides.index.name = "name"
    # then turn new index into column
    guides = guides.reset_index()
    # add 1 to each value in the new column
    guides.name = guides.name + 1
    # then prepend the string 'guides' to that column
    guides["name"] = "guide" + guides["name"].astype(str)
    # pairwise differences between distance from ATG
    differences = abs(
        guides["distance_from_ATG_bp"].values
        - guides["distance_from_ATG_bp"].values[:, None]
    )
    # turn back into dataframe
    differences = pd.DataFrame(
        differences, columns=guides.index, index=guides.index
    )

    # filter values less than 110 and more than 90
    df = guides
    results = []
    # created based on https://stackoverflow.com/questions/51409927/python-selecting-pairs-of-objects-from-a-data-frame

    def dist(a, b):
        """
        Calculates the Euclidean distance between two 3D-vectors.
        """
        diff = np.array(a) - np.array(b)
        d = np.sqrt(np.dot(diff, diff))
        return d

    for i, row1 in df.iterrows():

        # calculate distance between current coordinate and all original rows in the data
        df["distance"] = df.apply(
            lambda row2: dist(
                row1["distance_from_ATG_bp"], row2["distance_from_ATG_bp"]
            ),
            axis=1,
        )
        # filter only those within a specific distance and drop rows with same name as current coordinate
        df_tmp = df[
            (df["distance"] > lower_bound_distance)
            & (df["distance"] < upper_bound_distance)
            & (df["name"] != row1["name"])
        ].copy()
        # prepare final data
        df_tmp["name2"] = row1["name"]
        df_tmp["guide2"] = df_tmp["name"]
        df_tmp["guide1"] = df_tmp["name2"]
        df_tmp["pairs"] = list(zip(df_tmp["name"], df_tmp["name2"]))
        # remember data
        results.append(df_tmp)

    # combine all into one dataframe
    df = pd.concat(results)
    # select columns of interest
    df = df[["pairs", "guide1", "guide2", "distance"]]
    # split each tuple value into just the guide number so they can be sorted later

    df["number1"] = df.guide1.str.split("e", expand=True)[1]
    df["number2"] = df.guide2.str.split("e", expand=True)[1]
    # turn into integars
    df = df.astype({"number1": "int", "number2": "int"})
    # make tuple column
    df["pair_numbers"] = list(zip(df.number2, df.number1))

    def sort_tuple(idx):
        x, y = idx
        if y < x:
            return y, x
        return x, y

    # sort values of each tuple from low to high
    df["pair_numbers"] = df["pair_numbers"].apply(sort_tuple)

    # drop duplicates
    df.drop_duplicates(subset=["pair_numbers"], inplace=True)
    # add functionality so that it reduces guide pairs so that all guides are included at least once,
    # and if a guide appears in more than one pair, rank by distance between the individuals in each pair and
    # choose the pair with the shortest distance in between them as long as it doesn't remove the last remaining individual of a different guide
    # reset index so I can compare pair rows
    df = df.reset_index(drop=True)
    # append guide 2 to guide 1
    df["all_guides"] = df["number1"]
    df2 = df["all_guides"].append(df["number2"])
    # convert to dataframe
    df2 = pd.DataFrame(df2)
    # add column name
    df2.columns = ["guides_present_in_pairs"]
    # add other columns (merge on index)
    merged = pd.merge(
        df[
            [
                "pairs",
                "pair_numbers",
                "number1",
                "number2",
                "guide1",
                "guide2",
                "distance",
            ]
        ],
        df2,
        left_index=True,
        right_index=True,
    )

    # sort merged by distance
    merged.sort_values("distance", inplace=True)
    # drop duplicates from all_guides column, keeping shortest distance

    merged = merged.drop_duplicates(subset=["guides_present_in_pairs"])
    # sort values of each tuple from low to high
    merged.sort_values(["number1", "number2"], inplace=True)

    # make guides present in at least 1 pair
    guides_present_in_pairs = merged.sort_values("guides_present_in_pairs")
    guides_present_in_pairs[
        "guides_present_in_pairs"
    ] = "guide" + guides_present_in_pairs["guides_present_in_pairs"].astype(
        str
    )
    # merge old guides df with guides_present to get the old names
    guides_present_in_pairs = pd.merge(
        guides_present_in_pairs,
        guides,
        left_on="guides_present_in_pairs",
        right_on="name",
        how="left",
    )
    # drop duplicates in merged
    merged.drop_duplicates(subset=["pair_numbers"], inplace=True)
    # reset index
    merged = merged.reset_index(drop=True)
    # swap guide column values to layout if possible where all duplicates of a certain guide are in the same column
    # (so that later then all go into the same level1 Goldengate construct)
    # make dictionary for guide1 column
    dict1 = {}
    # make dictionary for guide 2 column
    dict2 = {}
    # iterrate over rows in merged df
    merged_copy = merged.copy()
    for i, row in merged_copy.iterrows():
        # record guidenames in column 1 and 2 for that row
        guide1 = merged.loc[i, "guide1"]
        guide2 = merged.loc[i, "guide2"]
        # if guide 1 is in column 1 or guide 2 is already in column 2, pass
        if any(guide1 == v for v in dict1.values()) or any(
            guide2 == v for v in dict2.values()
        ):
            pass
        # if guide 1 is in column 2 and guide 2 is not in column1, swap places with guide 2
        elif any(guide1 == v for v in dict2.values()) and any(
            guide2 != v for v in dict1.values()
        ):
            merged.loc[i, "guide1"] = guide2
            merged.loc[i, "guide2"] = guide1
        # record guidenames again incase they switched
        guide1 = merged.loc[i, "guide1"]
        guide2 = merged.loc[i, "guide2"]
        # do the same again with the potential switched columns
        # if guide 1 is in column 1 or guide 2 is already in column 2, pass
        if any(guide1 == v for v in dict1.values()) or any(
            guide2 == v for v in dict2.values()
        ):
            pass
        # if guide 1 is in column 2 and guide 2 is not in column1, swap places with guide 2
        elif any(guide1 == v for v in dict2.values()) and any(
            guide2 != v for v in dict1.values()
        ):
            merged.loc[i, "guide1"] = guide2
            merged.loc[i, "guide2"] = guide1
        # record guidenames again incase they switched
        guide1 = merged.loc[i, "guide1"]
        guide2 = merged.loc[i, "guide2"]
        # if they got swapped back then there was no way to sort them into the same column. Now do the same for guide2 in case nothing got switched
        # if guide 2 is in column 2 or guide1 in column 1, pass
        if any(guide2 == v for v in dict2.values()) or any(
            guide1 == v for v in dict1.values()
        ):
            pass
        # if guide 2 is in column 1 and guide 1 is not in column 2, swap guide 1 and guide 2
        elif any(guide2 == v for v in dict1.values()) and any(
            guide1 != v for v in dict2.values()
        ):
            merged.loc[i, "guide1"] = guide2
            merged.loc[i, "guide2"] = guide1
        # record guidenames again incase they switched
        guide1 = merged.loc[i, "guide1"]
        guide2 = merged.loc[i, "guide2"]
        # do the same again for the potentially switched columns
        # if guide 2 is in column 2 or guide1 in column 1, pass
        if any(guide2 == v for v in dict2.values()) or any(
            guide1 == v for v in dict1.values()
        ):
            pass
        # if guide 2 is in column 1 and guide 1 is not in column 2, swap guide 1 and guide 2
        elif any(guide2 == v for v in dict1.values()) and any(
            guide1 != v for v in dict2.values()
        ):
            merged.loc[i, "guide1"] = guide2
            merged.loc[i, "guide2"] = guide1
        # now add the columns values to the dictionary
        dict1[i] = merged.loc[i, "guide1"]
        dict2[i] = merged.loc[i, "guide2"]
    # df.reset_index(drop=True)
    return (
        guides_present_in_pairs[["guides_present_in_pairs", "oldname"]],
        merged[["pairs", "guide1", "guide2", "distance"]],
        guides[
            [
                "name",
                "oldname",
                "Benchling_position",
                "distance_from_ATG_bp",
                "strand",
                "sequence",
                "PAM",
                "specificity_score_2013",
                "specificity_score_2016",
                "efficiency_score",
            ]
        ].reset_index(drop=True),
    )


# %%
def list_used_guides(guide_pairs, guide_info):
    """function to take in the guide pairs and guide information dfs and return a list of guides which are in at least 1 pair"""
    # return list of primers which are in one of the selected pairs
    usedguides = guide_info[
        guide_info.name.isin(guide_pairs.guide1)
        | guide_info.name.isin(guide_pairs.guide2)
    ]
    print(
        f"number of guides in region after filtering by efficiency = {len(guide_info)}"
    )
    print(f"number of guides in pairs = {len(usedguides)}")
    return usedguides


# %%
def adapt_primers(gene_info_df, gene_name, output_file):
    """function to add the correct start and end sequences to the sgRNA (add the guide RNA scaffold)"""
    # add a G base to the start (5' end) of the primer if it doesn't already have one
    # iterate over rows
    for i, data in gene_info_df.iterrows():
        # if the spacer sequence doesn't start with a g, add a g to the start of the new sequence
        if (
            gene_info_df.loc[i, "sequence"][0] != "G"
            and gene_info_df.loc[i, "sequence"][0] != "g"
        ):
            gene_info_df.loc[i, "new_sequence"] = (
                "g" + gene_info_df.loc[i, "sequence"]
            )
        else:
            # else use the original sequence
            gene_info_df.loc[i, "new_sequence"] = gene_info_df.loc[
                i, "sequence"
            ]
        # add the BsaI to the start of the sequence
        gene_info_df.loc[i, "new_sequence"] = (
            "tgtGGTCTCtatt" + gene_info_df.loc[i, "new_sequence"]
        )
        # add the sequence modifications extending the stem of the secondary structure as in Chen et al 2013. Cell 155 (7):1479-1491.
        # (only add the first part of the extension, the rest is in the sgRNA scaffold) gtttaagagctatgctggaaac
        gene_info_df.loc[i, "new_sequence"] = (
            gene_info_df.loc[i, "new_sequence"] + "gtttaagagctatgctggaaac"
        )
        # change strand to F or R
        if gene_info_df.loc[i, "strand"] == 1:
            gene_info_df.loc[i, "strand"] = "F"
        elif gene_info_df.loc[i, "strand"] == -1:
            gene_info_df.loc[i, "strand"] = "R"
        # add new_sequence length column
        gene_info_df.loc[i, "length"] = len(
            gene_info_df.loc[i, "new_sequence"]
        )
        # add new column containing gene name and the new guide number
        gene_info_df.loc[i, "gene_guidename"] = (
            gene_name + "_" + gene_info_df.loc[i, "name"]
        )
    # change column order
    gene_info_df = gene_info_df[
        [
            "gene_guidename",
            "strand",
            "length",
            "new_sequence",
            "name",
            "oldname",
            "Benchling_position",
            "distance_from_ATG_bp",
            "PAM",
            "specificity_score_2013",
            "specificity_score_2016",
            "efficiency_score",
        ]
    ]
    # save the df
    gene_info_df.to_csv(output_file, header=1, index=None)
    return gene_info_df


# %%
def rename_guidepairs(df, gene_name, output_file):
    """function to rename the guides to include the gene name too"""
    df.guide1 = gene_name + "_" + df.guide1
    df.guide2 = gene_name + "_" + df.guide2
    # write to file
    df.to_csv(output_file, header=1, index=None)


# %%
# return df with list of primer pairs, and a df containing guides information
ARF9_guides_present, ARF9pairs, ARF9info = select_guides(
    ARF9_guides_file,
    region_start=12450889,
    region_end=12451388,
    lower_bound_distance=90,
    upper_bound_distance=140,
)
ARF18_guides_present, ARF18pairs, ARF18info = select_guides(
    ARF18_guides_file, 22887610, 22888109, 90, 140
)
DREB26_guides_present, DREB26pairs, DREB26info = select_guides(
    DREB26_guides_file, 7696155, 7696654, 90, 140
)
NLP7_guides_present, NLP7pairs, NLP7info = select_guides(
    NLP7_guides_file, 12479404, 12479903, 90, 140
)

# %%
# return list of primers which are in one of the selected pairs
ARF9_usedguides = list_used_guides(ARF9pairs, ARF9info)

# %%
ARF9_guides_present

# %%
# list the ARF9 pairs - note there are several guides present in more than 1 pair
ARF9pairs

# %%
# return list of primers which are in one of the selected pairs
ARF18_usedguides = list_used_guides(ARF18pairs, ARF18info)

# %%
ARF18_guides_present

# %%
# ARF18 pairs
ARF18pairs

# %%
# return list of primers which are in one of the selected pairs
DREB26_usedguides = list_used_guides(DREB26pairs, DREB26info)

# %%
DREB26info

# %%
DREB26_guides_present


# %%
DREB26pairs

# %%
DREB26info

# %%
# return list of primers which are in one of the selected pairs
NLP7_usedguides = list_used_guides(NLP7pairs, NLP7info)

# %%
NLP7pairs

# %%
NLP7_guides_present

# %%
# arf9 extending over part 2 open chromatin region ending at 12450824
ARF9_guides_present2, ARF9pairs2, ARF9info2 = select_guides(
    ARF9_guides_file,
    region_start=12450310,
    region_end=12451388,
    lower_bound_distance=90,
    upper_bound_distance=157,
)

# %%
ARF2_usedguides2 = list_used_guides(ARF9pairs2, ARF9info2)

# %%
ARF9_guides_present2

# %%
ARF9pairs2

# %%
ARF9info2

# %%
# arf18 extending until 30 PAMs
ARF18_guides_present2, ARF18pairs2, ARF18info2 = select_guides(
    ARF18_guides_file,
    region_start=22887459,
    region_end=22888109,
    lower_bound_distance=90,
    upper_bound_distance=149,
)

# %%
ARF18_guides_present2

# %%
ARF18pairs2

# %%
ARF18_usedguides2 = list_used_guides(ARF18pairs2, ARF18info2)

# %%
DREB26_guides_present2, DREB26pairs2, DREB26info2 = select_guides(
    DREB26_guides_file, 7695780, 7696654, 90, 168
)

# %%
DREB26_usedguides2 = list_used_guides(DREB26pairs2, DREB26info2)

# %%
DREB26_guides_present2

# %%
DREB26pairs2

# %%
NLP7_guides_present2, NLP7pairs2, NLP7info2 = select_guides(
    NLP7_guides_file, 12478898, 12479903, 90, 256
)

# %%
NLP7_usedguides2 = list_used_guides(NLP7pairs2, NLP7info2)

# %%
NLP7_guides_present2

# %%
NLP7pairs2

# %%
# output the csv files with new sequences

# %%
ARF9adapted_primers = adapt_primers(
    ARF9info2, "ARF9", "../../data/CRISPR_library/sgRNAs-ARF9_new.csv"
)

# %%
ARF18adapted_primers = adapt_primers(
    ARF18info2, "ARF18", "../../data/CRISPR_library/sgRNAs-ARF18_new.csv"
)

# %%
DREB26adapted_primers = adapt_primers(
    DREB26info2, "DREB26", "../../data/CRISPR_library/sgRNAs-DREB26_new.csv"
)

# %%
NLP7adapted_primers = adapt_primers(
    NLP7info2, "NLP7", "../../data/CRISPR_library/sgRNAs-NLP7_new.csv"
)

# %%
# rename guide pairs to include the gene name and then save as csv

# %%
# output the csv files with renamed guide pairs (which include the genename)
rename_guidepairs(
    ARF9pairs2, "ARF9", "../../data/CRISPR_library/ARF9guidepairs.csv"
)
rename_guidepairs(
    ARF18pairs2, "ARF18", "../../data/CRISPR_library/ARF18guidepairs.csv"
)
rename_guidepairs(
    DREB26pairs2, "DREB26", "../../data/CRISPR_library/DREB26guidepairs.csv"
)
rename_guidepairs(
    NLP7pairs2, "NLP7", "../../data/CRISPR_library/NLP7guidepairs.csv"
)

# %%
# output all_guides.csv
# concatenate the primer dfs for each gene
cat = pd.concat(
    [
        ARF9adapted_primers,
        ARF18adapted_primers,
        DREB26adapted_primers,
        NLP7adapted_primers,
    ]
)
# reset index
cat.reset_index(drop=True, inplace=True)
# reset index again this time making a new column with the inex numbers
cat.reset_index(inplace=True)
# add IDTname column for easy ordering
cat["IDT_name"] = ""
# rename columns
cat = cat.rename(
    columns={
        "index": "Name",
        "gene_guidename": "Description",
        "strand": "Strand",
        "new_sequence": "Sequence",
    }
)
# filter and order columns
cat = cat[["Name", "Description", "Strand", "length", "IDT_name", "Sequence"]]
# add 107 to Name column to match current primer number (change this number to your current number)
cat.Name = cat.Name + 107
# change name column to string
cat = cat.astype({"Name": "str"})
# prepend Name column with SW
cat.Name = "SW" + cat.Name
# save as csv
cat.to_csv("../../data/CRISPR_library/all_guides.csv", header=1, index=False)

# %%
