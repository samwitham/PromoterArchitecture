# Sort the protoplast luminescence data from the xlsx output from the Glariostar platereader. Use 2 input excels at a time (one firefly, one nanoluc)

import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


def xlsx_2_csv(xlsx):
    """ Function to read and convert xlsx file to csv file. Also return the data (name of the folder the xlsx is in)"""

    # Read in the xlsx file, second sheet
    file = pd.read_excel(xlsx, "End point", index_col=None)

    filename = os.path.basename(xlsx)
    removed_extension = os.path.splitext(filename)[0]
    path = Path(
        xlsx
    ).parent  # find parent directory to the one the xlsx fields are in
    # date = Path(xlsx)

    file.to_csv(
        f"{path}/{removed_extension}.csv", encoding="utf-8", index=False
    )


def combine_csvs(
    input_fluc,
    input_nluc,
    layout_csv,
    date,
    output_file_means,
    output_file_raw,
):
    """Function to combine two csv files containing luminescence data, and label values using layout csv file (plate layout)"""
    # read in files
    fluc = pd.read_csv(input_fluc, header=None)
    nluc = pd.read_csv(input_nluc, header=None)

    # if experiment metadata present, take data from fow 10
    if "User" in fluc.iloc[0, 0]:
        fluc = pd.read_csv(input_fluc, header=9)
    elif "Well Row" in fluc.iloc[0, 0]:
        print(fluc.iloc[0, 0])
        fluc = pd.read_csv(input_fluc, header=0)
    if "User" in nluc.iloc[0, 0]:
        nluc = pd.read_csv(input_nluc, header=9)
    elif "Well Row" in nluc.iloc[0, 0]:
        nluc = pd.read_csv(input_nluc, header=0)

    layout = pd.read_csv(layout_csv, header=0)

    # make new df with correct column names, including both fluc and nluc data
    combined = fluc[
        [
            "Well\nRow",
            "Well\nCol",
            "Content",
            "Average over replicates based on Blank corrected (No filter)",
        ]
    ].copy()
    combined.rename(
        columns={
            "Well\nRow": "well_row",
            "Well\nCol": "well_col",
            "Content": "content",
            "Average over replicates based on Blank corrected (No filter)": "fluc_luminescence",
        },
        inplace=True,
    )
    combined["nluc_luminescence"] = nluc[
        "Average over replicates based on Blank corrected (No filter)"
    ].copy()
    # mask any values less than 400 (turn into NaNs)
    combined["fluc_luminescence"] = combined.fluc_luminescence.mask(
        combined.fluc_luminescence < 400
    )
    combined["nluc_luminescence"] = combined.nluc_luminescence.mask(
        combined.nluc_luminescence < 400
    )
    # merge layout with combined
    combined_named = pd.merge(combined, layout, on=["well_row", "well_col"])
    # convert well_col column data type to string so it is excluded from the next bit
    combined_named.well_col = combined_named.well_col.astype(np.str)
    # add new column, nluc/fluc
    combined_named["nluc/fluc"] = (
        combined_named.nluc_luminescence / combined_named.fluc_luminescence
    )
    # remove NaNs
    combined_named_no_null = combined_named[
        pd.notnull(combined_named["nluc/fluc"])
    ]
    # add date to the data
    combined_named_no_null_date = combined_named_no_null.copy()
    combined_named_no_null_date["date"] = date
    # make csv of raw data
    combined_named_no_null_date.to_csv(
        output_file_raw, encoding="utf-8", index=False
    )
    # make new df with mean luminescence
    mean = (
        combined_named_no_null[["name", "condition", "nluc/fluc"]]
        .groupby(["name", "condition"])
        .mean()
        .reset_index()
    )
    # mean = combined_named_no_null[['name', 'nluc/fluc']].groupby('name').mean().reset_index()
    mean.rename(columns={"nluc/fluc": "mean_luminescence"}, inplace=True)
    # add standard error
    standard_error = (
        combined_named_no_null[["name", "condition", "nluc/fluc"]]
        .groupby(["name", "condition"])
        .sem()
        .reset_index()
    )
    # standard_error = combined_named_no_null[['name','nluc/fluc']].groupby('name').sem().reset_index()
    standard_error.rename(
        columns={"nluc/fluc": "standard_error"}, inplace=True
    )
    mean_samples = pd.merge(mean, standard_error, on=["name", "condition"])
    # #mean_samples = pd.merge(mean, standard_error, on='name')
    # add date of experiment
    mean_samples["date"] = date
    # create output file
    mean_samples.to_csv(output_file_means, encoding="utf-8", index=False)


# find all xlsx files recursively in the 'to_be_sorted' folder
xlsx_filenames = glob.glob(
    "../../data/luminescence/to_be_sorted/**/*.xlsx", recursive=True
)

# run the xlsx_2_csv function across all xlsx file in to_be_sorted folder
list(map(xlsx_2_csv, xlsx_filenames))

# #variables - change path/name of these to files you're interested in. Also change the date
# input_fluc = '../../data/luminescence/to_be_sorted/30.8.19/fluc.csv'
# input_nluc = '../../data/luminescence/to_be_sorted/30.8.19/nluc1.csv'
# layout = '../../data/luminescence/to_be_sorted/30.8.19/layout.csv'
# output = '../../data/luminescence/to_be_sorted/30.8.19/output.csv'
# date = '30.8.19'
date = "27.11.20"
input_fluc = (
    f"../../data/luminescence/to_be_sorted/{date}/leaf_coexpression_fluc.csv"
)
input_nluc = (
    f"../../data/luminescence/to_be_sorted/{date}/leaf_coexpression_nluc.csv"
)
layout = f"../../data/luminescence/to_be_sorted/{date}/layout.csv"
output = f"../../data/luminescence/to_be_sorted/{date}/output_means.csv"
output_raw = f"../../data/luminescence/to_be_sorted/{date}/output_raw.csv"


# input_fluc = f"../../data/luminescence/to_be_sorted/{date}/Fluc_repleteLEAF_PROTOPLASTs_TFcoexpression16.csv"
# input_nluc = f"../../data/luminescence/to_be_sorted/{date}/Nluc_repleteLEAF_PROTOPLASTs_TFcoexpression16.csv"
# layout = f"../../data/luminescence/to_be_sorted/{date}/layout_16.10.19.csv"
# output = (
#     f"../../data/luminescence/to_be_sorted/{date}/output_means_16.10.19.csv"
# )
# output_raw = (
#     f"../../data/luminescence/to_be_sorted/{date}/output_raw_16.10.19.csv"
# )
# input_fluc = '../../data/luminescence/to_be_sorted/26.11.19/nitrate_leaf_fluc.csv'
# input_nluc = '../../data/luminescence/to_be_sorted/26.11.19/nitrate_leaf_nluc.csv'
# layout = '../../data/luminescence/to_be_sorted/26.11.19/layout.csv'
# output = '../../data/luminescence/to_be_sorted/26.11.19/output_means.csv'
# output_raw = '../../data/luminescence/to_be_sorted/26.11.19/output_raw.csv'
# date = '26.11.19'


# run combine_csvs function
combine_csvs(input_fluc, input_nluc, layout, date, output, output_raw)
