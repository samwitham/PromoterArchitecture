# use env pybedtools
# parse arguments
import argparse
import io
import math
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

# chunks
from more_itertools import sliced
from pybedtools import BedTool


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="pacbio_analyse_overlapping_TFBSs"
    )
    parser.add_argument(
        "mapped_motifs_bed",
        type=str,
        help="Input location FIMO mapped motifs bed",
    )
    parser.add_argument(
        "mutations_tsv",
        type=str,
        help="Input location of tsv file containing mutations (no duplicates) after running pacbio_analyse_variantcall.py",
    )
    parser.add_argument(
        "gene_name",
        type=str,
        help="Name of gene being analysed",
    )
    parser.add_argument(
        "file_number",
        type=str,
        help="File number of file part",
    )

    parser.add_argument(
        "output_folder",
        type=str,
        help="Output folder location ending with a forward slash",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def read_in_files(mutations_file, mapped_motifs_bed):
    """read in the files"""
    # read in mapped motifs bed file
    mapped_motifs = pd.read_table(mapped_motifs_bed, sep="\t", header=None)
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
            "motif_chr",
            "motif_start",
            "motif_stop",
            "name_rep",
            "score",
            "motif_strand",
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
            "motif_chr",
            "motif_start",
            "motif_stop",
            "name_rep",
            "score",
            "motif_strand",
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
            "motif_chr",
            "motif_start",
            "motif_stop",
            "name_rep",
            "score",
            "motif_strand",
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

    mutations_df = pd.read_table(mutations_file, sep="\t", header=None)
    mutations_df_cols = [
        "chr",
        "plant_ID",
        "platename",
        "library",
        "first_reaction_primers",
        "second_reaction_primers",
        "guide",
        "guide_number",
        "aligned_sequence",
        "reference_sequence",
        "mutation_type",
        "read_number",
        "read_percentage",
        "insertion_positions",
        "deletion_positions",
        "substitution_positions",
        "insertion_cut_site_distance",
        "deletion_cut_site_distance",
        "substitution_cut_site_distance",
        "cut_site_promoter_position",
        "insertion_positions_relative_to_TSS",
        "insertion_genomic_positions",
        "deletion_positions_relative_to_TSS",
        "deletion_genomic_positions",
        "substitution_positions_relative_to_TSS",
        "substitution_genomic_positions",
    ]
    mutations_df.columns = mutations_df_cols
    # guide_pairs_df = pd.read_csv(guide_pairs,header=0)
    # only keep first 2 columns
    #   guide_cols = ['guide1','guide2']
    #  guide_pairs_df = guide_pairs_df[guide_cols]
    return mutations_df, mapped_motifs


def merge_bedfiles(bedfile, mapped_motifs_bed, output_buffer):
    """perform bedtools intersect on the two dfs"""
    df = BedTool(bedfile)
    motifs = BedTool(mapped_motifs_bed)
    # -wao =Write the original A and B entries plus the number of base pairs of overlap between the two features.
    # However, A features w/o overlap are also reported with a NULL B feature and overlap = 0
    intersect = df.intersect(motifs, wao=True)
    # Write to output_file
    # Each line in the file contains bed entry a and bed entry b that it overlaps plus the number of bp in the overlap so 19 columns
    output_buffer.write(str(intersect))
    # go back to beginning of buffer
    output_buffer.seek(0)
    mapped_motifs_bed.seek(0)
    return output_buffer


def find_overlapping_TFBSs(
    mutations_df_chunk, mapped_motifs_bed, mapped_motifs_bed_columns
):
    """function to find any overlapping TFBSs from FIMO mapped motif file"""
    # go back to start of buffer
    mapped_motifs_bed.seek(0)

    # turn mutations_df_chunk into list of dicts
    mutations_df_chunk_dict = mutations_df_chunk.to_dict(orient="records")

    for mutations_df_row_dict in mutations_df_chunk_dict:
        # mutations_df_index = mutations_df_row.index

        if mutations_df_row_dict["mutation_type"] == "None":
            pass
        else:
            # create temporary df in bed format

            cols = ["chr", "start", "stop", "mutation_type", "mutation_count"]
            temp_df = pd.DataFrame(columns=cols)
            chr = mutations_df_row_dict["chr"]
            # if not NaN
            if (
                pd.notna(mutations_df_row_dict["insertion_genomic_positions"])
                and mutations_df_row_dict["insertion_genomic_positions"]
                != "nan"
            ):
                # convert genomic positions from string to list
                # Convert string to list if more than one
                insertion_genomic_positions = (
                    mutations_df_row_dict["insertion_genomic_positions"]
                    .strip("][")
                    .split(", ")
                )
                # count which mutation number currently on to be added to the temporary bed file
                count = 0
                for gen_pos in insertion_genomic_positions:
                    count += 1
                    start = int(gen_pos)
                    stop = start + 1
                    mutation_type = "insertion"
                    # add to temp_df
                    temp_list = [chr, start, stop, mutation_type, count]
                    temp_df.loc[len(temp_df)] = temp_list
            if (
                pd.notna(mutations_df_row_dict["deletion_genomic_positions"])
                and mutations_df_row_dict["deletion_genomic_positions"]
                != "nan"
            ):
                # Convert string to list
                deletion_genomic_positions = (
                    mutations_df_row_dict["deletion_genomic_positions"]
                    .strip("][")
                    .split(", ")
                )
                # count which mutation number currently on to be added to the temporary bed file
                count = 0
                for gen_pos in deletion_genomic_positions:
                    count += 1
                    start = int(gen_pos)
                    stop = start + 1
                    mutation_type = "deletion"
                    # add to temp_df
                    temp_list = [chr, start, stop, mutation_type, count]
                    temp_df.loc[len(temp_df)] = temp_list
            if (
                pd.notna(
                    mutations_df_row_dict["substitution_genomic_positions"]
                )
                and mutations_df_row_dict["substitution_genomic_positions"]
                != "nan"
            ):
                # Convert string to list
                substitution_genomic_positions = (
                    mutations_df_row_dict["substitution_genomic_positions"]
                    .strip("][")
                    .split(", ")
                )
                # count which mutation number currently on to be added to the temporary bed file
                count = 0
                for gen_pos in substitution_genomic_positions:
                    count += 1
                    start = int(gen_pos)
                    stop = start + 1
                    mutation_type = "substitution"
                    # add to temp_df
                    temp_list = [chr, start, stop, mutation_type, count]
                    temp_df.loc[len(temp_df)] = temp_list
            # now do bedtools intersect to find which TFBSs overlap with which mutations
            # sort by chr then start
            temp_df = temp_df.sort_values(["chr", "start"]).reset_index(
                drop=True
            )
            # write to buffer
            temp_df_buffer = io.StringIO()
            temp_df.to_csv(temp_df_buffer, sep="\t", index=False, header=None)
            temp_df_buffer.seek(0)

            output_buffer = io.StringIO()

            output_buffer = merge_bedfiles(
                temp_df_buffer, mapped_motifs_bed, output_buffer
            )
            # read in output buffer as df
            output_df = pd.read_table(output_buffer, sep="\t")

            # get column names and rename columns
            output_df_cols = cols + mapped_motifs_bed_columns + ["bp_overlap"]
            output_df.columns = output_df_cols
            # for each mutation type get list of overlapping TFBSs. Add these to a dictionary column in the mutations_df_chunk
            # create defaultdicts with lists as values so that non-existing keys can be added to in one go
            insertion_overlapping_TFBS_family = defaultdict(list)
            insertion_overlapping_TFBS_AGI = defaultdict(list)

            deletion_overlapping_TFBS_family = defaultdict(list)
            deletion_overlapping_TFBS_AGI = defaultdict(list)

            substitution_overlapping_TFBS_family = defaultdict(list)
            substitution_overlapping_TFBS_AGI = defaultdict(list)

            # convert output_df into dict
            output_df_dict = output_df.to_dict(orient="records")
            for row_dict in output_df_dict:
                # if mutation and has a 1bp overlap
                # if no TF_AGI then pass
                if row_dict["TF_AGI"] == ".":
                    pass
                elif (
                    row_dict["mutation_type"] == "insertion"
                    and row_dict["bp_overlap"] > 0
                ):
                    # add insertion TFBS family information to dictionary for the correct mutation number
                    mut_count = row_dict["mutation_count"]
                    insertion_overlapping_TFBS_family[
                        f"insertion{mut_count}"
                    ] += [row_dict["TF_family"]]
                    # add insertion TFBS AGI information to dictionary for the correct mutation number
                    insertion_overlapping_TFBS_AGI[
                        f"insertion{mut_count}"
                    ] += [row_dict["TF_AGI"]]

                elif (
                    row_dict["mutation_type"] == "deletion"
                    and row_dict["bp_overlap"] > 0
                ):
                    # add deletion TFBS family information to dictionary for the correct mutation number
                    mut_count = row_dict["mutation_count"]
                    deletion_overlapping_TFBS_family[
                        f"deletion{mut_count}"
                    ] += [row_dict["TF_family"]]
                    # add deletion TFBS AGI information to dictionary for the correct mutation number
                    deletion_overlapping_TFBS_AGI[f"deletion{mut_count}"] += [
                        row_dict["TF_AGI"]
                    ]

                elif (
                    row_dict["mutation_type"] == "substitution"
                    and row_dict["bp_overlap"] > 0
                ):
                    # add substitution TFBS family information to dictionary for the correct mutation number
                    mut_count = row_dict["mutation_count"]
                    substitution_overlapping_TFBS_family[
                        f"substitution{mut_count}"
                    ] += [row_dict["TF_family"]]
                    # add substitution TFBS AGI information to dictionary for the correct mutation number
                    substitution_overlapping_TFBS_AGI[
                        f"substitution{mut_count}"
                    ] += [row_dict["TF_AGI"]]

            # add values to mutations_df row
            # first make overlapping TFBS families and AGIs unique
            insertion_overlapping_TFBS_family = dict(
                insertion_overlapping_TFBS_family
            )
            insertion_overlapping_TFBS_AGI = dict(
                insertion_overlapping_TFBS_AGI
            )
            deletion_overlapping_TFBS_family = dict(
                deletion_overlapping_TFBS_family
            )
            deletion_overlapping_TFBS_AGI = dict(deletion_overlapping_TFBS_AGI)
            substitution_overlapping_TFBS_family = dict(
                substitution_overlapping_TFBS_family
            )
            substitution_overlapping_TFBS_AGI = dict(
                substitution_overlapping_TFBS_AGI
            )
            # if empty dictionary, change to nan
            if insertion_overlapping_TFBS_family == {}:
                insertion_overlapping_TFBS_family = np.nan
            else:
                # keep only unique TFBS families
                for k, v in insertion_overlapping_TFBS_family.items():
                    # print(np.unique(v).astype(list))
                    insertion_overlapping_TFBS_family[k] = np.unique(
                        v
                    ).tolist()

            if insertion_overlapping_TFBS_AGI == {}:
                insertion_overlapping_TFBS_AGI = np.nan
            else:
                # keep only unique TFBS AGIs
                for k, v in insertion_overlapping_TFBS_AGI.items():
                    # print(np.unique(v).astype(list))
                    insertion_overlapping_TFBS_AGI[k] = np.unique(v).tolist()

            if deletion_overlapping_TFBS_family == {}:
                deletion_overlapping_TFBS_family = "nan"
            else:
                # keep only unique TFBS families
                for k, v in deletion_overlapping_TFBS_family.items():
                    # print(np.unique(v).astype(list))
                    deletion_overlapping_TFBS_family[k] = np.unique(v).tolist()

            if deletion_overlapping_TFBS_AGI == {}:
                deletion_overlapping_TFBS_AGI = np.nan
            else:
                # keep only unique TFBS AGIs
                for k, v in deletion_overlapping_TFBS_AGI.items():
                    # print(np.unique(v).astype(list))
                    deletion_overlapping_TFBS_AGI[k] = np.unique(v).tolist()

            if substitution_overlapping_TFBS_family == {}:
                substitution_overlapping_TFBS_family = np.nan
            else:
                # keep only unique TFBS families
                for k, v in substitution_overlapping_TFBS_family.items():
                    # print(np.unique(v).astype(list))
                    substitution_overlapping_TFBS_family[k] = np.unique(
                        v
                    ).tolist()

            if substitution_overlapping_TFBS_AGI == {}:
                substitution_overlapping_TFBS_AGI = np.nan
            else:
                # keep only unique TFBS AGIs
                for k, v in substitution_overlapping_TFBS_AGI.items():
                    # print(np.unique(v).astype(list))
                    substitution_overlapping_TFBS_AGI[k] = np.unique(
                        v
                    ).tolist()
            # append values to mutations_df_row_dict
            mutations_df_row_dict["insertion_overlapping_TFBS_family"] = str(
                insertion_overlapping_TFBS_family
            )
            mutations_df_row_dict["insertion_overlapping_TFBS_AGI"] = str(
                insertion_overlapping_TFBS_AGI
            )
            mutations_df_row_dict["deletion_overlapping_TFBS_family"] = str(
                deletion_overlapping_TFBS_family
            )
            mutations_df_row_dict["deletion_overlapping_TFBS_AGI"] = str(
                deletion_overlapping_TFBS_AGI
            )
            mutations_df_row_dict[
                "substitution_overlapping_TFBS_family"
            ] = str(substitution_overlapping_TFBS_family)
            mutations_df_row_dict["substitution_overlapping_TFBS_AGI"] = str(
                substitution_overlapping_TFBS_AGI
            )

            # go back to start of buffer
            mapped_motifs_bed.seek(0)
    return pd.DataFrame.from_dict(mutations_df_chunk_dict)


def chunkify(mutations_df, mapped_motifs_df, output_folder, gene, file_number):
    """function to prepare dfs and slice mutations_df into chunks to reduce memory before running find_overlapping_TFBSs()"""

    # for each guide containing mutations, create a temporary bed file containing each mutation and then do bedtools intersect to find which overlap TFBs
    # then add the TFBS names into a new column for that row

    # get column names from mapped_motifs_df
    mapped_motifs_bed_columns = list(mapped_motifs_df.columns)

    # turn mapped_motifs_df into a buffer
    mapped_motifs_bed = io.StringIO()
    mapped_motifs_df.to_csv(
        mapped_motifs_bed, sep="\t", index=False, header=None
    )
    # go back to start of buffer
    mapped_motifs_bed.seek(0)

    # add columns to mutations_df
    new_columns = [
        "insertion_overlapping_TFBS_family",
        "insertion_overlapping_TFBS_AGI",
        "deletion_overlapping_TFBS_family",
        "deletion_overlapping_TFBS_AGI",
        "substitution_overlapping_TFBS_family",
        "substitution_overlapping_TFBS_AGI",
    ]
    # first make a new df that will merge into mutations_df
    temp_new_df = pd.DataFrame(columns=new_columns)
    mutations_df = pd.concat([mutations_df, temp_new_df], axis=1)

    # first make certain columns string
    # make columns containing lists string for now so can use groupby
    to_string = [
        "insertion_positions",
        "deletion_positions",
        "substitution_positions",
        "insertion_cut_site_distance",
        "deletion_cut_site_distance",
        "substitution_cut_site_distance",
        "insertion_positions_relative_to_TSS",
        "insertion_genomic_positions",
        "deletion_positions_relative_to_TSS",
        "deletion_genomic_positions",
        "substitution_positions_relative_to_TSS",
        "substitution_genomic_positions",
    ]
    mutations_df[to_string] = mutations_df[to_string].astype(str)

    # convert mutations_df into chunks to reduce memory load
    CHUNK_SIZE = 5

    index_slices = sliced(range(len(mutations_df)), CHUNK_SIZE)
    # create list of chunks
    chunks = []
    for index_slice in index_slices:
        # go back to start of buffer
        mapped_motifs_bed.seek(0)
        chunk = mutations_df.iloc[
            index_slice
        ]  # your dataframe chunk ready for use
        new_chunk = find_overlapping_TFBSs(
            chunk, mapped_motifs_bed, mapped_motifs_bed_columns
        )
        chunks.append(new_chunk)

    # concatenate chunks into mutations_df
    mutations_df = pd.concat(chunks)
    # write out mutations_df
    mutations_df.to_csv(
        f"{output_folder}{gene}/{gene}_TFBSoverlapping.part{file_number}.part",
        sep="\t",
        index=False,
        header=None,
    )


def main(args):
    # parse arguments
    args = parse_args(args)

    # read in files
    mutations_df, mapped_motifs_df = read_in_files(
        args.mutations_tsv, args.mapped_motifs_bed
    )

    # find overlapping TFBSs by splitting the mutations_ARF9_df into chunks and iterating over each chunk, reducing memory usage
    mutations_df_overlapping_TFBS = chunkify(
        mutations_df,
        mapped_motifs_df,
        args.output_folder,
        args.gene_name,
        args.file_number,
    )

    # genotype the plant lines, producing an output tsv
    # mutations_df_genotyped = genotype_plant_lines(mutations_df_overlapping_TFBS,args.output_folder,args.gene_name)

    # produce an output tsv with one row per plant line with all guide site mutations found in that line
    # categorise_guide_pairs(mutations_df_genotyped,guide_pairs_df,args.output_folder,args.gene_name)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
