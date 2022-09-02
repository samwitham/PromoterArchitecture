# parse arguments
import argparse

# convert from string to dict
import ast
import itertools
import math
from collections import defaultdict

import numpy as np
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="pacbio_analyse_flattened_mutations"
    )
    parser.add_argument(
        "mutations_tsv",
        type=str,
        help="Input location of tsv file containing mutations (no duplicates) after running pacbio_analyse_overlapping_TFBSs_part2.py",
    )
    parser.add_argument(
        "gene_name",
        type=str,
        help="Name of gene being analysed",
    )
    parser.add_argument(
        "output_folder",
        type=str,
        help="Output folder location ending with a forward slash",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def read_in_files(mutations_file):
    """read in the files"""
    mutations_df = pd.read_table(mutations_file, sep="\t", header=0)
    # mutations_df_cols = ['chr','plant_ID','platename','library','first_reaction_primers','second_reaction_primers','guide','guide_number','aligned_sequence','reference_sequence','mutation_type','read_number','read_percentage','insertion_positions','deletion_positions','substitution_positions',
    # 'insertion_cut_site_distance','deletion_cut_site_distance','substitution_cut_site_distance','cut_site_promoter_position','insertion_positions_relative_to_TSS'
    # ,'insertion_genomic_positions','deletion_positions_relative_to_TSS','deletion_genomic_positions','substitution_positions_relative_to_TSS','substitution_genomic_positions']
    # mutations_df.columns = mutations_df_cols
    return mutations_df


def find_insertions_lengths(list_of_insertion_positions):
    # make copy:
    list_of_insertion_positions_copy = list_of_insertion_positions.copy()
    if len(list_of_insertion_positions_copy) == 1:
        return [1]

    elif len(list_of_insertion_positions_copy) > 1:
        # print(list_of_insertion_positions)
        ind = -1
        previous_mutation = []
        mutation_length = 1
        mutation_length_list = []
        # get length of list
        list_length = len(list_of_insertion_positions_copy)
        for mutation in list_of_insertion_positions_copy:
            # print(previous_mutation)
            # get index
            ind += 1
            # if index is 0, add it to previous mutation variable and move on to next mutation
            if ind == 0:
                previous_mutation = mutation
            # if index is one less than the length of the list
            elif ind == list_length - 1:
                # print(f'{list_length},{ind}')
                # if mutation is adjacent to previous mutation
                if mutation[0] == previous_mutation[1]:
                    # add 1 to mutation count
                    mutation_length += 1
                    # append mutation length to list
                    mutation_length_list.append(mutation_length)

                elif mutation[0] != previous_mutation[1]:
                    # append mutation length to list
                    mutation_length_list.append(mutation_length)
                    # reset mutation length
                    mutation_length = 1
                    # append mutation length to list
                    mutation_length_list.append(mutation_length)

            elif mutation[0] == previous_mutation[1]:
                mutation_length += 1
                previous_mutation = mutation
            elif mutation[0] != previous_mutation[1]:
                # append mutation length to mutation_length_list
                mutation_length_list.append(mutation_length)
                # reset mutation length
                mutation_length = 1
                previous_mutation = mutation

        return mutation_length_list


def find_mutation_lengths(list_of_mutation_positions):
    # make copy of list
    list_of_mutation_positions_copy = list_of_mutation_positions.copy()
    if len(list_of_mutation_positions_copy) == 1:
        return [1]
    elif len(list_of_mutation_positions_copy) > 1:
        ind = -1
        previous_mutation = int()
        mutation_length = 1
        mutation_length_list = []
        list_length = len(list_of_mutation_positions_copy)
        for mutation in list_of_mutation_positions_copy:
            # get index
            ind += 1
            # if index is 0, add it to previous mutation variable and move on to next mutation
            if ind == 0:
                previous_mutation = mutation
            elif ind == list_length - 1:
                # print(ind)
                # print(f'{list_length},{ind}')
                # if mutation is adjacent to previous mutation
                if mutation == previous_mutation + 1:
                    # add 1 to mutation count
                    mutation_length += 1
                    # append mutation length to list
                    mutation_length_list.append(mutation_length)

                elif mutation != previous_mutation + 1:
                    mutation_length_list.append(mutation_length)
                    # reset mutation length
                    mutation_length = 1
                    # append mutation length to list
                    mutation_length_list.append(mutation_length)

            elif mutation == previous_mutation + 1:
                mutation_length += 1
                previous_mutation = mutation
            elif mutation != previous_mutation + 1:
                # append mutation length to mutation_length_list
                mutation_length_list.append(mutation_length)
                # reset mutation length
                mutation_length = 1
                previous_mutation = mutation

        return mutation_length_list


def make_overlapping_TFBS_dict(
    input_AGIs, mutation_genomic_positions, output_dict, mutation_type
):
    """function to generate a dictionary of genomic positions overlapping TFBS AGIs"""
    # if no TFBSs overlapping insertions, pass
    if type(input_AGIs) == float and np.isnan(input_AGIs):
        pass
    # else if it is not nan carry on
    else:
        # convert from string to dict
        mutation_overlapping_TFBS_AGI = ast.literal_eval(input_AGIs)
        # make dictionary with genomic positions as keys, TFBS AGI a values
        TFBS_AGI_dict = {}
        # iterate through genomic positions
        for pos in mutation_genomic_positions:
            # get index
            pos_index = mutation_genomic_positions.index(pos)
            # if insertion is in insertion_overlapping_TFBS_AGI dict, then add it to the TFBS_AGI_dict
            if (
                f"{mutation_type}{pos_index+1}"
                in mutation_overlapping_TFBS_AGI
            ):
                # create dict
                TFBS_AGI_single_dict = {}
                TFBS_AGI = mutation_overlapping_TFBS_AGI[
                    f"{mutation_type}{pos_index+1}"
                ]
                TFBS_AGI_dict[pos] = TFBS_AGI
        # append to output dict
        # print(TFBS_AGI_dict)
        output_dict.update(TFBS_AGI_dict.copy())


def add_overlapping_TFBS_dict_to_default_dict(
    overlapping_TFBS_dict, default_dict_key, default_dict
):
    if overlapping_TFBS_dict == {}:
        # print(overlapping_TFBS_dict)
        overlapping_TFBS_dict = np.nan
        default_dict[default_dict_key] = overlapping_TFBS_dict

    else:
        # sort overlapping_TFBS_dict by key (genomic position)
        overlapping_TFBS_dict = dict(sorted(overlapping_TFBS_dict.items()))
        default_dict[default_dict_key] = overlapping_TFBS_dict
        # print(overlapping_TFBS_dict)
    return default_dict


def flatten_mutations(mutations_df, output_folder, gene_name):
    """function to put all of the mutations in a plant line in one row per allele"""
    # create a new df with only one row per plant ID
    one_row_each = mutations_df.plant_ID.drop_duplicates(
        keep="first"
    ).to_frame()
    # turn into a dict
    one_row_each_dict = one_row_each.to_dict(orient="records")
    # create dictionary of the mutations_df
    # mutations_df_dict = mutations_df.to_dict(orient="records")
    # create empty list
    flattened_mutations_list = []
    # iterate over mutations_df_dict adding plant IDs as keys and mutated guide sites as values
    for row_dict in one_row_each_dict:
        # create default dict of row
        dd = defaultdict(list, row_dict)
        plant_ID = row_dict["plant_ID"]
        # filter mutations_df to only include that plant ID
        filtered = mutations_df[mutations_df.plant_ID == plant_ID]
        # turn into a dict
        filtered_dict = filtered.to_dict(orient="records")
        # get platename, library, primers, genotype
        chromosome = filtered_dict[0]["chr"]
        platename = filtered_dict[0]["platename"]
        library = filtered_dict[0]["library"]
        first_reaction_primers = filtered_dict[0]["first_reaction_primers"]
        second_reaction_primers = filtered_dict[0]["second_reaction_primers"]
        insertion_genomic_positions_combined = []
        insertion_positions_relative_to_TSS_combined = []
        insertion_overlapping_TFBS_AGI_combined = {}
        insertion_overlapping_TFBS_family_combined = {}

        deletion_genomic_positions_combined = []
        deletion_positions_relative_to_TSS_combined = []
        deletion_overlapping_TFBS_AGI_combined = {}
        deletion_overlapping_TFBS_family_combined = {}

        substitution_genomic_positions_combined = []
        substitution_positions_relative_to_TSS_combined = []
        substitution_overlapping_TFBS_AGI_combined = {}
        substitution_overlapping_TFBS_family_combined = {}

        guides = []
        mutation_types = []
        genotypes = []
        genotype = ""

        # iterate through dict
        for row in filtered_dict:
            guide = row["guide"]
            # add to guide list
            guides += [guide]
            mutation_type = row["mutation_type"]
            # add to list
            mutation_types += [mutation_type]
            # add to overall mutation type
            # add to list
            genotype = row["genotype"]
            genotypes += [genotype]

            if "insertion" in mutation_type:
                # convert genomic positions from string to list
                insertion_genomic_positions = (
                    row["insertion_genomic_positions"].strip("][").split(", ")
                )

                # append to insertion_genomic_positions_combined
                insertion_genomic_positions_combined.append(
                    list(insertion_genomic_positions)
                )

                # for pos in insertion_genomic_positions:
                #     #add insertion genomic positions
                #     #print(pos)
                #     insertion_genomic_positions_combined += [pos]
                # convert genomic positions from string to list
                insertion_positions_relative_to_TSS = (
                    row["insertion_positions_relative_to_TSS"]
                    .strip("][")
                    .split(", ")
                )
                for pos in insertion_positions_relative_to_TSS:
                    # add insertion genomic positions
                    insertion_positions_relative_to_TSS_combined += [pos]
                # generate dictionary of genomic positions with overlapping AGIs
                make_overlapping_TFBS_dict(
                    row["insertion_overlapping_TFBS_AGI"],
                    insertion_genomic_positions,
                    insertion_overlapping_TFBS_AGI_combined,
                    "insertion",
                )
                # generate dictionary of genomic positions with overlapping TFBS families
                make_overlapping_TFBS_dict(
                    row["insertion_overlapping_TFBS_family"],
                    insertion_genomic_positions,
                    insertion_overlapping_TFBS_family_combined,
                    "insertion",
                )

            if "deletion" in mutation_type:
                # convert genomic positions from string to list
                deletion_genomic_positions = (
                    row["deletion_genomic_positions"].strip("][").split(", ")
                )
                for pos in deletion_genomic_positions:
                    # add insertion genomic positions
                    deletion_genomic_positions_combined += [pos]
                # convert genomic positions from string to list
                deletion_positions_relative_to_TSS = (
                    row["deletion_positions_relative_to_TSS"]
                    .strip("][")
                    .split(", ")
                )
                for pos in deletion_positions_relative_to_TSS:
                    # add insertion genomic positions
                    deletion_positions_relative_to_TSS_combined += [pos]
                # generate dictionary of genomic positions with overlapping AGIs
                make_overlapping_TFBS_dict(
                    row["deletion_overlapping_TFBS_AGI"],
                    deletion_genomic_positions,
                    deletion_overlapping_TFBS_AGI_combined,
                    "deletion",
                )
                # generate dictionary of genomic positions with overlapping TFBS families
                make_overlapping_TFBS_dict(
                    row["deletion_overlapping_TFBS_family"],
                    deletion_genomic_positions,
                    deletion_overlapping_TFBS_family_combined,
                    "deletion",
                )

            if "substitution" in mutation_type:
                # convert genomic positions from string to list
                substitution_genomic_positions = (
                    row["substitution_genomic_positions"]
                    .strip("][")
                    .split(", ")
                )
                for pos in substitution_genomic_positions:
                    # add insertion genomic positions
                    substitution_genomic_positions_combined += [pos]
                # convert genomic positions from string to list
                substitution_positions_relative_to_TSS = (
                    row["substitution_positions_relative_to_TSS"]
                    .strip("][")
                    .split(", ")
                )
                for pos in substitution_positions_relative_to_TSS:
                    # add insertion genomic positions
                    substitution_positions_relative_to_TSS_combined += [pos]
                # generate dictionary of genomic positions with overlapping AGIs
                make_overlapping_TFBS_dict(
                    row["substitution_overlapping_TFBS_AGI"],
                    substitution_genomic_positions,
                    substitution_overlapping_TFBS_AGI_combined,
                    "substitution",
                )
                # generate dictionary of genomic positions with overlapping TFBS families
                make_overlapping_TFBS_dict(
                    row["substitution_overlapping_TFBS_family"],
                    substitution_genomic_positions,
                    substitution_overlapping_TFBS_family_combined,
                    "substitution",
                )

        # if empty list, change to nan
        if insertion_genomic_positions_combined == []:
            insertion_genomic_positions_combined = np.nan
            insertion_sizes = np.nan
            dd[
                "insertion_genomic_positions_combined"
            ] = insertion_genomic_positions_combined
            dd["insertion_sizes"] = insertion_sizes
        else:
            # insertions = [int(l) for l in list(np.unique(insertion_genomic_positions_combined))]
            insertions = list(insertion_genomic_positions_combined)
            # remove duplicates
            insertions.sort()
            insertions_no_dups = list(
                insertions for insertions, _ in itertools.groupby(insertions)
            )
            # make integers [[int(float(j)) for j in i] for i in x]
            insertions_no_dups = [
                [int(float(j)) for j in i] for i in insertions_no_dups
            ]

            # print(insertions)
            dd["insertion_genomic_positions_combined"] = insertions_no_dups
            # insertions no
            # record length of mutations
            insertion_sizes = find_insertions_lengths(insertions_no_dups)
            dd["insertion_sizes"] = insertion_sizes
        # print(insertion_genomic_positions_combined)
        if insertion_positions_relative_to_TSS_combined == []:
            insertion_positions_relative_to_TSS_combined = np.nan
            dd[
                "insertion_positions_relative_to_TSS_combined"
            ] = insertion_positions_relative_to_TSS_combined
        else:
            # make integers
            dd["insertion_positions_relative_to_TSS_combined"] = [
                int(l)
                for l in list(
                    np.unique(insertion_positions_relative_to_TSS_combined)
                )
            ]
        # add_overlapping_TFBS_dict_to_default_dict
        dd = add_overlapping_TFBS_dict_to_default_dict(
            insertion_overlapping_TFBS_AGI_combined,
            "insertion_overlapping_TFBS_AGI_combined",
            dd,
        )
        dd = add_overlapping_TFBS_dict_to_default_dict(
            insertion_overlapping_TFBS_family_combined,
            "insertion_overlapping_TFBS_family_combined",
            dd,
        )

        if deletion_genomic_positions_combined == []:
            deletion_genomic_positions_combined = np.nan
            deletion_sizes = np.nan
            dd[
                "deletion_genomic_positions_combined"
            ] = deletion_genomic_positions_combined
            dd["deletion_sizes"] = deletion_sizes
        else:
            # make integers
            deletions = [
                int(l)
                for l in list(np.unique(deletion_genomic_positions_combined))
            ]

            dd["deletion_genomic_positions_combined"] = deletions
            # record length of mutations
            deletion_sizes = find_mutation_lengths(deletions)
            dd["deletion_sizes"] = deletion_sizes

        if deletion_positions_relative_to_TSS_combined == []:
            deletion_positions_relative_to_TSS_combined = np.nan
            dd[
                "deletion_positions_relative_to_TSS_combined"
            ] = deletion_positions_relative_to_TSS_combined
        else:
            # make integers
            dd["deletion_positions_relative_to_TSS_combined"] = [
                int(l)
                for l in list(
                    np.unique(deletion_positions_relative_to_TSS_combined)
                )
            ]
        # add_overlapping_TFBS_dict_to_default_dict
        dd = add_overlapping_TFBS_dict_to_default_dict(
            deletion_overlapping_TFBS_AGI_combined,
            "deletion_overlapping_TFBS_AGI_combined",
            dd,
        )
        dd = add_overlapping_TFBS_dict_to_default_dict(
            deletion_overlapping_TFBS_family_combined,
            "deletion_overlapping_TFBS_family_combined",
            dd,
        )

        if substitution_genomic_positions_combined == []:
            substitution_genomic_positions_combined = np.nan
            substitution_sizes = np.nan
            dd[
                "substitution_genomic_positions_combined"
            ] = substitution_genomic_positions_combined
            dd["substitution_sizes"] = substitution_sizes
        else:
            # make integers
            substitutions = [
                int(l)
                for l in list(
                    np.unique(substitution_genomic_positions_combined)
                )
            ]

            dd["substitution_genomic_positions_combined"] = substitutions
            # record length of mutations
            substitution_sizes = find_mutation_lengths(substitutions)
            dd["substitution_sizes"] = substitution_sizes

        if substitution_positions_relative_to_TSS_combined == []:
            substitution_positions_relative_to_TSS_combined = np.nan
            dd[
                "substitution_positions_relative_to_TSS_combined"
            ] = substitution_positions_relative_to_TSS_combined
        else:
            # make integers
            dd["substitution_positions_relative_to_TSS_combined"] = [
                int(l)
                for l in list(
                    np.unique(substitution_positions_relative_to_TSS_combined)
                )
            ]
        # add_overlapping_TFBS_dict_to_default_dict
        dd = add_overlapping_TFBS_dict_to_default_dict(
            substitution_overlapping_TFBS_AGI_combined,
            "substitution_overlapping_TFBS_AGI_combined",
            dd,
        )
        dd = add_overlapping_TFBS_dict_to_default_dict(
            substitution_overlapping_TFBS_family_combined,
            "substitution_overlapping_TFBS_family_combined",
            dd,
        )
        # make mutation types string
        mutation_types_unique = list(np.unique(mutation_types))
        # if contains '+', separate out
        # if length is 1, pass
        # make empty list
        mutation_types_unique_new = []
        if len(mutation_types_unique) == 1:
            mutation_types_unique_new = mutation_types_unique

        elif len(mutation_types_unique) > 1:
            for thing in mutation_types_unique:
                # iff contains '+' then split it on the +
                if "+" in thing:
                    split = thing.split("+")
                    # print(split)
                    # add split to list
                    mutation_types_unique_new += split
                else:
                    mutation_types_unique_new += [thing]
            # make unique
            mutation_types_unique_new = list(
                np.unique(mutation_types_unique_new)
            )

        # add overall genotype
        mutation_genotypes = list(np.unique(genotypes))
        # if length is one, use that genotype
        if len(mutation_genotypes) == 1:
            genotype = mutation_genotypes[0]
        elif len(mutation_genotypes) > 1:
            # if contains chimeric then chimeric
            if "chimeric" in mutation_genotypes:
                genotype = "chimeric"
            # else if contains biallelic then biallelic
            elif "biallelic" in mutation_genotypes:
                genotype = "biallelic"
            # else if contains heterzygous then heterozygous
            elif "heterozygous" in mutation_genotypes:
                genotype = "heterozygous"
            else:
                print(f"error, unknown genotype {mutation_genotypes}")

        # add new values to ddict
        dd["chr"] = chromosome
        dd["platename"] = platename
        dd["library"] = library
        dd["first_reaction_primers"] = first_reaction_primers
        dd["second_reaction_primers"] = second_reaction_primers
        dd["guides"] = guides
        dd["mutation_type"] = mutation_types_unique_new
        # unique values in lists
        dd["mutation_types"] = list(np.unique(mutation_types))
        dd["genotype"] = genotype
        dd["mutation_genotypes"] = list(np.unique(genotypes))
        # dd['insertion_genomic_positions_combined'] = list(np.unique(insertion_genomic_positions_combined))
        # dd['insertion_positions_relative_to_TSS_combined'] = list(np.unique(insertion_positions_relative_to_TSS_combined))
        # dd['deletion_genomic_positions_combined'] = list(np.unique(deletion_genomic_positions_combined))
        # dd['deletion_positions_relative_to_TSS_combined'] = list(np.unique(deletion_positions_relative_to_TSS_combined))
        # dd['substitution_genomic_positions_combined'] = list(np.unique(substitution_genomic_positions_combined))
        # dd['substitution_positions_relative_to_TSS_combined'] = list(np.unique(substitution_positions_relative_to_TSS_combined))

        # add ddict to list of ddicts
        flattened_mutations_list.append(dd)

    # turn ddict into df
    df_single_rows = pd.DataFrame.from_dict(flattened_mutations_list)

    # add plant_ID ID and number column that's an integer for sorting on
    df_single_rows["ID"] = df_single_rows.plant_ID.str.split("-").str[0]
    df_single_rows["ID_number"] = (
        df_single_rows.plant_ID.str.split("-").str[1].astype(int)
    )

    df_single_rows.sort_values(
        ["ID", "ID_number"], ascending=True, inplace=True
    )

    # change column order
    df_single_rows = df_single_rows[
        [
            "plant_ID",
            "chr",
            "platename",
            "library",
            "first_reaction_primers",
            "second_reaction_primers",
            "genotype",
            "mutation_type",
            "guides",
            "mutation_genotypes",
            "insertion_genomic_positions_combined",
            "insertion_positions_relative_to_TSS_combined",
            "insertion_sizes",
            "insertion_overlapping_TFBS_AGI_combined",
            "insertion_overlapping_TFBS_family_combined",
            "deletion_genomic_positions_combined",
            "deletion_positions_relative_to_TSS_combined",
            "deletion_sizes",
            "deletion_overlapping_TFBS_AGI_combined",
            "deletion_overlapping_TFBS_family_combined",
            "substitution_genomic_positions_combined",
            "substitution_positions_relative_to_TSS_combined",
            "substitution_sizes",
            "substitution_overlapping_TFBS_AGI_combined",
            "substitution_overlapping_TFBS_family_combined",
        ]
    ]

    # save this df
    df_single_rows.to_csv(
        f"{output_folder}/{gene_name}_TFBSoverlapping_genotyped_only_mutated_flattened.tsv",
        sep="\t",
        index=False,
        header=1,
    )


def main(args):
    # parse arguments
    args = parse_args(args)

    # read in files
    mutations_df = read_in_files(args.mutations_tsv)
    # read_in_files(args.input_folder,args.guide_pairs,args.gene_name)

    # produce a tsv with all mutations for one plant line in one row with no duplicate mutations
    flatten_mutations(mutations_df, args.output_folder, args.gene_name)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
