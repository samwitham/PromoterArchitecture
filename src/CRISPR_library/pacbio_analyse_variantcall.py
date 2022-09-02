import argparse
import fnmatch
import io
import os
import re
import subprocess

import numpy as np
import pandas as pd
from pyfaidx import Fasta


def parse_args(args):
    parser = argparse.ArgumentParser(description="pacbio_analyse_variantcall")
    parser.add_argument(
        "gene_folder",
        type=str,
        help="Input location of CRISPRESSO2 analysis",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output folder location ending with a forward slash",
    )
    parser.add_argument(
        "reference_gene_fasta",
        type=str,
        help="Input location of reference gene fasta file",
    )
    # parser.add_argument(
    #     "all_genes_longest_region_fasta",#reference_fasta
    #     type=str,
    #     help="Input location of fasta file containing sequences of the longest amplified region of all genes of interest",
    # )
    parser.add_argument(
        "all_genes_reference_promoter_bed",
        type=str,
        help="Input location of reference promoter bed file containing the longest amplified region of all genes of interest",
    )
    parser.add_argument(
        "all_promoters_bed",
        type=str,
        help="Input location of bed file containing all promoters in Arabidopsis genome",
    )
    # parser.add_argument(
    #     "mapped_motifs_bed",
    #     type=str,
    #     help="Input location of bed file containing all TFBSs scanned using FIMO in all Arabidopsis promoters",
    # )
    parser.add_argument(
        "plant_IDs",
        type=str,
        help="Input location of tsv file containing all plant IDs that were sequenced",
    )
    parser.add_argument(
        "gene_name",
        type=str,
        help="Name of gene being analysed",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


# code from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoCOREResources.pyx
# I converted it to pure python
def find_indels_substitutions(read_seq_al, ref_seq_al, _include_indx):
    re_find_indels = re.compile("(-*-)")

    # ref_positions holds the indices for which positions map back to the original reference
    # for example,
    #     1 2 3 4 5 6 7 8
    # ref A A T T G G C C
    #
    # and for a given alignment
    # ref A A T T - G G C C
    # aln A - T T T G G C C
    #     1 2 3 4-4 5 6 7 8 <ref positions. Note that the negative values/indices represent places that don't map back to the original reference
    ref_positions = []
    all_substitution_positions = []
    substitution_positions = []
    all_substitution_values = []
    substitution_values = []

    all_deletion_positions = []
    deletion_positions = []
    deletion_coordinates = []
    deletion_sizes = []

    start_deletion = (
        -1
    )  # the -1 value indicates that there currently isn't a deletion

    all_insertion_positions = []
    all_insertion_left_positions = []
    insertion_positions = []
    insertion_coordinates = []
    insertion_sizes = []

    start_insertion = (
        -1
    )  # the -1 value indicates that there currently isn't an insertion

    seq_len = len(ref_seq_al)
    include_indx_set = set(_include_indx)
    nucSet = set(["A", "T", "C", "G", "N"])
    idx = 0

    current_insertion_size = 0
    for idx_c, c in enumerate(ref_seq_al):
        # print(idx_c)
        if c != "-":
            ref_positions.append(idx)
            if (
                ref_seq_al[idx_c] != read_seq_al[idx_c]
                and read_seq_al[idx_c] != "-"
                and read_seq_al[idx_c] != "N"
            ):
                all_substitution_positions.append(idx)
                all_substitution_values.append(read_seq_al[idx_c])
                if idx in _include_indx:
                    substitution_positions.append(idx)
                    substitution_values.append(read_seq_al[idx_c])
            if start_insertion != -1:  # this is the end of an insertion
                all_insertion_left_positions.append(start_insertion)
                all_insertion_positions.append(start_insertion)
                all_insertion_positions.append(idx)
                if (
                    start_insertion in include_indx_set
                    and idx in include_indx_set
                ):
                    insertion_coordinates.append((start_insertion, idx))
                    insertion_positions.append(start_insertion)
                    insertion_positions.append(idx)
                    insertion_sizes.append(current_insertion_size)
                start_insertion = -1
            current_insertion_size = 0
            idx += 1
        else:  # the current ref position is -
            if idx == 0:
                ref_positions.append(-1)
            else:
                ref_positions.append(-idx)
            if (
                idx > 0 and start_insertion == -1
            ):  # this is the first index of an insertion
                start_insertion = idx - 1
            current_insertion_size += 1

        if (
            read_seq_al[idx_c] == "-" and start_deletion == -1
        ):  # this is the first part of a deletion
            if idx_c - 1 > 0:
                start_deletion = ref_positions[idx_c]
            else:
                start_deletion = 0
        elif (
            read_seq_al[idx_c] != "-" and start_deletion != -1
        ):  # this is the end of a deletion
            end_deletion = ref_positions[idx_c]
            all_deletion_positions.extend(range(start_deletion, end_deletion))
            if include_indx_set.intersection(
                range(start_deletion, end_deletion)
            ):
                deletion_positions.extend(range(start_deletion, end_deletion))
                deletion_coordinates.append((start_deletion, end_deletion))
                deletion_sizes.append(end_deletion - start_deletion)
            start_deletion = -1

    if start_deletion != -1:
        end_deletion = ref_positions[seq_len - 1]
        all_deletion_positions.extend(range(start_deletion, end_deletion))
        if include_indx_set.intersection(range(start_deletion, end_deletion)):
            deletion_positions.extend(range(start_deletion, end_deletion))
            deletion_coordinates.append((start_deletion, end_deletion))
            deletion_sizes.append(end_deletion - start_deletion)

    substitution_n = len(all_substitution_positions)
    deletion_n = len(all_deletion_positions)
    insertion_n = len(all_insertion_positions)

    return {
        "all_insertion_positions": all_insertion_positions,
        "all_insertion_left_positions": all_insertion_left_positions,
        "insertion_positions": insertion_positions,
        "insertion_coordinates": insertion_coordinates,
        "insertion_sizes": insertion_sizes,
        "insertion_n": insertion_n,
        "all_deletion_positions": all_deletion_positions,
        "deletion_positions": deletion_positions,
        "deletion_coordinates": deletion_coordinates,
        "deletion_sizes": deletion_sizes,
        "deletion_n": deletion_n,
        "all_substitution_positions": all_substitution_positions,
        "substitution_positions": substitution_positions,
        "all_substitution_values": np.array(all_substitution_values),
        "substitution_values": np.array(substitution_values),
        "substitution_n": substitution_n,
        "ref_positions": ref_positions,
    }


def find_guide_position_in_gene(
    gene, cut_site_reference_seq, guide_position, gene_fasta_location
):
    """function to find the relative promoter position of the current guide site position in the reference gene"""
    fasta_location = gene_fasta_location
    # read in fasta file
    fasta = Fasta(fasta_location)
    # get promoter sequence
    prom_seq = str(fasta[f"{gene}_promoter"])
    # find position of substring in string
    # remove 'N's from sequence
    cut_site_reference_seq = cut_site_reference_seq.replace("N", "")
    start_location = prom_seq.index(cut_site_reference_seq)
    site_promoter_position = start_location + guide_position
    return site_promoter_position


def get_reference_promoter_genomic_positions(
    reference_promoter_bed, all_promoters_bed
):
    """function to get the reference promoter genomic positions and to get the tSS genomic position"""
    # read in reference_promoter_bed
    reference_promoter_df = pd.read_table(
        reference_promoter_bed, sep="\t", header=None
    )
    cols = [
        "chr",
        "start",
        "stop",
        "promoter_name",
        "score",
        "strand",
    ]
    reference_promoter_df.columns = cols

    # read in all_promoters_bed
    all_promoters_df = pd.read_table(all_promoters_bed, sep="\t", header=None)
    cols2 = [
        "chr",
        "start",
        "stop",
        "AGI",
        "dot",
        "strand",
        "source",
        "type",
        "dot2",
        "attributes",
    ]
    all_promoters_df.columns = cols2

    # AGIs dictionary of the four genes
    AGI_dict = {
        "DREB26": "AT1G21910",
        "NLP7": "AT4G24020",
        "ARF18": "AT3G61830",
        "ARF9": "AT4G23980",
    }
    # make empty dict
    reference_genomic_positions = {}
    # get promoter genomic positions
    for gene, AGI in AGI_dict.items():
        promoter_genomic_pos_chromosome = int(
            reference_promoter_df[
                reference_promoter_df.promoter_name == f"{gene}_promoter"
            ].chr
        )
        promoter_genomic_pos_start = int(
            reference_promoter_df[
                reference_promoter_df.promoter_name == f"{gene}_promoter"
            ].start
        )
        promoter_genomic_pos_stop = int(
            reference_promoter_df[
                reference_promoter_df.promoter_name == f"{gene}_promoter"
            ].stop
        )

        # get TSS position
        TSS_pos = int(all_promoters_df[all_promoters_df.AGI == AGI].stop)
        # write to dict
        reference_genomic_positions[gene] = [
            AGI,
            promoter_genomic_pos_chromosome,
            promoter_genomic_pos_start,
            promoter_genomic_pos_stop,
            TSS_pos,
        ]
    return reference_genomic_positions


def find_genomic_position(
    reference_genomic_positions, gene, promoter_position
):
    """get the genomic position of the cut site"""
    positions_list = reference_genomic_positions[gene]
    # first get genomic position of the cut site
    promoter_genomic_pos_chromosome = positions_list[1]
    promoter_genomic_pos_start = positions_list[2]
    promoter_genomic_pos_stop = positions_list[3]
    # get the genomic location of the input promoter position
    genomic_position = promoter_position + promoter_genomic_pos_start
    # get position of cut site relative to TSS position
    # get TSS position
    TSS_pos = positions_list[4]
    position_relative_to_TSS = genomic_position - TSS_pos

    return (
        position_relative_to_TSS,
        genomic_position,
        promoter_genomic_pos_chromosome,
    )


def check_guide(
    root_dir,
    output,
    gene,
    reference_genomic_positions,
    plantIDs,
    gene_fasta_location,
):
    """read in the Alleles_frequency_table txt files"""
    # create the output file
    cols = [
        "chr",
        "platename",
        "library",
        "first_reaction_primers",
        "second_reaction_primers",
        "guide",
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
    ]  #'mutation_size',
    output_df = pd.DataFrame(columns=cols)

    # read in plant IDs table
    plant_ID_cols = [
        "plant_ID",
        "F_primer",
        "R_primer",
        "barcode_no",
        "platename",
    ]
    plantID_df = pd.read_table(plantIDs, sep="\t", header=0)
    plantID_df.columns = plant_ID_cols

    for subdir, dirs, files in os.walk(root_dir):
        for filename in fnmatch.filter(
            files, f"*Alleles_frequency_table_around_{gene}_*.txt"
        ):
            # read in as df
            df = pd.read_table(
                os.path.join(subdir, filename), sep="\t", header=0
            )
            # filter out rows with reads less than 8
            filtered_reads = df[df["#Reads"] >= 8]
            # save each row to df
            for index, row in filtered_reads.iterrows():

                # check if aligned and reference are different
                Aligned_Sequence = row.Aligned_Sequence
                Reference_Sequence = row.Reference_Sequence

                # split filename on .
                partitioned_string = filename.partition(".")
                # get second PCR reaction primers
                second_reaction_primers = partitioned_string[0]
                # get guide name
                guide = re.search(
                    f"[_]+({gene}(.*))",
                    partitioned_string[2].partition(".")[0],
                )[1]
                # get first reaction PCR primers and library number
                library_primers = re.search("_bc(.*)", subdir)[0]
                both_primers = re.findall("SW(\d*)", library_primers)
                # if 2 found, create first_reaction_primers value
                if len(both_primers) == 2:
                    first_reaction_primers = (
                        f"SW{both_primers[0]}_SW{both_primers[1]}"
                    )
                else:
                    first_reaction_primers = "NA"

                library = re.search("[^_bc]+(\d*)", library_primers)[0]
                read_number = row["#Reads"]
                read_percentage = row["%Reads"]

                # convert library number to different number (from 1017 to 1, 1018 to 2 etc)
                if int(library) == 1017:
                    new_library = 1
                elif int(library) == 1018:
                    new_library = 2
                elif int(library) == 1019:
                    new_library = 3
                elif int(library) == 1020:
                    new_library = 4
                elif int(library) == 1021:
                    new_library = 5
                elif int(library) == 1022:
                    new_library = 6
                # print(new_library)
                platename = f"p{new_library}{gene}"

                # remove dashes from string
                Aligned_Sequence_no_dashes = Aligned_Sequence.replace("-", "")
                Reference_Sequence_no_dashes = Reference_Sequence.replace(
                    "-", ""
                )
                if Aligned_Sequence == Reference_Sequence:
                    mutation_type = "None"
                    mutation_size = "NA"
                    insertion_positions = "NA"
                    deletion_positions = "NA"
                    substitution_positions = "NA"
                elif Aligned_Sequence != Reference_Sequence:
                    # find insertions, mutations, deletions
                    # get length of reference sequence for the index

                    # add 'N's to each end of the sequences so that edge case deletions/insertions aren't ignored
                    Aligned_Sequence_Ns = "N" + Aligned_Sequence + "N"
                    Reference_Sequence_Ns = "N" + Reference_Sequence + "N"
                    ref_index = [range(0, len(Reference_Sequence_Ns), 1)]
                    indels = find_indels_substitutions(
                        Aligned_Sequence_Ns, Reference_Sequence_Ns, ref_index
                    )
                    # NEED TO SUBTRACT ONE TO ALL POSITIONS BECAUSE USED AN ADDED 'N' BEFORE REFERENCE SEQUENCE WHEN FINDING INDELS

                    # print(indels)
                    # as well as labelling mutation types, get relative position of the mutation/mutations in the 40bp guide window
                    if indels["insertion_n"] != 0:
                        # NEED TO SUBTRACT ONE TO ALL POSITIONS BECAUSE USED AN ADDED 'N' BEFORE REFERENCE SEQUENCE WHEN FINDING INDELS
                        insertion_positions = [
                            x - 1 for x in indels["all_insertion_positions"]
                        ]
                        if indels["deletion_n"] != 0:
                            deletion_positions = [
                                x - 1 for x in indels["all_deletion_positions"]
                            ]
                            if indels["substitution_n"] != 0:
                                substitution_positions = [
                                    x - 1
                                    for x in indels[
                                        "all_substitution_positions"
                                    ]
                                ]
                                mutation_type = (
                                    "insertion+deletion+substitution"
                                )

                            elif indels["substitution_n"] == 0:
                                substitution_positions = "NA"
                                mutation_type = "insertion+deletion"
                        elif indels["deletion_n"] == 0:
                            deletion_positions = "NA"
                            if indels["substitution_n"] != 0:
                                substitution_positions = [
                                    x - 1
                                    for x in indels[
                                        "all_substitution_positions"
                                    ]
                                ]
                                mutation_type = "insertion+substitution"
                            elif indels["substitution_n"] == 0:
                                substitution_positions = "NA"
                                mutation_type = "insertion"
                    elif indels["insertion_n"] == 0:
                        insertion_positions = "NA"
                        if indels["deletion_n"] != 0:
                            deletion_positions = [
                                x - 1 for x in indels["all_deletion_positions"]
                            ]
                            if indels["substitution_n"] != 0:
                                substitution_positions = [
                                    x - 1
                                    for x in indels[
                                        "all_substitution_positions"
                                    ]
                                ]
                                mutation_type = "deletion+substitution"
                            elif indels["substitution_n"] == 0:
                                substitution_positions = "NA"
                                mutation_type = "deletion"
                        elif indels["deletion_n"] == 0:
                            deletion_positions = "NA"
                            if indels["substitution_n"] != 0:
                                substitution_positions = [
                                    x - 1
                                    for x in indels[
                                        "all_substitution_positions"
                                    ]
                                ]
                                mutation_type = "substitution"
                            elif indels["substitution_n"] == 0:
                                substitution_positions = "NA"
                                mutation_type = "None"
                # get length of reference sequence
                ref_length = len(Reference_Sequence_no_dashes)

                # get cut site position in whole promoter
                cut_site_promoter_position = find_guide_position_in_gene(
                    gene,
                    Reference_Sequence_no_dashes,
                    ref_length // 2,
                    gene_fasta_location,
                )
                # get cut site genomic position in whole promoter
                (
                    cut_site_position_relative_to_TSS,
                    cut_site_genomic_position,
                    promoter_genomic_pos_chromosome,
                ) = find_genomic_position(
                    reference_genomic_positions,
                    gene,
                    cut_site_promoter_position,
                )

                # add distance from guide cut site column to df

                # print(ref_length)
                if insertion_positions == "NA":
                    insertion_cut_site_distance = "NA"
                    insertion_positions_relative_to_TSS = "NA"
                    insertion_genomic_positions = "NA"
                else:
                    # make list of insertion_cut_site_distances (distance of insertion from cut site)
                    insertion_cut_site_distance = [
                        i - ref_length // 2 for i in insertion_positions
                    ]
                    # make list of insertion cut site genomic position and also the position relative to the Araport TSS
                    insertion_positions_relative_to_TSS = []
                    insertion_genomic_positions = []
                    for i in insertion_positions:
                        site_promoter_position = find_guide_position_in_gene(
                            gene,
                            Reference_Sequence_no_dashes,
                            i,
                            gene_fasta_location,
                        )
                        (
                            position_relative_to_TSS,
                            genomic_position,
                            promoter_genomic_pos_chromosome,
                        ) = find_genomic_position(
                            reference_genomic_positions,
                            gene,
                            site_promoter_position,
                        )
                        insertion_positions_relative_to_TSS.append(
                            position_relative_to_TSS
                        )
                        insertion_genomic_positions.append(genomic_position)

                if deletion_positions == "NA":
                    deletion_cut_site_distance = "NA"
                    deletion_positions_relative_to_TSS = "NA"
                    deletion_genomic_positions = "NA"

                else:
                    # make list of deletion positions (distance of deletion from cut site)
                    deletion_cut_site_distance = [
                        i - ref_length // 2 for i in deletion_positions
                    ]
                    # make list of deletion cut site genomic position and also the position relative to the Araport TSS
                    deletion_positions_relative_to_TSS = []
                    deletion_genomic_positions = []
                    for i in deletion_positions:
                        site_promoter_position = find_guide_position_in_gene(
                            gene,
                            Reference_Sequence_no_dashes,
                            i,
                            gene_fasta_location,
                        )
                        (
                            position_relative_to_TSS,
                            genomic_position,
                            promoter_genomic_pos_chromosome,
                        ) = find_genomic_position(
                            reference_genomic_positions,
                            gene,
                            site_promoter_position,
                        )
                        deletion_positions_relative_to_TSS.append(
                            position_relative_to_TSS
                        )
                        deletion_genomic_positions.append(genomic_position)

                if substitution_positions == "NA":
                    substitution_cut_site_distance = "NA"
                    substitution_positions_relative_to_TSS = "NA"
                    substitution_genomic_positions = "NA"
                else:
                    # make list of substitution positions (distance of substitution from cut site)
                    substitution_cut_site_distance = [
                        i - ref_length // 2 for i in substitution_positions
                    ]
                    # make list of substitution cut site genomic position and also the position relative to the Araport TSS
                    substitution_positions_relative_to_TSS = []
                    substitution_genomic_positions = []
                    for i in substitution_positions:
                        site_promoter_position = find_guide_position_in_gene(
                            gene,
                            Reference_Sequence_no_dashes,
                            i,
                            gene_fasta_location,
                        )
                        (
                            position_relative_to_TSS,
                            genomic_position,
                            promoter_genomic_pos_chromosome,
                        ) = find_genomic_position(
                            reference_genomic_positions,
                            gene,
                            site_promoter_position,
                        )
                        substitution_positions_relative_to_TSS.append(
                            position_relative_to_TSS
                        )
                        substitution_genomic_positions.append(genomic_position)

                # append list of values to output_df
                list = [
                    promoter_genomic_pos_chromosome,
                    platename,
                    library,
                    first_reaction_primers,
                    second_reaction_primers,
                    guide,
                    Aligned_Sequence,
                    Reference_Sequence,
                    mutation_type,
                    read_number,
                    read_percentage,
                    insertion_positions,
                    deletion_positions,
                    substitution_positions,
                    insertion_cut_site_distance,
                    deletion_cut_site_distance,
                    substitution_cut_site_distance,
                    cut_site_promoter_position,
                    insertion_positions_relative_to_TSS,
                    insertion_genomic_positions,
                    deletion_positions_relative_to_TSS,
                    deletion_genomic_positions,
                    substitution_positions_relative_to_TSS,
                    substitution_genomic_positions,
                ]  # mutation_size
                output_df.loc[len(output_df)] = list

    # merge plant id forward and reverse first reaction primers into one column to match the output df format
    ## Keep only the numbers before the hyphen.
    plantID_df.F_primer = plantID_df.F_primer.str.split("-").str[0]
    plantID_df.R_primer = plantID_df.R_primer.str.split("-").str[0]
    # merge the two columns
    plantID_df["first_reaction_primers"] = (
        plantID_df.F_primer + "_" + plantID_df.R_primer
    )
    # print(plantID_df)
    plantID_df = plantID_df[
        ["plant_ID", "first_reaction_primers", "platename"]
    ]
    # print(plantID_df)

    # add plant IDs
    output_df_new = pd.merge(
        output_df,
        plantID_df,
        how="left",
        on=["first_reaction_primers", "platename"],
    )

    # add guide number column that's an integer for sorting on
    output_df_new["guide_number"] = (
        output_df_new.guide.str.split("guide").str[1].astype(int)
    )

    # write out the output_df
    output_df_new.to_csv(
        f"{output}{gene}/{gene}_includingduplicates.tsv",
        sep="\t",
        index=False,
        header=1,
    )
    # find duplicate values
    # columns to group by
    columns = [
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
    output_df_new[to_string] = output_df_new[to_string].astype(str)

    col_order = [
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

    # merge duplicate values, adding the sum of the read_number and read_percentage
    output_df_merged = output_df_new.groupby(
        columns, sort=False, as_index=False
    ).agg({"read_number": "sum", "read_percentage": "sum"})
    # put columns back in the original order
    output_df_merged = output_df_merged[col_order]

    # sort by platename, plant ID then by guide number
    output_df_merged.sort_values(
        ["platename", "plant_ID", "guide_number"], ascending=True, inplace=True
    )

    # write out the output_df
    output_df_merged.to_csv(
        f"{output}{gene}/{gene}_merged.tsv", sep="\t", index=False, header=1
    )


def main(args):
    # parse arguments
    args = parse_args(args)

    # make directory for the files to be exported to
    dirName = f"{args.output}{args.gene_name}/"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # get reference promoter genomic positions along with TSS position
    reference_genomic_positions = get_reference_promoter_genomic_positions(
        args.all_genes_reference_promoter_bed, args.all_promoters_bed
    )
    # print(reference_genomic_positions)
    # read in the Alleles_frequency_table txt files and check each guide location for mutations
    check_guide(
        args.gene_folder,
        args.output,
        args.gene_name,
        reference_genomic_positions,
        args.plant_IDs,
        args.reference_gene_fasta,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
