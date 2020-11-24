import argparse
import io
import os
from pprint import pformat

import pandas as pd
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from pybedtools import BedTool
from pyfaidx import Fasta


def parse_args(args):
    """define arguments"""
    # action="store_true" automatically creates a default value of False
    parser = argparse.ArgumentParser(description="extract_promoter")
    # parser.add_argument('directory_path', type=str, help='Location of base directory')
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "promoter_length",
        help="Max length of promoters to be extracted upstream of the Araport 11 TSS",
    )
    parser.add_argument(
        "--remove_bidirectional",
        help="Exclude potentially bidirectional promoters with upstream TSS >2000bp from TSS in opposite direction.",
        action="store_true",
    )
    parser.add_argument(
        "--prevent_overlapping_genes",
        help="Reduce size of promoters if needed until they are not overlapping other genes.",
        action="store_true",
    )
    parser.add_argument(
        "--fiveUTR",
        help="Extend promoters to the first start codon.",
        action="store_true",
    )

    # parser.add_argument('--only_open_chromatin', help='Reduce size of promoters if needed so that they only fall within open chromatin.', action="store_true")
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def fasta_chromsizes(genome, output_file):
    """extracts chromosome sizes in a form compatible with BedTools.flank"""

    the_genome = Fasta(
        genome
    )  # using pyfaidx Fasta, parse the genome file as a fasta

    chromsizes = {}  # make dictionary called chromsizes
    for key in the_genome.keys():

        chromsizes[
            f"{key}"
        ] = f"({len(the_genome[key])})"  # add the chromosome name and length to dictionary
    # create empty string
    chromsizes_string = ""
    # iterate over chromsizes dictionary key/values. Add key, tab, value, newline to the string iteratively
    for k, v in chromsizes.items():
        chromsizes_string = chromsizes_string + f"{k}\t{v}\n"

    # write output file, deleting the parentheses created from the dictionary. This file is suitable for use in BedTools.flank
    with open(output_file, "w") as output:
        output.write(chromsizes_string.replace("(", "").replace(")", ""))


# def extract_genes(gene_gff,output_file):
#     """This function extracts all whole genes from a gff3 file, ignoring gene features, and adds them to an output file"""
#     #limit dictionary to genes
#     limit_info = dict(gff_type = ['gene'])
#     #open temporary buffer
#     tempbuffer = io.StringIO()
#     #open gff file, parse it, limiting to genes only. Save the file.
#     with open(gene_gff, 'r') as in_handle:
#             GFF.write(GFF.parse(in_handle, limit_info=limit_info),tempbuffer)

#     #go back to beginning of the buffer
#     tempbuffer.seek(0)

#     #now remove the unwanted annotations that were added by GFF.write eg. Chr1	annotation	remark	1	30425192	.	.	.	gff-version=3
#     #remove lines beginning with ##
#     with open(output_file, 'w') as newfile:
#         for line in tempbuffer:
#             line = line.strip() # removes hidden characters/spaces
#             if line[0] == "#":
#                 pass
#             else:
#                 #don't include lines beginning with ##
#                 newfile.write(line + '\n') #output to new file
#     #remove temporary buffer
#     tempbuffer.close()
#     #read in gff file
#     genes = pd.read_table(output_file, sep='\t', header=None)
#     cols2 = ['chr', 'source', 'type', 'start','stop','dot1','strand','dot2','attributes']
#     genes.columns = cols2
#     #remove all lines where source is annotation
#     no_annotation = genes[~genes.source.str.contains("annotation")]
#     #remove all lines that are not protein coding
#     protein_coding = no_annotation[no_annotation.attributes.str.contains('biotype=protein_coding')]
#     #create gff file containing no lines with annotation and no lines starting with ##
#     protein_coding.to_csv(output_file,index=False,sep='\t',header=None)


def extract_genes(gene_gff, output_overlapping, output_file):
    """This function extracts all whole genes from a gff3 file, ignoring gene features, and adds them to an output file. It also exports a file containing overlapping genes, and removes any overlapping genes from the final output (final output is protalso returns all genes that are protein coding regardless of overlap"""
    # limit dictionary to genes
    limit_info = dict(gff_type=["gene"])

    # open temporary buffer
    # output = open(temp, 'w')
    tempbuffer = io.StringIO()

    # open gff file, parse it, limiting to genes only. Save the file.
    with open(gene_gff, "r") as in_handle:
        GFF.write(GFF.parse(in_handle, limit_info=limit_info), tempbuffer)

    # go back to beginning of the buffer
    tempbuffer.seek(0)

    # now remove the unwanted annotations that were added by GFF.write eg. Chr1	annotation	remark	1	30425192	.	.	.	gff-version=3
    # remove lines beginning with ##
    with open(output_file, "w") as newfile:
        for line in tempbuffer:
            line = line.strip()  # removes hidden characters/spaces
            if line[0] == "#":
                pass
            else:
                # don't include lines beginning with ##
                newfile.write(line + "\n")  # output to new file
    # remove temporary buffer
    tempbuffer.close()

    # read in gff file
    genes = pd.read_table(output_file, sep="\t", header=None)
    cols2 = [
        "chr",
        "source",
        "type",
        "start",
        "stop",
        "dot1",
        "strand",
        "dot2",
        "attributes",
    ]
    genes.columns = cols2

    # remove all lines where source is annotation
    no_annotation = genes[~genes.source.str.contains("annotation")]

    # remove all lines that are not protein coding
    protein_coding = no_annotation[
        no_annotation.attributes.str.contains("biotype=protein_coding")
    ]

    # add AGI column
    gene_agi = protein_coding.assign(
        AGI=protein_coding.attributes.str.extract(r"ID=gene:(.*?)\;")
    )

    # remove all features which overlap:
    # make buffer of genes
    genes_buffer = io.StringIO()
    gene_agi.to_csv(genes_buffer, index=False, sep="\t", header=None)
    genes_buffer.seek(0)

    # create bedtools object of genes
    genes_bed = BedTool(genes_buffer)

    # merge protein coding genes with protein coding genes that they overlap
    # c = columns to apply function to
    # o = count number of merged promoters, name the first and last promoter that were merged
    merged = genes_bed.merge(c=10, o=["count_distinct", "first", "last"])

    # write to buffer
    merged_buffer = io.StringIO()
    merged_buffer.write(str(merged))
    merged_buffer.seek(0)

    # read in again as a dataframe
    overlapping = pd.read_table(merged_buffer, sep="\t", header=None)

    # name the columns
    cols3 = [
        "chr",
        "start",
        "stop",
        "number_of_overlaps",
        "first_overlap",
        "second_overlap",
    ]
    overlapping.columns = cols3

    # select only features made of more than one promoter that were merged as overlapping
    overlapping_only = overlapping[overlapping.number_of_overlaps >= 2]

    # save these overlapping genes to file
    overlapping_only.to_csv(
        output_overlapping, index=False, sep="\t", header=None
    )

    # keep only non-overlapping genes (they are labelled as having 1 overlap as they overlap themselves)
    final_genes_list = overlapping[overlapping.number_of_overlaps == 1]
    final_genes = gene_agi[gene_agi.AGI.isin(final_genes_list.first_overlap)]
    final_genes_noAGI = final_genes[cols2]

    # create gff file containing no lines with annotation and no lines starting with ##
    final_genes_noAGI.to_csv(output_file, index=False, sep="\t", header=None)

    return protein_coding


def add_promoter(genes_gff, chromsize, promoter_length):
    """This function adds a promoter of a certain length to each gene in the input file and exports an output pyBedTools object"""
    # output = open(output_location, 'w') #make output file with write capability
    # parse gff file containing only genes.
    genes = BedTool(genes_gff)
    # extract promoters upsteam using chromsize file and specified promoter length. r, no. of bp to add to end coordinate. s, based on strand.
    promoters = genes.flank(g=chromsize, l=promoter_length, r=0, s=True)

    return promoters


def remove_promoter_overlap(promoter_gff, all_genes_df, output_file):
    """function to create file containing promoters which overlap other genome features.
    Then create a df with only promoters which overlap other genes. Then shorten them so they are no longer overlapping.
    Then merge back with all the extracted promoters, keeping only the shortened ones"""
    # create buffer to save the protein_coding genes
    all_genes_gff = io.StringIO()
    all_genes_df.to_csv(all_genes_gff, sep="\t", header=None, index=False)
    all_genes_gff.seek(0)

    all_proms = BedTool(promoter_gff)  # read in files using BedTools
    features = BedTool(all_genes_gff)
    # report chromosome position of overlapping feature, along with the promoter which overlaps it (only reports the overlapping nucleotides, not the whole promoter length. Can use u=True to get whole promoter length)
    # f, the minimum overlap as fraction of A. F, nucleotide fraction of B (genes) that need to be overlapping with A (promoters)
    # wa, Write the original entry in A for each overlap.
    # wo,  Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported.
    # u, write original A entry only once even if more than one overlap
    intersect = all_proms.intersect(
        features, wo=True
    )  # could add u=True which indicates we want to see the promoters that overlap features in the genome
    # Write to output_file
    with open(output_file, "w") as output:
        # Each line in the file contains gff entry a and gff entry b that it overlaps plus the number of bp in the overlap so 19 columns
        output.write(str(intersect))
    # read in gff file
    overlapping_proms = pd.read_table(output_file, sep="\t", header=None)
    cols = [
        "chrA",
        "sourceA",
        "typeA",
        "startA",
        "stopA",
        "dot1A",
        "strandA",
        "dot2A",
        "attributesA",
        "chrB",
        "sourceB",
        "typeB",
        "startB",
        "stopB",
        "dot1B",
        "strandB",
        "dot2B",
        "attributesB",
        "bp_overlap",
    ]
    overlapping_proms.columns = cols
    # make new columns for the new start and stop for feature A
    # overlapping_proms['new_startA'] = overlapping_proms.startA
    # overlapping_proms['new_stopA'] = overlapping_proms.stopA
    # create dictionaary for AGI as key and starts and stops as values
    temp_dict = {}

    # iterate over rows
    for i, data in overlapping_proms.iterrows():
        # create AGI key
        key = (
            overlapping_proms.loc[i, "attributesA"]
            .split("ID=gene:")[1]
            .split(";")[0]
        )
        # if positive strand feature A, reduce length to the stop position of feature B + 1
        if overlapping_proms.loc[i, "strandA"] == "+":

            if key not in temp_dict:

                if (
                    overlapping_proms.loc[i, "stopA"]
                    <= overlapping_proms.loc[i, "stopB"]
                ):
                    start = overlapping_proms.loc[i, "stopB"] + 1
                    stop = overlapping_proms.loc[i, "stopB"] + 1
                    temp_dict[key] = [start, stop]
                else:
                    start = overlapping_proms.loc[i, "stopB"] + 1
                    # update dictionary
                    temp_dict[key] = [start, overlapping_proms.loc[i, "stopA"]]

            # if key is in the dictionary
            else:
                # if overlapping feature is downstream of the promoter or equal to promoter end, reduce both start and stop to the stop of the feature + 1
                if temp_dict[key][1] <= overlapping_proms.loc[i, "stopB"]:
                    start = overlapping_proms.loc[i, "stopB"] + 1
                    stop = overlapping_proms.loc[i, "stopB"] + 1
                    # overlapping_proms.drop(i, inplace=True)
                    temp_dict[key] = [start, stop]

                # else  reduce length of promoter to stop position of feature +1
                else:
                    # if promoter has already been shortened by another overlapping feature, pass
                    if temp_dict[key][0] > overlapping_proms.loc[i, "stopB"]:
                        pass
                    # else shorten the promoter
                    else:
                        start = overlapping_proms.loc[i, "stopB"] + 1
                        # update dictionary
                        temp_dict[key][0] = start

        # if negative strand feature A, reduce length of promoter to start position of feature -1
        elif overlapping_proms.loc[i, "strandA"] == "-":
            if key not in temp_dict:
                # if overlapping feature is downstream of the promoter, change promoter start and stop to the start position of the feature -1
                if (
                    overlapping_proms.loc[i, "startA"]
                    >= overlapping_proms.loc[i, "startB"]
                ):
                    stop = overlapping_proms.loc[i, "startB"] - 1
                    start = overlapping_proms.loc[i, "startB"] - 1
                    temp_dict[key] = [start, stop]
                else:
                    # if overlapping feature is downstream of the promoter, change promtoer start and stop to the start position of the feature -1
                    stop = overlapping_proms.loc[i, "startB"] - 1
                    temp_dict[key] = [overlapping_proms.loc[i, "startA"], stop]

            # if key is in dictionary already
            else:
                # if overlapping feature is downstream of the promoter, change promoter start and stop to the start position of the feature -1
                if temp_dict[key][0] >= overlapping_proms.loc[i, "startB"]:
                    stop = overlapping_proms.loc[i, "startB"] - 1
                    start = overlapping_proms.loc[i, "startB"] - 1
                    temp_dict[key] = [start, stop]

                    # overlapping_proms.drop(i, inplace=True)
                else:
                    # if promoter has already been shortened by another overlapping feature, pass
                    if temp_dict[key][1] < overlapping_proms.loc[i, "startB"]:
                        pass
                    else:
                        stop = overlapping_proms.loc[i, "startB"] - 1
                        temp_dict[key][1] = stop

    # create df with just feature A containing new start and stop locations
    new_feature_A = overlapping_proms[
        [
            "chrA",
            "sourceA",
            "typeA",
            "startA",
            "stopA",
            "dot1A",
            "strandA",
            "dot2A",
            "attributesA",
        ]
    ]

    # rename columns
    cols2 = [
        "chr",
        "source",
        "type",
        "start",
        "stop",
        "dot1",
        "strand",
        "dot2",
        "attributes",
    ]
    new_feature_A.columns = cols2
    # remove duplicates based on attributes
    no_dups = new_feature_A.drop_duplicates("attributes")
    # replace start and stop with those from dictionary
    no_dups_replaced = no_dups.copy()

    for i, data in no_dups.iterrows():
        key = (
            no_dups_replaced.loc[i, "attributes"]
            .split("ID=gene:")[1]
            .split(";")[0]
        )
        no_dups_replaced.loc[i, "start"] = temp_dict[key][0]
        no_dups_replaced.loc[i, "stop"] = temp_dict[key][1]

    # create a buffer of promoter_gff
    promoter_gff_buffer = io.StringIO()
    promoter_gff_buffer.write(str(promoter_gff))
    # go back to beginning of the buffer
    promoter_gff_buffer.seek(0)
    # read in all_proms as df
    all_proms_df = pd.read_table(promoter_gff_buffer, sep="\t", header=None)
    all_proms_df.columns = cols2
    # remove overlapping promoters from all_proms_df
    # first add AGI columns
    # add AGI column to promoters
    all_proms_df_agi = all_proms_df.assign(
        AGI=all_proms_df.attributes.str.extract(r"ID=gene:(.*?)\;")
    )
    new_feature_A_agi = no_dups_replaced.assign(
        AGI=no_dups_replaced.attributes.str.extract(r"ID=gene:(.*?)\;")
    )
    overlapping_proms_agi = overlapping_proms.assign(
        AGI=overlapping_proms.attributesA.str.extract(r"ID=gene:(.*?)\;")
    )
    # remove any promoters in all_prom_df that overlap other genes
    all_prom_removed_overlaps = all_proms_df_agi[
        ~all_proms_df_agi.AGI.isin(overlapping_proms_agi.AGI)
    ]

    # merge all_prom_removed_overlaps with new_feature_A_agi
    # outer merge adding suffix to new_feature_A
    merged = pd.merge(
        all_prom_removed_overlaps, new_feature_A_agi, how="outer"
    )
    # remove unwanted columns
    merged = merged[cols2]
    all_genes_gff.close()
    promoter_gff_buffer.close()
    # sort on chr and start
    merged.sort_values(["chr", "start"], inplace=True, ignore_index=True)

    return merged


def add_5UTR(promoter_gff, all_features_gff):
    """Function to extend the promoters to include the 5'UTR region until the start codon of the first CDS feature of the same gene.
    Also remove mitochondria and chloroplast features."""
    promoters = pd.read_table(
        promoter_gff, sep="\t", header=None, low_memory=False
    )

    cols = [
        "chr",
        "source",
        "type",
        "start",
        "stop",
        "dot1",
        "strand",
        "dot2",
        "attributes",
    ]
    promoters.columns = cols
    # make chr dtype string

    # add AGI column to promoters
    promoters = promoters.assign(
        AGI=promoters.attributes.str.extract(r"ID=gene:(.*?)\;")
    )

    # remove lines beginning with ##
    # create a buffer for feature output
    feature_buffer = io.StringIO()
    with open(all_features_gff, "r") as f:
        for line in f:
            line = line.strip()  # removes hidden characters/spaces
            if line[0] == "#":
                pass
            else:
                # don't include lines beginning with ##
                feature_buffer.write(line + "\n")  # output to buffer
    # go back to beginning of the buffer
    feature_buffer.seek(0)

    # read in feature buffer to df
    features = pd.read_table(
        feature_buffer, sep="\t", header=None, low_memory=False
    )

    features.columns = cols
    # filter features to contain only cds
    cds = features[features.type == "CDS"]
    # filter to only contain cds 1 in parent .1

    # remove all features whose parent isnt Parent=AGI.1
    # cds_1 = cds_1[cds_1.attributes.str.contains('Parent=transcript[A-Z]{2}\d{1}[A-Z]{1}\d{5}.1')]

    # keep only CDS:1;
    # cds_1 = cds_1[cds_1.attributes.str.contains(':CDS:1;')]
    # make chromosome column string
    cds = cds.astype({"chr": "str"})
    cds_1 = cds.copy()
    # remove mitochondria and chloroplast features
    cds_1 = cds[cds["chr"].isin(["1", "2", "3", "4", "5"])]

    # add AGI column to CDSs
    cds2 = cds_1.assign(AGI=cds_1.attributes.str.extract(r"ID=CDS:(.*?)\."))

    # Sort based on chromosome then start
    cds2 = cds2.sort_values(["chr", "start"]).reset_index(drop=True)
    # remove duplicates keeping the CDS closest to the promoter
    # if positive strand, remove the first AGI duplicate
    no_dups_pos = cds2[cds2.strand == "+"].drop_duplicates("AGI", keep="first")
    # if negative strand, remove the last AGI duplicate
    no_dups_neg = cds2[cds2.strand == "-"].drop_duplicates("AGI", keep="last")
    no_dups = pd.merge(no_dups_pos, no_dups_neg, how="outer")
    # sort by chr and start again
    no_dups = no_dups.sort_values(["chr", "start"]).reset_index(drop=True)

    # remove mitochondria and chloroplast features from promoters
    promoters = promoters[promoters["chr"].isin(["1", "2", "3", "4", "5"])]

    # merge on AGI
    merged = pd.merge(
        promoters, no_dups, on="AGI", how="left", suffixes=("", "_cds")
    )

    for i, v in merged.iterrows():
        # if positive strand feature A, reduce length to the stop position of feature B + 1
        if merged.loc[i, "strand"] == "+":
            # if overlapping feature is upstream or equal to start of the promoter, drop row
            if merged.loc[i, "start"] >= merged.loc[i, "start_cds"]:
                merged.drop(i, inplace=True)
            else:
                # extend promoter to include 5'UTR region up until the start of the cds
                merged.loc[i, "stop"] = merged.loc[i, "start_cds"] - 1
        elif merged.loc[i, "strand"] == "-":
            # if overlapping feature is upstream or equal to start of the promoter, drop row
            if merged.loc[i, "stop"] <= merged.loc[i, "stop_cds"]:
                merged.drop(i, inplace=True)
            else:
                # extend promoter to include 5'UTR region up until the start of the cds
                merged.loc[i, "start"] = merged.loc[i, "stop_cds"] + 1

    # remove unwanted columns
    merged_new = merged[
        [
            "chr",
            "source",
            "type",
            "start",
            "stop",
            "dot1",
            "strand",
            "dot2",
            "attributes",
        ]
    ]
    return merged_new


def filter_bad_proms(proms):
    """function to remove promoters whose start is higher or equal to the stop.
    Then calculate lengths and sort by length. Drop duplicates keeping the shortest promoter."""
    #     #Remove promoters with start higher than stop
    #     for i, v in proms.iterrows():
    #         if proms.loc[i,'start'] >=  proms.loc[i, 'stop']:
    #                 proms.drop(i, inplace=True)
    # sort based on length then keep the shortest, and remove lengths that are 0
    # Make promoter length column
    proms_length = proms.assign(length=(proms.stop - proms.start))
    # remove promoters with length <100
    removed = proms_length.loc[proms_length.length >= 100]
    # remove duplicates keeping the shortest promoter
    no_dups = removed.sort_values("length", ascending=True).drop_duplicates(
        "attributes", keep="first"
    )
    proms_dropped = no_dups[
        [
            "chr",
            "source",
            "type",
            "start",
            "stop",
            "dot1",
            "strand",
            "dot2",
            "attributes",
        ]
    ]
    # sort by chr then start
    proms_dropped = proms_dropped.sort_values(["chr", "start"]).reset_index(
        drop=True
    )
    return proms_dropped


def remove_characters_linestart(
    input_location, output_location, oldcharacters, newcharacters, linestart
):
    """this function removes characters from the start of each line in the input file and sends modified lines to output"""
    output = open(
        output_location, "w"
    )  # make output file with write capability
    # open input file
    with open(input_location, "r") as infile:
        # iterate over lines in file
        for line in infile:
            line = line.strip()  # removes hidden characters/spaces
            if line[0] == linestart:

                line = line.replace(
                    oldcharacters, newcharacters
                )  # remove characters from start of line, replace with new characters
            output.write(line + "\n")  # output to new file
    output.close()


def count_promoters(in_file, out_file):
    """this function creates a text file detailing the number of promoters in an input GFF file """
    examiner = GFFExaminer()
    # open input GFF file
    in_handle = open(in_file, "r")
    # output a text file, giving information such as no. of promoters in the file
    with open(out_file, "w") as fout:
        fout.write(pformat(examiner.available_limits(in_handle)))
    in_handle.close()


def remove_empty_lines(input_file):
    """function to strip empty line from end of file"""
    # make temporary buffer file
    temp_buffer = io.StringIO()
    with open(input_file) as f_input:
        data = f_input.read().rstrip("\n")
    temp_buffer.write(data)
    temp_buffer.seek(0)
    data2 = temp_buffer.read().rstrip("\n")

    with open(input_file, "w") as f_output:
        f_output.write(data2)


def bidirectional_proms(in_file, out_file):
    """this function creates a file containing all promoters with an upstream gene going in the other direction ie. potential bidirectional promoters"""
    # read in gff file
    promoters = pd.read_table(in_file, sep="\t", header=2)
    cols2 = [
        "chr",
        "source",
        "type",
        "start",
        "stop",
        "dot1",
        "strand",
        "dot2",
        "attributes",
    ]
    promoters.columns = cols2
    # make sure lines are sorted

    promoters = promoters.sort_values(["chr", "start"]).reset_index(drop=True)
    # if bidirectional
    promoters["bidirectional"] = "no"

    for i, data in promoters.iterrows():
        if i - 1 >= 0:
            # print(i,data)
            if (
                promoters.loc[i, "strand"] == "+"
                and promoters.loc[i - 1, "strand"] == "-"
                and promoters.loc[i, "start"] - promoters.loc[i - 1, "stop"]
                < 2000
            ):
                promoters.loc[i, "bidirectional"] = "yes"
                promoters.loc[i - 1, "bidirectional"] = "yes"

    # write to file
    promoters[promoters.bidirectional == "no"][
        [
            "chr",
            "source",
            "type",
            "start",
            "stop",
            "dot1",
            "strand",
            "dot2",
            "attributes",
        ]
    ].to_csv(out_file, index=False, sep="\t", header=None)


def main(args):
    # parse arguments
    args = parse_args(args)

    # make directory for the output files to be exported to
    dirName = f"../../data/output/{args.file_names}"
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " created")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

    # TSS_raw = f'../../data/TSS_data/AnnotatedPEATPeaks.txt'
    # TSS_renamedChr = '../..//data/TSS_data/AnnotatedPEATPeaks_renamedChr.gff'
    # genome2 = "../..//data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel_sl_simpleChr.fasta"
    genome = "../../data/genomes/TAIR10_chr_all.fas"

    # genes = f"../../data/genomes/Araport11_GFF3_genes_transposons.201606.gff"
    genes = "../../data/genomes/Arabidopsis_thaliana/annotation/Arabidopsis_thaliana.TAIR10.47.gff3"
    # genes_renamedChr = "../..//data/genomes/Araport11_GFF3_genes_transposons.201606_renamedChr.gff"
    # test_genes = f"../../data/genomes/test_genes.gff3"
    # testgenesonly_gff = f"../../data/genomes/testgenes_only.gff3"
    # temp = f"../../data/TSS_data/temp.txt"
    # TSS = "../..//data/TSS_data/AnnotatedPEATPeaks_renamedcol.gff"
    # TSS = "../..//data/TSS_data/TSStest.txt"
    # find_closest_TSS(genes,output,temp)
    genesoverlapping_gff = (
        f"../../data/output/{args.file_names}/genesoverlapping.gff3"
    )
    genesonly_gff = f"../../data/output/{args.file_names}/genesonly.gff3"
    promoters = f"../../data/output/{args.file_names}/promoters.gff3"
    promoters_5UTR = f"../../data/output/{args.file_names}/promoters_5UTR.gff3"
    # genes_bed = "../../data/genomes/genes.bed"
    # overlapping_promoters = f'../../data/output/{args.file_names}/promoters_overlapping.gff3'
    promoterandgenes_only_overlap = f"../../data/output/{args.file_names}/promoterandgenes_only_overlap.gff3"

    chromsizes_file = f"../../data/output/{args.file_names}/chromsizes.chr"
    # need to change temporary files to scratch directory
    chromsizes_file_renamedChr_temp = (
        f"../../data/output/{args.file_names}/chromsizes_renamedChr_temp.chr"
    )
    chromsizes_file_renamedChr = (
        f"../../data/output/{args.file_names}/chromsizes_renamedChr.chr"
    )
    # chromsizes_file2 = '../..//data/genomes/chromsizes2.chr'
    # promoters_renamedChr_temp = f'../../data/output/{args.file_names}/promoters_renamedChr_temp.gff3'
    # promoters_renamedChr_temp2 = f'../../data/genomes/{args.file_names}/promoters_renamedChr_temp2.gff3'
    # promoters_renamedChr = f'../../data/output/{args.file_names}/promoters_renamedChr.gff3'
    # promoters_5UTR_renamedChr = f'../../data/output/{args.file_names}/promoters_5UTR_renamedChr.gff3'
    nonbidirectional_promoters = (
        f"../../data/output/{args.file_names}/nonbidirectional_proms.gff3"
    )

    fasta_chromsizes(genome, chromsizes_file)

    # rename mitochondria and chloroplast to M and C
    remove_characters_linestart(
        chromsizes_file,
        chromsizes_file_renamedChr_temp,
        "mitochondria",
        "Mt",
        "m",
    )
    remove_characters_linestart(
        chromsizes_file_renamedChr_temp,
        chromsizes_file_renamedChr,
        "chloroplast",
        "Pt",
        "c",
    )
    os.remove(chromsizes_file_renamedChr_temp)

    # extract_genes(genes,genesonly_gff)
    # note - this changes chromosome no. to 1 rather than Chr1?? NOT SURE IT STILL DOES
    protein_coding = extract_genes(genes, genesoverlapping_gff, genesonly_gff)

    # createfile containing all nonbidirectional genes (bidirectional = genes with an upstream gene in the other direction ie. potential overlapping promoters)
    if args.remove_bidirectional:
        bidirectional_proms(genesonly_gff, nonbidirectional_promoters)
        selected_genes = nonbidirectional_promoters
    else:
        selected_genes = genesonly_gff

    # add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff
    # add_promoter(selected_genes,chromsizes_file_renamedChr,1000)

    # create file containing only promoters which overlap other genome features
    # promoter_overlap(promoters,genes,overlapping_promoters)

    # count no. of promoters in overlapping promoters file
    # count_promoters(overlapping_promoters, f'../../data/genomes/overlapping_promoters.txt')

    # all promoters
    # in_file = promoters
    # examiner = GFFExaminer()
    # in_handle = open(#in_file)
    # pprint(examiner.available_limits(in_handle))
    # in_handle.close()

    # promoters overlapping only genes
    if args.prevent_overlapping_genes:
        # add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff
        # note, this only removes overlap for protein coding genes
        promoters_incl_overlap = add_promoter(
            selected_genes, chromsizes_file_renamedChr, args.promoter_length
        )
        subtracted = remove_promoter_overlap(
            promoters_incl_overlap,
            protein_coding,
            promoterandgenes_only_overlap,
        )
        with open(promoters, "w") as f:
            subtracted.to_csv(f, index=False, sep="\t", header=None)

    else:
        # add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff
        promoters_gff = add_promoter(
            selected_genes, chromsizes_file_renamedChr, args.promoter_length
        )
        with open(promoters, "w") as f:
            f.write(str(promoters_gff))

    # add 5'utr?
    if args.fiveUTR:
        # extend promoter up until first CDS
        prom_UTR = add_5UTR(promoters, genes)
        with open(promoters_5UTR, "w") as f:
            prom_UTR.to_csv(f, index=False, sep="\t", header=None)
        # sort proms and proms+5'UTR, and remove AGI column
        proms_UTR = pd.read_table(promoters_5UTR, sep="\t", header=None)
        cols = [
            "chr",
            "source",
            "type",
            "start",
            "stop",
            "dot1",
            "strand",
            "dot2",
            "attributes",
        ]
        proms_UTR.columns = cols
        proms_UTR = proms_UTR.sort_values(["chr", "start"]).reset_index(
            drop=True
        )
        cleaned_proms_UTR = filter_bad_proms(proms_UTR)
        cleaned_proms_UTR.type = "promoter"
        cleaned_proms_UTR.to_csv(
            promoters_5UTR, index=False, sep="\t", header=None
        )
    else:
        pass
    # sort proms and proms+5'UTR, and remove AGI column
    proms = pd.read_table(promoters, sep="\t", header=None)

    cols = [
        "chr",
        "source",
        "type",
        "start",
        "stop",
        "dot1",
        "strand",
        "dot2",
        "attributes",
    ]
    # cols2 = ['chr', 'source', 'type', 'start','stop','dot1','strand','dot2','attributes']
    proms.columns = cols
    proms = proms.sort_values(["chr", "start"]).reset_index(drop=True)
    proms = proms[
        [
            "chr",
            "source",
            "type",
            "start",
            "stop",
            "dot1",
            "strand",
            "dot2",
            "attributes",
        ]
    ]

    cleaned_proms = filter_bad_proms(proms)

    # rename gene to promoter
    cleaned_proms.type = "promoter"

    cleaned_proms.to_csv(promoters, index=False, sep="\t", header=None)

    # #remove empty lines from end of the files. DO NOT DO THIS AS GFF2BED REMOVES THE LAST LINE - BED FILES SEEM TO NEED AN EMPTY LINE
    # remove_empty_lines(promoters)
    # remove_empty_lines(promoters_5UTR)

    # count no. of promoters in overlapping promoters file
    # count_promoters(promoterandgenes_only_overlap, f'../../data/genomes/promoterandgenes_only_overlap.txt')

    # examiner = GFFExaminer()
    # in_handle = open(in_file)
    # pprint.pprint(examiner.available_limits(in_handle))
    # in_handle.close()

    # remove the Chr from the chromosome names in promoters.gff3. #Replace M with mitochondria and C with chloroplast

    # remove_characters_linestart(promoters, promoters_renamedChr, 'Chr', '','C')
    # remove_characters_linestart(promoters_5UTR, promoters_5UTR_renamedChr, 'Chr', '','C')

    # remove_characters_linestart(promoters_renamedChr_temp, promoters_renamedChr_temp2, 'M', 'mitochondria','M')
    # remove_characters_linestart(promoters_renamedChr_temp2, promoters_renamedChr, 'C', 'chloroplast','C')
    # remove empty lines from end of the files
    # remove_empty_lines(promoters_renamedChr)
    # remove_empty_lines(promoters_5UTR_renamedChr)

    # os.remove(promoters_renamedChr_temp)
    # os.remove(promoters_renamedChr_temp2)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
