# create a bed file with the promoters upstream of the eukaryotic promoter database TSS
import argparse

import pandas as pd
from pybedtools import BedTool


def parse_args(args):
    """define arguments"""
    parser = argparse.ArgumentParser(description="create_EPD_promoters")
    parser.add_argument(
        "file_names",
        type=str,
        help="Name of folder and filenames for the promoters extracted",
    )
    parser.add_argument(
        "promoter_5UTR_bedfile",
        type=str,
        help="Input location of promoter-5UTR bedfile",
    )
    parser.add_argument(
        "EPD_TSS_bed",
        type=str,
        help="Input location of Eukaryotic Promoter Database TSS bed file",
    )
    parser.add_argument(
        "EPD_promoters",
        type=str,
        help="Output location of the EPD_promoters bedfile",
    )
    parser.add_argument(
        "EPD_promoters_5UTR",
        type=str,
        help="Output location of the EPD_promoters5UTR bedfile",
    )
    parser.add_argument(
        "flagged_proms_not_in_EPD",
        type=str,
        help="Output location of flagged promoters which are not in the Eukaryotic promoter database",
    )
    parser.add_argument(
        "flagged_EPD_TSS_in_CDS",
        type=str,
        help="Output location of flagged EPD_promoters which have TSSs in coding regions",
    )
    parser.add_argument(
        "flagged_EPD_overlappingprom_so_only5UTRs",
        type=str,
        help="Output location of flagged EPD_promoters which are overlapping other genes so they are only a shortened 5'UTR",
    )
    parser.add_argument(
        "promoter_length",
        type=int,
        help="Length of promoter to use upstream of EPD TSS in second promoter-5UTR bedfile",
    )
    parser.add_argument(
        "fiveUTR_length",
        type=int,
        help="Length of 5UTR to use downstream of EPD TSS in second promoter-5UTR bedfile",
    )
    parser.add_argument(
        "flagged_EPD_overlapping_specified_length",
        type=str,
        help="Output location of flagged EPD_promoters_5UTR which extend into coding regions",
    )

    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def create_EPD_proms(
    promoter_5UTR_bedfile,
    EPD_TSS_bed,
    EPD_promoters,
    flagged_proms_not_in_EPD,
    flagged_EPD_TSS_in_CDS,
    flagged_EPD_overlappingprom_so_only5UTRs,
    EPD_promoters_5UTR,
    promoter_length,
    fiveUTR_length,
    flagged_EPD_overlapping_specified_length,
):
    """function to split the promoter_5UTR bed file at the EPD_TSS and save an EPD_promoter file"""
    # read in files:
    # read in promoters with 5'UTRs
    promoter_5UTRs = pd.read_table(
        promoter_5UTR_bedfile, sep="\t", header=None
    )
    cols = [
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
    promoter_5UTRs.columns = cols
    # read in EPD_TSS_bed
    EPD_TSS_df = pd.read_table(
        EPD_TSS_bed, delim_whitespace=True, header=None, skiprows=4
    )
    cols = [
        "chr",
        "start",
        "stop",
        "transcript_EPD",
        "score_EPD",
        "strand_EPD",
        "thickstart_EPD",
        "thickend_EPD",
    ]
    EPD_TSS_df.columns = cols
    # add AGI column
    EPD_TSS_df["AGI"] = EPD_TSS_df.transcript_EPD.str.split("_", expand=True)[
        0
    ]
    # add TSS location column
    EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == "+", "TSS_EPD"] = EPD_TSS_df.loc[
        EPD_TSS_df.strand_EPD == "+", "thickstart_EPD"
    ]
    EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == "-", "TSS_EPD"] = (
        EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == "-", "thickend_EPD"] - 1
    )
    # merged with promoter_5UTRs
    merged = pd.merge(
        promoter_5UTRs, EPD_TSS_df, on="AGI", how="left", suffixes=("", "_EPD")
    )
    # remove NaN (promoters in EPD but not in promoters_5UTR)
    merged = merged[merged.source.notnull()]
    # flag genes which are in promoters_5UTR but are not in EPD
    flagged = merged[~merged.TSS_EPD.notnull()]
    # save flagged file genes which are not in EPD database
    BedTool.from_dataframe(flagged).saveas(flagged_proms_not_in_EPD)

    # remove NaN (promoters_5UTR but not in EPD)
    merged = merged[merged.TSS_EPD.notnull()]
    # make a copy of the df
    merged_copy = merged.copy()
    # make new columns with start and stop for EPD_promoter

    # iterate over rows
    for i, data in merged_copy.iterrows():
        # if positive strand gene, make the start of the EPD promoter equal to the the start position of the promoter_5UTR,
        # and the EPDprom stop the same as the EPD_TSS stop
        if merged_copy.loc[i, "strand"] == "+":
            merged.loc[i, "start_EPDprom"] = merged_copy.loc[i, "start"]
            merged.loc[i, "stop_EPDprom"] = merged_copy.loc[i, "TSS_EPD"]
        # if negative strand gene, make the start of the EPD promoter the same as EPD_TSS, and make the EPD promoter stop equal to the promoters stop
        if merged_copy.loc[i, "strand"] == "-":
            merged.loc[i, "start_EPDprom"] = merged_copy.loc[i, "TSS_EPD"]
            merged.loc[i, "stop_EPDprom"] = merged_copy.loc[i, "stop"]
    # make start and stop integars
    merged = merged.astype({"start_EPDprom": "int", "stop_EPDprom": "int"})
    # drop promoters whose EPD TSS falls within a coding region
    # make positive and negative strand masks
    pos_mask = merged[merged.strand == "+"]
    neg_mask = merged[merged.strand == "-"]
    # make dfs for pos and neg strands where the EPD TSS is in a coding region
    pos_in_CDS = pos_mask[pos_mask.stop < pos_mask.stop_EPDprom]
    neg_in_CDS = neg_mask[neg_mask.start > neg_mask.start_EPDprom]
    # concatenate those dfs
    hascodingTSS = pd.concat([pos_in_CDS, neg_in_CDS], ignore_index=True)
    # filter columns
    hascodingTSS = hascodingTSS[
        [
            "chr",
            "start_EPDprom",
            "stop_EPDprom",
            "AGI",
            "dot",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
        ]
    ]
    # sort on chr and start
    hascodingTSS.sort_values(
        ["chr", "start_EPDprom"], inplace=True, ignore_index=True
    )
    # save file with no UTRs
    BedTool.from_dataframe(hascodingTSS).saveas(flagged_EPD_TSS_in_CDS)
    # filter EPD_promoters which have their TSSs in coding regions
    merged = merged[~merged.AGI.isin(hascodingTSS.AGI)]

    # flag genes where the promoter is overlapping another gene so they consist of only a shortened 5'UTR
    only5UTR = merged[merged.start_EPDprom >= merged.stop_EPDprom]
    # filter columns
    only5UTR = only5UTR[
        [
            "chr",
            "start_EPDprom",
            "stop_EPDprom",
            "AGI",
            "dot",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
        ]
    ]
    # sort on chr and start
    only5UTR.sort_values(
        ["chr", "start_EPDprom"], inplace=True, ignore_index=True
    )
    # save file with no UTRs
    BedTool.from_dataframe(only5UTR).saveas(
        flagged_EPD_overlappingprom_so_only5UTRs
    )

    # filter genes where the promoter is overlapping another gene so they consist of only a shortened 5'UTR
    merged = merged[~(merged.start_EPDprom >= merged.stop_EPDprom)]

    # make a copy of the df for later
    merged_orig = merged.copy()

    # filter columns
    merged = merged[
        [
            "chr",
            "start_EPDprom",
            "stop_EPDprom",
            "AGI",
            "dot",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
        ]
    ]

    # sort on chr and start
    merged.sort_values(
        ["chr", "start_EPDprom"], inplace=True, ignore_index=True
    )

    # save file with no UTRs
    BedTool.from_dataframe(merged).saveas(EPD_promoters)

    # Create n-length promoter & m-length 5UTR file using EPD TSS
    # make copy of merged
    EPD_proms = merged.copy()

    # iterate over rows
    for i, data in EPD_proms.iterrows():
        # if positive strand gene, make the start of the EPD promoter equal to the the stop position minus promoter_length,
        # and the EPDprom stop extend to the 5UTR length
        if EPD_proms.loc[i, "strand"] == "+":
            merged.loc[i, "newstart_EPDprom"] = (
                EPD_proms.loc[i, "stop_EPDprom"] - promoter_length
            )
            merged.loc[i, "newstop_EPDprom"] = (
                EPD_proms.loc[i, "stop_EPDprom"] + fiveUTR_length
            )
        # if negative strand gene, make the start of the EPD promoter extend to the 5UTR length, and make the EPD promoter stop equal to the EPD TSS + the promoter length
        if EPD_proms.loc[i, "strand"] == "-":
            merged.loc[i, "newstart_EPDprom"] = (
                EPD_proms.loc[i, "start_EPDprom"] - fiveUTR_length
            )
            merged.loc[i, "newstop_EPDprom"] = (
                EPD_proms.loc[i, "start_EPDprom"] + promoter_length
            )
    # make start and stop integars
    merged = merged.astype(
        {"newstart_EPDprom": "int", "newstop_EPDprom": "int"}
    )

    # sort on chr and start
    merged.sort_values(
        ["chr", "newstart_EPDprom"], inplace=True, ignore_index=True
    )

    # select correct columns
    EPD_proms_shortened = merged[
        [
            "chr",
            "newstart_EPDprom",
            "newstop_EPDprom",
            "AGI",
            "dot",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
        ]
    ].copy()

    # flag and remove genes where the promoter start or 5UTR end are beyond the original Araport11 promoters_5'UTR file
    merged2 = pd.merge(
        EPD_proms_shortened,
        promoter_5UTRs,
        on="AGI",
        how="left",
        suffixes=("", "_ara"),
    )

    # make positive and negative strand masks
    pos_mask = merged2[merged2.strand == "+"]
    neg_mask = merged2[merged2.strand == "-"]
    # make dfs for pos and neg strands where the EPD TSS is in a coding region
    pos_in_CDS = pos_mask[
        (pos_mask.stop < pos_mask.newstop_EPDprom)
        | (pos_mask.start > pos_mask.newstart_EPDprom)
    ]
    neg_in_CDS = neg_mask[
        (neg_mask.start > neg_mask.newstart_EPDprom)
        | (neg_mask.stop < neg_mask.newstop_EPDprom)
    ]
    # concatenate those dfs
    extends_coding = pd.concat([pos_in_CDS, neg_in_CDS], ignore_index=True)
    # filter columns
    extends_coding = extends_coding[
        [
            "chr",
            "newstart_EPDprom",
            "newstop_EPDprom",
            "AGI",
            "dot",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
        ]
    ]
    # sort on chr and start
    extends_coding.sort_values(
        ["chr", "newstart_EPDprom"], inplace=True, ignore_index=True
    )
    # save file with flagged genes then overlap a coding region
    BedTool.from_dataframe(extends_coding).saveas(
        flagged_EPD_overlapping_specified_length
    )

    # shorten the EPD promoter/5UTRs if they are flagged
    merged_again = pd.merge(
        EPD_proms_shortened,
        merged_orig,
        on="AGI",
        how="left",
        suffixes=("", "_ara"),
    )
    # set newnewstop_EPDprom and newnewstart_EPDprom to newstart_EPDprom and newstop_EPDprom values
    merged_again["newnewstart_EPDprom"] = merged_again["newstart_EPDprom"]
    merged_again["newnewstop_EPDprom"] = merged_again["newstop_EPDprom"]
    EPD_proms_again = merged_again.copy()

    # iterate over rows
    for i, data in EPD_proms_again.iterrows():
        # if the start location of the araport promoter is greater than the start location of the EPD promoter, shorten the EPD promoter to the length of the araport promoter so it doesnt overlap a coding region
        if (
            EPD_proms_again.loc[i, "start"]
            > EPD_proms_again.loc[i, "newstart_EPDprom"]
        ):
            merged_again.loc[i, "newnewstart_EPDprom"] = EPD_proms_again.loc[
                i, "start"
            ]
        # if the stop location of the araport promoter is less than the stop location of the EPD promoter, shorten the EPD 5UTR to the length of the araport promoter so it doesnt overlap a coding region
        elif (
            EPD_proms_again.loc[i, "stop"]
            < EPD_proms_again.loc[i, "newstop_EPDprom"]
        ):
            merged_again.loc[i, "newnewstop_EPDprom"] = EPD_proms_again.loc[
                i, "stop"
            ]
        # if the EPD promoter/5UTR does not overlap a coding region then keep the start/stop location the same
    #         elif EPD_proms_again.loc[i, 'start'] > EPD_proms_again.loc[i, 'newstart_EPDprom']:
    #             merged_again.loc[i,'newnewstart_EPDprom'] = EPD_proms_again.loc[i, 'newstart_EPDprom']
    #         elif EPD_proms_again.loc[i, 'stop'] > EPD_proms_again.loc[i, 'newstop_EPDprom']:
    #             merged_again.loc[i,'newnewstop_EPDprom'] = EPD_proms_again.loc[i, 'newstop_EPDprom']

    # make start and stop integars
    merged_again = merged_again.astype(
        {"newnewstart_EPDprom": "int", "newnewstop_EPDprom": "int"}
    )
    # sort on chr and start
    merged_again.sort_values(
        ["chr", "newnewstart_EPDprom"], inplace=True, ignore_index=True
    )
    # find lengths of promoters/5UTRs
    merged_again["length"] = (
        merged_again.newnewstop_EPDprom - merged_again.newnewstart_EPDprom
    )
    print(merged_again.length.min())
    print(merged_again.length.max())

    # merged_again.to_csv(EPD_promoters_5UTR, sep='\t',index=False, header=True)
    # select columns
    EPD_proms_shortened = merged_again[
        [
            "chr",
            "newnewstart_EPDprom",
            "newnewstop_EPDprom",
            "AGI",
            "dot",
            "strand",
            "source",
            "type",
            "dot2",
            "attributes",
        ]
    ]

    # filter EPD_promoters which have their TSSs in coding regions
    # EPD_proms_shortened = EPD_proms_shortened[~EPD_proms_shortened.AGI.isin(extends_coding.AGI)]

    # save file with both specificed promoter length and 5UTR length
    BedTool.from_dataframe(EPD_proms_shortened).saveas(EPD_promoters_5UTR)


def main(args):
    # parse arguments
    args = parse_args(args)

    # create the EPD_promoters and various flagged gene files
    create_EPD_proms(
        args.promoter_5UTR_bedfile,
        args.EPD_TSS_bed,
        args.EPD_promoters,
        args.flagged_proms_not_in_EPD,
        args.flagged_EPD_TSS_in_CDS,
        args.flagged_EPD_overlappingprom_so_only5UTRs,
        args.EPD_promoters_5UTR,
        args.promoter_length,
        args.fiveUTR_length,
        args.flagged_EPD_overlapping_specified_length,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
