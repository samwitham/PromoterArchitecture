#create a bed file with the promoters upstream of the eukaryotic promoter database TSS
import pandas as pd
import argparse 
import io
from pybedtools import BedTool

# #define arguments
parser = argparse.ArgumentParser(description='create_EPD_promoters')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('promoter_5UTR_bedfile', type=str, help='Input location of promoter-5UTR bedfile')
parser.add_argument('EPD_TSS_bed', type=str, help='Input location of Eukaryotic Promoter Database TSS bed file')
parser.add_argument('EPD_promoters', type=str, help='Output location of the EPD_promoters bedfile')
parser.add_argument('flagged_proms_not_in_EPD', type=str, help='Output location of flagged promoters which are not in the Eukaryotic promoter database')
parser.add_argument('flagged_EPD_TSS_in_CDS', type=str, help='Output location of flagged EPD_promoters which have TSSs in coding regions')
parser.add_argument('flagged_EPD_overlappingprom_so_only5UTRs', type=str, help='Output location of flagged EPD_promoters which are overlapping other genes so they are only a shortened 5\'UTR')
args = parser.parse_args()

def create_EPD_proms(promoter_5UTR_bedfile,EPD_TSS_bed,EPD_promoters,flagged_proms_not_in_EPD,flagged_EPD_TSS_in_CDS,flagged_EPD_overlappingprom_so_only5UTRs):
    """function to split the promoter_5UTR bed file at the EPD_TSS and save an EPD_promoter file"""
    #read in files:
    #read in promoters with 5'UTRs
    promoter_5UTRs = pd.read_table(promoter_5UTR_bedfile,sep='\t',header=None)
    cols = ['chr','start','stop','AGI','dot','strand','source','type','dot2','attributes']
    promoter_5UTRs.columns = cols
    #read in EPD_TSS_bed
    EPD_TSS_df = pd.read_table(EPD_TSS_bed, delim_whitespace=True, header=None, skiprows=4)
    cols = ['chr','start','stop','transcript_EPD','score_EPD','strand_EPD','thickstart_EPD','thickend_EPD']
    EPD_TSS_df.columns = cols
    #add AGI column
    EPD_TSS_df['AGI'] = EPD_TSS_df.transcript_EPD.str.split('_',expand=True)[0]
    #add TSS location column
    EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '+', 'TSS_EPD'] = EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '+', 'thickstart_EPD']
    EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '-', 'TSS_EPD'] = EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '-', 'thickend_EPD'] -1
    #merged with promoter_5UTRs
    merged = pd.merge(promoter_5UTRs,EPD_TSS_df, on='AGI', how='left', suffixes=('','_EPD'))
    #remove NaN (promoters in EPD but not in promoters_5UTR)
    merged = merged[merged.source.notnull()]
    #flag genes which are in promoters_5UTR but are not in EPD
    flagged = merged[~merged.TSS_EPD.notnull()]
    #save flagged file genes which are not in EPD database
    BedTool.from_dataframe(flagged).saveas(flagged_proms_not_in_EPD)
    
    #remove NaN (promoters_5UTR but not in EPD)
    merged = merged[merged.TSS_EPD.notnull()]
    #make a copy of the df
    merged_copy = merged.copy()
    #make new columns with start and stop for EPD_promoter
    
    #iterate over rows
    for i,data in merged_copy.iterrows():
        #if positive strand gene, make the start of the EPD promoter equal to the the start position of the promoter_5UTR, 
        #and the EPDprom stop the same as the EPD_TSS stop
        if merged_copy.loc[i,'strand'] == '+':            
            merged.loc[i,'start_EPDprom'] = merged_copy.loc[i,'start']
            merged.loc[i,'stop_EPDprom'] = merged_copy.loc[i,'TSS_EPD']
        #if negative strand gene, make the start of the EPD promoter the same as EPD_TSS, and make the EPD promoter stop equal to the promoters stop
        if merged_copy.loc[i,'strand'] == '-':
            merged.loc[i,'start_EPDprom'] = merged_copy.loc[i,'TSS_EPD']
            merged.loc[i,'stop_EPDprom'] = merged_copy.loc[i,'stop']
    #make start and stop integars
    merged = merged.astype({'start_EPDprom': 'int','stop_EPDprom':'int'})
    #drop promoters whose EPD TSS falls within a coding region
    #make positive and negative strand masks
    pos_mask = merged[merged.strand=='+']
    neg_mask = merged[merged.strand=='-']
    #make dfs for pos and neg strands where the EPD TSS is in a cokding region
    pos_in_CDS = pos_mask[pos_mask.stop < pos_mask.stop_EPDprom]
    neg_in_CDS = neg_mask[neg_mask.start > neg_mask.start_EPDprom]
    #concatenate those dfs
    hascodingTSS = pd.concat([pos_in_CDS,neg_in_CDS],ignore_index=True)
    #filter columns
    hascodingTSS = hascodingTSS[['chr','start_EPDprom','stop_EPDprom','AGI','dot','strand','source','type','dot2','attributes']]
    #sort on chr and start
    hascodingTSS.sort_values(['chr','start_EPDprom'], inplace=True, ignore_index=True)
    #save file with no UTRs
    BedTool.from_dataframe(hascodingTSS).saveas(flagged_EPD_TSS_in_CDS)   
    #filter EPD_promoters which have their TSSs in coding regions
    merged = merged[~merged.AGI.isin(hascodingTSS.AGI)]
    
    #flag genes where the promoter is overlapping another gene so they consist of only a shortened 5'UTR
    only5UTR = merged[merged.start_EPDprom >= merged.stop_EPDprom]
    #filter columns
    only5UTR = only5UTR[['chr','start_EPDprom','stop_EPDprom','AGI','dot','strand','source','type','dot2','attributes']]
    #sort on chr and start
    only5UTR.sort_values(['chr','start_EPDprom'], inplace=True, ignore_index=True)
    #save file with no UTRs
    BedTool.from_dataframe(only5UTR).saveas(flagged_EPD_overlappingprom_so_only5UTRs)   
    
    
    #filter genes where the promoter is overlapping another gene so they consist of only a shortened 5'UTR
    merged = merged[~(merged.start_EPDprom >= merged.stop_EPDprom)]

    #filter columns
    merged = merged[['chr','start_EPDprom','stop_EPDprom','AGI','dot','strand','source','type','dot2','attributes']]
 
    #sort on chr and start
    merged.sort_values(['chr','start_EPDprom'], inplace=True, ignore_index=True)

    #save file with no UTRs
    BedTool.from_dataframe(merged).saveas(EPD_promoters)
    
#create the EPD_promoters and various flagged gene files
create_EPD_proms(args.promoter_5UTR_bedfile,args.EPD_TSS_bed,args.EPD_promoters,args.flagged_proms_not_in_EPD,args.flagged_EPD_TSS_in_CDS,args.flagged_EPD_overlappingprom_so_only5UTRs)