import pandas as pd
import argparse
import io
from pybedtools import BedTool
#define arguments
parser = argparse.ArgumentParser(description='create_5UTRs')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('promoter_bedfile', type=str, help='Input location of promoter bedfile')
parser.add_argument('promoter_5UTR_bedfile', type=str, help='Input location of promoter-5UTR bedfile')
parser.add_argument('five_UTR', type=str, help='Output location of the 5UTR bedfile')
parser.add_argument('five_UTR_swapped_strands', type=str, help='Output location of the 5UTR bedfile with artificially switched strands for sliding window analysis')
parser.add_argument('genes_no_5UTR', type=str, help='Output location of genes which have no 5UTR')
args = parser.parse_args()

def make_5UTRs(promoter_5UTR_bedfile, promoter_bedfile, five_UTR, five_UTR_swapped_strands, genes_no_5UTR):
    """function to create 5'UTR bed files (one normal, one with swapped strands for sliding window analysis) using input promoters.bed and promtoers_5UTR.bed"""
    #read in promoters with 5'UTRs
    promoter_5UTRs = pd.read_table(promoter_5UTR_bedfile,sep='\t',header=None)
    cols = ['chr','start','stop','promoter_AGI','dot','strand','source','type','dot2','attributes']
    promoter_5UTRs.columns = cols
    #read in promoters not 5'UTRs
    promoters = pd.read_table(promoter_bedfile,sep='\t',header=None)
    promoters.columns = cols
    #merged the dfs
    merged = pd.merge(promoter_5UTRs,promoters, on='promoter_AGI',how='left',suffixes=('','_promoters'))
    #remove NaN values (where only part of the 5UTR is in promoters_5UTR and the promoter is missing as it overlaps another gene)
    merged = merged[merged.start_promoters.notnull()]
    #create df copy
    merged_copy = merged.copy()
    #iterate over rows
    for i,data in merged_copy.iterrows():
        #if positive strand gene, make the start of the 5UTR equal to the the stop position of the promoter, and the 5UTR stop the same as the promoters_5UTR stop
        if merged_copy.loc[i,'strand'] == '+':            
            merged.loc[i,'start_5UTR'] = merged_copy.loc[i,'stop_promoters']
            merged.loc[i,'stop_5UTR'] = merged_copy.loc[i,'stop']
        #if negative strand gene, make the start of the 5'UTR the same as the start of the promoters_5UTR, and make the 5UTR stop equal to the promoters start
        if merged_copy.loc[i,'strand'] == '-':
            merged.loc[i,'start_5UTR'] = merged_copy.loc[i,'start']
            merged.loc[i,'stop_5UTR'] = merged_copy.loc[i,'start_promoters']
    #make start and stop integars
    merged = merged.astype({'start_promoters': 'int','stop_promoters':'int'})
 
    
    no_5UTR = merged[(merged.start == merged.start_promoters) & (merged.stop == merged.stop_promoters)]
    #filter columns
    no_5UTR = no_5UTR[['chr','start','stop','promoter_AGI','dot','strand','source','type','dot2','attributes']]
    #sort on chr and start
    no_5UTR.sort_values(['chr','start'], inplace=True, ignore_index=True)
    #save file with no UTRs
    BedTool.from_dataframe(no_5UTR).saveas(genes_no_5UTR)    
    #filter genes which don't have 5'UTRs (the promoter start and stop are identical to the promoter_5UTR start and stop)
    merged = merged[~((merged.start == merged.start_promoters) & (merged.stop == merged.stop_promoters))]
      
    #filter columns
    merged = merged[['chr','start_5UTR','stop_5UTR','promoter_AGI','dot','strand','source','type','dot2','attributes']]
    #change type to 5UTR
    merged['type'] = '5UTR'
    #sort on chr and start
    merged.sort_values(['chr','start_5UTR'], inplace=True, ignore_index=True)
    #make start_5UTR and stop_5UTR integars
    merged = merged.astype({'start_5UTR': 'int','stop_5UTR':'int'})

    #save as bed file
    BedTool.from_dataframe(merged).saveas(five_UTR)
    #make a df copy    
    merged_copy = merged.copy()
    #swap the strands by iterating over rows
    for i,data in merged_copy.iterrows():
        #if strand is positive, make it negative
        if merged_copy.loc[i,'strand'] == '+':
            merged.loc[i,'strand'] = '-'
        #if strand is negative, make it positive
        if merged_copy.loc[i,'strand'] == '-':
            merged.loc[i,'strand'] = '+'
    
    #save as bed file
    BedTool.from_dataframe(merged).saveas(five_UTR_swapped_strands)
    
    
make_5UTRs(args.promoter_5UTR_bedfile, args.promoter_bedfile, args.five_UTR, args.five_UTR_swapped_strands,args.genes_no_5UTR)