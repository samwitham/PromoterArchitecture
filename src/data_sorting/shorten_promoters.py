#Shorten promtoters to certain length upstream of ATG after doing FIMO analysis
from pybedtools import BedTool
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='shorten_promoter')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('promoterpref', type=str, help='Name of the prefix name for the promoters extracted')
parser.add_argument('promoter_5UTR_bed', type=str, help='Location of the extracted promoters/5UTRs')
parser.add_argument('promoter_length', type=int, help='Length to shorten promoters to')
args = parser.parse_args()

def shorten_promoters(promoter_5UTR_bed,promoter_length,shortened_promoter_5UTR_bed):
    """shorten promoters to certain length upstream of ATG and output to a file."""
    promoters = pd.read_table(promoter_5UTR_bed, sep='\t', header=None)
    cols = ['chr','start','stop','AGI','dot1','strand','source','type','dot2','attributes']
    promoters.columns = cols
    #shorten promoters to the promoter_length upstream of ATG
    promoters.loc[promoters.strand == '+', 'start'] = promoters.loc[promoters.strand == '+', 'stop'] - promoter_length
    promoters.loc[promoters.strand == '-', 'stop'] = promoters.loc[promoters.strand == '-', 'start'] + promoter_length
    #write to output file
    promoters.to_csv(shortened_promoter_5UTR_bed, sep='\t',header=None,index=False)
#set output location of shortened promoters
shortened_promoter_5UTR_bed = f'../../data/output/{args.file_names}/FIMO/{args.promoterpref}_{str(args.promoter_length)}bp.bed'    
#shorten promoters
shorten_promoters(args.promoter_5UTR_bed,args.promoter_length,shortened_promoter_5UTR_bed)