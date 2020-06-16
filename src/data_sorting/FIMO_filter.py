import pandas as pd
from pybedtools import BedTool
import numpy as np
import argparse
#create argparse variables
parser = argparse.ArgumentParser(description='FIMO_filter')
parser.add_argument('fimo_file', type=str, help='Input location of FIMO.tsv file')
parser.add_argument('promoter_bedfile', type=str, help='Input location of promoter bedfile')
parser.add_argument('motifs_bed', type=str, help='Output location of motifs bed file')
parser.add_argument('q_value', type=float, help='q_value threshold for filtering')
args = parser.parse_args()
#functions
def fimo_qfilter(fimo_file, q_value):
    """this uses a meme-suite version 5 fimo.tsv file, filters by a q-value, and returns a pandas df"""
    #read in fimo.tsv file
    fimo = pd.read_table(fimo_file, sep='\t')
    #rename sequence column to just the AGI code
    fimo.sequence_name = fimo.sequence_name.str.extract(r'(.*?)\:')
    #filter q_values to specified threshold
    fimo_qfilter = fimo[fimo['q-value'] <= q_value]    
    return fimo_qfilter
    
def fimo2bed(filtered_fimo_df, promoters_bed, output_bed):
    """This function creates a bed file using fimo.tsv motif file, and the promoter.bed file (chromosome number is used from this). It sorts the bedfile by chromosome then by start"""
    promoters = pd.read_table(promoters_bed, sep='\t', header=None) #read in promoter bed file
    #add column names
    cols = ['chr', 'start', 'stop', 'gene', 'dot', 'strand', 'source', 'type', 'dot2', 'details'] 
    promoters.columns = cols
    #merge promoters.bed with the fimo motif file
    merged = pd.merge(filtered_fimo_df, promoters, left_on='sequence_name', right_on='gene')
    #add motif start position to the promoter start position, and minus 1 (bed file is 0-based), and motif end to promoter start pos (bed file end coord is non-inclusive). Then format with no decimal place
    merged['correct_start'] = (merged['start_x'] + merged['start_y'] -1).astype(np.int64)
    merged['correct_stop'] = (merged['stop_x'] + merged['start_y']).astype(np.int64)
    #create motifs df in bed file column order
    motifs_df = merged.loc[:, ['chr','correct_start','correct_stop', 'motif_id', 'score', 'strand_x', 'sequence_name', 'p-value', 'q-value', 'matched_sequence']]
    #sort the df by chromosome then by motif start position
    motifs_df_sorted = motifs_df.sort_values(['chr','correct_start'])
    #create motifs bed file
    motifs = BedTool.from_dataframe(motifs_df_sorted).saveas(output_bed)
    
#same again using the functions
filtered_fimo = fimo_qfilter(args.fimo_file, args.q_value)
fimo2bed(filtered_fimo, args.promoter_bedfile, args.motifs_bed)