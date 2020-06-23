#Create a shortened promoter motifs_mapped file in correct format for the TF_diversity_plots_{shortenedpromoter}.py script
from pybedtools import BedTool
import pandas as pd
import argparse
import io

parser = argparse.ArgumentParser(description='motifs_mapped_intersect')
parser.add_argument('promoter_bed', type=str, help='Input location of shortened promoter file')
parser.add_argument('mapped_motif_bed', type=str, help='Input location of promoters mapped motif bed')
parser.add_argument('output_motifs_bed', type=str, help='Output location of mapped_motifs bed')
args = parser.parse_args()

def merge_windows_motifs(promoter_bed,mapped_motif_bed,output_file):
    """perform bedtools intersect on the two dfs"""
    promoters = BedTool(promoter_bed)
    motifs = mapped_motif_bed
    #-wao =Write the original A and B entries plus the number of base pairs of overlap between the two features.
    #However, A features w/o overlap are also reported with a NULL B feature and overlap = 0
    intersect = promoters.intersect(motifs, wao=True)
    #Write to output_buffer
    tempbuffer = io.StringIO()              
    tempbuffer.write(str(intersect))    
    #go back to beginning of the buffer
    tempbuffer.seek(0)
    #read buffer as pandas df
    df = pd.read_table(tempbuffer, sep='\t', header=None)
    cols = ['chr', 'start', 'stop', 'promoter_AGI','dot1','strand','source','type','dot2','attributes',
            'motif_chr','motif_start','motif_stop','name_rep', 'score', 'motif_strand',
            'promoter_AGI2', 'p-value', 'q-value', 'matched_sequence', 'TF_name', 'TF_family', 'TF_AGI','bp_overlap']
    df.columns = cols
#     #filter out unwanted columns
#     df=df[['motif_chr','motif_start','motif_stop','name_rep', 'score', 'motif_strand',
#             'promoter_AGI2', 'p-value', 'q-value', 'matched_sequence', 'TF_name', 'TF_family', 'TF_AGI']]
    #write file
    df.to_csv(output_file, sep='\t', index=False, header=None)
#merge windows with motifs and create output bed file
merge_windows_motifs(args.promoter_bed,args.mapped_motif_bed,args.output_motifs_bed)