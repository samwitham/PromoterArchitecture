import os
import argparse
import pandas as pd
import numpy as np
import io
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='rolling_window')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('promoter_bedfile', type=str, help='Input location of promoter bedfile')
parser.add_argument('window_bed', type=str, help='Output location of rolling window bed file')
parser.add_argument('window_size', type=int, help='Size of the rolling window in bp')
parser.add_argument('step_size', type=int, help='Size of the window offset in bp')
parser.add_argument('overlapping_proms', type=str, help='Output location of the overlapping promoter bed file')
args = parser.parse_args()


def flag_overlapping_proms(promoter_bed,output_bed):
    """function to take an input promoter bed file and output a bedfile containing a list of promoters which overlap"""
    #read in promoters
    proms_df = pd.read_table(promoter_bed, sep='\t', header=None)
    cols = ['chr','start','stop','AGI','dot1','strand','source','type','dot2','attributes']
    proms_df.columns = cols
    #create bedtools object of promoters
    proms_bed = BedTool(promoter_bed)
    #c = columns to apply function to
    #o = count number of merged promoters, name the first and last promoter that were merged
    merged = proms_bed.merge(c=4, o=['count_distinct','first', 'last'])
    #write to bufer
    merged_buffer = io.StringIO()
    merged_buffer.write(str(merged))
    merged_buffer.seek(0)
    #read as dataframe
    overlapping = pd.read_table(merged_buffer, sep='\t', header=None)
    cols2 = ['chr','start','stop', 'number_of_overlaps', 'first_overlap','second_overlap']
    overlapping.columns = cols2
    #select only features made of more than one promtoer that were merged as overlapping
    overlapping_only = overlapping[overlapping.number_of_overlaps >= 2]
    overlapping_only.to_csv(output_bed,index=False,sep='\t',header=False)  

def window_split(promoter_bed, output_bed, window_size, step_size):
    """function to split promoters into rolling windows"""
    #separate promoters into dfs by strand
    proms_df = pd.read_table(promoter_bed, sep='\t', header=None)
    cols1 = ['chr','start','stop','AGI','dot1','strand','source','type','dot2','attributes']
    proms_df.columns = cols1
    proms_pos = proms_df[proms_df.strand == '+']
    proms_neg = proms_df[proms_df.strand == '-']
    #fool bedtools makewindow so that the first window is actually made from the ATG for the positive strand
    proms_pos_copy = proms_pos.copy()
    proms_pos_copy['length'] = proms_pos.stop-proms_pos.start
    proms_pos_copy['altered_start'] = proms_pos.stop
    proms_pos_copy['altered_stop'] = proms_pos.start + 2*proms_pos_copy.length
    proms_changed = proms_pos_copy[['chr','altered_start','altered_stop','AGI','dot1','strand','source','type','dot2','attributes']]
   
    #write to temporary bed buffers
    pos_buffer = io.StringIO()
    neg_buffer = io.StringIO()
    proms_changed.to_csv(pos_buffer,index=False,sep='\t',header=None)
    proms_neg.to_csv(neg_buffer,index=False,sep='\t',header=None)
    pos_buffer.seek(0)
    neg_buffer.seek(0)
    
    pos_proms = BedTool(pos_buffer)
    neg_proms = BedTool(neg_buffer)
    #make the sliding windows
    #w = window size
    #s = step size
    #n = no. of windows - note there seems to be a bug, window size is one below what you put
    #i = srcwinnum (use source interval name with the window number)
    #reverse = reverse numbering of windows in output - used this for negative strand promoters so window 1 starts at the ATG
    #note - the windows in the reverse strand get cut off if I use a step_size bigger than 1, so I will manually remove windows of the incorrect step size
    windows_pos = BedTool().window_maker(b=pos_proms, w=window_size,s=50, i='srcwinnum')
    windows_neg = BedTool().window_maker(b=neg_proms, w=window_size,s=50, i='srcwinnum')#,reverse=True
    #make buffer
    window_pos_buffer = io.StringIO()
    window_neg_buffer = io.StringIO()
    #write bedfile_like buffer
    window_pos_buffer.write(str(windows_pos))
    window_pos_buffer.seek(0)
    window_neg_buffer.write(str(windows_neg))
    window_neg_buffer.seek(0)
    #create df
    window_pos_df = pd.read_table(window_pos_buffer, sep='\t', header=None)
    window_neg_df = pd.read_table(window_neg_buffer, sep='\t', header=None)
    cols = ['chr', 'start','stop','window']
    window_pos_df.columns = cols
    window_neg_df.columns = cols
    window_neg_df = window_neg_df.astype({'chr':'int','start': 'int', 'stop':'int', 'window':'str'})
    window_pos_df = window_pos_df.astype({'chr':'int','start': 'int', 'stop':'int', 'window':'str'})
    
    
    #reverse the start/stop changes that fooled bedtools makewindow
    pos_df_corrected = window_pos_df.copy()
    ### need to merge this df with proms_pos using AGI code. then make distance the original stop - 1 to the chunk stop
    #Make AGI column
    pos_df_corrected = pos_df_corrected.assign(AGI=pos_df_corrected.window.str.extract(r'(.*?)\_'))

    #Merge with proms_pos using AGI code
    merged = pd.merge(pos_df_corrected, proms_pos, on='AGI', how='left',suffixes=('','_wholeprom'))    
    
    merged['distance'] = merged.stop-merged.stop_wholeprom
    #create window length column
    merged['window_length'] = merged.stop-merged.start
    
    merged['correct_start'] = merged.stop - 2*merged.distance
    merged['correct_stop'] = merged.correct_start + merged.window_length
    #filter by correct_stop if it is outside the range of the promoter

    merged = merged[['chr','correct_start','correct_stop','window']]
    #rename columns
    merged.rename(columns={'correct_start':'start', 'correct_stop':'stop'}, inplace=True)
          
    
    #Merge positive and negative strand windows
    merged_all = pd.merge(window_neg_df,merged, how='outer')
    #merged_all.to_csv(output_bed,index=False,sep='\t',header=None)
    
    #filter lengths so they are only = 100bp
    window_lengths = merged_all.copy()
    window_lengths =  merged_all.assign(length=((merged_all.stop) - merged_all.start))
    removed= window_lengths.loc[window_lengths.length == 100]
    
    #sort by chr, start
    sorted_removed = removed.sort_values(['chr','start']).reset_index(drop=True)
    #remove length column
    sorted_removed = sorted_removed[['chr','start','stop','window']]
    #write to bed file
    sorted_removed.to_csv(output_bed,index=False,sep='\t',header=None)
    pos_buffer.close()
    neg_buffer.close()
    window_pos_buffer.close()
    window_neg_buffer.close()
    
directory_path = '../..'



#make directory for the output files to be exported to
#dirName = f'{args.directory_path}/data/output/{args.file_names}'
dirName = f'{directory_path}/data/output/{args.file_names}'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
#make directory for the output files to be exported to
#dirName = f'{args.directory_path}/data/output/{args.file_names}'
dirName = f'{directory_path}/data/output/{args.file_names}/rolling_window'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
#make directory for the output files to be exported to
#dirName = f'{args.directory_path}/data/output/{args.file_names}'
dirName = f'{directory_path}/data/output/{args.file_names}/rolling_window/TFBS_coverage_rw'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
flag_overlapping_proms(args.promoter_bedfile,args.overlapping_proms)
proms = window_split(args.promoter_bedfile, args.window_bed, args.window_size, args.step_size)
