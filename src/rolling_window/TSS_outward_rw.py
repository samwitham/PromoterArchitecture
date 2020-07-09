#merge the promoters and 5UTR rolling windows, 
#with window numbers going head to head around the Araport TSS (so window 1 of promoters becomes -1, while 5UTR window numbers stay the same)
import pandas as pd
import argparse
#define arguments
parser = argparse.ArgumentParser(description='create_5UTRs')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('five_UTR_windows', type=str, help='Input location of 5UTR sliding windows')
parser.add_argument('promoter_5UTR_bed', type=str, help='Input location of promoter 5UTR bed file to get strand information')
parser.add_argument('promoter_windows', type=str, help='Input location of promoter sliding windows')
parser.add_argument('TSS_outward_windows', type=str, help='Output location of the promoter + 5UTR sliding window numbers going outwards from the TSS')
parser.add_argument('window_size', type=int, help='Size of the rolling window in bp')
parser.add_argument('step_size', type=int, help='Size of the window offset in bp')
args = parser.parse_args()

def merge_window_files(promoter_windows,five_UTR_windows,promoter_5UTR_bed,TSS_outward_windows,window_size,step_size):
    """function to merge the promoters and 5UTR rolling windows
    with window numbers going head to head around the Araport TSS 
    (so window 1 of promoters becomes -1, while 5UTR window numbers stay the same). Also creates a window 0 across the TSS"""
    #read in the files
    fiveUTR=pd.read_table(five_UTR_windows, sep='\t',header=None)
    proms = pd.read_table(promoter_windows, sep='\t',header=None)
    promoters_5UTR = pd.read_table(promoter_5UTR_bed, sep='\t',header=None)    
    cols = ['chr','start','stop','name']
    cols2 = ['chr','start','stop','AGI','dot','strand','source','type','dot2','attributes']
    #rename columns
    fiveUTR.columns = cols
    proms.columns = cols
    promoters_5UTR.columns = cols2
    #split name into AGI and window number
    fiveUTR['AGI'] = fiveUTR.name.str.split('_',expand=True)[0]
    fiveUTR['window_number'] = fiveUTR.name.str.split('_',expand=True)[1]
    proms['AGI'] = proms.name.str.split('_',expand=True)[0]
    proms['window_number'] = proms.name.str.split('_',expand=True)[1]
    
    #first turn window number column into integars
    proms = proms.astype({'window_number': 'int'})
    fiveUTR = fiveUTR.astype({'window_number': 'int'})
    #make proms window numbers negative
    proms.window_number = -proms.window_number
    #merge fiveUTR and promoter_5UTR_bed on AGI
    fiveUTR_merged = pd.merge(fiveUTR,promoters_5UTR,on='AGI', how='left',suffixes=('','_proms'))
    #make an extra window upstream of the fiverUTR window 1 called window 0 (using window size and offset variables)
    #filter only window 1s
    fiveUTR_merged = fiveUTR_merged[fiveUTR_merged.window_number == 1]
    #create an empty list
    new_rows = []
    #iterate through fiveUTR_copy window number ones and append a row to fiveUTR with window 0
    for i,data in fiveUTR_merged.iterrows():
        #if strand is positive
        if fiveUTR_merged.loc[i,'strand'] == '+':
            #new window start is the window 1 start minus the step size
            start = fiveUTR_merged.loc[i,'start'] - step_size 
            #new window stop is the start plus the window size
            stop = start + window_size
            dict1 = {'chr':fiveUTR_merged.loc[i,'chr'],'start':start,'stop':stop, 'window_number':0,'AGI':fiveUTR_merged.loc[i,'AGI']}
            #append dict to new_rows list
            new_rows.append(dict1)
        #else of strand is negative
        elif fiveUTR_merged.loc[i,'strand'] == '-':
            #new window start is the window 1 stop minus the step size 
            start = fiveUTR_merged.loc[i,'stop'] - step_size
            #new window stop is start plus window size
            stop = start + window_size
            dict1 = {'chr':fiveUTR_merged.loc[i,'chr'],'start':start,'stop':stop, 'window_number':0,'AGI':fiveUTR_merged.loc[i,'AGI']}
            #append dict to new_rows list
            new_rows.append(dict1)
    #turn new_rows into a df
    rows_df = pd.DataFrame(new_rows)
    #remove name columns from fiveUTR
    fiveUTR = fiveUTR[['chr','start','stop','window_number','AGI']]
    #filter name column fromproms
    proms = proms[['chr','start','stop','window_number','AGI']]
    #concatenate rows_df and fiveUTR and proms
    all_windows = pd.concat([fiveUTR,rows_df,proms],ignore_index=True)
    #sort by chr and start
    all_windows.sort_values(['chr','start'], inplace=True, ignore_index=True)
    #merge AGI and window number into single column
    #turn window_number into string
    all_windows = all_windows.astype({'window_number': 'str'})
    all_windows['name'] = all_windows.AGI + '_' + all_windows.window_number
    #filter columns
    all_windows = all_windows[['chr','start','stop','name']]
    #save file
    all_windows.to_csv(TSS_outward_windows,sep='\t',index=None,header=None)
    
merge_window_files(args.promoter_windows,args.five_UTR_windows,args.promoter_5UTR_bed,args.TSS_outward_windows,args.window_size,args.step_size)