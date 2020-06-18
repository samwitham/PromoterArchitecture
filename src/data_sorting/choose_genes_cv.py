import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='choose_genes_cv')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('promoter_bedfile', type=str, help='Input location of promoter bedfile')
parser.add_argument('Czechowski_rankedcv', type=str, help='Input location of Czechowski et al 2005 ranked cv dataset reanalysed by Will Nash')
parser.add_argument('no_of_genes', type=int, help='Number of genes in each category to subset')
parser.add_argument('Czechowski_gene_categories', type=str, help='Output location of gene category subsets')
parser.add_argument('promoter_mapped_motifs', type=str, help='Input location of promoter mapped motifs bed file')
parser.add_argument('promoters_filtered_contain_motifs', type=str, help='output location of the promoter bed file filtered so that each promoter contains at least one TFBS')

args = parser.parse_args()

def remove_proms_no_TFBS(promoter_bedfile, promoter_mapped_motifs,promoters_filtered_contain_motifs):
    """remove promoters which had no TFBSs found within them after filtering the FIMO output. Create output file of these"""
    promoters = pd.read_table(promoter_bedfile, sep='\t', header=None)
    col = ['chr','start','stop','promoter_AGI','dot1', 'strand','source','type','dot2','attributes']
    promoters.columns = col
    mapped_motifs = pd.read_table(promoter_mapped_motifs, sep='\t', header=None)
    col2 = ['chr', 'start', 'stop', 'name_rep', 'score', 'strand', 'promoter_AGI', 'p-value', 'q-value', 'matched_sequence', 'TF_name', 'TF_family', 'TF_AGI']
    mapped_motifs.columns = col2

    merged = pd.merge(promoters,mapped_motifs, on='promoter_AGI', how='left',  suffixes=['', '_y'])
    #remove NaNs in TF_AGI column
    filtered_df = merged[merged.TF_AGI.notnull()]
    #reduce columns
    filtered_df = filtered_df[col]
    #filter duplicates
    idx = filtered_df.promoter_AGI.drop_duplicates().index
    #this will return filtered df
    no_dups = filtered_df.loc[idx,:]
    no_dups.to_csv(promoters_filtered_contain_motifs, sep='\t', header=None, index=False)    

def filter_genes(promoter_bed, select_genes_file):
    """filter out genes from the microarray data which aren't in the promoter_bed"""
    select_genes = pd.read_table(select_genes_file, sep='\t', header=None)
    cols = ['rank','probe_id','AGI','expression_mean','expression_SD','expression_CV','proportion_of_values_present_in_mas5','presence_in_araport11','constitutive_in_araport11']
    select_genes.columns = cols
    
    promoters = pd.read_table(promoter_bed, sep='\t', header=None)
    col = ['chr','start','stop','AGI','dot1', 'strand','source','type','dot2','attributes']
    promoters.columns = col

    merged = pd.merge(promoters, select_genes, on='AGI', how='left')
    #remove NaNs in expression_CV column
    filtered_df = merged[merged.expression_CV.notnull()]
    
    return filtered_df   

def subSet_onCV(in_df, out_dir, no_of_genes):
    '''
    Extract the constitutive, variable, and control subsets based on CV values
    '''
    #fliltering based on presence in the Araport11 annotation, define the first
    #n rows as the constitutive set and add label
    constitutive          = in_df[in_df.presence_in_araport11 == 1][0:no_of_genes]
    constitutive['state'] = 'constitutive'

    #define the last n rows as the variable set and add label
    variable          = in_df[in_df.presence_in_araport11 == 1][-no_of_genes:]
    variable['state'] = 'variable'

    #extract the rest of the rows as the control search space
    mid_range    = in_df[in_df.presence_in_araport11 == 1][(no_of_genes+1):-(no_of_genes+1)]

    #create 10 labelled bins
    mid_range['bins'] = pd.Series(pd.qcut(mid_range['expression_CV'], q = 10, precision = 2))

    #extract 10 random rows from these bins and label as the control set
    samples_from_bins = mid_range.groupby('bins').apply(pd.DataFrame.sample, 10, random_state = 2)
    samples_from_bins['state'] = 'control'

    #concatenate and write as output
    output_set = pd.concat([constitutive[['AGI', 'state']], variable[['AGI', 'state']], samples_from_bins[['AGI', 'state']]], ignore_index = True)
    output_set.to_csv(out_dir, sep = '\t', index = False, header=False)
    
    #function from expressionVar_subsets_plot.py
    #__author__ = "Will Nash"
    # __copyright__ = "Copyright 2020, The Earlham Institute"
    # __credits__ = ["Will Nash", "Wilfried Haerty"]
    # __license__ = "GPL"
    # __version__ = "1.0"
    # __maintainer__ = "Will Nash"
    # __email__ = "will.nash@earlham.ac.uk"
    # __status__ = "Testing"
    #__modified_by__ "Sam Witham"


#make directory for the output files to be exported to
#dirName = f'{args.directory_path}/data/output/{args.file_names}'
dirName = f'../../data/output/{args.file_names}/genes'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
    
remove_proms_no_TFBS(args.promoter_bedfile,args.promoter_mapped_motifs,args.promoters_filtered_contain_motifs)    
filtered = filter_genes(args.promoters_filtered_contain_motifs,args.Czechowski_rankedcv)
subSet_onCV(filtered,args.Czechowski_gene_categories,args.no_of_genes)