import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='choose_genes_cv')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('promoter_bedfile', type=str, help='Input location of promoter bedfile')
parser.add_argument('Czechowski_rankedcv', type=str, help='Input location of Czechowski et al 2005 ranked cv dataset reanalysed by Will Nash')
parser.add_argument('Mergner_rankedcv', type=str, help='Input location of Mergner et al 2020 ranked cv dataset')
parser.add_argument('no_of_genes', type=int, help='Number of genes in each category to subset')
parser.add_argument('Czechowski_gene_categories', type=str, help='Output location of microarray gene category subsets')
parser.add_argument('Mergner_gene_categories', type=str, help='Output location of RNAseq gene category subsets')
parser.add_argument('promoter_mapped_motifs', type=str, help='Input location of promoter mapped motifs bed file')
parser.add_argument('promoters_filtered_contain_motifs', type=str, help='output location of the promoter bed file filtered so that each promoter contains at least one TFBS')
parser.add_argument('Czechowski_allgenes', type=str, help='Output location of all filtered microarray genes')
parser.add_argument('Mergner_allgenes', type=str, help='Output location of all filtered RNAseq genes')
parser.add_argument('promoters_gff3', type=str, help='Input location of promoters gff3 file')
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
    
def remove_only5UTR(promoter5UTR_bedfile, promoters_gff3):
    """remove genes where only the Araport11 5'UTR is present due to the promoter overlapping other genes"""
    #read in df
    promoter5UTR = pd.read_table(promoter5UTR_bedfile,sep='\t', header=None)
    col = ['chr','start','stop','AGI','dot1', 'strand','source','type','dot2','attributes']
    promoter5UTR.columns = col
    
    # read in promoter only gff3
    promoter_no_5UTR_df = pd.read_table(promoters_gff3, sep='\t', header=None)
    col = ['chr', 'source', 'type', 'start','stop', 'dot1','strand','dot2','attributes']
    promoter_no_5UTR_df.columns = col
    #add AGI column
    promoter_no_5UTR_df_agi = promoter_no_5UTR_df.assign(AGI=promoter_no_5UTR_df.attributes.str.extract(r'ID=gene:(.*?)\;'))
    
    #filter promoters in promoter5UTR but not in promoter_no_5UTR_df_agi
    filtered = promoter5UTR[promoter5UTR.AGI.isin(promoter_no_5UTR_df_agi.AGI)]    
    
    #rename promoter5UTR_bedfile as including genes with only non-overlapping 5'UTRs
    
    oldextension = os.path.splitext(promoter5UTR_bedfile)[1]
    oldname = os.path.splitext(promoter5UTR_bedfile)[0]
    os.rename(promoter5UTR_bedfile, oldname + '_incl_only5UTR' + oldextension)
    
    #make a new file called the same name as promoter5UTR_bedfile
    filtered.to_csv(promoter5UTR_bedfile,sep='\t', header=None, index=False) 
    


def filter_genes_czechowski(promoter_bed, select_genes_file):
    """filter out genes from the microarray data which aren't in the promoter_bed"""
    select_genes = pd.read_table(select_genes_file, sep='\t', header=None)
    cols = ['rank','probe_id','AGI','expression_mean','expression_SD','expression_CV','proportion_of_values_present_in_mas5','presence_in_araport11','constitutive_in_araport11']
    select_genes.columns = cols
    #make AGI uppercase
    select_genes.AGI = select_genes.AGI.str.upper()
    promoters = pd.read_table(promoter_bed, sep='\t', header=None)
    col = ['chr','start','stop','AGI','dot1', 'strand','source','type','dot2','attributes']
    promoters.columns = col

    merged = pd.merge(promoters, select_genes, on='AGI', how='left')
    #remove NaNs in expression_CV column
    filtered_df = merged[merged.expression_CV.notnull()].copy()
    
    #sort by chr and start 
    filtered_df.sort_values(['chr','start'], inplace=True, ignore_index=True) 
    
    #save df
    filtered_df.to_csv(args.Czechowski_allgenes,sep='\t',columns=filtered_df.columns, index=False)
    
    #sort by CV value
    filtered_df.sort_values('expression_CV', inplace=True, ignore_index=True) 
    
    return filtered_df


def filter_genes_mergner(promoter_bed, select_genes_file):
    """filter out genes from the RNA-seq data which aren't in the promoter_bed"""
    select_genes = pd.read_csv(select_genes_file, header=0)
    cols = ['AGI','transcription_class','transcription_family','expression_CV']
    select_genes.columns = cols
    #all present in Araport 11 column
    select_genes['presence_in_araport11'] = 1
    
    promoters = pd.read_table(promoter_bed, sep='\t', header=None)
    col = ['chr','start','stop','AGI','dot1', 'strand','source','type','dot2','attributes']
    promoters.columns = col

    merged = pd.merge(promoters, select_genes, on='AGI', how='left')
    #remove NaNs in expression_CV column
    filtered_df = merged[merged.expression_CV.notnull()].copy()
    
    #sort by chr and start
    filtered_df.sort_values(['chr','start'], inplace=True, ignore_index=True)
    
    #save df
    filtered_df.to_csv(args.Mergner_allgenes,sep='\t',columns=filtered_df.columns,index=False)
    
    #sort by CV value
    filtered_df.sort_values('expression_CV', inplace=True, ignore_index=True)
    

    
    return filtered_df   

def subSet_onCV(in_df, out_dir, no_of_genes):
    '''
    Extract the constitutive, variable, and control subsets based on CV values
    '''
    #filtering based on presence in the Araport11 annotation, define the first
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
    sample = no_of_genes/10
    sample_integar = int(str(sample).replace('.0', '')) #convert sample to an integar
    samples_from_bins = mid_range.groupby('bins').apply(pd.DataFrame.sample, sample_integar, random_state = 2)
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
    
remove_only5UTR(args.promoter_bedfile, args.promoters_gff3)    
remove_proms_no_TFBS(args.promoter_bedfile,args.promoter_mapped_motifs,args.promoters_filtered_contain_motifs)
filtered_czechowski = filter_genes_czechowski(args.promoters_filtered_contain_motifs,args.Czechowski_rankedcv)
filtered_mergner = filter_genes_mergner(args.promoters_filtered_contain_motifs,args.Mergner_rankedcv)
#czechowksi subset
subSet_onCV(filtered_czechowski,args.Czechowski_gene_categories,args.no_of_genes)
#mergner subset
subSet_onCV(filtered_mergner,args.Mergner_gene_categories,args.no_of_genes)