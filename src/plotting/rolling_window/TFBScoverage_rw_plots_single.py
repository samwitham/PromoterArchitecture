import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.formula.api import ols
import researchpy as rp
import os
import argparse

parser = argparse.ArgumentParser(description='TFBScoverage_plots_rw')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('Czechowski_gene_categories', type=str, help='Input location of Czechowski gene categories text file')
parser.add_argument('bp_covered', type=str, help='Input location of TFBS base pairs covered file')
parser.add_argument('EPD_TSS_bed', type=str, help='Input location of eukaryotic promoter database transcription start site bed file')
parser.add_argument('promoter_bed', type=str, help='Input location of promoter bed file')       
parser.add_argument('promoter_no_5UTR', type=str, help='Input location of promoter no 5UTR bed file')  
parser.add_argument('foldername_prefix', type=str, help='Output folder name prefix to use')
parser.add_argument('root_chrom_bp_covered', type=str, help='Input location of root chromatin bed file')
parser.add_argument('shoot_chrom_bp_covered', type=str, help='Input location of shoot chromatin bed file')
parser.add_argument('rootshootintersect_chrom_bp_covered', type=str, help='Input location of rootshootintersect chromatin bed file')
parser.add_argument('palette', type=str, help='Optional replacement colour palette for plots',default = None, nargs="?")
parser.add_argument('author_name', type=str, help='Optional author name to add to output file names',default = 'Czechowski', nargs="?")
parser.add_argument('variable1_name', type=str, help='Optional replacement name for 2nd variable eg. non-specific',default = 'constitutive', nargs="?")
parser.add_argument('variable2_name', type=str, help='Optional variable 2 name eg. tissue_specific',default = 'variable', nargs="?")
args = parser.parse_args()


#make directory for the plots to be exported to
dirName = f'../../data/output/{args.file_names}/rolling_window/'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")


#make directory for the plots to be exported to
dirName = f'../../data/output/{args.file_names}/rolling_window/{args.foldername_prefix}'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")


#make directory for the plots to be exported to
dirName = f'../../data/output/{args.file_names}/rolling_window/{args.foldername_prefix}/plots'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
    
def process_input_files(bp_covered):
    """process and merge the input files into a df"""
    #read in coverage bedfile
    coverage_df = pd.read_table(bp_covered, sep='\t', header=None)
    col = ['chr','start','stop','name','no._of_overlaps', 'no._of_bases_covered','window_length','fraction_bases_covered']
    coverage_df.columns = col
    #add % bases covered column
    coverage_df['percentage_bases_covered'] = coverage_df.fraction_bases_covered * 100
    
    #Add AGI and window_number column to coverage_df
    coverage_df = coverage_df.assign(AGI=coverage_df.name.str.extract(r'(.*?)\_'))
    coverage_df['window_number'] = coverage_df.name.str.split('_', n=2, expand = True)[1]
    #coverage_df = coverage_df.assign(window_number=coverage_df.name.str[-1:])
    #remove windows with fewer than 100 promoters extending to that location
    coverage_df = coverage_df[coverage_df['window_number'].map(coverage_df['window_number'].value_counts()) > 99] 
    
    #allow colour codes in seaborn
    sns.set(color_codes=True)
    sns.set_style("ticks")
    sns.set_palette(args.palette)
    
    return coverage_df

def add_coverage(df,coverage_bed,suffix):
    """add % bp covered data from a bed file to the df. Prefix is a name added to any new columns"""
    #read in bed file
    coverage_df = pd.read_table(coverage_bed, sep='\t', header=None)
    cols = ['chr','start','stop','name','no._of_overlaps', 'no._of_bases_covered','window_length','fraction_bases_covered']
    coverage_df.columns = cols
    #add % bases covered column
    coverage_df['percentage_bases_covered'] = coverage_df['fraction_bases_covered'] * 100
    #filter columns
    coverage_df = coverage_df[['chr','start','stop','name','percentage_bases_covered']]
    #make df columns integars
    df = df.astype({'stop':'int','start':'int','chr':'int'})
    coverage_df = coverage_df.astype({'stop':'int','start':'int','chr':'int'})
    #merge the dfs
    merged = pd.merge(df,coverage_df, how='left',on=['chr','start','stop'],suffixes=('', f'_{suffix}'))
    #remove NaN
    #merged = merged[merged['name'].notnull()] 
    return merged

def rep_sample(df, col, n, random_state):
    """function to return a df with equal sample sizes
    taken from here: https://stackoverflow.com/questions/39457762/python-pandas-conditionally-select-a-uniform-sample-from-a-dataframe"""
    #identify number of categories
    nu = df[col].nunique()
    # find number of rows
    m = len(df)
    # integar divide total sample size by number of categories
    mpb = n // nu
    # multiply this by the number of categories and subtract from the number of samples to find the remainder
    mku = n - mpb * nu
    # make an array fileld with zeros corresponding to each category
    fills = np.zeros(nu)

    # make values in the array 1s up until the remainder
    fills[:mku] = 1

    # calculate sample sizes for each category
    sample_sizes = (np.ones(nu) * mpb + fills).astype(int)

    #group the df by categories
    gb = df.groupby(col)
    #define sample size function
    sample = lambda sub_df, i: sub_df.sample(sample_sizes[i], random_state = random_state)
    #run sample size function on each category
    subs = [sample(sub_df, i) for i, (_, sub_df) in enumerate(gb)]
    #return concatenated sub dfs
    return pd.concat(subs)


def windows_coords(output_prefix,variable_of_interest,variable_of_interest_name,variable_of_interest_df,promoter_bed,promoter_no_5UTR,
                   window_offset,EPD_TSS_bed,includeEPDTSS=False,includechromatin=False,chromatin_tissue_variable='percentage_bases_covered_rootshootintersect_chrom',
                   chromatin_tissue_variable_name='% open chromatin root and shoot intersect',x_range=False,estimator='median',ci=95, n_boot=10000, 
                   genetype=False, genetype2=False, genetype3=False):
    """function to add the centre of each window corresponding to each window no. and return a lineplot. Also add promoter length distribution, Araport TSS distribution,
    EPD TSS distribution (add the most common TSS as documented on eukaryotic promoter database Arabidopsis last modified on EPD 06/06/2018)"""
    #read in bed file
    promoter_df = pd.read_table(promoter_bed, sep='\t', header=None)
    col = ['chr','start','stop','AGI','dot1', 'strand','source','type','dot2','attributes']
    promoter_df.columns = col
    #add promoter length column
    promoter_df['length'] = promoter_df.stop-promoter_df.start  
    #merge promoter_bed with variable_of_interest_df on AGI
    merged = pd.merge(variable_of_interest_df, promoter_df, on='AGI',  how='left',suffixes=('','_wholeprom'))
    #remove NaN
    merged = merged[merged[variable_of_interest].notnull()]
    #make columns integars
    merged = merged.astype({'stop_wholeprom':'int','start_wholeprom':'int','start':'int','stop':'int'})
    #split merged into 2 dfs by strand
    pos = merged[merged.strand == '+'].copy()
    neg = merged[merged.strand == '-'].copy()
    
    #add variable of interest position column where position is the middle of the window, with whole length of the longest promoter
    pos['position'] = (pos.stop_wholeprom) - (pos.start+0.5*(pos.stop-pos.start))
    neg['position'] = (neg.start+0.5*(neg.stop-neg.start)) - neg.start_wholeprom
    merged2 = pd.merge(pos,neg,how='outer')

    merged2 = merged2.astype({'position': 'int64'})
    #get longest transcript TSS distribution (if Araport 11 definition used)
    promoter_no_5UTR_df = pd.read_table(promoter_no_5UTR, sep='\t', header=None)
    col = ['chr', 'source', 'type', 'start','stop', 'dot1','strand','dot2','attributes']
    promoter_no_5UTR_df.columns = col
    #add AGI column
    promoter_no_5UTR_df_agi = promoter_no_5UTR_df.assign(AGI=promoter_no_5UTR_df.attributes.str.extract(r'ID=gene:(.*?)\;'))
    #merged with windows    
    merged2 = pd.merge(merged2,promoter_no_5UTR_df_agi, on='AGI', how='left', suffixes=('','_no_UTR'))
    #remove NaN (promoters in promoters.gff but not in promoters_5UTR)
    merged2 = merged2[merged2.length.notnull()]
    #remove NaN (promoters in promoters_5UTR but not in promoters.gff - ie. only 5'UTRs)
    merged2 = merged2[merged2.chr_no_UTR.notnull()]
    #Get most common transcript TSS location from eukaryotic promoter database (last modified on EPD 06/06/2018)
    #Note - not all promoters have a TSS on EPD
    if includeEPDTSS==True:
        EPD_TSS_df = pd.read_table(EPD_TSS_bed, delim_whitespace=True, header=None, skiprows=4)
        cols = ['chr','start','stop','transcript_EPD','score_EPD','strand_EPD','thickstart_EPD','thickend_EPD']
        EPD_TSS_df.columns = cols
        #add AGI column
        EPD_TSS_df['AGI'] = EPD_TSS_df.transcript_EPD.str.split('_',expand=True)[0]
        #add TSS location column
        EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '+', 'TSS_EPD'] = EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '+', 'thickstart_EPD']
        EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '-', 'TSS_EPD'] = EPD_TSS_df.loc[EPD_TSS_df.strand_EPD == '-', 'thickend_EPD'] -1
        #merged with windows
        merged2 = pd.merge(merged2,EPD_TSS_df, on='AGI', how='left', suffixes=('','_EPD'))

        #remove NaN (promoters in EPD but not in promoters_5UTR)
        merged2 = merged2[merged2.length.notnull()]
        #transfrom EPD TSS column in the same way as the position column
        merged2.loc[merged2.strand == '-', 'TSS_transformed_EPD'] = merged2.loc[merged2.strand == '-', 'TSS_EPD'] - merged2.loc[merged2.strand == '-','start_wholeprom']
        merged2.loc[merged2.strand == '+', 'TSS_transformed_EPD'] = merged2.loc[merged2.strand == '+', 'stop_wholeprom'] - merged2.loc[merged2.strand == '+', 'TSS_EPD']
        #make integars
        merged2 = merged2.astype({'TSS_transformed_EPD':'float64'})
    
    #calculate longest promoter length based on window cutoff
    number_of_windows = len(variable_of_interest_df.window_number.unique())
    window_length = variable_of_interest_df.window_length.max()    
    length_of_longest_promoter = number_of_windows*(window_length-window_offset)

    #make 0 start like in bed files
    merged2.start_no_UTR = merged2.start_no_UTR -1
    #add Araport TSS location column
    #merged2['TSS'] = int()
    merged2.loc[merged2.strand == '+', 'TSS'] = merged2.loc[merged2.strand == '+', 'stop_no_UTR']
    merged2.loc[merged2.strand == '-', 'TSS'] = merged2.loc[merged2.strand == '-', 'start_no_UTR'] -1
    #transform TSS location in the same way as the position column
    merged2.loc[merged2.strand == '-', 'TSS_transformed_Araport11'] = merged2.loc[merged2.strand == '-', 'TSS'] - merged2.loc[merged2.strand == '-','start_wholeprom']
    merged2.loc[merged2.strand == '+', 'TSS_transformed_Araport11'] = merged2.loc[merged2.strand == '+', 'stop_wholeprom'] - merged2.loc[merged2.strand == '+', 'TSS']
   
    
    #make integars
    merged2 = merged2.astype({'start_no_UTR':'float64','stop_no_UTR':'float64','TSS':'float64','TSS_transformed_Araport11':'float64',f'{variable_of_interest}':'float64',f'{chromatin_tissue_variable}':'float64'})
    #return merged2[['AGI','strand','start','stop','start_wholeprom','stop_wholeprom','start_no_UTR','stop_no_UTR','TSS','TSS_transformed','position','chr_no_UTR','window_number']]


        
    #change estimator
    if estimator == 'mean':
        new_estimator=estimator
        
    if estimator == 'median':
        new_estimator = np.median
    #set number of subplots so can easily change all output possibilities, where subplotA is the top
    subplots = 2
    
    #make subplots
#     if includeEPDTSS == True:
#         subplots = subplots + 1
#         f, axes = plt.subplots(subplots, figsize=(10,10))
#         OpenChromplot = axes[subplots-subplots]
#         Araport11TSSplot = axes[subplots-(subplots-1)]
#         EPDTSSplot = axes[subplots-(subplots-2)]   
#         #promlengthsplot = axes[subplots-(subplots-3)]
#         variableofinterestplot = axes[subplots-(subplots-3)]
#     else:
#         f, axes = plt.subplots(subplots, figsize=(10,8))
#         OpenChromplot = axes[subplots-subplots]
#         Araport11TSSplot = axes[subplots-(subplots-1)]
#         #promlengthsplot = axes[subplots-(subplots-2)]
#         variableofinterestplot = axes[subplots-(subplots-2)]


   
    
    #check the plot axes variables are there. If they are not, assign None to them
    try:
        OpenChromplot
    except NameError:
        OpenChromplot = None
    try:
        Araport11TSSplot
    except NameError:
        Araport11TSSplot = None
    try:
        EPDTSSplot
    except NameError:
        EPDTSSplot = None
#     try:
#         promlengthsplot
#     except NameError:
#         promlengthsplot = None
    try:
        variableofinterestplot
    except NameError:
        variableofinterestplot = None
    
    #If EPD TSS plot is present, filter promoters which aren't in EPD to remove NaNs
    if EPDTSSplot != None:
        #remove NaN (promoters in promoters_5UTR but not in promoters.gff - ie. only 5'UTRs)
        merged2 = merged2[merged2.TSS_transformed_EPD.notnull()]
    
    if genetype!=False:
        #filter so only genetype subset present
        merged2 = merged2[merged2.gene_type.notnull()]
        #remove windows with fewer than 50 promoters extending to that location if looking at specific genetypes
        merged2 = merged2[merged2['window_number'].map(merged2['window_number'].value_counts()) > 49]
        #redefine longest promoter length based on window cutoff
        number_of_windows = len(merged2.window_number.unique())
        window_length = merged2.window_length.max()    
        length_of_longest_promoter = number_of_windows*(window_length-window_offset)
        
    #make all values of interest negative as upstream from ATG
    merged_positive = merged2.copy()
    if includeEPDTSS==True:
        merged2[['length','TSS_transformed_Araport11','position','TSS_transformed_EPD']] = -merged2[['length','TSS_transformed_Araport11','position','TSS_transformed_EPD']]
    else:
        merged2[['length','TSS_transformed_Araport11','position']] = -merged2[['length','TSS_transformed_Araport11','position']]

    if genetype==False:
        
        #length_of_longest_promoter = merged_positive.length.max()
        #if openchromplot variable present, add that plot
       
        #if variableofinterestplot variable present, add that plot

        #variable of interest lineplot
        plot = sns.lineplot(y=merged2[variable_of_interest], x=merged2.position, estimator=new_estimator,ci=ci, n_boot=n_boot).get_figure()
        #set titles and axes labels
        #plt.set_title(f'All promoters sliding windows {variable_of_interest_name}', weight='bold')
        plt.ylabel(f'{estimator} {variable_of_interest_name}')
        plt.xlabel('position upstream of ATG')

               
    elif genetype2==False:
        #filter so only genetype subset present
        merged2 = merged2[merged2.gene_type.notnull()]
        
          #if openchromplot variable present, add that plot

        

        #if variableofinterestplot variable present, add that plot
        #variable of interest lineplot
        plot = sns.lineplot(y=merged2[merged2.gene_type == genetype][variable_of_interest], 
                        x=merged2[merged2.gene_type == genetype].position, estimator=new_estimator,ci=ci, n_boot=n_boot).get_figure() 
        #set titles and axes labels
        #plt.set_title(f'{genetype} {variable_of_interest_name}', weight='bold')
        plt.ylabel(f'{estimator} {variable_of_interest_name}')
        plt.xlabel('position upstream of ATG')
            
        
        #set y axis as maximum mean window % bp covered of all genetype subset        
        #variableofinterestplot.set_ylim([0,merged2.groupby('window_number')[variable_of_interest].median().max()+20])
        #set x axis range if specified 
#         if x_range==False:
#             pass
#         else:
#             length_of_longest_promoter = x_range
        
#         #for all subplots:
#         for n in axes:
#             #remove grids
#             n.grid(False)
#             n.set_xlim([-length_of_longest_promoter,0])        
#         f.tight_layout()
        
    elif genetype3==False:
        #filter so only genetype subset present
        merged2 = merged2[merged2.gene_type.notnull()]
        
        ## make a subselection of categories so all sample sizes are equal
        # first select only the relevant genetypes
        merged2 = merged2[merged2.gene_type.isin([genetype,genetype2])]
        # make each promoter unique
        merged2_unique = merged2.drop_duplicates('AGI')
        # identify sample size of the minimum category
        minimum_sample_size = merged2_unique.gene_type.value_counts().min()
        # print this
        print(f'sample size in each category = {minimum_sample_size}')
        # multiply this by the number of categories
        total_sample_size = minimum_sample_size * len(merged2_unique.gene_type.unique())
        #select equal sample sizes of each category with a random state of 1 so it's reproducible
        equal_samplesizes = rep_sample(merged2_unique, 'gene_type',total_sample_size,random_state = 1)
        # now filter out genes which were not selected using the minimum sample size
        to_remove = merged2_unique[~merged2_unique.AGI.isin(equal_samplesizes.AGI)]
        merged2 = merged2[~merged2.AGI.isin(to_remove.AGI)]
        
        
        #if variableofinterestplot variable present, add that plot

        #lineplot variable of interest
        plot=sns.lineplot(y=merged2[merged2.gene_type == genetype][variable_of_interest], x=merged2[merged2.gene_type == genetype].position,
                        label=genetype,estimator=new_estimator,ci=ci, n_boot=n_boot,).get_figure()
        sns.lineplot(y=merged2[merged2.gene_type == genetype2][variable_of_interest], x=merged2[merged2.gene_type == genetype2].position,
                        label=genetype2,estimator=new_estimator,ci=ci, n_boot=n_boot) 
            #set titles & axes names
        #plt.set_title(f'{variable_of_interest_name}', weight='bold')
        plt.ylabel(f'{estimator} {variable_of_interest_name}')
        plt.xlabel('position upstream of ATG')
        #set y axis as maximum mean window % bp covered of all genetype subset        
        #axes[2].set_ylim([0,merged2.groupby('window_number').percentage_bases_covered.median().max()+20])
        #gene_type labels
        gene_type_labels = [genetype, genetype2]
        # Create the legend
        plt.legend()

#         if x_range==False:
#             pass
#         else:
#             length_of_longest_promoter = x_range

#         #for all subplots:
#         for n in axes:
#             #remove grids
#             n.grid(False)
#             n.set_xlim([-length_of_longest_promoter,0])
#         f.tight_layout()
    else:
        #filter so only genetype subset present
        merged2 = merged2[merged2.gene_type.notnull()]
        
        ## make a subselection of categories so all sample sizes are equal    
        # make each promoter unique
        merged2_unique = merged2.drop_duplicates('AGI')
        # identify sample size of the minimum category
        minimum_sample_size = merged2_unique.gene_type.value_counts().min()
        # print this
        print(f'sample size in each category = {minimum_sample_size}')
        # multiply this by the number of categories
        total_sample_size = minimum_sample_size * len(merged2_unique.gene_type.unique())
        #select equal sample sizes of each category with a random state of 1 so it's reproducible
        equal_samplesizes = rep_sample(merged2_unique, 'gene_type',total_sample_size,random_state = 1)
        # now filter out genes which were not selected using the minimum sample size
        to_remove = merged2_unique[~merged2_unique.AGI.isin(equal_samplesizes.AGI)]
        merged2 = merged2[~merged2.AGI.isin(to_remove.AGI)]
        
        #if variableofinterestplot variable present, add that plot
        #lineplot
        plot=sns.lineplot(y=merged2[merged2.gene_type == genetype][variable_of_interest], x=merged2[merged2.gene_type == genetype].position,
                        label=genetype,estimator=new_estimator,ci=ci, n_boot=n_boot).get_figure()
        sns.lineplot(y=merged2[merged2.gene_type == genetype2][variable_of_interest], x=merged2[merged2.gene_type == genetype2].position,
                        label=genetype2,estimator=new_estimator,ci=ci, n_boot=n_boot)
        sns.lineplot(y=merged2[merged2.gene_type == genetype3][variable_of_interest], x=merged2[merged2.gene_type == genetype3].position,
                        label=genetype3,estimator=new_estimator,ci=ci, n_boot=n_boot)
            #set titles & axes names      
        #plt.set_title(f'{variable_of_interest_name}', weight='bold')
        plt.ylabel(f'{estimator} {variable_of_interest_name}')
        plt.xlabel('position upstream of ATG')
        #set y axis as maximum mean window % bp covered of all genetype subset        
        #axes[2].set_ylim([0,merged2.groupby('window_number').percentage_bases_covered.median().max()+20])
        #gene_type labels
        gene_type_labels = [genetype, genetype2, genetype3]
        # Create the legend
        plt.legend()
        #set x axis length
#         if x_range==False:
#             pass
#         else:
#             length_of_longest_promoter = x_range

        
#         #for all subplots:
#         for n in axes:
#             #remove grids
#             n.grid(False)
#             n.set_xlim([-length_of_longest_promoter,0])
#             leg = n.legend()
#         f.tight_layout()
    #set x axis range if specified    
    if x_range==False:
        pass
    else:
        length_of_longest_promoter = x_range
        

    
    #remove grids
    plt.grid(False)
    #plt.set_xlim([(-length_of_longest_promoter-50),0])
    #set a tight layout
    plt.tight_layout()
    #save figure
    plt.savefig(f'../../data/output/{args.file_names}/rolling_window/{args.foldername_prefix}/plots/{output_prefix}_{variable_of_interest}_{estimator}_sliding_window.pdf', format='pdf')
    #remove plot
    plt.clf()  
    return merged2

def plot_length(df, genetype=False, genetype2=False, genetype3=False):
    """function to plot length distribution of promoter"""
    #make lengths positive by squaring and then square rooting
    df.length = (df.length**2)**(1/2)
    
    if genetype==False:
        dist_plot = df['length']
        #create figure with no transparency
        dist_plot_fig = sns.distplot(dist_plot).get_figure()

    elif genetype2==False:
        sns.distplot(df[df.gene_type == genetype].length, label=genetype)
        
    elif genetype3==False:
        plot = sns.distplot(df[df.gene_type == genetype].length,hist=None,label=genetype)
        sns.distplot(df[df.gene_type == genetype2].length,hist=None,label=genetype2).get_figure()
        plt.legend()
    else:
        plot = sns.distplot(df[df.gene_type == genetype].length,hist=None,label=genetype)
        sns.distplot(df[df.gene_type == genetype2].length,hist=None,label=genetype2)
        sns.distplot(df[df.gene_type == genetype3].length,hist=None,label=genetype3)
        plt.legend()
        
    #remove plot
    plt.clf() 

def add_genetype(coverage,gene_categories):
    """function to add gene type to the df, and remove random genes"""
    select_genes = pd.read_table(gene_categories, sep='\t', header=None)
    cols = ['AGI', 'gene_type']
    select_genes.columns = cols

    merged = pd.merge(select_genes, coverage, on='AGI', how='left')   
   # no_random = merged_renamed[merged_renamed.gene_type != 'random']
  #  no_random.reset_index(drop=True, inplace=True)
    
    return merged

coverage_df = process_input_files(args.bp_covered)
    
#add root chromatin coverage data
coverage_df = add_coverage(coverage_df, args.root_chrom_bp_covered,'root_chrom')
#add shoot chromatin coverage data
coverage_df = add_coverage(coverage_df, args.shoot_chrom_bp_covered,'shoot_chrom')
#add rootshootintersect chromatin coverage data
coverage_df = add_coverage(coverage_df, args.rootshootintersect_chrom_bp_covered,'rootshootintersect_chrom')

#all promoters in genome
merged = windows_coords('all_proms','percentage_bases_covered','% bp covered',coverage_df,args.promoter_bed,args.promoter_no_5UTR, 50,args.EPD_TSS_bed)

#add gene type categories
coverage_prom_types = add_genetype(coverage_df, args.Czechowski_gene_categories)

rolling_rootshootintersect = windows_coords(f'{args.author_name}_genetypenocontrol','percentage_bases_covered','% bp covered',coverage_prom_types,args.promoter_bed,args.promoter_no_5UTR,
                         50,args.EPD_TSS_bed, includeEPDTSS=False,estimator='median',  genetype=args.variable1_name, genetype2=args.variable2_name)

#plot including control variable
rolling_incl_control = windows_coords(f'{args.author_name}_genetype','percentage_bases_covered','% bp covered',coverage_prom_types,args.promoter_bed,args.promoter_no_5UTR, 50,args.EPD_TSS_bed,estimator='median',x_range=1500,  genetype=args.variable1_name, genetype2=args.variable2_name, genetype3='control')