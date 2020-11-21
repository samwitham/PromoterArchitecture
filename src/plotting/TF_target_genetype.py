import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from pingouin import kruskal
import scikit_posthocs as sp
import os
from statannot import add_stat_annotation
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='TF_target_genetype')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('tf_gene_categories', type=str, help='Input location of transcription factor gene categories')
parser.add_argument('all_genes_ranked', type=str, help='Input location of all ranked genes')
parser.add_argument('promoter_TF', type=str, help='Input location of mapped TF-promoter target file')
parser.add_argument('dependent_variable', type=str, help='Dependent variable name [tau or expression_cv]')
parser.add_argument('output_folder_name', type=str, help='Output folder name')
parser.add_argument('variable1_name', type=str, help='Optional replacement name for 1st variable eg. non-specific',default = 'constitutive', nargs="?")
parser.add_argument('variable2_name', type=str, help='Optional replacement name for 2nd variable eg. tissue&condition_specific',default = 'variable', nargs="?")
parser.add_argument('author_name', type=str, help='Optional replacement name for author in reference to the geneset',default = 'Czechowski', nargs="?")
parser.add_argument('palette', type=str, help='Optional replacement colour palette for plots',default = None, nargs="?")

args = parser.parse_args()

#make directory for the plots to be exported to
dirName = f'../../data/output/{args.file_names}/{args.output_folder_name}'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
#make directory for the plots to be exported to
dirName = f'../../data/output/{args.file_names}/{args.output_folder_name}plots'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")
    
def read_files(tf_gene_categories,promoter_TF):
    """read in the files as dfs"""    
    #all_gene_categories_dfs = pd.read_table(all_gene_categories, sep='\t',header=0)
    promoter_TF_df = pd.read_table(promoter_TF,sep='\t',header=0)
    tf_gene_categories_df = pd.read_table(tf_gene_categories, sep='\t',header=None)
    cols = ['TF_AGI','gene_type']
    tf_gene_categories_df.columns = cols
    #merge the dfs
    merged = pd.merge(tf_gene_categories_df, promoter_TF_df, on = 'TF_AGI', how='left')
    #remove NaN
    if 'TF_AGI_allgenes' in merged.columns:
        merged_filtered = merged[merged.TF_AGI_allgenes.notna()]
    elif 'AGI code' in merged.columns:
        merged_filtered = merged[merged['AGI code'].notna()]
        
    
    

    return tf_gene_categories_df,promoter_TF_df,merged_filtered

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

def dunn_posthoc_test(df,dependent_variable, between):
    """dunn_posthoc tests with bonferroni multiple correction"""
    return sp.posthoc_dunn(df, val_col=dependent_variable, group_col=between, p_adjust='bonferroni')

def calculate_mean_CV(df,dependent_variable):
    """calculate the mean coefficient of variation of the promoters the TF binds to"""
    #group by TF and calculate mean for each TF
    means = df.groupby('TF_AGI')[dependent_variable].mean()
    #turn into a dataframe
    means_df = pd.DataFrame(means)
    #turn the index into a new column
    means_df.reset_index(level=0, inplace=True)
    #name columns
    cols = ['TF_AGI', f'mean_{dependent_variable}']
    means_df.columns = cols
    
        
    #group by promoter and calculate SD (standard deviation) for each promoter
    #sd = df.groupby('promoter_AGI')['expression_CV'].std()
    #turn into a dataframe
    #sd_df = pd.DataFrame(sd)
    #turn the index into a new column
    #sd_df.reset_index(level=0, inplace=True)
    #name columns
    #cols = ['promoter_AGI', 'sd']
    #sd_df.columns = cols
    
    #merge the dfs
    #merged = pd.merge(means_df,sd_df)
    return means_df

def merge_genetype(df, gene_categories):
    """merge df with gene_categories file adding the genetype of the promoters (if in top 100 constitutive or top 100 variable promoters)"""
    gene_cats = pd.read_table(gene_categories, sep='\t', header=None)
    cols = ['TF_AGI','gene_type']
    gene_cats.columns = cols
    merged = pd.merge(gene_cats,df, on='TF_AGI', how='left')
    #drop NaN
    merged_filtered = merged.dropna()
    #reset index
    merged_filtered_index = merged_filtered.reset_index(drop=True)
    
    return merged_filtered_index

def make_plot(df,x_variable, y_variable,x_label, y_label, output_prefix, plot_kind,palette,mean=False):
    """function to make and save plot"""
    #allow colour codes in seaborn
    sns.set(color_codes=True)
    sns.set_style("whitegrid")
    #plot
    x=x_variable
    y=y_variable
    order=[args.variable1_name, args.variable2_name, "control"]
    
    #set colour palette
    colours = sns.color_palette(palette)
    
   #make copy of df
    if mean == False:
        merged2 = df.copy()
   
        merged2_unique = merged2.drop_duplicates(['promoter_AGI'],keep='last')
    elif mean == True:
        merged2_unique = df.copy()
    #make sample sizes equal for comparison
    # identify sample size of the minimum category
    minimum_sample_size = merged2_unique.gene_type.value_counts().min()
    # print this
    print(f'sample size in each category = {minimum_sample_size}')
    #save sample size as file
    with open(f'../../data/output/{args.file_names}/{args.output_folder_name}plots/number_of_genes_in_each_category_{args.dependent_variable}.txt','w') as file:
        file.write('number_of_genes_in_each_category='+str(minimum_sample_size))
    
    
    # multiply this by the number of categories
    total_sample_size = minimum_sample_size * len(merged2_unique.gene_type.unique())
    #select equal sample sizes of each category with a random state of 1 so it's reproducible
    equal_samplesizes = rep_sample(merged2_unique, 'gene_type',total_sample_size,random_state = 1)
    # now filter out genes which were not selected using the minimum sample size
    to_remove = merged2_unique[~merged2_unique.TF_AGI.isin(equal_samplesizes.TF_AGI)]
    df = df[~df.TF_AGI.isin(to_remove.TF_AGI)]
    
    #if violin plot don't extend past datapoints
    if plot_kind == 'violin': 
        plot = sns.catplot(x=x, y=y, data=df, kind=plot_kind,order=order,palette=colours, cut=0)
    else:
        plot = sns.catplot(x=x, y=y, data=df, kind=plot_kind,order=order,palette=colours)
    #plot points
    ax = sns.swarmplot(x=x, y=y, data=df, color=".25",order=order)
    #add significance if necessary - dunn's posthocs with multiple Bonferroni correction
    stat = dunn_posthoc_test(df,y_variable,x_variable)
    #label box pairs
    box_pairs=[(args.variable1_name, args.variable2_name),(args.variable1_name, "control"),(args.variable2_name, "control")]
    #make empty list of p_values
    p_values = []
    #populate the list of p_values according to the box_pairs
    for pair in box_pairs:
        print(pair)
        #select p value for each pair
        p = stat.loc[pair[0],pair[1]]
        p_values.append(p)


    
    #add stats annotation to the plot
    test_results = add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                                      box_pairs=box_pairs,
                                      text_format='star',
                                      loc='outside',verbose=2,
                                      perform_stat_test=False,
                                       pvalues=p_values, test_short_name='Dunn')
    
    #change axes labels
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    #tight layout
    plt.tight_layout()  
    #save figure
    ax.get_figure().savefig(f'../../data/output/{args.file_names}/{args.output_folder_name}plots/{output_prefix}_{plot_kind}.pdf', format='pdf')      
    
tf_gene_categories_dfs,promoter_TF_df,merged = read_files(args.tf_gene_categories,args.promoter_TF)

#plot the CV for each promoter gene_type - whole promoter individual TF CVs
#make_plot(merged,'gene_type', 'expression_CV','Gene type', 'Promoter target expression CV', 'Czechowski_CV', 'box',palette)

#calculate mean CV of targets per TF
Czechowski_means = calculate_mean_CV(merged,args.dependent_variable)

#add TF gene_type categories
Czechowski_means_genetype = merge_genetype(Czechowski_means, args.tf_gene_categories)

#plot the mean CV for each promoter gene_type - whole promoter mean TF CVs
make_plot(Czechowski_means_genetype,'gene_type', f'mean_{args.dependent_variable}','Gene type', f'Promoter target mean {args.dependent_variable}', f'{args.author_name}_{args.dependent_variable}_mean', 'box',args.palette,mean=True)