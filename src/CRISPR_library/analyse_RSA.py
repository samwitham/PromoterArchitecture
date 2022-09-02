# %%
#script to read in .csv output files from SmartRoot analysis, concatenate them and then analyse and make plots
#use qpcr conda environment

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import os
import glob 
import sys
#import argparse
import statsmodels.api as sm
#sats annotations
from statannotations.Annotator import Annotator

import statsmodels.formula.api as smf
from bioinfokit.analys import stat
from statsmodels.stats.multicomp import pairwise_tukeyhsd

#estimated marginal mean contrasts
# from rpy2.robjects import Formula, pandas2ri
# import rpy2.robjects as ro  
# from rpy2.robjects.packages import importr
# stats = importr('stats')
# emmeans = importr('emmeans')
# from rpy2.robjects.conversion import localconverter
from patsy import dmatrix
# from pymer4.utils import get_resource_path
# from pymer4.models import Lm
import pingouin as pg


# %%
# #define arguments
# parser = argparse.ArgumentParser(description='Analyse SmartRoot output')
# parser.add_argument('-i', '--input', help='input directory', required=True)
# parser.add_argument('-o', '--output', help='output directory', required=True)





# %%
#function to recursively find all .csv files in a directory and concatenate them into a single dataframe
def concat_csv_recursive(PATH, EXT):
    #find all .csv files in the directory
    csv_files = [file for path, subdir, fname in os.walk(PATH) 
                for file in glob.glob(os.path.join(path, EXT))]
        #glob.glob(f'{directory}/{EXT}', recursive=True)
    #print(csv_files)
    #initialise empty dataframe
    df = pd.DataFrame()
    #loop through all files and concatenate them into a single dataframe
    for file in csv_files:
        df = pd.concat([df, pd.read_csv(file)], ignore_index=True)
    return df

# %%
def sort_data(df,output_location):
    #sort dataframe by sample name
    df = df.sort_values(by=['image'])
    #remove duplicate rows
    df = df.drop_duplicates(keep='first')
    #make nitrate concentration column using image column
    df['nitrate_concentration'] = df['image'].str.split('_').str[1]
    #make sample name column using image column
    df['sample_name'] = df['image'].str.split('_').str[0]
    #make plate column
    df['plate'] = df['sample_name']+'_'+df['image'].str.split('_').str[2]
    #remove spaces from column names
    df.columns = df.columns.str.replace(' ', '')
   
    
    #make several new columns
    #first make new df which will contain one row per plant
    df_plant = df[df.root_order == 0]
    #remove all lines which have no length or which are NaN
    df_plant = df_plant[df_plant.length.notnull()]
    #df_plant = df_plant[df_plant.length != 0]
    # df_plant = df.groupby(['sample_name', 'plate', 'nitrate_concentration',root_ontology]).agg({'image':'count', 'nitrate_concentration':'first', 'sample_name':'first', 'plate':'first'})
    #print(df_plant)
    ## PR = primary root length (cm)
    #change length column to PR
    df_plant['PR'] = df_plant['length']
    # LR = lateral root number (visible from scan)
    #for each root in df_plant, count the number of rows whose parent root in df is the same as the root id in df_plant
    df_plant['LR'] = df_plant.apply(lambda row: df[(df.parent == row.root) & (df.root_order == 1)].shape[0], axis=1)
    #make list of first order lateral root ids
    df_plant['LR_ids'] = df_plant.apply(lambda row: df[(df.parent == row.root) & (df.root_order == 1)].root.tolist(), axis=1)
    #for each id in LR_ids, count the number of rows whose parent root in df is the same as the root id 
    df_plant['LR_2nd_order'] = df_plant.apply(lambda row: df[(df.parent.isin(row.LR_ids)) & (df.root_order == 2)].shape[0], axis=1)
    # LRL = total lateral root length (all LRs added together - cm). Have separate column for 2nd order lateral roots
    
    df_plant['LRL_1st_order'] = df_plant.apply(lambda row: df[(df.parent == row.root) & (df.root_order == 1)].length.sum(), axis=1)
    df_plant['LRL_2nd_order'] = df_plant.apply(lambda row: df[(df.parent.isin(row.LR_ids)) & (df.root_order == 2)].length.sum(), axis=1)
    #add LRL and 2nd order LRL to get total LRL
    df_plant['LRL'] = df_plant['LRL_1st_order'] + df_plant['LRL_2nd_order']
    # ALRL = average lateral root length (LRL/LR - cm)
    df_plant['ALRL'] = (df_plant.LRL) / (df_plant.LR)
   # df_plant['ALRL'] = df_plant.apply(lambda row: (row.LRL / row.LR, axis=1)
    # TRL = total root length (PR + LRL)
    df_plant['TRL'] = df_plant.PR + df_plant.LRL
    # LRD = lateral root density (LR/PR)
    df_plant['LRD'] = df_plant.LR / df_plant.PR
    # LRL_div_TRL = percentage of LRL contributing to TRL (LRL/TRL)
    df_plant['LRL_div_TRL'] = (df_plant.LRL) / df_plant.TRL

    #add genotype column
    df_plant['genotype'] = df_plant['root_name'].str.split('_').str[0]
    #remove spaces from genotype
    df_plant['genotype'] = df_plant['genotype'].str.replace(' ', '')

    #add log columns for PR, LR, LR_2nd_order, LRL, LRL_2nd_order. ALRL, TRL, LRD, LRL_div_TRL
    df_plant['log_PR'] = np.log(df_plant.PR)
    df_plant['log_LR'] = np.log(df_plant.LR)
    df_plant['log_LR_2nd_order'] = np.log(df_plant.LR_2nd_order)
    df_plant['log_LRL'] = np.log(df_plant.LRL)
    df_plant['LRL_1st_order'] = df_plant.LRL_1st_order
    df_plant['log_LRL_2nd_order'] = np.log(df_plant.LRL_2nd_order)
    df_plant['log_ALRL'] = np.log(df_plant.ALRL)
    df_plant['log_TRL'] = np.log(df_plant.TRL)
    df_plant['log_LRD'] = np.log(df_plant.LRD)
    df_plant['log_LRL_div_TRL'] = np.log(df_plant.LRL_div_TRL)
    






    #save df as tsv file
    df_plant.to_csv(f'{output_location}/single_plant_data.tsv', sep='\t', index=False)
    #count number of plants for each plant line

    #partition variation across mutants relative to wild type using principal component analysis of all RSA traits
    #do stats: Using a two-way ANOVA, three phenotypic categories: genotype effects in both nitrogen conditions (genotype-dependent), genotype effects in only one condition (nitrogen-condition-dependent) or genotype by nitrogen condition-dependent effects 
    

    #print(len(df))
    return df, df_plant

# %%
#set matplotlib rc parameters
def set_rc_params():
    #set matplotlib default parameters
    rcParams['xtick.major.width'] = 2
    rcParams['ytick.major.width'] = 2
    rcParams['axes.linewidth'] = 2
    rcParams['lines.linewidth'] = 2
    #remove top and right lines
    rcParams['axes.spines.top'] = False
    rcParams['axes.spines.right'] = False
    #font size
    fontsize = 14
    rcParams['font.size'] = fontsize
    #for getting the microsoft font Arial working, please follow this guide: https://alexanderlabwhoi.github.io/post/2021-03-missingfont/
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    #allow font to be edited later in pdf editor
    #make svg text editable
    rcParams['svg.fonttype'] = 'none'
    rcParams ['pdf.fonttype'] = 42 
    #align y-axis top most tick with end of axis
    rcParams['axes.autolimit_mode'] = 'round_numbers'
    #set margins to ensure any error bars fit
    rcParams['axes.xmargin'] = 0.2
    rcParams['axes.ymargin'] = 0.2
    #define bar width
    #bar_width = 0.65
    #allow math text to be displayed
    #rcParams['mathtext.default'] = 'regular'
    

# %%
def qqplots(df,variables,sample_name, output_location):
        """function to make qq plots"""
        #make qq plots for each variable
        for var in variables:
            #run anova
            #only run if not empty array
            try:

                anova = smf.ols(f'{var} ~ genotype*nitrate_concentration + plate', data=df).fit()
                #make qq  of residuals
                _ = sm.qqplot(anova.resid, line='s')

                
                #save figure
                plt.savefig(f'{output_location}/qqplot_{var}_{sample_name}_residuals.svg',format="svg",
                                    bbox_inches="tight",transparent=True)
                plt.cla()
                plt.close('all')
        
            except ValueError:
                print(f'{sample_name}_{var} is empty array, skipping')
                pass

# %%
def boxplot(df,var,y_label,sample_name,box_pair_p_values, output_location):
    """function to make box plots"""
    #print(df.genotype.unique())
    #print(df.nitrate_concentration.unique())
    #make box plot
    #print(df.genotype.unique())
    #print(box_pair_p_values)
    #filter dict by significance and put in a new dictionary
    # box_pairs_significant = {}
    # for k,v in box_pair_p_values.items():
    #     if v <0.05:
    #         box_pairs_significant[k] = v

    # def convert_pvalue_to_asterisks(pvalue):
    #     if pvalue <= 0.001:
    #         return "***"
    #     elif pvalue <= 0.01:
    #         return "**"
    #     elif pvalue <= 0.05:
    #         return "*"
    #     return "ns"
    
    
    order = ['1mM','10mM']
    fig_args = {'x':'nitrate_concentration', 'y':var,'data':df, 'order':order, 'dodge':True,'hue':'genotype','hue_order':['col0',sample_name]}
    configuration = {'test':None, 'text_format':'star', 'pvalue_thresholds':[[1e-3, "***"],[1e-2, "**"],[0.05, "*"],[1, "ns"]]}#"pairs":list(box_pairs_significant.keys()),"pvalues":list(box_pairs_significant.values()), 'loc':'inside'
    _ = plt.figure(figsize=(5,5))

    fig = sns.boxplot(**fig_args, linewidth=2, palette=["white", "grey"])
    fig = sns.swarmplot(**fig_args, color='black', palette=["black", "black"],size=4)
    #get pairs and pvalues
    #print(box_pairs_significant)
    pairs=list(box_pair_p_values.keys())
    #print(f'pairs={pairs}')
    
    pvalues=list(box_pair_p_values.values())
    #print(f'pvalues={pvalues}')
    # #add statsannotator = Annotator(fig, pairs, **fig_args,verbose=False)
    annotator = Annotator(fig, pairs, **fig_args,verbose=False, show_non_significant=False)#show_non_significant=False will be added in the next version of statsannotator
    #annotator.set_pvalues(pvalues)
    annotator.configure(**configuration)
    
    annotator.set_pvalues_and_annotate(pvalues)

    _ = fig.set(xlabel='$KNO_{3}$ concentration', ylabel=y_label)
    #set y axis limit to start at 0
    _ = plt.ylim(0,None)

    ##plot legend, excluding legend from swarm plot
    h,l = fig.get_legend_handles_labels()
    #change name of label
    l[0] = "Col-0"
    l[1] = sample_name
    #l[2] = "1 mM nitrate"     
    leg = plt.legend(h[0:2],l[0:2],fontsize=14,frameon=False,)#.set_linewidth(2)#,bbox_to_anchor=(0,0.85), loc='best',
    #set linewith of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2)
    # change axes labels
    #_ = plt.ylabel('Relative expression (a.u.)')

    #save plot to file
    plt.savefig(
                    f'{output_location}/{var}_{sample_name}_boxplot.pdf',
                    format="pdf",
                    bbox_inches="tight",transparent=True)
    plt.savefig(
                    f'{output_location}/{var}_{sample_name}_boxplot.svg',
                    format="svg",
                    bbox_inches="tight",transparent=True)
    plt.cla()  # clear axis              
    plt.close('all')   
        


    


# %%
def marginal_effects(df,variables,sample_name, output_location):
        """function to calculate marginal means for interaction genotype*nitrate_concentration"""
        for var, y_label in variables.items():
            #remove NaN values from dataframe
            df = df.dropna(subset=[var]).copy()
            
            #only run if not empty array
            #try:

            #run anova
            anova = smf.ols(f'{var} ~ genotype*nitrate_concentration + plate', data=df).fit()
            #get marginal means, save as txt file            
            #anova.summary_frame().to_csv(f'{output_location}/marginal_means_{var}.tsv', sep='\t')
            #save anova summary to tsv
            #use type 1 anova to test interaction term. If not significant, refit without the interaction term and use Type-II to test the main effects
            table = sm.stats.anova_lm(anova, type=1)
            #print table columns
            #print(f'cols={table.columns}')
            
            #check if interaction term genotype:nitrate_concentration is significant
            if table.loc['genotype:nitrate_concentration']['PR(>F)'] >= 0.05:
                #print(table.loc['genotype:nitrate_concentration']['PR(>F)'])
                #if not significant, refit without the interaction term and use Type-II to test the main effects
                
                
                anova = smf.ols(f'{var} ~ genotype+nitrate_concentration + plate', data=df).fit()
                #filter columns of df
                df_filtered = df[['genotype','nitrate_concentration',var]].copy()
                #drop na
                df_filtered = df_filtered.dropna(subset=['genotype','nitrate_concentration',var]).copy()
                #
                
                # res = stat()
                # # for main effect Genotype, tukey posthocs
                # res.tukey_hsd(df=d_melt, res_var='value', xfac_var='genotype', anova_model=f'{var}~C(genotype)+C(years)+C(Genotype):C(years)')



                # #create copy of df
                # df_copy = df.copy()
                
                #table_type2 = sm.stats.anova_lm(anova, type=2)
                #run tukey posthocs 
                #first make new column with interactions
                #df_copy = df.copy()
                #filter columns
                #df_copy = df_copy[['genotype','nitrate_concentration',var,]]
                
                # df.loc[:,'combination'] = df_copy['genotype'] +'/'+ df_copy['nitrate_concentration']
                # #remove nan
                # clean_df = df.filter(items=[var, 'combination']).dropna()
                # #remove unwanted combinations, only keep var/10mM vs col0/10mM and var/1mM vs col0/1mM
                # clean_df = clean_df[clean_df['combination'].isin(['col0/10mM',f'var/10mM','col0/1mM','var/1mM'])]
                #print(f'cleandf = {clean_df}')

               

                # perform multiple pairwise comparison (Tukey HSD)
                #posthoc = sm.stats.multipletests(table_type2['PR(>F)'], alpha=0.05, method='fdr_bh')


                
                #m_comp = pairwise_tukeyhsd(endog=clean_df[var], groups=clean_df['combination'],alpha=0.05)
                #convert to df
                #m_comp_df = pd.DataFrame(data=m_comp.results_table.data[1:], columns=m_comp.results_table.data[0])
                #m_comp_df = pd.DataFrame()
                #write stats to file
                # with open(f'{output_location}/stats/marginal_means_{var}_{sample_name}.txt', 'w') as f:
                
                #     f.write(f'anova_type_1:\n{table}\ngenotype:nitrate_concentration is not significant so use type 2 anova excluding interaction term\nanova_type_2:\n{table_type2}\n{anova.summary()}\nTukey_post_hocs:\n{m_comp}')


                #print(m_comp)
                #get box pairs and p values for adding stats annotations
                #p_values = pd.DataFrame(data=m_comp._results_table.data[1:] , columns=m_comp._results_table.data[0])
                # print(f'pvalues = {p_values}')
                # print(p_values.columns)
                #filter df columns

                
                genotypes_unique = df['genotype'].unique()
                length_samples = len(genotypes_unique)
                
                box_pair_p_values = {}
                for x in range (0, (length_samples)):                        
                    if genotypes_unique[x] != 'col0':
                        #perform estimated marginal means contrasts between genotype*nitrate interaction
                        #create model, dropping the plate effects because it seems non-significant in most/all cases
                        #test assumption that variances are all equal
                        #split df by nitrate concentration
                        df_10mM = df_filtered[df_filtered['nitrate_concentration'] == '10mM']
                        df_1mM = df_filtered[df_filtered['nitrate_concentration'] == '1mM']
                        var_10mM = pg.homoscedasticity(df_10mM, group='genotype', dv=f'{var}')#Levene test
                        var_1mM = pg.homoscedasticity(df_1mM, group='genotype', dv=f'{var}')

                          #print(f'df_1mM_genotype variance {var} {sample_name}: {var_1mM}')
                       # print(f'df_10mM_genotype variance {var} {sample_name}: {var_10mM}')
                       #followed this https://www.reneshbedre.com/blog/anova.html
                        res = stat()
                        res.anova_stat(df=df_filtered, res_var=var, anova_model=f'{var}~C(genotype)+C(nitrate_concentration)')
                        #print(res.anova_summary)
                        # for main effect genotype
                        # res.tukey_hsd(df=df_filtered, res_var=var, xfac_var='genotype', anova_model=f'{var}~C(genotype)+C(nitrate_concentration)+C(genotype):C(nitrate_concentration)')
                        # genotype_tukey = res.tukey_summary
                        # #for main effect nitrate_concentration
                        # res.tukey_hsd(df=df_filtered, res_var=var, xfac_var='nitrate_concentration', anova_model=f'{var}~C(genotype)+C(nitrate_concentration)')
                        # nitrate_concentration_tukey = res.tukey_summary
                        # for interaction effect between genotype and nitrate_concentration
                        res.tukey_hsd(df=df_filtered, res_var=var, xfac_var=['genotype','nitrate_concentration'], anova_model=f'{var}~C(genotype)+C(nitrate_concentration)')
                        genotype_nitrate_concentration_tukey = res.tukey_summary

                        #get p values
                        p_values = genotype_nitrate_concentration_tukey.copy()
                        #turn group1 and group2 column into string
                        # p_values['group1'] = p_values['group1'].astype(str)
                        # p_values['group2'] = p_values['group2'].astype(str)
                        #print(genotype_nitrate_concentration_tukey['group1'][0])
                        




                        # aov = pg.anova(dv=f'{var}', between=['genotype','nitrate_concentration'], data=df,ss_type=2)


                       

                        #write stats to file
                        with open(f'{output_location}/stats/marginal_means_{var}_{sample_name}.txt', 'w') as f:
                
                            f.write(f'first test assumption that variances are all equal\ndf_1mM_genotype variance {var} {sample_name}:\n {var_1mM}\ndf_10mM_genotype variance {var} {sample_name}:\n {var_10mM}\nanova_type_1:\n{table}\ngenotype:nitrate_concentration is not significant so use type 2 anova excluding interaction term\nanova_type_2:\n{res.anova_summary}\nTukey_post_hocs:\ngenotype_nitrate_concentration_tukey\n{p_values}')
                        

                        #add to box_pair_p_values dictionary the box pair as the key and the p value as the value
                        box_pair_p_values[(('1mM','col0'),('1mM',genotypes_unique[x]))] = p_values.loc[((p_values['group1'] == (f'col0', '1mM')) & (p_values['group2'] == (f'{genotypes_unique[x]}', '1mM')))|((p_values['group1'] == (f'{genotypes_unique[x]}', '1mM')) & (p_values['group2'] == (f'col0', '1mM'))),'p-value'].values[0]
                        box_pair_p_values[(('10mM','col0'),('10mM',genotypes_unique[x]))] = p_values.loc[((p_values['group1'] == (f'col0', '10mM')) & (p_values['group2'] == (f'{genotypes_unique[x]}', '10mM')))|((p_values['group1'] == (f'{genotypes_unique[x]}', '10mM')) & (p_values['group2'] == (f'col0', '10mM'))),'p-value'].values[0]
                        # box_pair_p_values[(('1mM','col0'),('1mM',genotypes_unique[x]))] = p_values.loc[((p_values['group1'] == f'(col0, 1mM)') & (p_values['group2'] == f'({genotypes_unique[x]}, 1mM)'))|((p_values['group1'] == f'({genotypes_unique[x]}, 1mM)') & (p_values['group2'] == f'(col0, 1mM)')),'p-value'].values[0]
                        # box_pair_p_values[(('10mM','col0'),('10mM',genotypes_unique[x]))] = p_values.loc[((p_values['group1'] == f'(col0, 10mM)') & (p_values['group2'] == f'({genotypes_unique[x]}, 10mM)'))|((p_values['group1'] == f'({genotypes_unique[x]}, 10mM)') & (p_values['group2'] == f'(col0, 10mM)')),'p-value'].values[0]
                      

                        #pairwise tests (contrasts between genotype*nitrate interaction)

                        # interaction_contrasts = pg.pairwise_corr(dv=var, between='genotype',within='nitrate_concentration', data=df, padjust='fdr_bh', effsize='hedges', correction = 'auto', alpha =0.05, parametric=True, marginal=True,interaction=False, subject='subject') #Benjamini/Hochberg FDR correction, Hedges g effect size method, correction = 'auto' (it will automatically use Welch T-test when the sample sizes are unequal)
                        # print(interaction_contrasts)



                        #I used https://glennwilliams.me/blog/2021/09/07/estimating-marginal-means-and-pairwise-tests-by-hand-in-python/ for help with the code

                        # model = smf.ols(f'{var} ~ genotype*nitrate_concentration', data=df_copy).fit()
                        # #First, we need to set up a grid allowing us to see the unique combinations for the levels of each factor.
                        # grid = np.array(np.meshgrid(
                        #     df_copy["genotype"].unique(),
                        #     df_copy["nitrate_concentration"].unique()
                        #     )).reshape(2, 4).T
                        # #make into a dataframe
                        # grid = pd.DataFrame(grid, columns=['genotype', 'nitrate_concentration'])
                        # #print(grid)
                        # mat = dmatrix(
                        #     "C(genotype, Treatment(1))*C(nitrate_concentration, Treatment(1))", 
                        #     grid, 
                        #     return_type = "matrix"
                        #     )
                        # print(mat)
                        # #Now we have a design matrix we can get the betas from our fixed effects like so.
                        # # grab beta coefficients (means)
                        # #print(model)
                        # betas = model.params
                        # #print(betas)
                        # #calculate the estimated marginal means
                        # emmeans = grid
                        # emmeans["means"] = np.dot(mat, betas)
                        # #get the standard errors for our estimated marginal means by first getting the variance covariance matrix using the cov_params() 
                        # #reduce the variance covariance matrix to exclude our random effects. We do that here by filtering out any terms including the strings “Var” or “Cor” which by default statsmodels uses to indicate random effects terms.
                        # #in this case there are no random effects so nothing changes
                        # vcov = model.cov_params()
                        # print(vcov)
                        # # vcov2= vcov[~vcov.index.str.contains('Var|Cor')]
                        # # vcov3 = vcov2.loc[:,~vcov2.columns.str.contains('Var|Cor')]
                        # # if vcov.equals(vcov3):
                        # #     print('vcov is the same as vcov3')
                        # #Next add standard errors to the means by getting the square root of the diagonal for the design matrix multiplied by the variance covariance matrix multiplied by the the transpose of the design matrix.
                        # emmeans['SE'] = np.sqrt(np.diagonal(mat @ vcov) @ mat.T)
                        # # emmeans["se"] = np.sqrt(np.diag(np.dot(np.dot(mat.T, vcov), mat)))
                        # # emmeans["se"] = np.sqrt(np.dot(np.diagonal(np.dot(mat, vcov)), mat.T))
                        # print(emmeans)



                        
                            
                            
                            






                        #make r dataframe
                        # with localconverter(ro.default_converter + pandas2ri.converter):
                        #     r_df = ro.conversion.py2rpy(df_copy)

                        
                        # #make formula in R language
                        # fmla = Formula(f'{var} ~ genotype+nitrate_concentration')
                        # env = fmla.environment
                        # env[var] = df_copy[var]
                        # env['genotype'] = df_copy['genotype']
                        # env['nitrate_concentration'] = df_copy['nitrate_concentration']
                                                
                        # model = stats.lm(formula=fmla, data=r_df)
                        # #print(model)
                        # fmla2 = Formula('pairwise ~ genotype*nitrate_concentration')
                        # contrast(Simple.Effects.By.Type,Set1,adjust='none')
                        # interaction_emmeans = emmeans(data=model, formula=fmla2)
                        # print(interaction_emmeans)
                        # # #run posthocs
                        # # sp.posthoc_ttest(df, val_col='genotype', group_col='nitrate_concentration', p_adjust='holm')
                        # interaction_emmeans = emmeans(model, 'pairwise ~ Genotype+NO3_Level')
                        # print(interaction_emmeans)
                        # Compute post-hoc tests - Compare each level of IV3 to each other level of IV3, within each level of IV4. Use default Tukey HSD p-values.
                        # marginal_estimates, comparisons = model.post_hoc(
                        # marginal_vars="genotype", grouping_vars="nitrate_concentration"
                        # )
                        # print(marginal_estimates)
                        



                        
                
                
                
                
            if table.loc['genotype:nitrate_concentration']['PR(>F)'] < 0.05:
                #if significant interaction effect, analyse nitrate concentrations separately using one-way ANOVA
                #first split dataframe into separate dataframes for each nitrate concentration
                print(f'{var}{sample_name} is significant for genotype*nitrate_concentration interaction')
                df_low = df[df.nitrate_concentration == '1mM']
                df_high = df[df.nitrate_concentration == '10mM']
                anova_low_nitrate = smf.ols(f'{var} ~ genotype + plate', data=df_low).fit()
                anova_high_nitrate = smf.ols(f'{var} ~ genotype + plate', data=df_high).fit()
                table_low = sm.stats.anova_lm(anova_low_nitrate, type=2)
                table_high = sm.stats.anova_lm(anova_high_nitrate, type=2)
                #print(table_high)
                #write stats to file
                with open(f'{output_location}/stats/marginal_means_{var}_{sample_name}.txt', 'w') as f:
                
                    f.write(f'anova_type_1:\n{table}\ngenotype*nitrate_concentration is significant so analyse each nitrate concentration separately\nanova_type_2_1mM_nitrate:\n{table_low}\nanova_type_2_10mM_nitrate:\n{table_high} \n{anova.summary()}')

                

            
                #get p values
                p_value_low_nitrate_df = pd.DataFrame(data=table_low)
                p_value_high_nitrate_df = pd.DataFrame(data=table_high)
                
                #get box pairs and p values for adding stats annotations
                genotypes_unique = df['genotype'].unique()
                length_samples = len(genotypes_unique)
                box_pair_p_values = {}
                for x in range (0, (length_samples)):                        
                    if genotypes_unique[x] != 'col0':                            
                        box_pair_p_values[(('1mM','col0'),('1mM',genotypes_unique[x]))] = p_value_low_nitrate_df.loc['genotype','PR(>F)']
                        box_pair_p_values[(('10mM','col0'),('10mM',genotypes_unique[x]))] = p_value_high_nitrate_df.loc['genotype','PR(>F)']
                        

            #PR(>F)
            #make boxplots
            #remove all string before the first underscore in the variable name, and return all subsequent string     

            #split var string on _
            no_log_var = var.split('_')[1:]
            no_log_var = '_'.join(no_log_var)
            #print(no_log_var)
            #make boxplots
            #first filter df
            boxplot_df = df.filter(items=[no_log_var, 'nitrate_concentration','genotype']).dropna().copy()
            #get column types
            #print(boxplot_df.dtypes)
            
            boxplot(boxplot_df,no_log_var,y_label,sample_name,box_pair_p_values, f'{output_location}/boxplots')

            # except ValueError:
            #     print(f'{sample_name}_{var} is empty array, skipping')
            #     print(df[df.sample_name == sample_name][var])
            #     pass

# %%
#function to analyse data and make plots
def analyse_data(df_plant,sample_name, output_location):
    """function to run anovas and make boxplots"""
    #anova_PR <- lm(logPR ~ Genotype*NO3_Level + Plate, data = Roots1)
    #change -inf values to NaN using .loc
    df_plant.loc[df_plant['log_PR'] == -np.inf, 'log_PR'] = np.nan
    df_plant.loc[df_plant['log_LR'] == -np.inf, 'log_LR'] = np.nan
    df_plant.loc[df_plant['log_LRL'] == -np.inf, 'log_LRL'] = np.nan
    df_plant.loc[df_plant['log_ALRL'] == -np.inf, 'log_ALRL'] = np.nan
    df_plant.loc[df_plant['log_TRL'] == -np.inf, 'log_TRL'] = np.nan
    df_plant.loc[df_plant['log_LRD'] == -np.inf, 'log_LRD'] = np.nan
    df_plant.loc[df_plant['log_LRL_div_TRL'] == -np.inf, 'log_LRL_div_TRL'] = np.nan
    df_plant.loc[df_plant['log_LR_2nd_order'] == -np.inf, 'log_LR_2nd_order'] = np.nan
    df_plant.loc[df_plant['log_LRL_2nd_order'] == -np.inf, 'log_LRL_2nd_order'] = np.nan
   

    # anova_PR = smf.ols('PR ~ genotype*nitrate_concentration + plate', data=df_plant).fit()
    #check anova assumptions
    # print(anova_PR.summary())
    # fig = sm.qqplot(anova_PR.resid, line='s')
    #save figure
    #make directory for the plots to be exported to
    output_dir = f'{output_location}/qqplots'
    
        
    # fig.savefig(f'{output_location}/qqplots/qqplot_PR.png')
    #log_PR residuals look mainly normal from the qqplot, (points at the extreme ends can be discounted)
    variables = ['PR','log_PR','LR','log_LR','LRL','log_LRL','ALRL','log_ALRL','TRL','log_TRL','LRD','log_LRD','LRL_div_TRL','log_LRL_div_TRL','LR_2nd_order','log_LR_2nd_order','LRL_2nd_order','log_LRL_2nd_order']
    #variables_logs = ['log_PR','log_LR','log_LRL','log_ALRL','log_TRL','log_LRD','log_LRL_div_TRL','log_LR_2nd_order','log_LRL_2nd_order']
    variables_logs_dict = {'log_PR':'Primary root length (cm)','log_LR':'Number of lateral roots','log_LRL':'Total lateral root length (cm)','log_ALRL':'Average lateral root length (cm)','log_TRL':'Total root length (cm)','log_LRD':'Lateral root density','log_LRL_div_TRL':'Ratio of lateral root length to\ntotal root length (LRL/TRL)',}#'log_LR_2nd_order':'Number of second order lateral roots','log_LRL_2nd_order':'Second order lateral root length (cm)'
    qqplots(df_plant,variables,sample_name, output_dir)
    #I will only use log transformed data
    #run anovas and calculate marginal effects for interaction genotype*nitrate_concentration
    
    marginal_effects(df_plant,variables_logs_dict, sample_name, output_location)
    return df_plant



    # ANOVA table using bioinfokit v1.0.3 or later (it uses wrapper script for anova_lm)

    # res = stat()
    # res.anova_stat(df=df_plant, res_var='PR', anova_model='PR ~ genotype*nitrate_concentration + plate')
    # res.anova_summary
    # #generate QQ-plot from standardized residuals
    # # res.anova_std_residuals are standardized residuals obtained from ANOVA (check above)
    # # sm.qqplot(res.anova_std_residuals, line='45')
    # # plt.xlabel("Theoretical Quantiles")
    # # plt.ylabel("Standardized Residuals")
    # # plt.show()
    # res.qq_plot(df=df_plant, res_var='PR', anova_model='PR ~ genotype*nitrate_concentration + plate')
        





# %%
#main function
def main(args):
    #read in arguments
    #input_dir = args.input
    input_dir = f'../../data/CRISPR_library/images/rsa_output'
    #output_dir = args.output
    output_dir = f'../../data/CRISPR_library'
    #make directory for the plots to be exported to
    output_dir = f'{output_dir}/smartroot_plots'
    try:
        # Create target Directory
        os.mkdir(output_dir)
        print("Directory " , output_dir ,  " created") 
    except FileExistsError:
        print("Directory " , output_dir ,  " already exists")



    output_dir2 = f'{output_dir}/stats'
    try:
        # Create target Directory
        os.mkdir(output_dir2)
        print("Directory " , output_dir2 ,  " created") 
    except FileExistsError:
        print("Directory " , output_dir2 ,  " already exists")
    output_dir2 = f'{output_dir}/boxplots'
    try:
        # Create target Directory
        os.mkdir(output_dir2)
        print("Directory " , output_dir2 ,  " created") 
    except FileExistsError:
        print("Directory " , output_dir2 ,  " already exists")

    output_dir2 = f'{output_dir}/qqplots'
    try:
        # Create target Directory
        os.mkdir(output_dir2)
        print("Directory " , output_dir2 ,  " created") 
    except FileExistsError:
        print("Directory " , output_dir2 ,  " already exists")

    #read in and concatenate .csv files
    df = concat_csv_recursive(input_dir, '*.csv')
    #print(df.head())

    
    #sort data
    df,df_plant = sort_data(df,output_dir)
    #analyse dataframe and make plots
    #analyse_data(output_dir)
    #set matplotlib rc parameters
    set_rc_params()
    #first split into separate dataframes for each sample_name
    #then analyse each dataframe and make plots


    
    
    for sample_name in df_plant['sample_name'].unique():
        #get dataframe for each sample_name
        df_sample = df_plant[df_plant['sample_name'] == sample_name].copy()
        #analyse dataframe and make plots
        df_plants = analyse_data(df_sample,sample_name,output_dir)





    #analyse_data(df_plant, output_dir)
    
    #save dataframe to csv file
    df.to_csv(f'{output_dir}/all_smartroot_data.csv', index=False)

# %%
if __name__ == '__main__':
    main(sys.argv)


