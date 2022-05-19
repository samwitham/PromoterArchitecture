import pandas as pd
import numpy as np
import io
from pybedtools import BedTool
from collections import defaultdict
import math
#chunks
from more_itertools import sliced
from itertools import combinations
#parse arguments
import argparse
import shutil
import glob


def parse_args(args):
    parser = argparse.ArgumentParser(description="pacbio_analyse_overlapping_TFBSs")

    parser.add_argument(
        "input_folder",
        type=str,
        help="Input location of tsv part files from pacbio_analyse_overlapping_TFBSs_part1.py",
    )
    parser.add_argument(
        "guide_pairs",
        type=str,
        help="Input location of guide_pairs .csv (only first 2 cols will be used)",
    )
    parser.add_argument(
        "gene_name",
        type=str,
        help="Name of gene being analysed",
    )
    parser.add_argument(
        "output_folder",
        type=str,
        help="Output folder location ending with a forward slash",
    )
    return parser.parse_args(
        args
    )  # let argparse grab args from sys.argv itself to allow for testing in module import


def read_in_files(input_folder,guide_pairs,gene):
    """read in the files"""
    #read in guide pairs file
    guide_pairs_df = pd.read_csv(guide_pairs,header=0)
    #only keep first 2 columns
    guide_cols = ['guide1','guide2']
    guide_pairs_df = guide_pairs_df[guide_cols]

    mutations_df_cols = ['chr','plant_ID','platename','library','first_reaction_primers','second_reaction_primers','guide','guide_number','aligned_sequence','reference_sequence','mutation_type','read_number','read_percentage','insertion_positions','deletion_positions','substitution_positions',
    'insertion_cut_site_distance','deletion_cut_site_distance','substitution_cut_site_distance','cut_site_promoter_position','insertion_positions_relative_to_TSS'
    ,'insertion_genomic_positions','deletion_positions_relative_to_TSS','deletion_genomic_positions','substitution_positions_relative_to_TSS','substitution_genomic_positions','insertion_overlapping_TFBS_family','insertion_overlapping_TFBS_AGI','deletion_overlapping_TFBS_family','deletion_overlapping_TFBS_AGI','substitution_overlapping_TFBS_family','substitution_overlapping_TFBS_AGI']




    # #now read in the mutations file using pandas
    mutations_df = pd.read_table(f'{input_folder}{gene}_TFBSoverlapping.tsv', sep='\t', header = None)
    #add columns
    mutations_df.columns = mutations_df_cols
    # #remove duplicates since multiple columns of same name were added
    # mutations_df = mutations_df.drop_duplicates(keep='first')
    # #save file
    mutations_df.to_csv(f'{input_folder}{gene}_TFBSoverlapping.tsv', sep="\t", index=False, header=True)
    
    return mutations_df,guide_pairs_df




# #written in dask format rather than pandas
def genotype_plant_lines(mutations_df,output_folder,gene):
    """function to decide whether each plant line is -homozygous 
    -biallelic - each mutated at same location/site twice
    -chimeric - different mutations in different cells or tissues - dont analyse if chimeric - 1 small region, if more than 2 mutations then probably chimeric - ie.    probably still has tDNA
    Multiple guide sites - eg. multiple agro strains in each plant. Check for paired guides and whether guides that aren’t meant to be paired are paired.
    """
    #homozygouslines are those with no duplicated guides in each plant line
    mutations_df.loc[~mutations_df.duplicated(['plant_ID','guide']),'genotype'] = 'homozygous'
    #if 70% of reads or more are the same then mark as homozygous
    mutations_df.loc[mutations_df.read_percentage >= 70, 'genotype'] = 'homozygous'
    #if 10% of reads or less are that mutation, mark genotype as 'nan
    mutations_df.loc[mutations_df.read_percentage <=10, 'genotype'] = 'nan'
    #if between 10 and 70% of reads, mark genotype as heterozygous
    mutations_df.loc[(mutations_df.read_percentage > 10) & (mutations_df.read_percentage < 70), 'genotype'] = 'heterozygous'

    #create a count column
    mutations_df['number_of_different_alleles'] = int()
    #create non wild type count column
    mutations_df['number_of_different_non_wt_alleles'] = int()
    
    #create a sum_of_count column
    #count how many duplicates there are for heterozygous lines.
    mutations_df.loc[mutations_df.genotype == 'heterozygous','number_of_different_alleles'] = mutations_df[mutations_df.genotype == 'heterozygous'].groupby(['plant_ID','guide'])['number_of_different_alleles'].transform('count')
    #make non wild type count
    mutations_df.loc[(mutations_df.genotype == 'heterozygous')&~(mutations_df.mutation_type == 'None'),'number_of_different_non_wt_alleles'] = mutations_df[(mutations_df.genotype == 'heterozygous')&~(mutations_df.mutation_type == 'None')].groupby(['plant_ID','guide'])['number_of_different_non_wt_alleles'].transform('count')

    #if heterozygous count is 1 then homozygous
    mutations_df.loc[(mutations_df.genotype == 'heterozygous') & (mutations_df.number_of_different_alleles == 1),'genotype'] = 'homozygous'

    #if homozygous then count is 1
    mutations_df.loc[mutations_df.genotype == 'homozygous','number_of_different_alleles'] = 1

    #if number_of_different_alleles is two and no wild type is present in either of the groups of reads then biallelic
    mutations_df.loc[(mutations_df['number_of_different_non_wt_alleles']==2)&~(mutations_df.mutation_type == 'None'),'genotype'] = 'biallelic'
    
    #if number_of_different_alleles more than two then chimeric
    mutations_df.loc[mutations_df['number_of_different_alleles']>2,'genotype'] = 'chimeric'

    #change column order
    mutations_df = mutations_df[['chr', 'plant_ID', 'platename', 'library', 'first_reaction_primers',
       'second_reaction_primers', 'guide', 'guide_number', 'aligned_sequence',
       'reference_sequence', 'mutation_type','genotype', 'read_number', 'read_percentage',
       'insertion_positions', 'deletion_positions', 'substitution_positions',
       'insertion_cut_site_distance', 'deletion_cut_site_distance',
       'substitution_cut_site_distance', 'cut_site_promoter_position',
       'insertion_positions_relative_to_TSS', 'insertion_genomic_positions',
       'deletion_positions_relative_to_TSS', 'deletion_genomic_positions',
       'substitution_positions_relative_to_TSS',
       'substitution_genomic_positions', 'insertion_overlapping_TFBS_family',
       'insertion_overlapping_TFBS_AGI', 'deletion_overlapping_TFBS_family',
       'deletion_overlapping_TFBS_AGI', 'substitution_overlapping_TFBS_family',
       'substitution_overlapping_TFBS_AGI','number_of_different_alleles','number_of_different_non_wt_alleles'
       ]]
    #remove genotype 'nan'
    mutations_df = mutations_df[~(mutations_df.genotype == 'nan')]
    
    #save df
    mutations_df.to_csv(f'{output_folder}{gene}/{gene}_TFBSoverlapping_genotyped.tsv', sep="\t", index=False, header=1)

    return mutations_df
    
def categorise_guide_pairs(mutations_df_genotyped, guide_pairs_df, output_folder,gene):
    """function to check which guide pairs where delivered to each plant line and to check whether mutations at more than one guide site are within the other guide pair"""
    #first get unique guides based on 2 columns in guide_pairs_df    
    unique_guides = pd.concat([guide_pairs_df['guide1'],guide_pairs_df['guide2']]).unique()
    #for guides in unique guides list, add dictionary value of all other potential guides that it is paired with
    #create defaultdict with lists as values so that non-existing keys can be added to in one go
    guide_dict = defaultdict(list)
    for guide in unique_guides:
        #check for instances of that guide in the first column
        filtered_col_1 = guide_pairs_df[guide_pairs_df.guide1==guide]['guide2'].to_list()
        #then do the same for the second column
        filtered_col_2 = guide_pairs_df[guide_pairs_df.guide2==guide]['guide1'].to_list()
        #create list of unique values
        list_of_unique = np.unique((filtered_col_1 + filtered_col_2)).astype(list)
        #append each value in list to dict
        for val in list_of_unique:
            guide_dict[guide] += [val]
        
    #turn into normal dict
    guide_dict = dict(guide_dict)

    #mutations_df_genotyped df, filter out non mutated guides. Then create a dictionary of guides which were mutated. Then compare the guide pairs in guide_dict to the mutated guides.
    #make a new df with one plant line per row, detailing how many guides sites were mutated, which guide sites were mutated and then if more than 1 guide site is mutated
    #finally, merge this df with the mutations_df_genotyped (applying the same value in each column within each plant line)
    #filter out non mutated guides
    only_mutated_guides = mutations_df_genotyped[~(mutations_df_genotyped.mutation_type=='None')]
    #create a dictionary of guides which were mutated with plant ID as the key
    #create dictionary of the only_mutated_guides df
    only_mutated_guides_dict = only_mutated_guides.to_dict(orient="records")
    #create another dict with plant IDs as the key
    plant_ID_guide_dict = defaultdict(list)
    #iterate over only_mutated_guides_dict adding plant IDs as keys and mutated guide sites as values
    for row_dict in only_mutated_guides_dict:
        plant_ID = row_dict['plant_ID']       
        guide = row_dict['guide']
        #add to plant_ID_guide_dict if not there already
        plant_ID_guide_dict[plant_ID] += [guide]

        
    #create a new df with only one row per plant ID and only containing plant IDs that aren't wild type
    one_row_each_no_wildtype = only_mutated_guides.plant_ID.drop_duplicates(keep='first').to_frame()

    #for each plant line check whether the mutated guides fall within a guide pair
    #first turn into dict
    one_row_each_no_wildtype_dict = one_row_each_no_wildtype.to_dict(orient="records")
    #print(one_row_each_no_wildtype_dict)
    #turn into defaultdict
    one_row_each_no_wildtype_ddict = []

    #iterate through dict add more values/columns names
    for row_dict in one_row_each_no_wildtype_dict:
        #create default dict of row
        dd = defaultdict(list,row_dict)

        plant_ID = row_dict['plant_ID']
        plant_ID_guides = plant_ID_guide_dict[plant_ID]
        dd['mutated_guide_sites'] += list(np.unique(plant_ID_guides))
        #add list of guides to dict

        plant_ID_guide_pair_dict = defaultdict(list)
        
        for guide in dd['mutated_guide_sites']:
            #create temporary list 
            temp_list = []
            #get approved guide pairs
            guide_pairs = guide_dict[guide]
            #print(guide_pairs)
            for g in guide_pairs:
                temp_list.append((guide,g))
            
            #append to approved guide pairs            
            plant_ID_guide_pair_dict['guide_pairs'] += set(temp_list)

        #check for approved guide pairs
        guidepair_combos = list(combinations(dd['mutated_guide_sites'],2))
        #print(plant_ID_guide_pair_dict['guide_pairs'])
        guidepair_intersection =  list(set(guidepair_combos).intersection(plant_ID_guide_pair_dict['guide_pairs']))
        non_approved_guidepairs = list(set(guidepair_combos).difference(plant_ID_guide_pair_dict['guide_pairs']))
        #find unique guides in non_approved_guidepairs
        unique_guidepair_intersection = list(np.unique(guidepair_intersection))
        unique_nonapproved_guides = list(np.unique(non_approved_guidepairs))
        #find guides in nonapproved guides list that aren't approved
        list_of_nonapproved = list(set(unique_nonapproved_guides).difference(unique_guidepair_intersection))
        #if guide in non_approved guidepairs is not in guidepair_intersection
        #append to ddict
        dd['approved_guide_pairs'] = guidepair_intersection
        dd['guides_not_in_pairs'] = list_of_nonapproved
    
        #estimate agro strain number
        #round up the division using math.ceil
        #first create list of mutated_guide_sites which aren't in approved pairs.
        another_list = []
        for pair in dd['approved_guide_pairs']:
            another_list.append(pair)
        another_unique_list = list(np.unique(another_list))
        temp_set = list(set(dd['mutated_guide_sites']).difference(another_unique_list))
        #divide mutated_guide_sites which aren't in approved pairs by 2, rounding up, and then add 1 for every approved guide pair
        number_of_strains = int(math.ceil(len(temp_set)/2)) + len(dd['approved_guide_pairs'])
        dd['estimated_agro_strain_number'] = number_of_strains
        #append defaultdict to list
        one_row_each_no_wildtype_ddict.append(dd)
    #turn ddict into df
    df_single_rows = pd.DataFrame.from_dict(one_row_each_no_wildtype_ddict)
    #merge this df with the only_mutated_guides df
    only_mutated_guides = pd.merge(only_mutated_guides,df_single_rows, on='plant_ID', how='left')
    #save this df
    only_mutated_guides.to_csv(f'{output_folder}{gene}/{gene}_TFBSoverlapping_genotyped_only_mutated.tsv', sep="\t", index=False, header=1)


def main(args):
    # parse arguments
    args = parse_args(args)

    #read in files
    mutations_df,guide_pairs_df =  read_in_files(args.input_folder,args.guide_pairs,args.gene_name)
    #read_in_files(args.input_folder,args.guide_pairs,args.gene_name)

    # #genotype the plant lines, producing an output tsv
    mutations_df_genotyped = genotype_plant_lines(mutations_df,args.output_folder,args.gene_name)

    #produce an output tsv with one row per plant line with all guide site mutations found in that line
    categorise_guide_pairs(mutations_df_genotyped,guide_pairs_df,args.output_folder,args.gene_name)

if __name__ == "__main__":
    import sys

    main(sys.argv[1:])