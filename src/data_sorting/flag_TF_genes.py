import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='flag_TF_genes')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('gene_categories', type=str, help='Input location of Czechowski gene categories text file')
parser.add_argument('TF_list', type=str, help='Input location of Arabidopsis transcription factor list')
parser.add_argument('output', type=str, help='Output location of flagged TF genes')
args = parser.parse_args()

def flagTFs(genes_of_interest, TF_list, output):
    """Function to export a file with all genes which are TFs from a list of genes of interest"""
    all_genes = pd.read_table(genes_of_interest, sep='\t', header=None, names=['Gene_ID','gene_type'])
    TFs = pd.read_table(TF_list, sep='\t', header=1)
    #merge the dfs so that only TFs are included    
    TFs_of_interest = pd.merge(all_genes, TFs, on='Gene_ID')
    #select only columns of interest

    TFs_of_interest = TFs_of_interest[['Gene_ID', 'gene_type','Family']]
    #drop duplicates
    TFs_of_interest = TFs_of_interest.drop_duplicates('Gene_ID')
    #make output file containing only genes from genes of interest which are transcription factors
    TFs_of_interest.to_csv(output, sep='\t', index=False)
    
flagTFs(args.gene_categories,args.TF_list,args.output)