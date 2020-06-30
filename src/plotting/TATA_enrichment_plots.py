import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='TATA_enrichment_plots')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('gat_TATA_constitutive_output', type=str, help='Location of constitutive promoter gat analysis output')
parser.add_argument('gat_TATA_variable_output', type=str, help='Location of variable promoter gat analysis output')
parser.add_argument('output_prefix', type=str, help='Output prefix to add to plot file name')

args = parser.parse_args()

def create_plot(gat_TATA_constitutive_output,gat_TATA_variable_output):
    """import and process the raw outputs after running gat (Genomic association tester). Then create barplot of constitutive and variable gene TATA enrichment"""
    #import gat output files as dfs
    constitutive = pd.read_table(gat_TATA_constitutive_output, sep='\t', header=0)
    variable = pd.read_table(gat_TATA_variable_output, sep='\t', header=0)
    #merge dfs
    merged = pd.concat([constitutive,variable], ignore_index=True)
    #set style to ticks
    sns.set(style="ticks", color_codes=True)
    #bar chart, 95% confidence intervals
    plot = sns.barplot(x="annotation", y="l2fold", data=merged,order=["constitutive", "variable"])
    plot.axhline(0, color='black')
    plt.xlabel("Gene type")
    plt.ylabel("Log2-fold enrichment over background").get_figure().savefig(f'../../data/output/{args.file_names}/TATA/plots/{args.output_prefix}_log2fold.pdf', format='pdf')
    
#Create barplot
create_plot(args.gat_TATA_constitutive_output,args.gat_TATA_variable_output)
