#this script maps the 4th column of the promoter_motifs.bed (from Dap-seq dataset only) to the Arabidopsis gene ID nomenclature (eg. AT4G3800). It also adds a TF name column
import pandas as pd
import argparse
from pybedtools import BedTool
import re
import numpy as np

parser = argparse.ArgumentParser(description='map_motif_ids')
parser.add_argument('motifs_bed', type=str, help='Input location of motif bed file')
parser.add_argument('geneIDtable', type=str, help='Input location of geneIDtable file derived from the .json (see .ipynb notebook)')
parser.add_argument('motifs_bed_mapped', type=str, help='Output location of motif bed file')
args = parser.parse_args()
#motif_mapping_table = '../../data/FIMO/motif_data/motif_map_IDs.txt'
#motifs_bed = '../../data/FIMO/promoters_renamedChr_motifs.bed'
#motifs_bed_mapped = '../../data/FIMO/promoters_renamedChr_motifs_mapped.bed'

def rename_motif(motifs_bed):
    """function to add a TF name column in the motifs.bed file (for DAP-seq cistrome motifs only)"""
    #read in motifs_bed_file
    motifs = pd.read_table(motifs_bed, sep='\t', header=None)
    cols = ['chr', 'start', 'stop', 'name_rep', 'score', 'strand', 'sequence_name', 'p-value', 'q-value', 'matched_sequence']
    motifs.columns = cols
    motifs['TF'] = motifs.name_rep    
    capitalise = lambda x: x.upper()
    motifs.TF = motifs.TF.apply(capitalise)
    #replace characters upto and including the '.' in the TF column
    motifs.TF = motifs.TF.replace('^[^\.]*\.', '', regex=True)
    #replace characters after the '_'
    motifs.TF = motifs.TF.replace('_.*$', '', regex=True)
    return motifs

def map_ID2(motifs, geneIDtable, output_bed):   
    """function to rename the TF column values in a motifs.bed file to the Arabidopsis gene ID nomenclature using geneIDtable. (for DAP-seq cistrome motifs only). Outputs a bed file."""
    #remove '_m' from end of name_rep value in motifs
    motifs.name_rep = motifs.name_rep.str.replace('_m1', '')
    merged = pd.merge(motifs,geneIDtable, on='name_rep')
    #print(merged.shape)
    #make bed file
    sorted_motifs = merged.sort_values(['chr','start'])
    bed = BedTool.from_dataframe(sorted_motifs).saveas(output_bed)
    
#make df of geneIDtable
geneIDtable = pd.read_table(args.geneIDtable, sep='\t', header=None)
cols = ['name_rep', 'at_id']
geneIDtable.columns = cols
    
motifs_renamed = rename_motif(args.motifs_bed)
map_ID2(motifs_renamed, geneIDtable, args.motifs_bed_mapped)
    
    
    #Created by Sam Witham 12/11/19