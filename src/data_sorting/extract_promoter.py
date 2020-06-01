from pybedtools import BedTool
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from pyfaidx import Fasta
import pandas as pd
import os
from pprint import pprint, pformat
import argparse
import io

parser = argparse.ArgumentParser(description='extract_promoter')
parser.add_argument('directory_path', type=str, help='Location of base directory')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('--remove_bidirectional', help='Exclude potentially bidirectional promoters with upstream TSS >2000bp from TSS in opposite direction.', action="store_true")
parser.add_argument('--prevent_overlapping_genes', help='Reduce size of promoters if needed until they are not overlapping other genes.', action="store_true")
parser.add_argument('--only_open_chromatin', help='Reduce size of promoters if needed so that they only fall within open chromatin.', action="store_true")
args = parser.parse_args()

def fasta_chromsizes(genome, output_file):
    """extracts chromosome sizes in a form compatible with BedTools.flank"""
    
    the_genome = Fasta(genome) #using pyfaidx Fasta, parse the genome file as a fasta
     
    chromsizes = {} #make dictionary called chromsizes
    for key in the_genome.keys():       
    
        chromsizes[f'Chr{key}'] = f'({len(the_genome[key])})' #add the chromosome name and length to dictionary        
    #create empty string
    chromsizes_string = ''
    #iterate over chromsizes dictionary key/values. Add key, tab, value, newline to the string iteratively
    for k,v in chromsizes.items():
        chromsizes_string = chromsizes_string + f'{k}\t{v}\n'
    
    #write output file, deleting the parentheses created from the dictionary. This file is suitable for use in BedTools.flank
    with open(output_file, 'w') as output:            
        output.write(chromsizes_string.replace('(','').replace(')',''))


def extract_genes(gene_gff,temp,output_file):
    """This function extracts all whole protein coding genes from a gff3 file, ignoring gene features, and adds them to an output file"""
    #limit dictionary to genes
    limit_info = dict(gff_type = ['gene'])
    #open temporary file
    output = open(temp, 'w')
    #open gff file, parse it, limiting to genes only. Save the file.
    with open(gene_gff, 'r') as in_handle:                    
            GFF.write(GFF.parse(in_handle, limit_info=limit_info),output)
    output.close()
    
    #now remove the unwanted annotations that were added by GFF.write eg. Chr1	annotation	remark	1	30425192	.	.	.	gff-version=3
    #remove lines beginning with ##
    with open(temp, 'r') as tempfile, open(output_file, 'w') as newfile:        
        for line in tempfile:
            line = line.strip() # removes hidden characters/spaces
            if line[0] == "#":
                pass
            else:                
                #don't include lines beginning with ##
                newfile.write(line + '\n') #output to new file
    #remove temporary file
    os.remove(temp)       
    #read in gff file
    genes = pd.read_table(output_file, sep='\t', header=0)
    cols2 = ['chr', 'source', 'type', 'start','stop','dot1','strand','dot2','attributes']
    genes.columns = cols2
    #remove all lines where source is annotation
    no_annotation = genes[~genes.source.str.contains("annotation")]
    #remove all lines that are not protein coding
    protein_coding = no_annotation[no_annotation.attributes.str.contains('locus_type=protein_coding')]
    #create gff file containing no lines with annotation and no lines starting with ##
    protein_coding.to_csv(output_file,index=False,sep='\t',header=0)


def add_promoter(genes_gff,chromsize,promoter_length):
    """This function adds a promoter of a certain length to each gene in the input file and exports an output pyBedTools object"""
    #output = open(output_location, 'w') #make output file with write capability
    #parse gff file containing only genes.
    genes = BedTool(genes_gff)
    #extract promoters upsteam using chromsize file and specified promoter length. r, no. of bp to add to end coordinate. s, based on strand.
    promoters = genes.flank(g=chromsize, l=promoter_length, r=0, s=True)
    
    return promoters

def remove_promoter_overlap(promoter_gff, all_genes_gff, output_file):
    """function to create file containing promoters which overlap other genome features. 
    Then create a df with only promoters which overlap other genes. Then shorten them so they are no longer overlapping. 
    Then merge back with all the extracted promoters, keeping only the shortened ones"""
    all_proms = BedTool(promoter_gff) #read in files using BedTools
    features = BedTool(all_genes_gff)
    #report chromosome position of overlapping feature, along with the promoter which overlaps it (only reports the overlapping nucleotides, not the whole promoter length. Can use u=True to get whole promoter length)
    #f, the minimum overlap as fraction of A. F, nucleotide fraction of B (genes) that need to be overlapping with A (promoters)
    #wa, Write the original entry in A for each overlap.
    #wo,  Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. 
    #u, write original A entry only once even if more than one overlap
    intersect = all_proms.intersect(features, wo=True) #could add u=True which indicates we want to see the promoters that overlap features in the genome
    #Write to output_file
    with open(output_file, 'w') as output:
        #Each line in the file contains gff entry a and gff entry b that it overlaps plus the number of bp in the overlap so 19 columns
        output.write(str(intersect))    
    #read in gff file
    overlapping_proms = pd.read_table(output_file, sep='\t', header=0)
    cols = ['chrA', 'sourceA', 'typeA', 'startA','stopA','dot1A','strandA','dot2A','attributesA','chrB', 'sourceB', 'typeB', 'startB','stopB','dot1B','strandB','dot2B','attributesB','bp_overlap']
    overlapping_proms.columns = cols
    #make new columns for the new start and stop for feature A
    overlapping_proms['new_startA'] = overlapping_proms.startA
    overlapping_proms['new_stopA'] = overlapping_proms.startA
   
    #iterate over rows
    for i,data in overlapping_proms.iterrows():
        #if positive strand feature A, reduce length to the stop position of feature B + 1
        if overlapping_proms.loc[i,'strandA'] == '+':        
            overlapping_proms.loc[i, 'new_startA'] = overlapping_proms.loc[i, 'stopB'] + 1
        #if negative strand feature A, reduce length of promoter to start position of feature -1
        elif overlapping_proms.loc[i,'strandA'] == '-': 
            overlapping_proms.loc[i, 'new_stopA'] = overlapping_proms.loc[i, 'startB'] - 1
    #create df with just feature a containing new start and stop locations
    new_feature_A = overlapping_proms[['chrA', 'sourceA', 'typeA', 'new_startA', 'new_stopA', 'dot1A', 'strandA', 'dot2A', 'attributesA']]
    #rename columns
    cols2 = ['chr', 'source', 'type', 'start','stop','dot1','strand','dot2','attributes']
    new_feature_A.columns = cols2
    #create a buffer of promoter_gff
    promoter_gff_buffer = io.StringIO()
    promoter_gff_buffer.write(str(promoter_gff))
    #go back to beginning of the buffer
    promoter_gff_buffer.seek(0)
    #read in all_proms as df
    all_proms_df = pd.read_table(promoter_gff_buffer, sep='\t', header=0)
    all_proms_df.columns = cols2
    #merge all_proms_df with new_feature_A keeping the promoters in new_feature_A and only adding promoters from all_proms_df that aren't in new_feature_A
    #outer merge adding suffix to new_feature_A
    merged = pd.merge(all_proms_df, new_feature_A, how='left', on='attributes',suffixes=('','_reduced'))
    #if start_reduced column from new_feature_A is NaN, pass
    #replace NaN in start_reduced with string 'NAN'
    merged.fillna('NAN', inplace=True)
    for index,row in overlapping_proms.iterrows():
        if merged.loc[index, 'start_reduced'] == 'NAN':
        #merged.isnull(row.start_reduced):
            pass
        #if not NaN, replace start and stop from all_proms_df with those from new_feature_A
        else:
            merged.loc[index, 'start'] = merged.loc[index, 'start_reduced']
            merged.loc[index, 'stop'] = merged.loc[index, 'stop_reduced']
    

#         merged.loc[merged.start] == merged.loc[merged.start_reduced]
#         merged.loc[merged.stop] == merged.loc[merged.stop_reduced]
   
    #merged.columns = cols2
    promoter_gff_buffer.close()
    new_merged = merged[cols2]
    #remove trailing decimal .0 from start and stop
    new_merged = new_merged.astype({'start': 'int'})
    new_merged = new_merged.astype({'stop': 'int'})
    return new_merged

def remove_characters_linestart(input_location,output_location,oldcharacters,newcharacters,linestart):
    """this function removes characters from the start of each line in the input file and sends modified lines to output"""
    output = open(output_location, 'w') #make output file with write capability
    #open input file
    with open(input_location, 'r') as infile:  
        #iterate over lines in fuile
        for line in infile:
            line = line.strip() # removes hidden characters/spaces
            if line[0] == linestart:
                                 
                line = line.replace(oldcharacters, newcharacters) #remove characters from start of line, replace with new characters        
            output.write(line + '\n') #output to new file
    output.close()

def count_promoters(in_file, out_file):
    """this function creates a text file detailing the number of promoters in an input GFF file """
    examiner = GFFExaminer()
    #open input GFF file
    in_handle = open(in_file,'r')
    #output a text file, giving information such as no. of promoters in the file
    with open(out_file, 'w') as fout:
        fout.write(pformat(examiner.available_limits(in_handle)))
    in_handle.close()

    
def bidirectional_proms(in_file, out_file):
    """this function creates a file containing all promoters with an upstream gene going in the other direction ie. potential bidirectional promoters"""
    #read in gff file
    promoters = pd.read_table(in_file, sep='\t', header=2)
    cols2 = ['chr', 'source', 'type', 'start','stop','dot1','strand','dot2','attributes']
    promoters.columns = cols2
    #make sure lines are sorted
  
    promoters = promoters.sort_values(['chr','start']).reset_index(drop=True)
    #if bidirectional
    promoters['bidirectional'] = 'no'
    
    for i,data in promoters.iterrows():
        if i-1 >= 0:
            #print(i,data)
            if promoters.loc[i, 'strand'] == '+' and promoters.loc[i-1, 'strand'] == '-' and promoters.loc[i, 'start'] - promoters.loc[i-1, 'stop'] < 2000:
                promoters.loc[i, 'bidirectional'] = 'yes'
                promoters.loc[i-1, 'bidirectional'] = 'yes'

    with open(out_file, 'w') as output:  
        promoters[promoters.bidirectional == 'no'][['chr', 'source', 'type', 'start','stop','dot1','strand','dot2','attributes']].to_csv(out_file,index=False,sep='\t',header=0)


#TSS_raw = f'{args.directory_path}/data/TSS_data/AnnotatedPEATPeaks.txt'
#TSS_renamedChr = '{args.directory_path}//data/TSS_data/AnnotatedPEATPeaks_renamedChr.gff'
#genome2 = "{args.directory_path}//data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel_sl_simpleChr.fasta"
genome = f'{args.directory_path}/data/genomes/TAIR10_chr_all.fas'



genes = f"{args.directory_path}/data/genomes/Araport11_GFF3_genes_transposons.201606.gff"
#genes_renamedChr = "{args.directory_path}//data/genomes/Araport11_GFF3_genes_transposons.201606_renamedChr.gff"
#test_genes = f"{args.directory_path}/data/genomes/test_genes.gff3"
#testgenesonly_gff = f"{args.directory_path}/data/genomes/testgenes_only.gff3"
#temp = f"{args.directory_path}/data/TSS_data/temp.txt"
#TSS = "{args.directory_path}//data/TSS_data/AnnotatedPEATPeaks_renamedcol.gff"
#TSS = "{args.directory_path}//data/TSS_data/TSStest.txt"
#find_closest_TSS(genes,output,temp)
genesonly_gff = f"{args.directory_path}/data/genomes/{args.file_names}/genesonly.gff3"
promoters = f"{args.directory_path}/data/genomes/{args.file_names}/promoters.gff3"
#genes_bed = "{args.directory_path}//data/genomes/genes.bed"
overlapping_promoters = f'{args.directory_path}/data/genomes/{args.file_names}/promoters_overlapping.gff3'
promoterandgenes_only_overlap = f'{args.directory_path}/data/genomes/{args.file_names}/promoterandgenes_only_overlap.gff3'


chromsizes_file = f'{args.directory_path}/data/genomes/{args.file_names}/chromsizes.chr'
#need to change temporary files to scratch directory
chromsizes_file_renamedChr_temp = f'{args.directory_path}/data/genomes/{args.file_names}/chromsizes_renamedChr_temp.chr'
chromsizes_file_renamedChr = f'{args.directory_path}/data/genomes/{args.file_names}/chromsizes_renamedChr.chr'
#chromsizes_file2 = '{args.directory_path}//data/genomes/chromsizes2.chr'
promoters_renamedChr_temp = f'{args.directory_path}/data/genomes/{args.file_names}/promoters_renamedChr_temp.gff3'
promoters_renamedChr_temp2 = f'{args.directory_path}/data/genomes/{args.file_names}/promoters_renamedChr_temp2.gff3'
promoters_renamedChr = f'{args.directory_path}/data/genomes/{args.file_names}/promoters_renamedChr.gff3'
nonbidirectional_promoters = '../../data/genomes/{args.file_names}/nonbidirectional_proms.gff3'
temp_gff = f"{args.directory_path}data/genomes/genesonly_no_annotation_temp.gff3"
#make directory for the output files to be exported to
dirName = f'{args.directory_path}/data/genomes/{args.file_names}'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")


fasta_chromsizes(genome, chromsizes_file)




#rename mitochondria and chloroplast to M and C
remove_characters_linestart(chromsizes_file, chromsizes_file_renamedChr_temp, 'mitochondria','M', 'C')
remove_characters_linestart(chromsizes_file_renamedChr_temp, chromsizes_file_renamedChr, 'chloroplast','C','C')
os.remove(chromsizes_file_renamedChr_temp)



#extract_genes(genes,genesonly_gff)
#note - this changes chromosome no. to 1 rather than Chr1?? NOT SURE IT STILL DOES
extract_genes(genes,temp_gff,genesonly_gff)

#createfile containing all nonbidirectional genes (bidirectional = genes with an upstream gene in the other direction ie. potential overlapping promoters)
if args.remove_bidirectional:
    bidirectional_proms(genesonly_gff, nonbidirectional_promoters)
    selected_genes = nonbidirectional_promoters
else:
    selected_genes = genesonly_gff



#add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff
# add_promoter(selected_genes,chromsizes_file_renamedChr,1000)

#create file containing only promoters which overlap other genome features
#promoter_overlap(promoters,genes,overlapping_promoters)

#count no. of promoters in overlapping promoters file
#count_promoters(overlapping_promoters, f'{args.directory_path}/data/genomes/overlapping_promoters.txt')

# all promoters
#in_file = promoters
#examiner = GFFExaminer()
#in_handle = open(#in_file)
#pprint(examiner.available_limits(in_handle))
#in_handle.close()

#promoters overlapping only genes
if args.prevent_overlapping_genes:
    #add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff
    promoters_incl_overlap = add_promoter(selected_genes,chromsizes_file_renamedChr,1000)
    subtracted = remove_promoter_overlap(promoters_incl_overlap,selected_genes,promoterandgenes_only_overlap)
    with open(promoters,'w') as f:
        subtracted.to_csv(f,index=False,sep='\t',header=0)
    
else:
    #add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff
    promoters_gff = add_promoter(selected_genes,chromsizes_file_renamedChr,1000)
    with open(promoters,'w') as f:
        f.write(str(promoters_gff))

#count no. of promoters in overlapping promoters file
#count_promoters(promoterandgenes_only_overlap, f'{args.directory_path}/data/genomes/promoterandgenes_only_overlap.txt')

#examiner = GFFExaminer()
#in_handle = open(in_file)
#pprint.pprint(examiner.available_limits(in_handle))
#in_handle.close()

#remove the Chr from the chromosome names in promoters.gff3. Replace M with mitochondria and C with chloroplast
remove_characters_linestart(promoters, promoters_renamedChr_temp, 'Chr', '','C')
remove_characters_linestart(promoters_renamedChr_temp, promoters_renamedChr_temp2, 'M', 'mitochondria','M')
remove_characters_linestart(promoters_renamedChr_temp2, promoters_renamedChr, 'C', 'chloroplast','C')


os.remove(promoters_renamedChr_temp)
os.remove(promoters_renamedChr_temp2)