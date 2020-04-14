from pybedtools import BedTool
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import pymysql.cursors
from pyfaidx import Fasta
import pandas as pd
import os
from pprint import pprint, pformat



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


def extract_genes(gene_gff,output_file):
    """This function extracts all whole genes from a gff3 file, ignoring gene features, and adds them to an output file"""
    #limit dictionary to genes
    limit_info = dict(gff_type = ['gene'])
    #open output file
    output = open(output_file, 'w')
    #open gff file, parse it, limiting to genes only. Save the file.
    with open(gene_gff, 'r') as in_handle:                    
            GFF.write(GFF.parse(in_handle, limit_info=limit_info),output)
    output.close()
    pass
   


def add_promoter(genes_gff,chromsize,promoter_length,output_file):
    """This function adds a promoter of a certain length to each gene in the input file and exports to an output file"""
    #output = open(output_location, 'w') #make output file with write capability
    #parse gff file containing only genes.
    genes = BedTool(genes_gff)
    #extract promoters upsteam using chromsize file and specified promoter length. r, no. of bp to add to end coordinate. s, based on strand.
    promoters = genes.flank(g=chromsize, l=promoter_length, r=0, s=True)
    #write to file         
    with open(output_file,'w') as f:
        f.write(str(promoters))



def promoter_overlap(promoter_gff,allfeatures_gff,output_file):
    """function to create file containing promoters which overlap other genome features"""
    #read in files using BedTools
    promoters = BedTool(promoter_gff) 
    features = BedTool(allfeatures_gff)
    #report chromosome position of overlapping feature, along with the promoter which overlaps it (only reports the overlapping nucleotides, not the whole promoter length. Can use u=True to get whole promoter length)
    #f, the minimum overlap as fraction of A. F, nucleotide fraction of B (genes) that need to be overlapping with A (promoters)
    #wa, Write the original entry in A for each overlap.
    #u, write original A entry only once even if more than one overlap
    intersect = promoters.intersect(features,f=0.001, F=0.001, u=True, wa=True) #could add u=True which indicates we want to see the promoters that overlap features in the genome
    #write to file
    with open(output_file, 'w') as output:
        output.write(str(intersect))

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
    """this function create a file containing all promoters with an upstream gene going in the other direction ie. potential bidirectional promoters"""
    promoters = pd.read_table(in_file, sep='\t', header=4)
    cols2 = ['promoter_AGI', 'gene_type']
    promoters.columns = cols2
    
def bidirectional_proms(in_file, out_file):
    """this function create a file containing all promoters with an upstream gene going in the other direction ie. potential bidirectional promoters"""
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

directory_path = '/ei/workarea/group-eg/project_PromoterArchitecturePipeline'

TSS_raw = f'{directory_path}/data/TSS_data/AnnotatedPEATPeaks.txt'
#TSS_renamedChr = '{directory_path}//data/TSS_data/AnnotatedPEATPeaks_renamedChr.gff'
#genome2 = "{directory_path}//data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel_sl_simpleChr.fasta"
genome = f'{directory_path}/data/genomes/TAIR10_chr_all.fas'



genes = f"{directory_path}/data/genomes/Araport11_GFF3_genes_transposons.201606.gff"
#genes_renamedChr = "{directory_path}//data/genomes/Araport11_GFF3_genes_transposons.201606_renamedChr.gff"
test_genes = f"{directory_path}/data/genomes/test_genes.gff3"
testgenesonly_gff = f"{directory_path}/data/genomes/testgenes_only.gff3"
temp = f"{directory_path}/data/TSS_data/temp.txt"
#TSS = "{directory_path}//data/TSS_data/AnnotatedPEATPeaks_renamedcol.gff"
#TSS = "{directory_path}//data/TSS_data/TSStest.txt"
#find_closest_TSS(genes,output,temp)
genesonly_gff = f"{directory_path}/data/genomes/genesonly.gff3"
promoters = f"{directory_path}/data/genomes/promoters.gff3"
#genes_bed = "{directory_path}//data/genomes/genes.bed"
overlapping_promoters = f'{directory_path}/data/genomes/promoters_overlapping.gff3'
promoterandgenes_only_overlap = f'{directory_path}/data/genomes/promoterandgenes_only_overlap.gff3'


chromsizes_file = f'{directory_path}/data/genomes/chromsizes.chr'
chromsizes_file_renamedChr_temp = f'{directory_path}/data/genomes/chromsizes_renamedChr_temp.chr'
chromsizes_file_renamedChr = f'{directory_path}/data/genomes/chromsizes_renamedChr.chr'
#chromsizes_file2 = '{directory_path}//data/genomes/chromsizes2.chr'
promoters_renamedChr_temp = f'{directory_path}/data/genomes/promoters_renamedChr_temp.gff3'
promoters_renamedChr_temp2 = f'{directory_path}/data/genomes/promoters_renamedChr_temp2.gff3'
promoters_renamedChr = f'{directory_path}/data/genomes/promoters_renamedChr.gff3'
nonbidirectional_promoters = '../../data/genomes/nonbidirectional_proms.gff3'

fasta_chromsizes(genome, chromsizes_file)




#rename mitochondria and chloroplast to M and C
remove_characters_linestart(chromsizes_file, chromsizes_file_renamedChr_temp, 'mitochondria','M', 'C')
remove_characters_linestart(chromsizes_file_renamedChr_temp, chromsizes_file_renamedChr, 'chloroplast','C','C')
os.remove(chromsizes_file_renamedChr_temp)



#extract_genes(genes,genesonly_gff)
#note - this changes chromosome no. to 1 rather than Chr1
extract_genes(genes,genesonly_gff)

#createfile containing all nonbidirectional genes (bidirectional = genes with an upstream gene in the other direction ie. potential overlapping promoters)
bidirectional_proms(genesonly_gff, nonbidirectional_promoters)

#add 1000 bp promoters upstream of genes, using chromsizes file, input gene annotation file (gff) and output promoters gff file
add_promoter(nonbidirectional_promoters,chromsizes_file_renamedChr,1000,promoters)

#create file containing only promoters which overlap other genome features
promoter_overlap(promoters,genes,overlapping_promoters)

#count no. of promoters in overlapping promoters file
count_promoters(overlapping_promoters, f'{directory_path}/data/genomes/overlapping_promoters.txt')

# all promoters
#in_file = promoters
#examiner = GFFExaminer()
#in_handle = open(#in_file)
#pprint(examiner.available_limits(in_handle))
#in_handle.close()


#promoters overlapping only genes
promoter_overlap(promoters,genesonly_gff,promoterandgenes_only_overlap)

#count no. of promoters in overlapping promoters file
count_promoters(promoterandgenes_only_overlap, f'{directory_path}/data/genomes/promoterandgenes_only_overlap.txt')

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