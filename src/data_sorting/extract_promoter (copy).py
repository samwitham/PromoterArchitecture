from pybedtools import BedTool
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import pymysql.cursors
from pyfaidx import Fasta
import pandas as pd
import os
from pprint import pprint, pformat
import argparse

#create argparse variables


def fasta_chromsizes(genome, output_file):
    """extracts chromosome sizes in a form compatible with BedTools.flank (as a dictionary)"""
    #output = open(output_file, 'w')
    the_genome = Fasta(genome) #using pyfaidx Fasta, parse the genome file as a fasta
     
    chromsizes = {} #make dictionary called chromsizes
    for key in the_genome.keys():
       
        #chromsizes[key] = f'({0}, {len(the_genome[key])})' #add the chromosome name and length to dictionary
        chromsizes[f'Chr{key}'] = f'({len(the_genome[key])})' #add the chromosome name and length to dictionary
        #chromsizes[key] = f'({len(the_genome[key])})' #add the chromosome name and length to dictionary
    #output_file.write(chrom)
    chromsizes_string = ''
    
    for k,v in chromsizes.items():
        chromsizes_string = chromsizes_string + f'{k}\t{v}\n'
    
    
    with open(output_file, 'w') as output:            
        output.write(chromsizes_string.replace('(','').replace(')',''))


def extract_genes(gene_gff,output_file):
    """This function extracts all whole genes from a gff3 file, ignoring gene features, and adds them to an output file"""
    limit_info = dict(gff_type = ['gene'])
    #matchedLine = ''
    output = open(output_file, 'w')
    
    with open(gene_gff, 'r') as in_handle:                    
            GFF.write(GFF.parse(in_handle, limit_info=limit_info),output)
    output.close()
   


def add_promoter(genes_gff,chromsize,promoter_length,output_file):
    """This function adds a promoter of a certain length to each gene in the input file and exports to an output file"""
    #output = open(output_location, 'w') #make output file with write capability
    genes = BedTool(genes_gff)
    promoters = genes.flank(g=chromsize, l=promoter_length, r=0, s=True)
      

   
         
    with open(output_file,'w') as f:
        f.write(str(promoters))



def promoter_overlap(promoter_gff,allfeatures_gff,output_file):
    """function to create file containing promoters which overlap other genome features"""
    promoters = BedTool(promoter_gff) #read in files using BedTools
    features = BedTool(allfeatures_gff)
    #report chromosome position of overlapping feature, along with the promoter which overlaps it (only reports the overlapping nucleotides, not the whole promoter length. Can use u=True to get whole promoter length)
    intersect = promoters.intersect(features, u=True) #could add u=True which indicates we want to see the promoters that overlap features in the genome
    with open(output_file, 'w') as output:
        output.write(str(intersect))




def find_closest_TSS(gene_gff,TSS_gff,output_location):
    """this reads in the genes gff file and TSS gff file, finds the closest gene each TSS belongs to"""




def remove_characters_linestart(input_location,output_location,oldcharacters,newcharacters,linestart):
    """this function removes characters from the start of each line in the input file and sends modified lines to output"""
    output = open(output_location, 'w') #make output file with write capability
    with open(input_location, 'r') as infile:  
        
        for line in infile:
            line = line.strip() # removes hidden characters/spaces
            if line[0] == linestart:
                                 
                line = line.replace(oldcharacters, newcharacters) #remove characters from start of line, replace with new characters        
            output.write(line + '\n') #output to new file
    output.close()

def count_promoters(in_file, out_file):
    """this function creates a text file detailing the number of promoters in an input GFF file """
    examiner = GFFExaminer()
    in_handle = open(in_file,'r')
    with open(out_file, 'w') as fout:
        fout.write(pformat(examiner.available_limits(in_handle)))
    in_handle.close()

# def gff2bed(gff_file,bed_file):
#     """This function renames the third column of a gff3 file with the ID=ATetc., and adds them to an output file"""
    
    
#     output = open(output_file, 'w')
    
#     with open(gene_gff, 'r') as in_handle:                    
#             GFF.write(GFF.parse(in_handle, limit_info=limit_info),output)
#     output.close()


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

fasta_chromsizes(genome, chromsizes_file)




#rename mitochondria and chloroplast to M and C
remove_characters_linestart(chromsizes_file, chromsizes_file_renamedChr_temp, 'mitochondria','M', 'C')
remove_characters_linestart(chromsizes_file_renamedChr_temp, chromsizes_file_renamedChr, 'chloroplast','C','C')
os.remove(chromsizes_file_renamedChr_temp)



#extract_genes(genes,genesonly_gff)
extract_genes(genes,genesonly_gff)
#this changes chromosome no. to 1 rather than Chr1



add_promoter(genesonly_gff,chromsizes_file_renamedChr,1000,promoters)



#create file containing only promoters which overlap other genome features
promoter_overlap(promoters,genes,overlapping_promoters)


#examine this gff3 promoter file and compare to all promoters
#in_file = overlapping_promoters
#examiner = GFFExaminer()
#in_handle = open(in_file)
#pprint(examiner.available_limits(in_handle))
#in_handle.close()

count_promoters(overlapping_promoters, f'{directory_path}/data/genomes/overlapping_promoters.txt')

# all promoters
#in_file = promoters
#examiner = GFFExaminer()
#in_handle = open(#in_file)
#pprint(examiner.available_limits(in_handle))
#in_handle.close()


#promoters overlapping only genes
promoter_overlap(promoters,genesonly_gff,promoterandgenes_only_overlap)


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