from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.cli.ncbi_gene_results_to_python import NCBIgeneToPythonCli
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
import os
from itertools import islice
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='go_term_enrichment')
parser.add_argument('file_names', type=str, help='Name of folder and filenames for the promoters extracted')
parser.add_argument('go_directory', type=str, help='Directory location of go term enrichment files')
parser.add_argument('background_gene_set', type=str, help='Location of background gene set')
parser.add_argument('NCBI_gene_list', type=str, help='Location of NCBI gene list')
parser.add_argument('genes_of_interest', type=str, help='Location of genes of interest')
parser.add_argument('variable1_name', type=str, help='Optional replacement name for 2nd variable eg. non-specific',default = 'constitutive', nargs="?")
parser.add_argument('variable2_name', type=str, help='Optional replacement name for 2nd variable eg. tissue_specific',default = 'variable', nargs="?")
parser.add_argument('author_name', type=str, help='Optional replacement name for author in reference to the geneset',default = 'Czechowski', nargs="?")
args = parser.parse_args()

NCBI_gene_list_filtered = 'gene_result_filtered.txt'

#make directory for the output files to be exported to
dirName = args.go_directory
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " created") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")

def dl_files(go_directory):
    """function to download latest ontologies and associations files from geneontology.org
    specify the directory to download the files to"""
    
    #change to go directory
    os.chdir(go_directory)
    
    # Get http://geneontology.org/ontology/go-basic.obo
    obo_fname = download_go_basic_obo()
    
    #print go file version:
    with open(obo_fname) as fin:
        for line in islice(fin, 1, 2):
            print(line)
    
    #download gene2go annotation file
    fin_gene2go = download_ncbi_associations()
    
    return obo_fname,fin_gene2go


def generate_background(background_gene_set, NCBI_gene_list,NCBI_gene_list_filtered):
    """function to filter the NCBI protein copding arabidopsis gene list to only include genes in my background set"""
    #load background genes
    background_df = pd.read_table(background_gene_set, sep='\t', header=0)
    
    #load NCBI gene list
    NCBI_list = pd.read_table(NCBI_gene_list, sep='\t', header=0)
    
    #Add AGI column to NCBI_list
    NCBI_list_AGI = NCBI_list.assign(AGI=NCBI_list.Aliases.str.extract(r'(.*?)\,'))
    #if NAN in the AGI column then copy the Alias column
    NCBI_list_AGI.AGI.fillna(NCBI_list.Aliases, inplace=True)
    
    #filter the NCBI gene list to the background genes of interest
    NCBI_list_filtered = NCBI_list_AGI[NCBI_list_AGI.AGI.isin(background_df.AGI)]
    #reset  index
    NCBI_list_filtered.reset_index(inplace=True, drop=True)
    
    #remove AGI column
    NCBI_list_filtered = NCBI_list_filtered[['tax_id','Org_name','GeneID','CurrentID','Status','Symbol','Aliases','description',
                                             'other_designations','map_location','chromosome','genomic_nucleotide_accession.version',
                                             'start_position_on_the_genomic_accession','end_position_on_the_genomic_accession','orientation','exon_count','OMIM']]
    
    #save file
    NCBI_list_filtered.to_csv(NCBI_gene_list_filtered,sep='\t',header=1,index=None)

    
    #convert NCBI gene tsv file to a Python module as in https://github.com/tanghaibao/goatools/blob/master/notebooks/backround_genes_ncbi.ipynb
    NCBIgeneToPythonCli().tsv_to_py(NCBI_gene_list_filtered,'genes_ncbi_3702_proteincoding.py')
       
    #Initialize a GOEA object
    #The GOEA object holds the Ontologies, Associations, and background.
    #Numerous studies can then be run withough needing to re-load the above items
    
def load_files(obo_fname,fin_gene2go):
    """function to load ontologies, associations and background gene set and then initialise a GOEA object"""
    
    #import the python module created in generate_background()
    #find specificy the current folder lcoation as the location of the module
    import sys
    sys.path.insert(1, '.')
    #import the module
    from genes_ncbi_3702_proteincoding import GENEID2NT as GeneID2nt_ara
    
    #load ontologies
    obodag = GODag(obo_fname)
    
    #load associations
    # Read NCBI's gene2go. Store Arabidopsis thaliana annotations in a list of named tuples
    objanno = Gene2GoReader(fin_gene2go, taxids=[3702])
    # Get namespace2association where:
    #    namespace is:
    #        BP: biological_process               
    #        MF: molecular_function
    #        CC: cellular_component
    #    assocation is a dict:
    #        key: NCBI GeneID
    #        value: A set of GO IDs associated with that gene
    ns2assoc = objanno.get_ns2assc()

    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated Arabidopsis genes".format(NS=nspc, N=len(id2gos)))
    

    
    goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_ara.keys(), # List of filtered Arabidopsis protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    return goeaobj,obodag,ns2assoc

def map_AGI2NCBI(NCBI_gene_list_filtered,genes_of_interest):
    """map AGI gene IDs to NCBI IDs"""
    #read in NCBI list
    NCBI_list = pd.read_table(NCBI_gene_list_filtered, sep='\t', header=0)
    
    #Add AGI column to NCBI_list
    NCBI_list_AGI = NCBI_list.assign(AGI=NCBI_list.Aliases.str.extract(r'(.*?)\,'))
    #if NAN in the AGI column then copy the Alias column
    NCBI_list_AGI.AGI.fillna(NCBI_list.Aliases, inplace=True)
    
    #read in genes of interest
    genes_of_interest_df = pd.read_table(genes_of_interest,sep='\t',header=None)
    cols = ['AGI','gene_type']
    genes_of_interest_df.columns = cols
    
    #merge the dfs
    merged = pd.merge(genes_of_interest_df, NCBI_list_AGI, how='left', on='AGI')
    
    #save the file
    merged.to_csv('mapped_NCBI.tsv',sep='\t',header=True, index=False)
    
    return merged

def run_GOEA(study_genes,gene_class):
    "run Gene Ontology Enrichment Analysis. Provide study genes in dictionary key format and a gene_class string"
    #extract geneIDs
    geneids_study = study_genes
    
    #run goea
    goea_results_all = goeaobj.run_study(geneids_study)
    #filter to only significant results
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    #write to txt file
    goeaobj.wr_txt(f"{gene_class}_sig_enrichment.txt", goea_results_sig)
    
    
    
    
    #plot results
    plot_results(f"{gene_class}" + "_sig_enrichment_{NS}.png", goea_results_sig)
    

    
    return goea_results_sig,goea_results_all
    
## Download Ontologies and Associations    
obo_fname,fin_gene2go = dl_files(args.go_directory)

## Load Ontologies, Associations and Background gene set then initialize a GOEA object
#Filter NCBI gene list and create python module as in https://github.com/tanghaibao/goatools/blob/master/notebooks/backround_genes_ncbi.ipynb
generate_background(args.background_gene_set, args.NCBI_gene_list,NCBI_gene_list_filtered)

#load files aand initialise a GOEA object
goeaobj,obodag,ns2assoc = load_files(obo_fname,fin_gene2go)

## Read study genes

#map AGI gene ID to NCBI gene ID
study_genes = map_AGI2NCBI(NCBI_gene_list_filtered,args.genes_of_interest)

#select gene_types of interest
constitutive_genes = study_genes[study_genes.gene_type == 'constitutive']
variable_genes = study_genes[study_genes.gene_type == args.variable2_name]

#make dictionary where key is geneid and value is gene symbol
constitutive_geneid2symbol = dict(zip(constitutive_genes.GeneID,constitutive_genes.Symbol))
variable_geneid2symbol = dict(zip(variable_genes.GeneID,variable_genes.Symbol))

#get NCBI gene ids in dictionary format
constitutive_NCBI_IDs = constitutive_geneid2symbol.keys()
variable_NCBI_IDs = variable_geneid2symbol.keys()

## Run Gene Ontology Enrichment Analysis (GOEA) using Benjamini/Hochberg FDR correction
#keep only significant results
constitutive_goea_sig,constitutive_goea_all = run_GOEA(constitutive_NCBI_IDs,f'constitutive_{args.author_name}')
variable_goea_sig,variable_goea_all = run_GOEA(variable_NCBI_IDs,f'{args.variable2_name}_{args.author_name}')
# #plot colouring by significance - constitutive
# # GO terms colored by P-value:

# #     pval < 0.005 (light red)
# #     pval < 0.01 (light orange)
# #     pval < 0.05 (yellow)
# #     pval > 0.05 (grey) Study terms that are not statistically significant
# plot_gos("constitutive_intracellular_300genes_gene_counts.png",
#          ['GO:0006886'], # Source GO ids
#          obodag,
#          goea_results=constitutive_goea_all) # Use pvals for coloring

# #plot colouring by significance - constitutive and add gene names
# # GO terms colored by P-value:

# #     pval < 0.005 (light red)
# #     pval < 0.01 (light orange)
# #     pval < 0.05 (yellow)
# #     pval > 0.05 (grey) Study terms that are not statistically significant
# plot_gos("constitutive_intracellular_300genes.png",
#          ['GO:0006886'], # Source GO ids
#          obodag,
#          goea_results=constitutive_goea_all,
#         id2symbol=constitutive_geneid2symbol,
#         study_items=8, # Only only 8 gene Symbols max on GO terms
#         items_p_line=4, # Print 3 genes per line) # Use pvals for coloring
#         )

# #plot colouring by significance - variable
# # GO terms colored by P-value:

# #     pval < 0.005 (light red)
# #     pval < 0.01 (light orange)
# #     pval < 0.05 (yellow)
# #     pval > 0.05 (grey) Study terms that are not statistically significant
# plot_gos("variable_intracellular_300genes_gene_counts.png",
#          ['GO:0006886'], # Source GO ids
#          obodag,
#          goea_results=variable_goea_all) # Use pvals for coloring


# #plot colouring by significance - variable and add gene names
# # GO terms colored by P-value:

# #     pval < 0.005 (light red)
# #     pval < 0.01 (light orange)
# #     pval < 0.05 (yellow)
# #     pval > 0.05 (grey) Study terms that are not statistically significant
# plot_gos("variable_intracellular_300genes.png",
#          ['GO:0006886'], # Source GO ids
#          obodag,
#          goea_results=variable_goea_all,
#         id2symbol=variable_geneid2symbol,
#         study_items=8, # Only only 8 gene Symbols max on GO terms
#         items_p_line=4, # Print 3 genes per line) # Use pvals for coloring
#         )