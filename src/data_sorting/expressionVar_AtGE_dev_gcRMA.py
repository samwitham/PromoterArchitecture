#!/usr/bin/python
#=====================================================================================================================================================
'''
################################################################################

expressionVar_bins.py

parser for the expression variation data from http://www.weigelworld.org/resourc
es/microarray/AtGenExpress/AtGE_dev_gcRMA.txt.zip/ (which is no longer hosted).

Also take output from Affymetrix MAS5 software (Bioconductor affy) for presence
/ abscence estimation.

Calculates CV as per 2005 paper

################################################################################
'''

import os
import sys
import csv
import time
import argparse

import numpy as np

from operator import itemgetter
from statistics import mean, stdev
from collections import OrderedDict

def loader (in_data):
    '''
    Load raw data into memory. For Expression data, inverse log2 is performed to
    recreate expression value.
    '''
    #create empty dictionary
    processed_dict = {}

    #open the input file
    with open(in_data) as dat:

        #create an iterator, splitting the raw data on the tab
        iter_raw = csv.reader(dat, delimiter = '\t')

        #extract the header line
        header_raw = next(iter_raw)

        #iterate through the above
        for line in iter_raw:

            #check for duplicates in the input
            if line[0] not in processed_dict:

                #handle the downloaded data differently as not a raw matrix
                if 'AtGE_dev_gcRMA' in in_data:
                    #probe_id and gene_id joined by :, used as key
                    processed_dict[':'.join([line[0],line[1]])] = {}
                    for count, item in enumerate(line[2:]):
                        #invert the log2 transform on the matrix
                        processed_dict[':'.join([line[0],line[1]])][header_raw[count + 1]] = 2**float(item)
                else:
                    #handling the mas5 calls
                    processed_dict[line[0]] = {}

                    for count, item in enumerate(line[1:]):
                        processed_dict[line[0]][header_raw[count + 1]] = item
            else:
                print("\nDuplicated probe id detected during loading of raw matrix\n")
                sys.exit()

    return processed_dict


def mas5Stats (expressionSet, mas5calls):
    '''
    Calculate mean, SD, CV, and %P for all probes. Filter values that appeear
    <80% Present based on mas5 calls
    '''
    #create balnk dict to store output
    statistics_dict = {}

    #iterate through expression data
    for full_probe_id, exp_value_dict in expressionSet.items():

        #Calculate stats on the expression values
        expression_mean  = mean([float(x) for x in expressionSet[full_probe_id].values()])
        expression_stdev = stdev([float(x) for x in expressionSet[full_probe_id].values()])
        expression_CV    = expression_stdev / expression_mean

        #extract probe_id only
        probe_id      = full_probe_id.split(':')[0]

        #count the number of these removed
        present_count = [x for x in mas5calls[probe_id].values()].count('P')

        #count the total number of samples
        n_samples     = len([x for x in expressionSet[full_probe_id].values()])

        #calculate the total proportion of samples Present
        prop_present  = int((present_count / n_samples)*100)

        #filter values that appear below 80% presence and collate output into result dict
        #also filter Affymetrix normalisation probes
        if prop_present >= 80 and not probe_id.startswith('AFFX'):
            statistics_dict[probe_id] = [full_probe_id.split(':')[1].upper(),expression_mean, expression_stdev, expression_CV, prop_present]

    return statistics_dict

def sortOutput(input_dict):
    '''
    Sort by CV value.
    '''
    sorted_stats = OrderedDict(sorted(input_dict.items(), key = lambda item: item[1][3]))

    return sorted_stats

def currentAnnotCheck(input_stats, input_annot, Araport_housekeeping_set):
    '''
    create flag based on genes present in current annotation of Athal
    '''
    #create blank list to store current gene ids
    gene_list = []
    hk_list   = []

    #load annotation
    with open(input_annot) as arabidopsis_annotation:
        for line in csv.reader(arabidopsis_annotation, delimiter = '\t'):
            if not line[0].startswith('#'):
                if line[2] == 'gene':
                    gene_list.append((line[8].split(';')[0].split(':')[1]))

    with open(Araport_housekeeping_set) as ara_hk_in:
        for line in ara_hk_in:
            hk_list.append(line.split('.')[0])

    print('Checking', str(len(gene_list)), 'gene ids loaded from current annotation.')
    print('Checking', str(len(hk_list)), ' housekeeping genes loaded from Araport annotation.\n')

    for stats in input_stats.values():
        if stats[0] not in gene_list:
            stats.append('0')
        else:
            stats.append('1')

        if stats[0] not in hk_list:
            stats.append('0')
        else:
            stats.append('1')





def main (raw_expr, mas5_pres, ara_annot, araport_houskeeeping, output_directory):
    print('\nexpressionVar_bins.py')
    print('\nReading AtGE_dev_gcRMA expression set from: ', raw_expr)
    print('Reading mas5 presence/absence calls from: ', mas5_pres, '\n')

    mas5_expSet = loader(raw_expr)
    print('Expression set loaded.')

    mas5_calls = loader(mas5_pres)
    print('mas5 presence/absense calls loaded.\n')

    print('Calculating stats...')
    expSet_statistics = sortOutput(mas5Stats(mas5_expSet, mas5_calls))
    print('Done.\n')

    currentAnnotCheck(expSet_statistics, ara_annot, araport_houskeeeping)

    with open(os.path.join(output_directory, 'AtGE_dev_gcRMA__all_probes__CV.tsv'), 'w') as all_out:
        all_out.write('\t'.join('rank,probe_id,gene_id,mean,SD,CV,propP,current,Araport_const'
        count = 1
        for k, v in expSet_statistics.items():
            all_out.write('\t'.join([str(count),k,'\t'.join([str(x) for x in v ])])+'\n')
            count = count + 1
            # time.sleep(.1)
    print('All probe output written to: ', os.path.join(output_directory, 'AtGE_dev_gcRMA__all_probes__CV.tsv'))
    # with open(os.join(output_directory, 'AtGE_dev_gcRMA__all_probes__CV.tsv'), 'w') as all_out:


    # print(expSet_statistics['253287_at'])
    # print(mas5_pres_dict['253287_at'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('raw_expression_in', type=str, help='input file: the final output of the bash based eqtl work on FIMO output')
    parser.add_argument('mas5_presence', type=str, help='affmetrix presence/abscence data')
    parser.add_argument('annotation',    type=str, help='current Athal annotation')
    parser.add_argument('housekeeping_set',  type=str, help='Araport housekeeping genes from Data S4 from Cheng et al. 2016')
    parser.add_argument('out_dir',       type=str, help='path to directory where the output should be written')

    args = parser.parse_args()

    main(args.raw_expression_in, args.mas5_presence, args.annotation, args.housekeeping_set,  args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
