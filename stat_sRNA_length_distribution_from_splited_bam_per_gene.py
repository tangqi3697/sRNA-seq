#! /usr/bin/python

import getopt, sys, datetime, os
import numpy as np
import pandas as pd
import pysam

################
# 获取命令行参数。
def usage():
    print("For example:\n\n   python stat_sRNA_length_distribution_from_splited_bam_per_gene.py -i input.bam -I gene.ids.txt -l 18,31 -f False -t uniq -o sRNA_length_distribution \n")
    print("Parameters:");
    print("  -i:  output directory of split_featurecounts_bam_by_feature_tags.py.");
    print("  -I:  IDs of gene/intron etc.");
    print("  -l:  length range (default: 18,30).");
    print("""  -f:  Assign fractional counts to features (Choose: True/False, default: False).
        When 'True' is specified, each reported alignment from a multi-mapping read will carry a fractional count of 1/x, instead of 1 (one),
        where x is the total number of alignments reported for the same read. .""");
    print("  -t:  mapping type (Choose: all,multi,uniq, default:uniq).");
    print("  -o:  output name.");
    print("  -h:  print this page.");
    sys.exit()

def obtainParameter():
    opts, args = getopt.getopt(sys.argv[1:], "hi:I:l:f:t:o:")
    
    parameters = {
        'len_range':[18, 30],
        'fraction': 'False',
        'mapping_type':'uniq',
    }
    if opts:
        for op, value in opts:
            if op == "-i": parameters['input_file'] = value; print(op, value.split('/')[-1])
            elif op == "-I": parameters['ID_file'] = value; print(op, value.split('/')[-1])
            elif op == "-l": parameters['len_range'] = [int(x) for x in value.split(',')]; print(op, value.split('/')[-1])
            elif op == "-f": parameters['fraction'] = value; print(op, value.split('/')[-1])
            elif op == "-t": parameters['mapping_type'] = value; print(op, value.split('/')[-1])
            elif op == "-o": parameters['output_name'] = value; print(op, value.split('/')[-1])
            elif op == "-h": usage()
    else: usage()
    if parameters['fraction'] == 'True': parameters['fraction'] = True
    if parameters['fraction'] == 'False': parameters['fraction'] = False
    return(parameters)


def get_len_dist_df(parameters, ID_list):
    ID_dict = {}
    for ID in ID_list:
        if ID not in ID_dict:
            ID_dict[ID] = ''
    len_list = range(parameters['len_range'][0], parameters['len_range'][1]+1)
    zeros_df = np.zeros((len(ID_list), len(len_list)+1))
    len_dist_df = pd.DataFrame(zeros_df, index=ID_list, columns=[ x for x in len_list]+['Total'])
    return(len_dist_df, ID_dict)

def stat_sRNA_length_by_fraction_with_all(parameters, len_dist_df, ID_dict):
    with pysam.AlignmentFile(parameters['input_file'], 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XT'):
                IDs = line.get_tag('XT').split(',')
                for ID in IDs:
                    if line.query_length >= parameters['len_range'][0] and line.query_length <= parameters['len_range'][1] and ID in ID_dict:
                        if line.has_tag('XY'):
                            len_dist_df.loc[ID, line.query_length] += line.get_tag('XZ')/len(IDs)
                            len_dist_df.loc[ID, 'Total'] += line.get_tag('XZ')/len(IDs)
    return(len_dist_df)

def stat_sRNA_length_by_fraction_with_multi(parameters, len_dist_df, ID_dict):
    with pysam.AlignmentFile(parameters['input_file'], 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XY') and line.get_tag('XY') != 'U':
                if line.has_tag('XT'):
                    IDs = line.get_tag('XT').split(',')
                    for ID in IDs:
                        if line.query_length >= parameters['len_range'][0] and line.query_length <= parameters['len_range'][1] and ID in ID_dict:
                            if line.has_tag('XY'):
                                len_dist_df.loc[ID, line.query_length] += line.get_tag('XZ')/len(IDs)
                                len_dist_df.loc[ID, 'Total'] += line.get_tag('XZ')/len(IDs)
    return(len_dist_df)

def stat_sRNA_length_by_fraction_with_uniq(parameters, len_dist_df, ID_dict):
    with pysam.AlignmentFile(parameters['input_file'], 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XY') and line.get_tag('XY') == 'U':
                if line.has_tag('XT'):
                    IDs = line.get_tag('XT').split(',')
                    for ID in IDs:
                        if line.query_length >= parameters['len_range'][0] and line.query_length <= parameters['len_range'][1] and ID in ID_dict:
                            if line.has_tag('XY'):
                                len_dist_df.loc[ID, line.query_length] += line.get_tag('XZ')/len(IDs)
                                len_dist_df.loc[ID, 'Total'] += line.get_tag('XZ')/len(IDs)
    return(len_dist_df)

def stat_sRNA_length_no_fraction_with_all(parameters, len_dist_df, ID_dict):
    with pysam.AlignmentFile(parameters['input_file'], 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XT'):
                IDs = line.get_tag('XT').split(',')
                for ID in IDs:
                    if line.query_length >= parameters['len_range'][0] and line.query_length <= parameters['len_range'][1] and ID in ID_dict:
                        if line.has_tag('XY'):
                            len_dist_df.loc[ID, line.query_length] += line.get_tag('XZ')
                            len_dist_df.loc[ID, 'Total'] += line.get_tag('XZ')
    return(len_dist_df)

def stat_sRNA_length_no_fraction_with_multi(parameters, len_dist_df, ID_dict):
    with pysam.AlignmentFile(parameters['input_file'], 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XY') and line.get_tag('XY') != 'U':
                if line.has_tag('XT'):
                    IDs = line.get_tag('XT').split(',')
                    for ID in IDs:
                        if line.query_length >= parameters['len_range'][0] and line.query_length <= parameters['len_range'][1] and ID in ID_dict:
                            if line.has_tag('XY'):
                                len_dist_df.loc[ID, line.query_length] += line.get_tag('XZ')
                                len_dist_df.loc[ID, 'Total'] += line.get_tag('XZ')
    return(len_dist_df)

def stat_sRNA_length_no_fraction_with_uniq(parameters, len_dist_df, ID_dict):
    with pysam.AlignmentFile(parameters['input_file'], 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XY') and line.get_tag('XY') == 'U':
                if line.has_tag('XT'):
                    IDs = line.get_tag('XT').split(',')
                    for ID in IDs:
                        if line.query_length >= parameters['len_range'][0] and line.query_length <= parameters['len_range'][1] and ID in ID_dict:
                            if line.has_tag('XY'):
                                len_dist_df.loc[ID, line.query_length] += line.get_tag('XZ')
                                len_dist_df.loc[ID, 'Total'] += line.get_tag('XZ')
    return(len_dist_df)


def stat_sRNA_length_distribution_from_splited_bam_per_gene(parameters):
    ID_list = pd.read_csv(parameters['ID_file'], header=None, index_col=None, sep=',')[0].tolist()
    len_dist_df, ID_dict = get_len_dist_df(parameters, ID_list)
    if parameters['fraction']:
        if parameters['mapping_type'] == 'all':
            len_dist_df = stat_sRNA_length_by_fraction_with_all(parameters, len_dist_df, ID_dict)
        elif parameters['mapping_type'] == 'multi':
            len_dist_df = stat_sRNA_length_by_fraction_with_multi(parameters, len_dist_df, ID_dict)
        elif parameters['mapping_type'] == 'uniq':
            len_dist_df = stat_sRNA_length_by_fraction_with_uniq(parameters, len_dist_df, ID_dict)
    else:
        if parameters['mapping_type'] == 'all':
            len_dist_df = stat_sRNA_length_no_fraction_with_all(parameters, len_dist_df, ID_dict)
        elif parameters['mapping_type'] == 'multi':
            len_dist_df = stat_sRNA_length_no_fraction_with_multi(parameters, len_dist_df, ID_dict)
        elif parameters['mapping_type'] == 'uniq':
            len_dist_df = stat_sRNA_length_no_fraction_with_uniq(parameters, len_dist_df, ID_dict)
    len_dist_df.index.name = 'ID'
    len_dist_df.to_csv(parameters['output_name'], header=True, index=True, sep='\t')

if __name__ == '__main__':

    parameters = obtainParameter()

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)
    
    ################
    # 导入注释文件中的exon区域：
    stat_sRNA_length_distribution_from_splited_bam_per_gene(parameters)

    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))
