#! /usr/bin/python

import getopt, sys, datetime, os
import pandas as pd
import pysam

################
# 获取命令行参数。
def usage():
    print("For example:\n\n   python stat_sRNA_length_distribution_from_splited_bam.py -i input.bam -l 18,30 -p priority.txt -o sRNA_length_distribution.tsv \n")
    print("Parameters:");
    print("  -i:  output directory of split_featurecounts_bam_by_feature_tags.py.");
    print("  -l:  length range (default: 18,30).");
    print("  -p:  priority file of features.");
    print("  -t:  mapping type (default: all,multi,uniq).");
    print("  -o:  output name.");
    print("  -h:  print this page.");
    sys.exit()

def obtainParameter():
    opts, args = getopt.getopt(sys.argv[1:], "hi:l:p:t:o:")
    
    parameters = {
        'len_range':[18, 30],
        'mapping_type':'all,multi,uniq'.split(','),
    }
    if opts:
        for op, value in opts:
            if op == "-i": parameters['input_file'] = value; print(op, value.split('/')[-1])
            elif op == "-l": parameters['len_range'] = [int(x) for x in value.split(',')]; print('-s:', value.split('/')[-1])
            elif op == "-p": parameters['priority_file'] = value; print(op, value.split('/')[-1])
            elif op == "-t": parameters['mapping_type'] = value.split(','); print(op, value.split('/')[-1])
            elif op == "-o": parameters['output_name'] = value; print(op, value.split('/')[-1])
            elif op == "-h": usage()
    else: usage()
    return(parameters)


def get_curr_file_sRNA_len_distribution(curr_file):
    len_dict = {'all':{}, 'multi':{}, 'uniq':{}}
    with pysam.AlignmentFile(curr_file, 'r') as curr_f:
        for line in curr_f:
            if line.has_tag('XZ'):
                if line.query_length not in len_dict['all']:
                    len_dict['all'][line.query_length] = 0
                    len_dict['multi'][line.query_length] = 0
                    len_dict['uniq'][line.query_length] = 0
                if line.has_tag('XY'):
                    len_dict['all'][line.query_length] += line.get_tag('XZ')
                    if line.get_tag('XY') == 'U':
                        len_dict['uniq'][line.query_length] += line.get_tag('XZ')
                    else:
                        len_dict['multi'][line.query_length] += line.get_tag('XZ')

    type_list, len_list, count_list = [], [], []
    for mapping_type in len_dict:
        for key,value in len_dict[mapping_type].items():
            type_list.append(mapping_type)
            len_list.append(key)
            count_list.append(value)

    len_df = pd.DataFrame({'type':type_list,'length':len_list, 'Count':count_list})
    return(len_df)

def stat_sRNA_length_distribution_from_splited_bam(parameters):
    priority_list = pd.read_csv(parameters['priority_file'], header=None, index_col=None, sep=',')[0].tolist()
    len_dist_df = pd.DataFrame(columns=['type', 'length', 'Count'])
    for priority in priority_list:
        curr_file = parameters['input_file'] + '.' + priority + '.bam'
        len_df = get_curr_file_sRNA_len_distribution(curr_file)
        len_df['feature'] = priority
        for mapping_type in parameters['mapping_type']:
            template_df = pd.DataFrame({'len':range(parameters['len_range'][0], (parameters['len_range'][1] + 1))})
            mapping_type_len_dist_df = len_df[(len_df['type']==mapping_type)]
            template_df = pd.merge(template_df, mapping_type_len_dist_df, left_on='len', right_on='length', how='left')
            del template_df['len']
            len_dist_df = pd.concat([len_dist_df, template_df])
    

    len_dist_RPM_df = pd.DataFrame(columns=['type', 'length', 'Count', 'RPM'])
    for mapping_type in parameters['mapping_type']:
        mapping_type_len_dist_df = len_dist_df[(len_dist_df['type']==mapping_type)]
        mapping_type_len_dist_df['RPM'] = round((10**6 * mapping_type_len_dist_df['Count']) / mapping_type_len_dist_df['Count'].sum(), 3)
        len_dist_RPM_df = pd.concat([len_dist_RPM_df, mapping_type_len_dist_df])

    len_dist_RPM_df['length'] = len_dist_RPM_df['length'].astype('int')
    len_dist_RPM_df.to_csv(parameters['output_name'], header=True, index=False, sep='\t')

if __name__ == '__main__':

    parameters = obtainParameter()

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)
    
    ################
    # 导入注释文件中的exon区域：
    stat_sRNA_length_distribution_from_splited_bam(parameters)

    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))
