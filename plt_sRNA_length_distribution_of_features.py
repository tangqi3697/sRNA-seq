#! /usr/bin/bash

import getopt, sys, datetime
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def usage():
    print("For example:\n\n    python plt_sRNA_length_distribution_of_features.py.py -i features.length_distribution.tsv -p sRNA_source_feature.priority.txt -o features.length_distribution.pdf \n")
    print("Parameters:");
    print("  -i:  length_distribution file from stat_sRNA_length_distribution_from_splited_bam.py.");
    print("  -p:  priority file of features.");
    print("  -o:  the figure name.");
    print("  -h:  print this page.");
    sys.exit()

def obtainParameter():
    opts, args = getopt.getopt(sys.argv[1:], "hi:p:o:")
    parameters = {}
    if opts:
        print('Input Parameters:')
        for op, value in opts:
            if op == "-i":   parameters['input_file'] = value; print(op, ":", value)
            elif op == "-p": parameters['priority_file'] = value; print(op, ":", value)
            elif op == "-o": parameters['output_file'] = value; print(op, ":", value)
            elif op == "-h": usage()
    else: usage()
    return(parameters)


def plt_sRNA_length_distribution_of_features(parameters):
    priority_list = pd.read_csv(parameters['priority_file'], header=None, index_col=None, sep='\t')[0].tolist()
    df = pd.read_csv(parameters['input_file'], header=0, index_col=None, sep='\t')

    sns.axes_style('white')
    for value in ['RPM', 'Count']:
        for mapping_type in df['type'].unique():
            plt.close()
            curr_df = df[df['type']==mapping_type][['length', 'feature', value]]
            sns.barplot(x='length', y=value, hue='feature', data=curr_df, hue_order=priority_list, palette='tab20c')
            plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0, fontsize='x-small')
            plt.title('{} of {} mapping sRNA'.format(value, mapping_type))
            output_file1 = '.'.join(parameters['output_file'].split('.')[:-1] + [mapping_type, value] +  [parameters['output_file'].split('.')[-1]])
            plt.savefig(output_file1, dpi=600, bbox_inches='tight')
            ############
            curr_df_stack = curr_df.pivot(index='length',columns='feature', values=value)[priority_list]
            plt.close()
            curr_df_stack.plot(kind='bar', stacked=True, colormap='tab20c')
            plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0, fontsize='x-small')
            plt.ylabel(value)
            plt.title('{} of {} mapping sRNA'.format(value, mapping_type))
            output_file1 = '.'.join(parameters['output_file'].split('.')[:-1] + [mapping_type, value, 'stack'] +  [parameters['output_file'].split('.')[-1]])
            plt.savefig(output_file1, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    parameters = obtainParameter()

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)

    plt_sRNA_length_distribution_of_features(parameters)

    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))
