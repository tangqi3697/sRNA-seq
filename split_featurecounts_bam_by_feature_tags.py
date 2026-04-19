#! /usr/bin/python

import getopt, sys, datetime, os
import pysam

################
# 获取命令行参数。
def usage():
    print("For example:\n\n   python split_featurecounts_bam_by_feature_tags.py -i input.bam -p priority.txt -o splited-featureCounts \n")
    print("Parameters:");
    print("  -i:  bam annotation with biotype from featureCounts by parameter '-R'.");
    print("  -p:  priority file of features.");
    print("  -o:  output name.");
    print("  -h:  print this page.");
    sys.exit()

def obtainParameter():
    opts, args = getopt.getopt(sys.argv[1:], "hi:p:o:")
    if opts:
        for op, value in opts:
            if op == "-i": input_file = value; print("-i:", value.split('/')[-1])
            elif op == "-p": priority_file = value; print('-p:', value.split('/')[-1])
            elif op == "-o": output_dir = value; print('-o:', value.split('/')[-1])
            elif op == "-h": usage()
    else: usage()
    return(input_file, priority_file, output_dir)

def build_priority_dict(priority_file):
    priority_dict, priority_no = {float('inf'): 'Others'}, 1
    with open(priority_file, 'r') as priority_f:
        for line in priority_f:
            priority = line.strip()
            if priority not in priority_dict:
                priority_dict[priority] = priority_no
                priority_dict[priority_no] = priority
                priority_no += 1
    return(priority_dict)


def write_bam_into_tmp_by_XT_tag(input_file, tmp_file, priority_dict):
    with pysam.AlignmentFile(input_file, 'r') as input_f:
        with pysam.AlignmentFile(tmp_file, 'wb', header=input_f.header) as tmp_f:
            for line in input_f:
                priority_no = float('inf')
                if line.has_tag('XT'):
                    features = line.get_tag('XT')
                    for feature in features.split(','):
                        if feature in priority_dict:
                            if priority_dict[feature] < priority_no:
                                priority_no = priority_dict[feature]
                        # else:
                        #     print('{} is not in priority file.'.format(feature))
                    line.set_tag('XT', priority_dict[priority_no], value_type='Z', replace=True)
                    tmp_f.write(line)

def split_tmp_file_by_XT_tag(priority_file, tmp_file, output_dir):
    with open(priority_file, 'r') as priority_f:
        for line in priority_f:
            feature_name = line.strip()
            feature_file = output_dir + '.' + feature_name + '.bam'
            with pysam.AlignmentFile(tmp_file, 'r') as tmp_f:
                with pysam.AlignmentFile(feature_file, 'wb', header=tmp_f.header) as feature_f:
                    for line in tmp_f:
                        if line.has_tag('XT'):
                            if feature_name == line.get_tag('XT'):
                                mapping_index = -1
                                for m in range(len(line.get_tags())):
                                    if line.get_tags()[m][0] == 'XS':
                                        mapping_index = m
                                line.set_tags(line.get_tags()[:mapping_index])
                                feature_f.write(line)

def split_featurecounts_bam_by_feature_tags(input_file, priority_file, output_dir):
    priority_dict = build_priority_dict(priority_file)
    tmp_file = input_file + '.tmp'
    write_bam_into_tmp_by_XT_tag(input_file, tmp_file, priority_dict)
    split_tmp_file_by_XT_tag(priority_file, tmp_file, output_dir)

    os.remove(tmp_file)


if __name__ == '__main__':

    input_file, priority_file, output_dir = obtainParameter()

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)
    
    ################
    # 导入注释文件中的exon区域：
    split_featurecounts_bam_by_feature_tags(input_file, priority_file, output_dir)

    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))
