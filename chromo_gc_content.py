# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:34:56 2016

@author: daniellefever
"""

#!/usr/bin/python

from __future__ import division



import argparse
import os
import sys
import csv

from Bio import SeqIO
from Bio.SeqUtils import GC
#from util_functions import *

parser = argparse.ArgumentParser()
req = parser.add_argument_group('required arguments')

req.add_argument('-i',
    help = 'big gbk file',
    metavar = '',
    required=True)

args = parser.parse_args()


def main():
    header = (
    'chromosome_id,' +
    'gene_name,' +
    'gene_id,' +
   # 'gi_num,' +
    'start_up_gc,' +
    'start_down_gc,' +
    'end_up_gc,' +
    'end_down_gc,' +
    'start_up_seq,' +
    'start_down_seq,' +
    'end_up_seq,' +
    'end_down_seq'
    )

    file_dir =  os.path.dirname(os.path.abspath(args.i))

    raw_dir = os.path.join(file_dir, 'raw_csvs')
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)

    unique_dir = os.path.join(file_dir, 'unique_csvs')
    if not os.path.exists(unique_dir):
        os.makedirs(unique_dir)

    #chrom_seqrecs = []
    for seqrec in SeqIO.parse(open(args.i), 'gb'):
        if seqrec.id.split('_')[0] == 'NC':
            gc_writer(seqrec, raw_dir, header)

    unique_cds_maker(raw_dir, unique_dir, header)
    concat_list = [i.split(',') for i in concat_list_maker(unique_dir, header)]
    with open('concat_file.csv', 'w+') as cf:
        cf_writer = csv.writer(cf)
        #cf_writer.writerows((concat_list.split(',')))
        cf_writer.writerows((concat_list))

def unique_cds_maker(raw_dir, unique_dir, header):
    header = header + '\n'
    for raw_cnt_csv in os.listdir(raw_dir):
        unique_cds_set = set()
        with open(os.path.join(raw_dir,raw_cnt_csv)) as cur_f:
            next(cur_f)
            for line in cur_f:
                unique_cds_set.add(line)


        unique_list = [(i.strip('\r\n').split(',')) for i in unique_cds_set]
        cur_file_handle = os.path.join(unique_dir, (raw_cnt_csv.split('.')[0] + '_unique_feats.csv'))

        with open(cur_file_handle, 'w+') as ind_unique_csv:
            ind_unique_csv.write(header)

            #ind_csv_writer = csv.writer(ind_unique_csv, lineterminator='\n')
            ind_csv_writer = csv.writer(ind_unique_csv)

            ind_csv_writer.writerows(unique_list)







def gc_writer(parsed_gb, out_dir, header, feature_string = 'CDS', consec_bases = 750):
    '''
    This function finds all the features delineated by the feature_string.
    The defualt is for CDS, and four regions with default size 750 bp
    flanking up/down stream of start and stop sequences are extracted.

    The output is a csv file with the gene name, id, the GC content
    of the 4 regions, and the actual sequences of those regions.

    '''
    chromo_id = parsed_gb.id
    cur_handle = os.path.join(out_dir, (str(chromo_id) + '.csv'))
    str_seq = str(parsed_gb.seq)
    with open(cur_handle, 'w+') as f:
        f.write(header)

        cds_seqrec_list = feature_getter(parsed_gb, feature_string)
        for feature_num, feature_seqrec in enumerate(cds_seqrec_list):
            if feature_num % 1000 == 0:
                sys.stdout.write((chromo_id + ' ') + '{0}\r'.format(feature_num),)
                sys.stdout.flush()
            gene_name = feature_seqrec.qualifiers['gene'][0]
            gene_id = gene_id_finder(feature_seqrec)
            start_upseq, start_downseq, end_upseq, end_downseq = \
            seq_streams(str_seq, feature_seqrec, consec_bases)
            start_up_gc = GC(start_upseq)
            start_down_gc = GC(start_downseq)
            end_up_gc = GC(end_upseq)
            end_down_gc = GC(end_downseq)

            f.write(
                str(chromo_id) + ',' +
                str(gene_name) + ',' +
                str(gene_id) + ',' +
                #str(gi_num) + ',' +
                str(start_up_gc) + ',' +
                str(start_down_gc) + ',' +
                str(end_up_gc) + ',' +
                str(end_down_gc) + ',' +
                str(start_upseq) + ',' +
                str(start_downseq) + ',' +
                str(end_upseq) + ',' +
                str(end_downseq) + '\n'  #+ ',' +
                )


def gene_id_finder(feature_seqrec):
    for db_xref_item in feature_seqrec.qualifiers['db_xref']:
        if 'GeneID' in db_xref_item:
            gene_id_num = db_xref_item.split(':')[1]
    return int(gene_id_num)

def full_fp_list(input_dir):
    file_list = []
    for file in os.listdir(input_dir):
        full_fp = os.path.join(input_dir, file)
        file_list.append(full_fp)
    return file_list

def seq_streams(parsed_gb, feature_seqrec, consec_bases = 750):
    '''
    This function takes returns the 4 flanking regions
    on the up & down stream regions flanking the start and
    stop positions of a given feature. The default len of
    flanking regions is 750 bases.

    '''
    start_chrome_ind = feature_seqrec.location.start.position
    end_chrome_ind = feature_seqrec.location.end.position

    start_upseq = parsed_gb[(start_chrome_ind - consec_bases):start_chrome_ind]
    start_downseq = parsed_gb[(start_chrome_ind + 3):(start_chrome_ind + consec_bases + 3)]
    end_upseq = parsed_gb[(end_chrome_ind - consec_bases):end_chrome_ind]
    end_downseq = parsed_gb[(end_chrome_ind + 3):(end_chrome_ind + consec_bases + 3)]

    return(start_upseq, start_downseq, end_upseq, end_downseq)

def feature_getter(parsed_gb, feature_string='CDS', spec_gene=None):
    cds_seqs = []
    for index, cur_feature in enumerate(parsed_gb.features):
        if cur_feature.type == feature_string:
            if spec_gene is None:
                extracted_feature = parsed_gb.features[index]
                feature_seqrec = extracted_feature.extract(parsed_gb.seq)
                #cds_seqs.append((extracted_feature,feature_seqrec))
                cds_seqs.append(extracted_feature)
            if spec_gene is not None:
                gene_id = cur_feature.qualifiers['gene']
                if spec_gene in gene_id[0]:
                    extracted_feature = parsed_gb.features[index]
                    feature_seqrec = extracted_feature.extract(parsed_gb.seq)
                    cds_seqs.append(extracted_feature)
                    #cds_seqs.append((extracted_feature,feature_seqrec))
    return(cds_seqs)

def unique_cds_getter(raw_count_path):
    with open('concat_unique_cds.csv', 'w+') as concat_file:
        concat_file.write(
            'chromosome_id,' +
            'gene_name,' +
            'gene_id,' +
           # 'gi_num,' +
            'start_up_gc,' +
            'start_down_gc,' +
            'end_up_gc,' +
            'end_down_gc,' +
            'start_up_seq,' +
            'start_down_seq,' +
            'end_up_seq,' +
            'end_down_seq,' +
            '\n'
        )

        for cur_csv in os.listdir(raw_count_path):
            #print cur_csv
            if cur_csv.split('.')[-1] == 'csv':
                cur_file_list = [line for line in open(os.path.join(raw_count_path,cur_csv), 'rU')]
                unique_cds = list(set(cur_file_list))
                #return unique_cds
                #unique_cds_with_chromo_id = [cur_csv.split('.')[0] + i for i in unique_cds]
                cur_file_handle = (cur_csv.split('.')[0] + '_unique_feats.csv')
                #return unique_cds
                with open(cur_file_handle, 'w+') as ind_unique_csv:
                    ind_csv_writer = csv.writer(ind_unique_csv)
                    ind_csv_writer.writerows((unique_cds.strip('\r\n').split(',')))


def concat_list_maker(unique_cds_dir, header):
    concat_list = []
    concat_list.append(header)

    for cur_csv in os.listdir(unique_cds_dir):
        if cur_csv.split('.')[-1] == 'csv':
            #cur_csv_list = [i.strip('\r') for i in open(os.path.join(unique_cds_dir,cur_csv))]
            cur_csv_list = [i.strip('\r\n') for i in open(os.path.join(unique_cds_dir,cur_csv))]
            #cur_csv_list = [i for i in open(os.path.join(unique_cds_dir,cur_csv))]
        try:
            for cur_line in cur_csv_list[1:]:
                concat_list.append(cur_line)
                #concat_list.append((cur_line.strip('\"\"')))
                #concat_list.append((cur_line.strip('\n')))
        except UnboundLocalError:
            continue

    return concat_list









if __name__ == '__main__':
    main()
