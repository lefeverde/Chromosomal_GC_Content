# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:34:56 2016

@author: daniellefever
"""

#!/usr/bin/python

from __future__ import division



import argparse
import os
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
    chrom_seqrecs = []
    for seqrec in SeqIO.parse(open(args.i), 'gb'):
        gc_writer(seqrec)


def gc_writer(parsed_gb, feature_string = 'CDS', consec_bases = 750):
    '''
    This function finds all the features delineated by the feature_string.
    The defualt is for CDS, and four regions with default size 750 bp
    flanking up/down stream of start and stop sequences are extracted.

    The output is a csv file with the gene name, id, the GC content
    of the 4 regions, and the actual sequences of those regions.

    '''
    chromo_id = parsed_gb.id
    cur_handle = str(chromo_id) + '.csv'

    with open(cur_handle, 'w+') as f:
        f.write(
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

        cds_seqrec_list = feature_getter(parsed_gb, feature_string)
        for feature_seqrec in cds_seqrec_list:
            gene_name = feature_seqrec.qualifiers['gene'][0]
            gene_id = gene_id_finder(feature_seqrec)

            '''
            try:
                gene_id = feature_seqrec.qualifiers['db_xref'][2].split(':')[1]
                if 'HGNC' in gene_id:
                    continue
            except IndexError:
                continue
            '''


            #gi_num = feature_seqrec.qualifiers['db_xref'][0].split(':')[1]
            start_upseq, start_downseq, end_upseq, end_downseq = \
            seq_streams(parsed_gb, feature_seqrec, consec_bases)
            start_up_gc = GC(start_upseq.seq)
            start_down_gc = GC(start_downseq.seq)
            end_up_gc = GC(end_upseq.seq)
            end_down_gc = GC(end_downseq.seq)

            f.write(
                str(gene_name) + ',' +
                str(gene_id) + ',' +
                #str(gi_num) + ',' +
                str(start_up_gc) + ',' +
                str(start_down_gc) + ',' +
                str(end_up_gc) + ',' +
                str(end_down_gc) + ',' +
                str(start_upseq.seq) + ',' +
                str(start_downseq.seq) + ',' +
                str(end_upseq.seq) + ',' +
                str(end_downseq.seq) + '\n'  #+ ',' +
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


if __name__ == '__main__':
    main()
