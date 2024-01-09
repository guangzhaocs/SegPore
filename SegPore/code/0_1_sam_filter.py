# -*- coding: UTF-8 -*-
"""
@Project : cuda_check_after_bonito 
@File    : mimimap2_test.py
@Author  : Guangzhao Cheng
@Date    : 2022/11/22 15:37 
@Github  : 
"""

import csv
import os
import argparse
import pysam
import pandas as pd


def generate_trans_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    trans2idx_dict = data.set_index('transcript')['index'].to_dict()
    idx2trans_dict = data.set_index('index')['transcript'].to_dict()
    return trans2idx_dict, idx2trans_dict


def write_to_csv(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def filter_sam(sam_file, keep_file_name, del_file_name, trans2idx_dict):

    with pysam.Samfile(sam_file, 'r') as sf:
        for read in sf:
            if read.flag == 0 and read.query_name and read.reference_name:
                write_to_csv(keep_file_name, [read.query_name, trans2idx_dict[read.reference_name], read.flag,
                                              read.reference_start, read.reference_end])
            elif read.flag != 0 and read.query_name and read.reference_name:
                write_to_csv(del_file_name, [read.query_name, trans2idx_dict[read.reference_name], read.flag,
                                             read.reference_start, read.reference_end])
            else:
                pass


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--root-dir', type=str, default='/scratch/cs/infantbiome/chengg1/base_calling/dataset')
    parser.add_argument('--dataset', type=str, default='dinopore_h9_KO')
    parser.add_argument('--sample', type=str, default='FAK06641')
    parser.add_argument('--ref-fa', type=str, default='Homo_sapiens.GRCh37.cdna.all.fa')
    args = parser.parse_args()

    trans2idx_dict, idx2trans_dict = generate_trans_dict(os.path.join(args.root_dir, args.dataset, "0_reference",
                                                                      args.ref_fa + ".t2idx.csv"))

    sam_file_name = os.path.join(args.root_dir, args.dataset, "1_fast5", args.sample, "reads-ref.sorted.sam")
    keep_file_name = os.path.join(args.root_dir, args.dataset, "1_fast5", args.sample, "reads-ref.filter.keep.csv")
    del_file_name = os.path.join(args.root_dir, args.dataset, "1_fast5", args.sample, "reads-ref.filter.del.csv")
    if os.path.exists(keep_file_name):
        os.remove(keep_file_name)
    if os.path.exists(del_file_name):
        os.remove(del_file_name)

    print('*', args.dataset, args.sample)
    write_to_csv(keep_file_name, ["read_name", "ref_idx", "flag", "reference_start", "reference_end"])
    write_to_csv(del_file_name, ["read_name", "ref_idx", "flag", "reference_start", "reference_end"])
    filter_sam(sam_file_name, keep_file_name, del_file_name, trans2idx_dict)