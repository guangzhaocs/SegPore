# -*- coding: utf-8 -*-
# @Time    : 2021/12/28 15:34
# @Author  : Guangzhao Cheng
# @FileName: combine_eventalign.py
import os
import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import h5py


def read_kmer_dict(kmer_file_name):
    kmer_pd = pd.read_csv(kmer_file_name)
    ref_idx_list = set(list(kmer_pd['ref_idx']))
    kmer_dict = dict()
    for _ref_idx in ref_idx_list:
        _kmer_pd = kmer_pd[kmer_pd['ref_idx'] == _ref_idx]
        _ref_idx_dict = dict()
        for index, row in _kmer_pd.iterrows():
            _ref_idx_dict[row['pos']] = {'kmer': row['kmer'], 'kmer_idx': row['kmer_idx']}
        kmer_dict[_ref_idx] = _ref_idx_dict
    return kmer_dict


def write_to_file(file_name, data):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(data)


def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])


def split_eventalign(eventalign_file_name, save_dir, kmer_dict):
    """
    eventalign_file:
               0         1        2      3       4          5           6
         ['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']
    """
    print('* Start to processing ... ')
    num_lines = wc_count(eventalign_file_name)

    with open(eventalign_file_name, 'r') as f:
        for i, line in enumerate(tqdm(f, total=num_lines)):
            line = line.strip().replace('\n', '').replace('\r', '').split('\t')
            _read_idx = int(line[0])
            _contig = int(line[1])
            _pos = int(line[2])
            _kmer = int(line[3])
            _mean = float(line[4])
            _start_idx = int(line[5])
            _end_idx = int(line[6])

            _kmer_dict = kmer_dict[_contig]
            if _pos in _kmer_dict:
                assert _kmer == _kmer_dict[_pos]['kmer_idx']
                _save_file = os.path.join(save_dir, str(_contig) + "_" + str(_pos) + "_" + _kmer_dict[_pos]['kmer'] +
                                          '.csv')
                write_to_file(_save_file, [_read_idx, _contig, _pos, _kmer, _mean, _start_idx, _end_idx])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Combine the same kmers in eventalign file and write to a new file.')
    parser.add_argument('--root-dir', type=str, default='/scratch/cs/nanopore/chengg1/base_calling/dataset')
    parser.add_argument('--dataset', type=str, default='cell_genomics')
    parser.add_argument('--sample', type=str, default='1623_native')
    parser.add_argument('--ref-fa', type=str, default='ec.fasta')
    args = parser.parse_args()

    ref_fa = os.path.join(args.root_dir, args.dataset, "0_reference", args.ref_fa)

    # read kmer list
    kmer_file_name = os.path.join(args.root_dir, args.dataset, "0_reference", args.ref_fa + '.kmer.csv')
    kmer_dict = read_kmer_dict(kmer_file_name)

    # read file name
    eventalign_combine_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                                args.sample + "-eventalign-combined.txt")

    # save folder name
    save_dir = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample, "kmer")

    print(' *', args.dataset, args.sample)

    # ----------------------------------------------------------------------------------------
    # step 2: split eventalign-combined.txt into kmer dir
    # ----------------------------------------------------------------------------------------
    if os.path.exists(save_dir):
        import shutil
        shutil.rmtree(save_dir)
    os.mkdir(save_dir)
    split_eventalign(eventalign_combine_file_name, save_dir, kmer_dict)
    print(" * split over !")