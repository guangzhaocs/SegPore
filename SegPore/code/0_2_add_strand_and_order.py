# -*- coding: utf-8 -*-
# @Time    : 2021/12/28 15:34
# @Author  : Guangzhao Cheng
# @FileName: combine_eventalign.py
import os
import pandas as pd
import numpy as np
import h5py
from tqdm import tqdm
import argparse
import csv
import shutil

base_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "X": "X"}


def reverse_com(kmer):
    kmer = kmer[::-1]
    com_kmer = "".join(base_dict[base] for base in kmer)
    return com_kmer


def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])


def generate_kmer_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    kmer2idx_dict = data.set_index('kmer')['index'].to_dict()
    idx2kmer_dict = data.set_index('index')['kmer'].to_dict()
    return kmer2idx_dict, idx2kmer_dict


def generate_trans_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    trans2idx_dict = data.set_index('transcript')['index'].to_dict()
    idx2trans_dict = data.set_index('index')['transcript'].to_dict()
    return trans2idx_dict, idx2trans_dict


def read_transcript_fa(transcript_fasta_file_name):
    fasta = open(transcript_fasta_file_name, "r")
    entries, separate_by_pipe="", False
    for ln in fasta:
        entries += ln
    entries = entries.split(">")
    dict = {}
    for entry in entries:
        entry = entry.split("\n")
        if len(entry[0].split()) > 0:
            id = entry[0].split(' ')[0]
            seq = "".join(entry[1:])
            dict[id] = [seq]
    return dict


def write_to_csv(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(row)


def write_event_combine_file(read_data, read_strand, read_order, eventalign_add_file):
    if read_strand == '+':
        assert read_order == 1

    if read_strand == '-':
        assert read_order == 0
        return
    
    read_data = np.array(read_data)
    df = pd.DataFrame(read_data, columns=['contig', 'position', 'reference_kmer', 'read_index', 'strand', 'event_index',
                                          'event_level_mean',  'event_stdv', 'event_length', 'model_kmer', 'model_mean',
                                          'model_stdv', 'standardized_level', 'start_idx', 'end_idx'])
    df['position'] = df['position'].apply(int)
    df['read_index'] = df['read_index'].apply(int)
    df['event_index'] = df['event_index'].apply(int)
    df['start_idx'] = df['start_idx'].apply(int)
    df['end_idx'] = df['end_idx'].apply(int)
    df['map_strand'] = read_strand
    df.to_csv(eventalign_add_file, index=False, header=False, sep='\t', mode='a+')


def add_strand_eventalign(eventalign_file_name, eventalign_add_strand_file_name):

    # read the header and the first line of the eventalign file.
    f = open(eventalign_file_name, "r")
    header = f.readline().strip().replace('\n', '').replace('\r', '').split('\t')
    column_num = len(header)
    assert header == ['contig', 'position', 'reference_kmer', 'read_index', 'strand', 'event_index', 'event_level_mean',
                      'event_stdv', 'event_length', 'model_kmer', 'model_mean', 'model_stdv', 'standardized_level',
                      'start_idx', 'end_idx']

    first_line = f.readline().strip().replace('\n', '').replace('\r', '').split('\t')
    assert len(first_line) == column_num

    # parameters for current read
    curent_read_index = first_line[3]
    curent_read_start_pos = int(first_line[-2])
    curent_read_order = -1
    if first_line[9] != "NNNNN":
        if first_line[2] == first_line[9]:
            current_read_strand = "+"
        else:
            current_read_strand = "-"
    else:
        current_read_strand = "*"
    curent_read_data = [first_line]
    f.close()

    print('* Start to processing ... ')
    num_lines = wc_count(eventalign_file_name)

    with open(eventalign_file_name, 'r') as f:
        for i, line in enumerate(tqdm(f, total=num_lines)):
            if i > 1:
                line = line.strip().replace('\n', '').replace('\r', '').split('\t')
                if len(line) != column_num:
                    print(f"line {i+1} error. ")
                    continue

                _index = line[3]
                _ref_kmer = line[2]
                _modle_kmer = line[9]
                _start_idx = int(line[-2])

                # hold current read
                if curent_read_index == line[3]:

                    # check strand
                    if _modle_kmer != "NNNNN":
                        if current_read_strand == "*":
                            if _ref_kmer == _modle_kmer:
                                current_read_strand = "+"
                            else:
                                current_read_strand = "-"
                        elif current_read_strand == "+":
                            assert _ref_kmer == _modle_kmer
                        else:
                            assert _ref_kmer == reverse_com(_modle_kmer)

                    # check order
                    if curent_read_order == -1:
                        if curent_read_start_pos > _start_idx:
                            curent_read_order = 1
                        else:
                            curent_read_order = 0
                    else:
                        if curent_read_order == 1:
                            assert curent_read_start_pos > _start_idx
                        else:
                            assert curent_read_start_pos < _start_idx
                    curent_read_start_pos = _start_idx
                    curent_read_data.append(line)

                # a new read
                else:
                    write_event_combine_file(curent_read_data, current_read_strand, curent_read_order,
                                             eventalign_add_strand_file_name)

                    curent_read_index = line[3]
                    curent_read_start_pos = int(line[-2])
                    curent_read_order = -1
                    if line[9] != "NNNNN":
                        if line[2] == line[9]:
                            current_read_strand = "+"
                        else:
                            current_read_strand = "-"
                    else:
                        current_read_strand = "*"
                    curent_read_data = [line]

                # the last line
                if i == num_lines - 1:
                    write_event_combine_file(curent_read_data, current_read_strand, curent_read_order,
                                             eventalign_add_strand_file_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Combine the same kmers in eventalign file and write to a new file.')
    parser.add_argument('--root-dir', type=str, default='/scratch/cs/infantbiome/chengg1/base_calling/dataset')
    parser.add_argument('--dataset', type=str, default='dinopore_h9_KO')
    parser.add_argument('--sample', type=str, default='FAK06641')
    args = parser.parse_args()

    print('*', args.dataset, args.sample)

    # input file name
    root_dir = os.path.join(args.root_dir, args.dataset, "1_fast5", args.sample)
    eventalign_file_name = os.path.join(root_dir,  args.sample + "-eventalign.txt")

    # output file name
    eventalign_add_strand_file_name = os.path.join(root_dir, args.sample + "-eventalign-add-strand.txt")
    if os.path.exists(eventalign_add_strand_file_name):
        os.remove(eventalign_add_strand_file_name)

    new_header = ['contig', 'position', 'reference_kmer', 'read_index', 'strand', 'event_index',
                  'event_level_mean', 'event_stdv', 'event_length', 'model_kmer', 'model_mean',
                  'model_stdv', 'standardized_level', 'start_idx', 'end_idx', 'map_strand']
    write_to_csv(eventalign_add_strand_file_name, new_header)

    add_strand_eventalign(eventalign_file_name, eventalign_add_strand_file_name)

    print("* original line = ", wc_count(eventalign_file_name))
    print("* add file line = ", wc_count(eventalign_add_strand_file_name))