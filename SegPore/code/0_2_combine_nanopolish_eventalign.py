# -*- coding: utf-8 -*-
# @Time    : 2021/12/28 15:34
# @Author  : https://github.com/guangzhaocs/SegPore
# @FileName: combine_nanopolish_eventalign.py
import os
import pandas as pd
import numpy as np
import h5py
from tqdm import tqdm
import argparse
import csv
import shutil


def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])


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


def read_summary(summary_file_name):
    """
    Read the eventalign summary.txt and return a dict.

    The summary.txt contains the follow attributes:
    ['read_index', 'read_name', 'fast5_path', 'model_name', 'strand', 'num_events', 'num_steps',
    'num_skips', 'num_stays', 'total_duration', 'shift', 'scale', 'drift', 'var']

    :param summary_file_name
    :return:
          eventalign_summary_dict:
          - key: read_index
          - value: {'read_index': '0', 'read_name': 'bd76fe67-db4b-4d72-95df-e23f07aea51a', ... ,  'num_stays': '144',
                    'total_duration': '0.36', 'shift': '-3.697', 'scale': '0.949', 'drift': '0.000', 'var': '1.523'}

     read_index     read_name      fast5_path     model_name      strand  num_events      num_steps       num_skips
      num_stays       total_duration  shift   scale   drift   var

    """
    summary_pd = pd.read_csv(summary_file_name, sep='\t')
    print("* summary_pd len: ", len(summary_pd))
    print("* summary_pd read_index len: ", len(list(set(list(summary_pd['read_index'])))))
    print("* summary_pd read_name len: ", len(list(set(list(summary_pd['read_name'])))))

    summary_pd['fast5_path'] = summary_pd['fast5_path'].apply(lambda x: x.split('/')[-1])
    idx2name_dict = summary_pd.set_index('read_index')['read_name'].to_dict()
    idx2path_dict = summary_pd.set_index('read_index')['fast5_path'].to_dict()

    return idx2name_dict, idx2path_dict


def write_to_csv(file_name, row, delimiter=','):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(row)


def write_event_combine_file(read_index, contig, read_data, combine_file_name):

    read_data = np.array(read_data)
    df = pd.DataFrame(read_data, columns=['pos', 'kmer', 'mean', 'start_idx', 'end_idx'])
    df['pos'] = df['pos'].apply(int)
    df['start_idx'] = df['start_idx'].apply(int)
    df['end_idx'] = df['end_idx'].apply(int)
    df['contig'] = contig
    df['read_idx'] = read_index
    df = df[['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']]
    df.to_csv(combine_file_name, index=False, header=False, sep='\t', mode='a+')


def combine_eventalign(eventalign_file_name, combine_file_name, summary_file_name, shift,
                       summary_idx2name_dict, summary_idx2path_dict, ref_dict):
    """
    Combine the same kmers in eventalign file and write to a new file.

    Original file name:  eventalign.txt
    Combined file name:  eventalign_combine.txt

    Example:
    Original file:
    ================================================================================================================
    contig	position	reference_kmer	read_index	strand	event_index	event_level_mean ...	start_idx	end_idx
    ----------------------------------------------------------------------------------------------------------------
    gi|545778205|gb|U00096.3|:c514859-514401	3	ATGGAG	0	t	16538	98.58	...	81407	81411
    gi|545778205|gb|U00096.3|:c514859-514401	3	ATGGAG	0	t	16537	97.60	...	81403	81407
    gi|545778205|gb|U00096.3|:c514859-514401	3	ATGGAG	0	t	16536	104.00	...	81398	81403
    gi|545778205|gb|U00096.3|:c514859-514401	3	ATGGAG	0	t	16535	89.95	...	81392	81398
    ================================================================================================================

    After combining:
    =======================================================================================
    contig	read_index	position	trans_position	kmer	mean	start_idx	end_idx
    ---------------------------------------------------------------------------------------
    gi|545778205|gb|U00096.3|:c514859-514401	0	3	5	ATGGAG	97.53	81392	81411
    =======================================================================================

    trans_position = position + shift(default: 2)
    You can ignore the trans_position if you do not need it.
    """

    # read the header and the first line of the eventalign file.
    f = open(eventalign_file_name, "r")
    header = f.readline().strip().replace('\n', '').replace('\r', '').split('\t')
    column_num = len(header)
    assert header == ['contig', 'position', 'reference_kmer', 'read_index', 'strand', 'event_index', 'event_level_mean',
                      'event_stdv', 'event_length', 'model_kmer', 'model_mean', 'model_stdv', 'standardized_level',
                      'start_idx', 'end_idx', 'map_strand']

    first_line = f.readline().strip().replace('\n', '').replace('\r', '').split('\t')
    assert len(first_line) == column_num

    # parameters for current read
    curent_contig = first_line[0]
    curent_read_index = first_line[3]
    curent_read_strand = first_line[15]
    curent_read_data = list()

    # parameters for current kmer
    curent_kmer = first_line[2]
    curent_position_id = int(first_line[1])
    curent_position_start_idx = int(first_line[13])
    curent_position_end_idx = int(first_line[14])
    curent_kmer_indice_len_list = [abs(curent_position_end_idx - curent_position_start_idx)]
    curent_kmer_mean_list = [float(first_line[6])]
    f.close()

    print('* Start to processing ... ')
    num_lines = wc_count(eventalign_file_name)

    with open(eventalign_file_name, 'r') as f:
        for i, line in enumerate(tqdm(f, total=num_lines)):
            if i > 1:
                line = line.strip().replace('\n', '').replace('\r', '').split('\t')
                if len(line) != column_num:
                    print(i, " -- ", line)
                    continue

                _position_start_idx = int(line[13])
                _position_end_idx = int(line[14])

                # hold current read
                if curent_contig == line[0] and curent_read_index == line[3]:

                    # hold current kmer
                    if curent_kmer == line[2] and curent_position_id == int(line[1]):

                        curent_kmer_indice_len_list.append(abs(_position_end_idx - _position_start_idx))
                        curent_kmer_mean_list.append(float(line[6]))
                        if curent_position_start_idx > _position_start_idx:
                            curent_position_start_idx = _position_start_idx
                        else:
                            curent_position_end_idx = _position_end_idx

                    # a new kmer in this read
                    else:
                        # first add previous kmer to curent_read_data
                        curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                        curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                        assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)

                        _sum_len = sum(curent_kmer_indice_len_list)

                        curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                        _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                        _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                        curent_read_data.append([curent_position_id + shift, curent_kmer,
                                                 _kmer_evnet_mean, curent_position_start_idx, curent_position_end_idx])

                        # then update paras for the new kmer
                        curent_kmer = line[2]
                        curent_position_id = int(line[1])
                        curent_position_start_idx = int(line[13])
                        curent_position_end_idx = int(line[14])
                        curent_kmer_indice_len_list = [abs(_position_end_idx - _position_start_idx)]
                        curent_kmer_mean_list = [float(line[6])]

                # a new read
                else:
                    # first add previous kmer data to the list
                    curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                    curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                    assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)

                    _sum_len = sum(curent_kmer_indice_len_list)

                    curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                    _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                    _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                    curent_read_data.append([curent_position_id + shift, curent_kmer,
                                             _kmer_evnet_mean, curent_position_start_idx, curent_position_end_idx])

                    # then write previous kmer data to a file
                    _end_idx = [curent_read_data[0][-1], curent_read_data[0][-2],
                                curent_read_data[-1][-1], curent_read_data[-1][-2]]
                    if int(curent_read_index) in summary_idx2name_dict:
                        _ref = ref_dict[curent_contig][0][curent_read_data[0][0] - shift: curent_read_data[-1][0] +
                                                                                          shift + 1]
                        combined_kmer = [curent_contig, curent_read_index, curent_read_strand,
                                         curent_read_data[0][0], curent_read_data[-1][0],
                                         min(_end_idx), max(_end_idx),
                                         summary_idx2name_dict[int(curent_read_index)],
                                         summary_idx2path_dict[int(curent_read_index)], _ref]
                        write_to_csv(summary_file_name, combined_kmer)

                        write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                                 combine_file_name)


                    # final update paras for the new kmer
                    curent_contig = line[0]
                    curent_read_index = line[3]
                    curent_read_strand = line[15]
                    curent_read_data = list()

                    curent_kmer = line[2]
                    curent_position_id = int(line[1])
                    curent_position_start_idx = int(line[13])
                    curent_position_end_idx = int(line[14])
                    curent_kmer_indice_len_list = [abs(_position_end_idx - _position_start_idx)]
                    curent_kmer_mean_list = [float(line[6])]

                # the last line
                if i == num_lines - 1:

                    # first for kmer
                    curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                    curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                    assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)
                    # assert sum(curent_kmer_indice_len_list) == (curent_position_end_idx - curent_position_start_idx)
                    _sum_len = sum(curent_kmer_indice_len_list)

                    curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                    _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                    _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                    curent_read_data.append([curent_position_id + shift, curent_kmer,
                                             _kmer_evnet_mean, curent_position_start_idx, curent_position_end_idx])

                    # then for read
                    _end_idx = [curent_read_data[0][-1], curent_read_data[0][-2],
                                curent_read_data[-1][-1], curent_read_data[-1][-2]]
                    if int(curent_read_index) in summary_idx2name_dict:
                        # _ref = ref_dict[curent_contig][0][trans_position_start - shift: trans_position_end + shift+1]
                        _ref = ref_dict[curent_contig][0][curent_read_data[0][0] - shift: curent_read_data[-1][0] +
                                                                                          shift + 1]
                        combined_kmer = [curent_contig, curent_read_index, curent_read_strand,
                                         curent_read_data[0][0], curent_read_data[-1][0],
                                         min(_end_idx), max(_end_idx),
                                         summary_idx2name_dict[int(curent_read_index)],
                                         summary_idx2path_dict[int(curent_read_index)], _ref]
                        write_to_csv(summary_file_name, combined_kmer)

                        write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                                 combine_file_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Combine the same kmers in eventalign file and write to a new file.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    parser.add_argument('--ref-fa', type=str, default='demo.fa')
    parser.add_argument('--shift', type=int, default=2)
    args = parser.parse_args()

    ref_dict = read_transcript_fa(os.path.join(args.root_dir, "0_reference", args.ref_fa))
    print('* Generate ref dict over! ')

    # input file name
    eventalign_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_eventalign_add_strand.txt")
    summary_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_summary.txt")

    # output file name
    eventalign_summary_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_eventalign_summary.csv")
    eventalign_combine_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_eventalign_combined.txt")

    idx2name_dict, idx2path_dict = read_summary(summary_file_name)
    print('* Read summary over! ')

    # --------------------------------------------------------------
    # check output files
    if os.path.exists(eventalign_combine_file_name):
        os.remove(eventalign_combine_file_name)
    combine_eventalign_header = ['read_idx', 'contig', 'pos',  'kmer', 'mean', 'start_idx', 'end_idx']
    write_to_csv(eventalign_combine_file_name, combine_eventalign_header, '\t')

    if os.path.exists(eventalign_summary_file_name):
        os.remove(eventalign_summary_file_name)
    summary_file_header = ['contig', 'read_index', 'strand',
                           'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
                           'read_name', 'fast5_path', 'ref']
    write_to_csv(eventalign_summary_file_name, summary_file_header)

    # --------------------------------------------------------------
    # combine eventalign file
    combine_eventalign(eventalign_file_name, eventalign_combine_file_name, eventalign_summary_file_name,
                       args.shift, idx2name_dict, idx2path_dict, ref_dict)
