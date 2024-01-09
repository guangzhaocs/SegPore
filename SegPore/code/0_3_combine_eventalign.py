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
    header = None
    column_num = 0
    eventalign_summary_dict = dict()
    error_reads_list = list()

    with open(summary_file_name, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip().replace('\n', '').replace('\r', '').split('\t')
            if i == 0:
                column_num = len(line)
                header = line[0:]
            else:
                if len(line) == column_num:
                    # assert len(line) == column_num, f'Error in line {i} !'
                    read = zip(header, line[0:])
                    eventalign_summary_dict[int(line[0])] = dict(read)
                else:
                    print(line)
                    error_reads_list.append(int(line[0]))

    return eventalign_summary_dict, error_reads_list


def write_to_csv(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


# def read_and_write_multi_fast5(summary_dict, read_index, contig, event_arr, original_fast5_dir):
#
#     read_index = int(read_index)
#     fast5_file_name = summary_dict[read_index]['fast5_path']
#     fast5_file_name = fast5_file_name.split('/')[-1]
#     fast5_file_name = os.path.join(original_fast5_dir, fast5_file_name)
#     read_name = summary_dict[read_index]['read_name']
#
#     with h5py.File(fast5_file_name, 'r+') as fast5_data:
#        for read_key in fast5_data:
#            read_name_ = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
#            if read_name_ == read_name:
#                if "Analyese" in fast5_data[read_key]:
#                    del fast5_data[read_key]["Analyese"]
#                event_data = fast5_data[read_key].create_group("Analyese")
#                event_data['Events'] = event_arr
#                event_data.attrs.create("read_index", read_index, shape=None, dtype=None)
#                event_data.attrs.create("contig", contig, shape=None, dtype=None)
#                return


def write_event_combine_file(read_index, contig, read_data, combine_file_name, trans2idx_dict):

    read_data = np.array(read_data)
    df = pd.DataFrame(read_data, columns=['pos', 'kmer', 'mean', 'start_idx', 'end_idx'])
    df['pos'] = df['pos'].apply(int)
    df['kmer'] = df['kmer'].apply(int)
    df['start_idx'] = df['start_idx'].apply(int)
    df['end_idx'] = df['end_idx'].apply(int)
    df['contig'] = trans2idx_dict[contig]
    df['read_idx'] = read_index
    df = df[['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']]
    df.to_csv(combine_file_name, index=False, header=False, sep='\t', mode='a+')

# combine_eventalign(eventalign_file_name, eventalign_combine_file_name, eventalign_summary_file_name,
#                        args.shift, summary_dict, kmer2idx_dict, trans2idx_dict, error_reads_list, ref_dict)
def combine_eventalign(eventalign_file_name, combine_file_name, summary_file_name, shift,
                       summary_dict, kmer2idx_dict, trans2idx_dict, ref_dict):
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
                      'start_idx', 'end_idx']

    first_line = f.readline().strip().replace('\n', '').replace('\r', '').split('\t')
    assert len(first_line) == column_num

    # parameters for current read
    curent_contig = first_line[0]
    curent_read_index = first_line[3]
    curent_read_strand = first_line[4]
    curent_read_data = list()
    curent_read_order = -1

    # parameters for current kmer
    curent_kmer = first_line[2]
    curent_position_id = int(first_line[1])
    curent_position_start_idx = int(first_line[-2])
    curent_position_end_idx = int(first_line[-1])
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

                _position_start_idx = int(line[-2])
                _position_end_idx = int(line[-1])

                # hold current read
                if curent_contig == line[0] and curent_read_index == line[3]:

                    # check order [order = 1: reverse]
                    if curent_read_order == -1:
                        if curent_position_start_idx > _position_start_idx:
                            curent_read_order = 1
                        else:
                            curent_read_order = 0

                    # hold current kmer
                    if curent_kmer == line[2] and curent_position_id == int(line[1]):

                        curent_kmer_indice_len_list.append(abs(_position_end_idx - _position_start_idx))
                        curent_kmer_mean_list.append(float(line[6]))
                        if curent_position_start_idx > _position_start_idx:
                            curent_position_start_idx = _position_start_idx
                            assert curent_read_order == 1
                        else:
                            curent_position_end_idx = _position_end_idx
                            assert curent_read_order == 0

                    # a new kmer in this read
                    else:
                        # first add previous kmer to curent_read_data
                        curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                        curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                        assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)
                        # assert sum(curent_kmer_indice_len_list) == (curent_position_end_idx - curent_position_start_idx)

                        _sum_len = sum(curent_kmer_indice_len_list)

                        curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                        _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                        _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                        curent_read_data.append([curent_position_id + shift, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                                 curent_position_start_idx, curent_position_end_idx])

                        # then update paras for the new kmer
                        curent_kmer = line[2]
                        curent_position_id = int(line[1])
                        curent_position_start_idx = int(line[-2])
                        curent_position_end_idx = int(line[-1])
                        curent_kmer_indice_len_list = [abs(_position_end_idx - _position_start_idx)]
                        curent_kmer_mean_list = [float(line[6])]

                # a new read
                else:
                    # first add previous kmer data to the list
                    curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                    curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                    assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)
                    # assert sum(curent_kmer_indice_len_list) == (curent_position_end_idx - curent_position_start_idx)

                    _sum_len = sum(curent_kmer_indice_len_list)

                    curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                    _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                    _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                    curent_read_data.append([curent_position_id + shift, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                             curent_position_start_idx, curent_position_end_idx])

                    # then write previous kmer data to a file
                    _end_idx = [curent_read_data[0][-1], curent_read_data[0][-2],
                                curent_read_data[-1][-1], curent_read_data[-1][-2]]
                    if int(curent_read_index) in summary_dict:
                        # _ref = ref_dict[curent_contig][0][trans_position_start - shift: trans_position_end + shift+1]
                        _ref = ref_dict[curent_contig][0][curent_read_data[0][0] - shift: curent_read_data[-1][0] +
                                                                                          shift + 1]
                        combined_kmer = [curent_read_order, curent_contig, curent_read_index, curent_read_strand,
                                         curent_read_data[0][0], curent_read_data[-1][0],
                                         min(_end_idx), max(_end_idx),
                                         summary_dict[int(curent_read_index)]["read_name"],
                                         summary_dict[int(curent_read_index)]["fast5_path"].split('/')[-1], _ref]
                        write_to_csv(summary_file_name, combined_kmer)

                        write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                                 combine_file_name, trans2idx_dict)


                    # final update paras for the new kmer
                    curent_contig = line[0]
                    curent_read_index = line[3]
                    curent_read_strand = line[4]
                    curent_read_data = list()
                    curent_read_order = -1

                    curent_kmer = line[2]
                    curent_position_id = int(line[1])
                    curent_position_start_idx = int(line[-2])
                    curent_position_end_idx = int(line[-1])
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

                    curent_read_data.append([curent_position_id + shift, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                             curent_position_start_idx, curent_position_end_idx])

                    # then for read
                    _end_idx = [curent_read_data[0][-1], curent_read_data[0][-2],
                                curent_read_data[-1][-1], curent_read_data[-1][-2]]
                    if int(curent_read_index) in summary_dict:
                        # _ref = ref_dict[curent_contig][0][trans_position_start - shift: trans_position_end + shift+1]
                        _ref = ref_dict[curent_contig][0][curent_read_data[0][0] - shift: curent_read_data[-1][0] +
                                                                                          shift + 1]
                        combined_kmer = [curent_read_order, curent_contig, curent_read_index, curent_read_strand,
                                         curent_read_data[0][0], curent_read_data[-1][0],
                                         min(_end_idx), max(_end_idx),
                                         summary_dict[int(curent_read_index)]["read_name"],
                                         summary_dict[int(curent_read_index)]["fast5_path"].split('/')[-1], _ref]
                        write_to_csv(summary_file_name, combined_kmer)

                        write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                                 combine_file_name, trans2idx_dict)


def add_ref(eventalign_summary_file, ref_dict, shift):

    data = pd.read_csv(eventalign_summary_file)
    data['ref'] = None

    # summary_file_header = ['order', 'contig', 'read_index', 'strand',
    #                        'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
    #                        'read_name', 'fast5_path', 'ref']

    for i in tqdm(range(len(data))):
        contig = data.iat[i, 1]
        trans_position_start = data.iat[i, 4]
        trans_position_end = data.iat[i, 5]
        data.iat[i, 10] = ref_dict[contig][0][trans_position_start - shift: trans_position_end + shift + 1]

    data.to_csv(eventalign_summary_file, index=False, header=True)

    all_read_name_list = list(data["read_name"])
    all_read_name_set = set(all_read_name_list)
    if len(all_read_name_set) == len(all_read_name_list):
        print(" *  All reads in eventalign summary file are unique !! ")
    else:
        print(" *  Error !!! All reads in eventalign summary file are not unique.")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Combine the same kmers in eventalign file and write to a new file.')
    parser.add_argument('--root-dir', type=str, default='/scratch/cs/nanopore/chengg1/base_calling/dataset')
    parser.add_argument('--dataset', type=str, default='mes')
    parser.add_argument('--sample', type=str, default='mES_KO')
    parser.add_argument('--shift', type=int, default=2)
    parser.add_argument('--ref-fa', type=str, default='gencode.vM18.transcripts.fa')
    args = parser.parse_args()
    ref_fa = os.path.join(args.root_dir, args.dataset, "0_reference", args.ref_fa)

    # read reference data
    kmer2idx_dict, idx2kmer_dict = generate_kmer_dict('/scratch/cs/nanopore/chengg1/base_calling/dataset/xpore/'
                                                      '0_reference/model_kmer_motif/model_idx_kmer.csv')

    trans2idx_dict, idx2trans_dict = generate_trans_dict(os.path.join(args.root_dir, args.dataset, "0_reference",
                                                                      args.ref_fa + ".t2idx.csv"))

    # input file name
    eventalign_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                        args.sample + "-eventalign.txt")
    summary_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                     args.sample + "-summary.txt")

    # output file name
    eventalign_summary_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                                args.sample + "-eventalign-summary.csv")
    eventalign_combine_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                                args.sample + "-eventalign-combined.txt")

    print('*', args.dataset, args.sample)

    summary_dict, error_reads_list = read_summary(summary_file_name)

    # ref_dict = read_transcript_fa(ref_fa)
    # print(' *  Generate ref dict ! ')
    #
    # # --------------------------------------------------------------
    # # check output files
    # if os.path.exists(eventalign_combine_file_name):
    #     os.remove(eventalign_combine_file_name)
    # # eventalign_combine_file_name does not Ahave header
    #
    # if os.path.exists(eventalign_summary_file_name):
    #     os.remove(eventalign_summary_file_name)
    # summary_file_header = ['order', 'contig', 'read_index', 'strand',
    #                        'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
    #                        'read_name', 'fast5_path', 'ref']
    # write_to_csv(eventalign_summary_file_name, summary_file_header)
    #
    # # --------------------------------------------------------------
    # # combine eventalign file
    # summary_dict, error_reads_list = read_summary(summary_file_name)
    # print(" * Error reads list : ", error_reads_list)
    # combine_eventalign(eventalign_file_name, eventalign_combine_file_name, eventalign_summary_file_name,
    #                    args.shift, summary_dict, kmer2idx_dict, trans2idx_dict, error_reads_list, ref_dict)
