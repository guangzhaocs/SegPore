# -*- coding: utf-8 -*-
# @Time    : 2021/12/28 15:34
# @FileName: combine_eventalign.py
# @Github  : https://github.com/guangzhaocs/SegPore
import os
import sys
from scipy.stats import norm
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


def read_model_kmer(file_name):
    """
    The file is without header.
    Read model kmer and return a dict.
    :param file_name:
    :return:
    """
    model_kmer_pd = pd.read_csv(file_name, header=None)
    model_kmer_pd.columns = ["model_kmer", "model_mean_un", "model_stdv_un", "model_mean_mod", "model_stdv_mod"]
    model_kmer_pd['model_mean_un'] = model_kmer_pd['model_mean_un'].apply(float)
    model_kmer_pd['model_stdv_un'] = model_kmer_pd['model_stdv_un'].apply(float)
    model_kmer_pd['model_mean_mod'] = model_kmer_pd['model_mean_mod'].apply(float)
    model_kmer_pd['model_stdv_mod'] = model_kmer_pd['model_stdv_mod'].apply(float)
    model_kmer_pd = model_kmer_pd.to_dict('records')
    model_kmer_dict = dict()
    for item in model_kmer_pd:
        model_kmer_dict[item['model_kmer']] = item
    return model_kmer_dict


def write_to_csv(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(row)


def get_mod_sta(_kmer, _mu):

    pdf_un = np.round(norm.pdf(_mu, loc=model_kmer_dict[_kmer]["model_mean_un"],
                               scale=model_kmer_dict[_kmer]["model_stdv_un"]), 5)
    pdf_mod = np.round(norm.pdf(_mu, loc=model_kmer_dict[_kmer]["model_mean_mod"],
                                scale=model_kmer_dict[_kmer]["model_stdv_mod"]), 5)

    if pdf_un >= pdf_mod:
        return 0
    else:
        return 1


def write_event_combine_file(read_index, contig, read_data, combine_file_name):

    if len(read_data) == 0:
        return

    read_data = np.array(read_data)
    df = pd.DataFrame(read_data, columns=['pos', 'kmer', 'kmer_idx', 'mean', 'start_idx', 'end_idx', 'event_len'])
    df['pos'] = df['pos'].apply(int)
    df['mean'] = df['mean'].apply(float)
    df['kmer_idx'] = df['kmer_idx'].apply(int)
    df['start_idx'] = df['start_idx'].apply(int)
    df['end_idx'] = df['end_idx'].apply(int)
    df['event_len'] = df['event_len'].apply(int)
    df['contig'] = contig
    df['read_idx'] = read_index
    df['mod'] = df.apply(lambda x: get_mod_sta(x['kmer'], x['mean']), axis=1)
    df = df[['read_idx', 'contig', 'pos', 'kmer', 'kmer_idx', 'mean', 'start_idx', 'end_idx', 'event_len', 'mod']]
    df.to_csv(combine_file_name, index=False, header=False, sep='\t', mode='a+')


def combine_eventalign(input_file_name, combine_file_name, kmer2idx_dict):
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

    with open(input_file_name, 'r') as f:
        for idx, line in enumerate(f):
            first_line = line.strip().replace('\n', '').replace('\r', '').split('\t')
            if first_line[3] != "-" and first_line[8] != "-1":
                # parameters for current read
                curent_read_index = first_line[0]
                curent_contig = first_line[2]
                curent_read_data = list()
                curent_read_order = -1

                # parameters for current kmer
                curent_kmer = first_line[3]
                curent_position_id = int(first_line[4])
                curent_position_start_idx = int(first_line[7])
                curent_position_end_idx = int(first_line[8])
                curent_kmer_indice_len_list = [int(first_line[9])]
                curent_kmer_mean_list = [float(first_line[5])]
                break

    print("* start idx :", idx)

    print('* Start to processing ... ')
    num_lines = wc_count(input_file_name)

    with open(input_file_name, 'r') as f:
        for i, line in enumerate(tqdm(f, total=num_lines)):
            if i > idx:
                line = line.strip().replace('\n', '').replace('\r', '').split('\t')
                if line[3] == "-" or line[8] == "-1":
                    # the last line
                    if i == num_lines - 1:
                        # first for kmer
                        curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                        curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                        assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)
                        _sum_len = sum(curent_kmer_indice_len_list)
                        if _sum_len != 0:
                            curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                            _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                            _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                            curent_read_data.append(
                                [curent_position_id, curent_kmer, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                 curent_position_start_idx, curent_position_end_idx, _sum_len])

                        # then for read
                        write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                                 combine_file_name)

                    continue

                _position_start_idx = int(line[7])
                _position_end_idx = int(line[8])

                # hold current read
                if curent_contig == line[2] and curent_read_index == line[0]:

                    # check order [order = 1: reverse]
                    if curent_read_order == -1:
                        if curent_position_start_idx > _position_start_idx:
                            curent_read_order = 1
                        else:
                            curent_read_order = 0
                        # print(" * curent_position_start_idx ,", curent_position_start_idx)
                        # print(" * _position_start_idx ,", _position_start_idx)
                        # print(" * update order : ", curent_read_order)


                    # hold current kmer
                    if curent_kmer == line[3] and curent_position_id == int(line[4]):

                        curent_kmer_indice_len_list.append(int(line[9]))
                        curent_kmer_mean_list.append(float(line[5]))
                        if curent_position_start_idx > _position_start_idx:
                            curent_position_start_idx = _position_start_idx
                            assert curent_read_order == 1
                        else:
                            # print(line)
                            # print(curent_position_end_idx)
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
                        if _sum_len != 0:
                            curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                            _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                            _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                            curent_read_data.append([curent_position_id, curent_kmer, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                                     curent_position_start_idx, curent_position_end_idx, _sum_len])

                        # then update paras for the new kmer
                        curent_kmer = line[3]
                        curent_position_id = int(line[4])
                        curent_position_start_idx = int(line[7])
                        curent_position_end_idx = int(line[8])
                        curent_kmer_indice_len_list = [int(line[9])]
                        curent_kmer_mean_list = [float(line[5])]

                # a new read
                else:
                    # first add previous kmer data to the list
                    curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                    curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                    assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)
                    # assert sum(curent_kmer_indice_len_list) == (curent_position_end_idx - curent_position_start_idx)

                    _sum_len = sum(curent_kmer_indice_len_list)
                    if _sum_len != 0:
                        curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                        _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                        _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                        curent_read_data.append([curent_position_id, curent_kmer, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                                 curent_position_start_idx, curent_position_end_idx, _sum_len])

                    # then write previous read data to a file

                    write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                             combine_file_name)

                    # final update paras for the new kmer
                    curent_contig = line[2]
                    curent_read_index = line[0]
                    curent_read_data = list()
                    curent_read_order = -1

                    curent_kmer = line[3]
                    curent_position_id = int(line[4])
                    curent_position_start_idx = int(line[7])
                    curent_position_end_idx = int(line[8])
                    curent_kmer_indice_len_list = [int(line[9])]
                    curent_kmer_mean_list = [float(line[5])]
                    # print(30 * "-")
                    # print(" * start a new read ---- ", line)

                # the last line
                if i == num_lines - 1:

                    # first for kmer
                    curent_kmer_indice_len_list = np.array(curent_kmer_indice_len_list)
                    curent_kmer_mean_list = np.array(curent_kmer_mean_list)

                    assert len(curent_kmer_indice_len_list) == len(curent_kmer_mean_list)
                    # assert sum(curent_kmer_indice_len_list) == (curent_position_end_idx - curent_position_start_idx)
                    _sum_len = sum(curent_kmer_indice_len_list)
                    if _sum_len != 0:
                        curent_kmer_indice_len_list = curent_kmer_indice_len_list / _sum_len
                        _kmer_evnet_mean = np.sum(curent_kmer_indice_len_list * curent_kmer_mean_list)
                        _kmer_evnet_mean = np.around(_kmer_evnet_mean, 3)

                        curent_read_data.append([curent_position_id, curent_kmer, kmer2idx_dict[curent_kmer], _kmer_evnet_mean,
                                                 curent_position_start_idx, curent_position_end_idx, _sum_len])

                    # then for read
                    write_event_combine_file(curent_read_index, curent_contig, curent_read_data,
                                             combine_file_name)


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

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    parser.add_argument('--eventalign', type=str, default="segpore_eventalign_2D")
    args = parser.parse_args()

    # read reference data
    kmer2idx_dict, idx2kmer_dict = generate_kmer_dict(os.path.join(args.root_dir, "0_reference", "model_idx_kmer.csv"))

    # generate model kmer dict
    model_kmer_dict = read_model_kmer(os.path.join(args.root_dir, "0_reference", "model_kmer_m6A_without_header.csv"))

    # input file name
    eventalign_file_name = os.path.join(args.root_dir, "5_align", args.sample,  args.eventalign + ".txt")

    # output file name
    eventalign_combine_file_name = os.path.join(args.root_dir, "5_align", args.sample,  args.eventalign + "_combined.txt")
    if os.path.exists(eventalign_combine_file_name):
        os.remove(eventalign_combine_file_name)
    write_to_csv(eventalign_combine_file_name,
                 ['read_idx', 'contig', 'pos', 'kmer', 'kmer_idx', 'mean', 'start_idx', 'end_idx', 'event_len', 'mod'])

    # eventalign_combine_file_name does not have header

    combine_eventalign(eventalign_file_name, eventalign_combine_file_name, kmer2idx_dict)
    print('* Done !!!')
