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


def read_polya_results_tsv(polya_results_file_name):
    """
    Read PolyA file
    :param polya_results_file_name:
    :return:
    """
    polya_results_pd = pd.read_csv(polya_results_file_name, sep='\t').to_dict('records')
    polya_results_dict = dict()
    for item in polya_results_pd:
        polya_results_dict[item['readname']] = item
    return polya_results_dict


def read_event_summary_file(eventalign_summary_file, trans2idx_dict):
    """
   ['order', 'contig', 'read_index', 'strand',
    'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
    'read_name', 'fast5_path', 'ref']
    """

    data = pd.read_csv(eventalign_summary_file)
    data['read_index'] = data['read_index'].apply(lambda x: int(x))
    data['contig'] = data['contig'].apply(lambda x: trans2idx_dict[x])
    idx2fast5_dict = data.set_index('read_index')['fast5_path'].to_dict()
    idx2readName_dict = data.set_index('read_index')['read_name'].to_dict()
    idx2fullInfo_arr = data.to_dict("records")
    idx2fullInfo_dict = dict()
    for key in idx2fullInfo_arr:
        idx2fullInfo_dict[key['read_index']] = key
    return idx2fast5_dict, idx2readName_dict, idx2fullInfo_dict


def read_filter_file(filter_file_name):
    """
    ["read_name", "ref_idx", "flag", "reference_start", "reference_end"]
    """
    data = pd.read_csv(filter_file_name)
    data['ref_idx'] = data['ref_idx'].apply(lambda x: str(x))
    data['key'] = data['read_name'] + "#" + data['ref_idx']
    return list(data['key'])


def generate_trans_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    trans2idx_dict = data.set_index('transcript')['index'].to_dict()
    idx2trans_dict = data.set_index('index')['transcript'].to_dict()
    return trans2idx_dict, idx2trans_dict


def write_to_file(file_name, data):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)


def write_event_combine_file(summary_dict, read_index, contig, read_data, save_dir):

    curent_read_name = summary_dict[int(read_index)]["read_name"]
    _read_file_name = os.path.join(save_dir, read_index + '#' + curent_read_name + "#" + contig + '.txt')
    if os.path.exists(_read_file_name):
        return
        # os.remove(_read_file_name)
    with open(_read_file_name, 'w', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        for _data in read_data:
            writer.writerow(_data)


def split_eventalign(eventalign_file_name, save_dir, filter_arr, idx2fast5_dict, idx2readname_dict):
    """
    eventalign_file:
               0         1        2      3       4          5           6
         ['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']
    """
    print('* Start to processing ... ')
    num_lines = sum(1 for _ in open(eventalign_file_name, 'r'))

    with open(eventalign_file_name, 'r') as f:
        for i, line in enumerate(tqdm(f, total=num_lines)):
            line = line.strip().replace('\n', '').replace('\r', '').split('\t')
            if i == 0:
                last_read_idx = int(line[0])
                last_contig = int(line[1])
                last_pos = int(line[2])
                last_kmer = int(line[3])
                last_mean = float(line[4])
                last_start_idx = int(line[5])
                last_end_idx = int(line[6])
                curr_read_data = [[last_read_idx, last_contig, last_pos, last_kmer, last_mean, last_start_idx,
                                   last_end_idx]]
            else:
                if last_read_idx != int(line[0]):
                    fast5_name = idx2fast5_dict[last_read_idx]
                    read_name = idx2readname_dict[last_read_idx]
                    save_file_name = os.path.join(save_dir, fast5_name + '.csv')
                    if read_name + "#" + str(last_contig) not in filter_arr:
                        write_to_file(save_file_name, curr_read_data)
                    curr_read_data = list()

                last_read_idx = int(line[0])
                last_contig = int(line[1])
                last_pos = int(line[2])
                last_kmer = int(line[3])
                last_mean = float(line[4])
                last_start_idx = int(line[5])
                last_end_idx = int(line[6])
                curr_read_data.append([last_read_idx, last_contig, last_pos, last_kmer, last_mean, last_start_idx,
                                       last_end_idx])

            # the last line
            if i == num_lines - 1:
                fast5_name = idx2fast5_dict[last_read_idx]
                read_name = idx2readname_dict[last_read_idx]
                save_file_name = os.path.join(save_dir, fast5_name + '.csv')
                if read_name + "#" + str(last_contig) not in filter_arr:
                    write_to_file(save_file_name, curr_read_data)



def generate_event_data_dict(tmp_event_file_name, idx2fullInfo_dict):
    """
    :param idx2fullInfo_dict:  ['order', 'contig', 'read_index', 'strand',
                                'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
                                'read_name', 'fast5_path', 'ref']
    :return:
    """
    tmp_event_dict = dict()
    tmp_event_data = pd.read_csv(tmp_event_file_name, header=None)
    tmp_event_data.columns = ['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']
    tmp_event_data = tmp_event_data.groupby(['read_idx', 'contig'])
    for site, group in tmp_event_data:
        _read_index = site[0]
        _read_info_dict = idx2fullInfo_dict[_read_index]
        _read_name = _read_info_dict['read_name']
        group = group[['pos', 'kmer', 'mean', 'start_idx', 'end_idx']]
        _read_info_dict['event'] = group.values.tolist()
        tmp_event_dict[_read_name] = _read_info_dict
    return tmp_event_dict


def check_fast5(fast5_dir):
    all_files = os.listdir(fast5_dir)
    all_files = ['FAK01563_771a38df9d9b678a485923dc79a6a5d5d59c115c_17.fast5']
    for fast5_file_name in all_files:
        print(10 * "-", fast5_file_name, 10 * "-")
        fast5_file_full_name = os.path.join(fast5_dir, fast5_file_name)
        with h5py.File(fast5_file_full_name, 'r+') as fast5_data:
            for read_key in fast5_data:
                if "NanopolishEvent" in fast5_data[read_key]:
                    print("event : ", read_key)
                    del fast5_data[read_key]["NanopolishEvent"]

                if "Normalized" in fast5_data[read_key]:
                    print("normal : ", read_key)
                    del fast5_data[read_key]["Normalized"]


# def read_and_write_multi_fast5(fast5_dir, tmp_dir, idx2fullInfo_dict, idx2trans_dict):
#
#     all_files = os.listdir(fast5_dir)
#     for fast5_file_name in all_files:
#         print(" * ", fast5_file_name)
#         tmp_event_file_name = os.path.join(tmp_dir, fast5_file_name + ".csv")
#         if not os.path.exists(tmp_event_file_name):
#             continue
#         tmp_event_dict = generate_event_data_dict(tmp_event_file_name, idx2fullInfo_dict)
#         # key: read_name
#         # {order, contig, read_index, strand, trans_position_start, trans_position_end,
#         #  start_idx, end_idx, read_name, fast5_path, ref, event}
#
#         fast5_file_full_name = os.path.join(fast5_dir, fast5_file_name)
#         with h5py.File(fast5_file_full_name, 'r+') as fast5_data:
#             for read_key in fast5_data:
#                 read_name = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
#                 if read_name in tmp_event_dict:
#
#                     if "NanopolishEvent" in fast5_data[read_key]:
#                         # del fast5_data[read_key]["NanopolishEvent"]
#                         continue
#
#                     event_data = fast5_data[read_key].create_group("NanopolishEvent")
#                     event_data['Events'] = tmp_event_dict[read_name]['event']
#                     event_data['Reference'] = tmp_event_dict[read_name]['ref']
#                     event_data.attrs.create("read_index", tmp_event_dict[read_name]['read_index'], shape=None, dtype=None)
#                     event_data.attrs.create("contig_index", tmp_event_dict[read_name]['contig'], shape=None, dtype=None)
#                     # event_data.attrs.create("contig_name", idx2trans_dict[int(tmp_event_dict[read_name]['contig'])],
#                     #                         shape=None, dtype="S")
#                     event_data.attrs.create("order", tmp_event_dict[read_name]['order'], shape=None, dtype=None)
#                     # event_data.attrs.create("strand", tmp_event_dict[read_name]['strand'], shape=None, dtype="S")
#                     event_data.attrs.create("trans_position_start", tmp_event_dict[read_name]['trans_position_start'],
#                                             shape=None, dtype=None)
#                     event_data.attrs.create("trans_position_end", tmp_event_dict[read_name]['trans_position_end'],
#                                             shape=None, dtype=None)
#                     event_data.attrs.create("start_idx", tmp_event_dict[read_name]['start_idx'], shape=None, dtype=None)
#                     event_data.attrs.create("end_idx", tmp_event_dict[read_name]['end_idx'], shape=None, dtype=None)
#
#         del tmp_event_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Combine the same kmers in eventalign file and write to a new file.')
    parser.add_argument('--root-dir', type=str, default='/scratch/cs/nanopore/chengg1/base_calling/dataset')
    parser.add_argument('--dataset', type=str, default='mes')
    parser.add_argument('--sample', type=str, default='mES_KO')
    parser.add_argument('--shift', type=int, default=2)
    parser.add_argument('--ref-fa', type=str, default='gencode.vM18.transcripts.fa')
    args = parser.parse_args()
    ref_fa = os.path.join(args.root_dir, args.dataset, "0_reference", args.ref_fa)

    # read file name
    multi_fast5_dir = os.path.join(args.root_dir, args.dataset, "2_fast5", args.sample,  "fast5")
    filter_file_name = os.path.join(args.root_dir, args.dataset, "3_bamtx", args.sample, "filter.csv")
    eventalign_summary_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                                args.sample + "-eventalign-summary.csv")
    eventalign_combine_file_name = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample,
                                                args.sample + "-eventalign-combined.txt")

    # save file name
    save_dir = os.path.join(args.root_dir, args.dataset, "3_nanopolish_res", args.sample, "tmp_dir")

    print(' *', args.dataset, args.sample)

    # ----------------------------------------------------------------------------------------
    # step 1: read some files
    # ----------------------------------------------------------------------------------------
    trans2idx_dict, idx2trans_dict = generate_trans_dict(os.path.join(args.root_dir, args.dataset, "0_reference",
                                                                      args.ref_fa + ".t2idx.csv"))
    idx2fast5_dict, idx2readname_dict, idx2fullInfo_dict = read_event_summary_file(eventalign_summary_file_name,
                                                                                   trans2idx_dict)
    filter_arr = read_filter_file(filter_file_name)

    print(" * read all files over !")

    # ----------------------------------------------------------------------------------------
    # step 2: split eventalign-combined.txt into fast5 dir
    # ----------------------------------------------------------------------------------------
    if os.path.exists(save_dir):
        import shutil
        shutil.rmtree(save_dir)
    os.mkdir(save_dir)
    split_eventalign(eventalign_combine_file_name, save_dir, filter_arr, idx2fast5_dict, idx2readname_dict)
    print(" * split over !")