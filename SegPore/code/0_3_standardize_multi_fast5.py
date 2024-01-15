# -*- coding: utf-8 -*-
# @Time    : 2021/12/28 15:34
# @Author  : https://github.com/guangzhaocs/SegPore
# @FileName: standardize_multi_fast5.py
import os
import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import h5py
import sys
POLYA_STANDARD_MU = 108.9
POLYA_STANDARD_SIGMA = 1.67


def generate_kmer_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    kmer2idx_dict = data.set_index('kmer')['index'].to_dict()
    idx2kmer_dict = data.set_index('index')['kmer'].to_dict()
    return kmer2idx_dict, idx2kmer_dict


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def percentile(arr, q=5):
    """
    Remove < 5% and > 95%.
    :param arr:
    :param q:
    :return:
    """
    low_percen = np.percentile(arr, q)
    hign_percen = np.percentile(arr, 100 - q)
    filter_arr = []
    for x in arr:
        if low_percen <= x <= hign_percen:
            filter_arr.append(x)

    return np.array(filter_arr)


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


def read_event_summary_file(eventalign_summary_file):
    """
   ['order', 'contig', 'read_index', 'strand',
    'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
    'read_name', 'fast5_path', 'ref']
    """

    data = pd.read_csv(eventalign_summary_file)
    data['read_index'] = data['read_index'].apply(lambda x: int(x))
    idx2fullInfo_arr = data.to_dict("records")
    idx2fullInfo_dict = dict()
    for key in idx2fullInfo_arr:
        idx2fullInfo_dict[key['read_index']] = key
    return idx2fullInfo_dict


def write_to_file(file_name, data):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)


def generate_event_data_dict(tmp_event_file_name, idx2read_dict, kmer2idx_dict):
    tmp_event_dict = dict()
    tmp_event_data = pd.read_csv(tmp_event_file_name, sep='\t')
    tmp_event_data['kmer'] = tmp_event_data['kmer'].apply(lambda x: kmer2idx_dict[x])
    tmp_event_data = tmp_event_data.groupby(['read_idx', 'contig'])
    for site, group in tmp_event_data:
        _read_index = site[0]
        _read_info_dict = idx2read_dict[_read_index]
        _read_name = _read_info_dict['read_name']
        group = group[['pos', 'kmer', 'mean', 'start_idx', 'end_idx']]
        _read_info_dict['event'] = group.values.tolist()
        tmp_event_dict[_read_name] = _read_info_dict
    return tmp_event_dict


def convert_raw_signal_to_pA_value(raw_signal, offset, range, digitisation):
    """
    Transform the raw signal to pico-ampere current values.
      𝑠𝑖𝑔𝑛𝑎𝑙_𝑖𝑛_𝑝𝑖𝑐𝑜_𝑎𝑚𝑝𝑒𝑟𝑒 = (𝑟𝑎𝑤_𝑠𝑖𝑔𝑛𝑎𝑙_𝑣𝑎𝑙𝑢𝑒 + 𝑜𝑓𝑓𝑠𝑒𝑡) ∗ 𝑟𝑎𝑛𝑔𝑒 / 𝑑𝑖𝑔𝑖𝑡𝑖𝑠𝑎𝑡𝑖𝑜𝑛
    :param raw_signal:
    :param offset:
    :param range:
    :param digitisation:
    :return: signal_in_pico_ampere
    """
    return (raw_signal + offset) * range / digitisation


def normalize_fast5_with_polyA(signal_in_pico_ampere, polyA_start, polyA_end):

    signal_pts = 50

    polyA_len = polyA_end - polyA_start
    polyA_sig = signal_in_pico_ampere[polyA_start: polyA_end]

    polyA_extract_start = int(polyA_len * 0.25)
    polyA_extract_end = int(polyA_len * 0.75)
    polyA_extract_sig = polyA_sig[polyA_extract_start: polyA_extract_end]

    polya_extract_smooth_sig = smooth(polyA_sig, signal_pts)[polyA_extract_start: polyA_extract_end]

    polyA_mu = np.mean(polya_extract_smooth_sig)
    polya_res_arr = polyA_extract_sig - polya_extract_smooth_sig
    polyA_sigma = np.std(polya_res_arr)

    if polyA_sigma < 3.0:
        ployA_sigma_percentile = np.std(percentile(polya_res_arr))
        normalized_signal = ((signal_in_pico_ampere - polyA_mu) / ployA_sigma_percentile) * POLYA_STANDARD_SIGMA + \
                           POLYA_STANDARD_MU
        normalized_signal = np.around(np.array(normalized_signal), decimals=3)
        return normalized_signal, np.around(polyA_mu, decimals=5), np.around(ployA_sigma_percentile, decimals=5)
    else:
        return None, None, None


def write_and_normalize_fast5(fast5_file, split_event_file, idx2read_dict, polya_results_dict, kmer2idx_dict):

    tmp_event_dict = generate_event_data_dict(split_event_file, idx2read_dict, kmer2idx_dict)
    print(f" * generate event dict over ! dict length is {len(tmp_event_dict)}.")

    # da848fa9-1322-4fea-b550-7efb32b014b6
    # {'contig': 'ENST00000273480.3', 'read_index': 0, 'strand': '+', 'trans_position_start': 16,
    # 'trans_position_end': 902, 'start_idx': 9867, 'end_idx': 53405,
    # 'read_name': 'da848fa9-1322-4fea-b550-7efb32b014b6',
    # 'fast5_path': 'HEK293T_WT_rep1_FAK27249_demo_0.fast5', 'ref': 'TAGGC...',
    # 'event': [...] }

    with h5py.File(fast5_file, 'r+') as fast5_data:
        for read_key in fast5_data:
            read_name = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
            if read_name in tmp_event_dict:
                event_data = fast5_data[read_key].create_group("NanopolishEvent")
                event_data['Events'] = tmp_event_dict[read_name]['event']
                event_data['Reference'] = tmp_event_dict[read_name]['ref']
                event_data.attrs.create("read_index", tmp_event_dict[read_name]['read_index'], shape=None, dtype=None)
                event_data.attrs.create("contig", tmp_event_dict[read_name]['contig'], shape=None, dtype=None)
                event_data.attrs.create("trans_position_start", tmp_event_dict[read_name]['trans_position_start'],
                                        shape=None, dtype=None)
                event_data.attrs.create("trans_position_end", tmp_event_dict[read_name]['trans_position_end'],
                                        shape=None, dtype=None)
                event_data.attrs.create("start_idx", tmp_event_dict[read_name]['start_idx'], shape=None, dtype=None)
                event_data.attrs.create("end_idx", tmp_event_dict[read_name]['end_idx'], shape=None, dtype=None)

            if read_name in polya_results_dict:
                polyA_start = int(polya_results_dict[read_name]['polya_start'])
                polyA_end = int(polya_results_dict[read_name]['transcript_start'])
                raw_signal = fast5_data[read_key]["Raw"]['Signal'][:]

                # read channel attributes parameters
                channel_id = fast5_data[read_key]['channel_id']
                digitisation = channel_id.attrs['digitisation']
                offset = channel_id.attrs['offset']
                range = channel_id.attrs['range']

                # transform the raw signal to pico-ampere current values
                signal_in_pico_ampere = convert_raw_signal_to_pA_value(raw_signal, offset, range, digitisation)

                # normalized pico-ampere current values
                normalized_signal, polyA_mu, ployA_sigma_percentile = \
                    normalize_fast5_with_polyA(signal_in_pico_ampere, polyA_start, polyA_end)

                if normalized_signal is not None:
                    norm_data = fast5_data[read_key].create_group("Normalized")
                    norm_data['Signal'] = normalized_signal
                    norm_data.attrs.create("polyA_start", polyA_start, shape=None, dtype=None)
                    norm_data.attrs.create("polyA_end", polyA_end, shape=None, dtype=None)
                    norm_data.attrs.create("polyA_mu", polyA_mu, shape=None, dtype=None)
                    norm_data.attrs.create("polyA_sigma", ployA_sigma_percentile, shape=None, dtype=None)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Standardize the fast5 file.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    args = parser.parse_args()

    # read reference data
    kmer2idx_dict, idx2kmer_dict = generate_kmer_dict(os.path.join(args.root_dir, "0_reference", "model_idx_kmer.csv"))

    # file names
    fast5_file_name = os.path.join(args.root_dir, "1_fast5", args.sample + ".standardized.fast5")
    eventalign_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_eventalign_combined.txt")
    eventalign_summary_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_eventalign_summary.csv")
    polyA_file_name = os.path.join(args.root_dir, "3_nanopolish", args.sample + "_polya_pass.tsv")

    if (not os.path.exists(fast5_file_name)) or (not os.path.exists(eventalign_file_name)) \
            or (not os.path.exists(polyA_file_name)):
        print(" *  fast5 file or eventalign file or ployA file do not exist!")
        sys.exit()

    # ----------------------------------------------------------------------------------------
    # step 1: read summary dict and polyA dict
    # ----------------------------------------------------------------------------------------
    idx2fullInfo_dict = read_event_summary_file(eventalign_summary_file_name)
    # key: 0
    # {'contig': 'ENST00000273480.3', 'read_index': 0, 'strand': '+',
    # 'trans_position_start': 16, 'trans_position_end': 902,
    # 'start_idx': 9867, 'end_idx': 53405, 'read_name': 'da848fa9-1322-4fea-b550-7efb32b014b6',
    # 'fast5_path': 'HEK293T_WT_rep1_FAK27249_demo_0.fast5', 'ref': 'TAGGCAC...' }

    polya_results_dict = read_polya_results_tsv(polyA_file_name)
    # da848fa9-1322-4fea-b550-7efb32b014b6
    # {'readname': 'da848fa9-1322-4fea-b550-7efb32b014b6', 'contig': 'ENST00000273480.3', 'position': 14,
    # 'leader_start': 44.0, 'adapter_start': 1875.0, 'polya_start': 7772.0, 'transcript_start': 9939.0,
    # 'read_rate': 111.56, 'polya_length': 75.22, 'qc_tag': 'PASS'}

    # ----------------------------------------------------------------------------------------
    # step 2: split eventalign-combined.txt into fast5 dir
    # ----------------------------------------------------------------------------------------
    write_and_normalize_fast5(fast5_file_name, eventalign_file_name,
                              idx2fullInfo_dict, polya_results_dict, kmer2idx_dict)
    print(' * All over !!!')
