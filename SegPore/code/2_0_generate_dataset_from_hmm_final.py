# -*- coding: utf-8 -*-
# @Time    : 2022/7/7 17:20
# @FileName: generate_dataset_for_align.py
# @Software: PyCharm
# @Github  : https://github.com/guangzhaocs/SegPore

import os, csv
import numpy as np
import h5py
import sys
from csv import reader
import pandas as pd
import argparse
from bisect import bisect_left, bisect_right


Base2Id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}


def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])


def write(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def read_read_name_file(file_name):
    arr_list = list()
    with open(file_name, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            arr_list.append(row[0])
    return arr_list


def read_int_csv_file(file_name):
    arr_list = list()
    with open(file_name, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            arr = [int(r) for r in row]
            arr_list.append(arr)
    return arr_list


def read_float_csv_file(file_name):
    arr_list = list()
    with open(file_name, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            arr = [float(r) for r in row]
            arr_list.append(arr)
    return arr_list


def check_file(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def generate_sub_mu_sigma_arr(signal_start_idx, signal_end_idx, border_arr, mu_arr, sigma_arr, len_arr):

    b_st_sig_ind_arr = border_arr[::2]
    b_en_sig_ind_arr = border_arr[1::2]

    assert len(b_st_sig_ind_arr) == len(b_en_sig_ind_arr) == len(mu_arr) == len(sigma_arr) == len(len_arr), \
        f"b_st_sig_ind_arr_len = {len(b_st_sig_ind_arr)}, " \
        f"b_en_sig_ind_arr_len = {len(b_en_sig_ind_arr)}, " \
        f"mu_arr_len = {len(mu_arr)}, " \
        f"sigma_arr_len = {len(sigma_arr)}, len_arr = {len(len_arr)}"


    assert 0 <= signal_start_idx <= border_arr[-1]
    assert 0 <= signal_end_idx <= border_arr[-1]

    _start_idx = bisect_left(b_st_sig_ind_arr, signal_start_idx)
    _end_idx = bisect_right(b_st_sig_ind_arr, signal_end_idx)

    return mu_arr[_start_idx: _end_idx], sigma_arr[_start_idx: _end_idx], len_arr[_start_idx: _end_idx], \
           b_st_sig_ind_arr[_start_idx: _end_idx], b_en_sig_ind_arr[_start_idx: _end_idx]


def check_res(sub_len_arr, sub_start_arr, sub_end_arr):
    for _len, _start, _end in zip(sub_len_arr, sub_start_arr, sub_end_arr):
        assert _start + _len <= _end, f"start={_start}, end={_end}, len={_len}"


def porcess_one_dir(root_dir, sample):

    hmm_init_dir = os.path.join(args.root_dir, "4_hhmm", "hhmm_init")
    hmm_output_dir = os.path.join(args.root_dir, "4_hhmm", "hhmm_output")
    hmm_final_dir = os.path.join(args.root_dir, "4_hhmm", "hhmm_final")
    eventalign_dir = os.path.join(args.root_dir, "5_align")

    # read folder
    readname_file_name = os.path.join(hmm_init_dir, sample + "-readname.csv")
    border_file_name = os.path.join(hmm_output_dir, sample, "res_border.csv")
    mu_file_name = os.path.join(hmm_final_dir, sample, "curr_mu.csv")
    sigma_file_name = os.path.join(hmm_final_dir, sample, "curr_sigma.csv")
    len_file_name = os.path.join(hmm_final_dir, sample, "curr_len.csv")

    # save folder
    save_folder = os.path.join(eventalign_dir, sample)
    if os.path.exists(save_folder):
        import shutil
        shutil.rmtree(save_folder)
    os.mkdir(save_folder)

    save_mu_file = os.path.join(save_folder, "mu.csv")
    save_sigma_file = os.path.join(save_folder, "sigma.csv")
    save_label_file = os.path.join(save_folder, "label.csv")
    save_readname_file = os.path.join(save_folder, "readname.csv")
    save_start_file = os.path.join(save_folder, "start_pos.csv")
    save_end_file = os.path.join(save_folder, "end_pos.csv")
    save_len_file = os.path.join(save_folder, "len.csv")

    check_file(save_mu_file)
    check_file(save_sigma_file)
    check_file(save_label_file)
    check_file(save_readname_file)
    check_file(save_start_file)
    check_file(save_end_file)
    check_file(save_len_file)

    # step 1: generate label list
    # ['read_name', 'read_index', 'contig', 'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx', 'seq']
    read_pd = pd.read_csv(readname_file_name)
    label_dict = dict()
    for index, row in read_pd.iterrows():
        label_dict[index] = list(row)

    # step 2: read mu, sigma, border
    border_list = read_int_csv_file(border_file_name)
    mu_list = read_float_csv_file(mu_file_name)
    sigma_list = read_float_csv_file(sigma_file_name)
    len_list = read_int_csv_file(len_file_name)

    assert len(label_dict) == len(border_list) == len(mu_list) == len(sigma_list) == len(len_list), \
        f"error, lengths are not the same"

    # step 3: generate_dataset
    for i in label_dict:

        label_item = label_dict[i]
        # [0'read_name', 1'read_index', 2'contig', 3'trans_position_start', 4'trans_position_end',
        #  5'start_idx', 6'end_idx', 7'seq']

        signal_start_idx = label_item[5]
        signal_end_idx = label_item[6]
        seq = label_item[7]
        label_seq = [Base2Id[base] for base in seq]

        border_arr = border_list[i]
        border_arr = border_arr[1:]  # cuda output: the first element is idx
        mu_arr = mu_list[i]
        sigma_arr = sigma_list[i]
        len_arr = len_list[i]

        try:
            if len(mu_arr) == 1:
                assert mu_arr[0] == -1
                print("--- mu is -1 !!!!")
                continue

            sub_mu_arr, sub_sigma_arr, sub_len_arr, sub_start_arr, sub_end_arr = generate_sub_mu_sigma_arr(
                signal_start_idx,
                signal_end_idx,
                border_arr,
                mu_arr, sigma_arr, len_arr)

            check_res(sub_len_arr, sub_start_arr, sub_end_arr)
            assert len(sub_sigma_arr) == len(sub_mu_arr) == len(sub_start_arr) == len(sub_end_arr) == len(
                sub_len_arr)

            sub_mu_arr = np.around(sub_mu_arr, decimals=4)
            sub_sigma_arr = np.around(sub_sigma_arr, decimals=4)

            write(save_mu_file, sub_mu_arr)
            write(save_sigma_file, sub_sigma_arr)
            write(save_len_file, sub_len_arr)
            write(save_label_file, label_seq)
            write(save_start_file, sub_start_arr)
            write(save_end_file, sub_end_arr)
            write(save_readname_file, label_item[0:-1])

        except:
            print("--- fail !!!!", folder_name)

    # check
    start_lines = wc_count(save_start_file)
    end_lines = wc_count(save_end_file)
    label_lines = wc_count(save_label_file)
    mu_lines = wc_count(save_mu_file)
    sigma_lines = wc_count(save_sigma_file)
    len_lines = wc_count(save_len_file)
    readname_lines = wc_count(save_readname_file)

    if start_lines == end_lines == label_lines == mu_lines == sigma_lines == readname_lines == len_lines:
        print(f" * {sample} : success ！")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    args = parser.parse_args()
    porcess_one_dir(args.root_dir, args.sample)