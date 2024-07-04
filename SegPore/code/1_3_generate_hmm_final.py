# -*- coding: UTF-8 -*-
"""
@Project : base_calling_workflow 
@File    : generate_hmm_final.py
@Date    : 2022/11/24 10:07 
@Github  : https://github.com/guangzhaocs/SegPore
"""
import os
import csv
import sys
import numpy as np
from csv import reader
import argparse
import pandas as pd


def wc_count(filename):
    with open(filename) as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


def write_to_file(file_name, row):

    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def read_res_csv_file(file_name):
    """
    In cuda hmm res, the first element is idx.
    :param file_name:
    :return:
    """
    arr_dict = dict()
    with open(file_name, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            arr = [int(r) for r in row]
            idx = arr[0]
            arr = arr[1:]
            arr_dict[idx] = arr
    return arr_dict


def check_state_arr(state_arr):
    """
    state in state_arr must be in [0,1,2,3,4]
    :param state_arr:
    :return:  True or False
    """
    if not (state_arr[0] == 1 or state_arr[0] == 2 or state_arr[0] == 3 or state_arr[0] == 4):
        return False
    for i in state_arr[1:]:
        if not (i == 0 or i == 1 or i == 2 or i == 3 or i == 4):
            return False
    return True


def cal_mu_and_sigma_arr(signal_arr, state_arr, border_arr):

    b_st_sig_ind_arr = border_arr[::2]
    b_en_sig_ind_arr = border_arr[1::2]
    t_st_sig_ind_arr = border_arr[1:-1:2]
    t_en_sig_ind_arr = border_arr[2:-1:2]

    mu_arr = list()
    sigma_arr = list()
    len_arr = list()

    prev_mu_arr = list()
    prev_sigma_arr = list()
    prev_len_arr = list()

    next_mu_arr = list()
    next_sigma_arr = list()
    next_len_arr = list()

    for b_st_idx, b_ed_idx in zip(b_st_sig_ind_arr, b_en_sig_ind_arr):
        tmp_sig_arr = np.array(signal_arr[b_st_idx: b_ed_idx])
        tmp_state_arr = np.array(state_arr[b_st_idx: b_ed_idx])

        tmp_curr_idx = np.where(tmp_state_arr == 1)[0]
        tmp_prev_idx = np.where(tmp_state_arr == 2)[0]
        tmp_next_idx = np.where(tmp_state_arr == 3)[0]

        curr_len = len(tmp_curr_idx)
        prev_len = len(tmp_prev_idx)
        next_len = len(tmp_next_idx)

        # for current component
        if curr_len == len(tmp_sig_arr):
            mu_arr.append(np.mean(tmp_sig_arr))
            sigma_arr.append(np.std(tmp_sig_arr))
        elif curr_len == 0:
            mu_arr.append(0.0)
            sigma_arr.append(0.0)
        else:
            tmp_curr_sig_arr = tmp_sig_arr[tmp_curr_idx]
            mu_arr.append(np.mean(tmp_curr_sig_arr))
            sigma_arr.append(np.std(tmp_curr_sig_arr))

        # for prev component
        if prev_len == len(tmp_sig_arr):
            prev_mu_arr.append(np.mean(tmp_sig_arr))
            prev_sigma_arr.append(np.std(tmp_sig_arr))
        elif prev_len == 0:
            prev_mu_arr.append(0.0)
            prev_sigma_arr.append(0.0)
        else:
            tmp_prev_sig_arr = tmp_sig_arr[tmp_prev_idx]
            prev_mu_arr.append(np.mean(tmp_prev_sig_arr))
            prev_sigma_arr.append(np.std(tmp_prev_sig_arr))

        # for next component
        if next_len == len(tmp_sig_arr):
            next_mu_arr.append(np.mean(tmp_sig_arr))
            next_sigma_arr.append(np.std(tmp_sig_arr))
        elif next_len == 0:
            next_mu_arr.append(0.0)
            next_sigma_arr.append(0.0)
        else:
            tmp_next_sig_arr = tmp_sig_arr[tmp_next_idx]
            next_mu_arr.append(np.mean(tmp_next_sig_arr))
            next_sigma_arr.append(np.std(tmp_next_sig_arr))

        len_arr.append(curr_len)
        prev_len_arr.append(prev_len)
        next_len_arr.append(next_len)

    for t_st_idx, t_ed_idx in zip(t_st_sig_ind_arr, t_en_sig_ind_arr):
        assert (state_arr[t_st_idx: t_ed_idx]).all() == 0

    for b_st_idx, b_ed_idx in zip(b_st_sig_ind_arr, b_en_sig_ind_arr):
        assert (state_arr[b_st_idx: b_ed_idx]).all() != 0

    return mu_arr, sigma_arr, len_arr, prev_mu_arr, prev_sigma_arr, prev_len_arr, next_mu_arr, \
                next_sigma_arr, next_len_arr


def generate_all_mu_sigma(signal_file, res_border_file, res_state_file, read_info_file, save_dir):

    read_pd = pd.read_csv(read_info_file)
#     read_info_dict = dict()
#     for index, row in read_pd.iterrows():
#         read_info_dict[index] = list(row)

    border_arr_dict = read_res_csv_file(res_border_file)
    state_arr_dict = read_res_csv_file(res_state_file)
    assert len(border_arr_dict) == len(state_arr_dict) == len(read_pd)

    with open(signal_file, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for idx, row in enumerate(csv_reader):
            signal_arr = np.array([float(r) for r in row])

            if idx in border_arr_dict:
                assert idx in state_arr_dict

                border_arr = np.array(border_arr_dict[idx])
                state_arr = np.array(state_arr_dict[idx])
                assert len(state_arr) == len(signal_arr)

                if check_state_arr(state_arr):
                    mu_arr, sigma_arr, len_arr, prev_mu_arr, prev_sigma_arr, prev_len_arr, next_mu_arr, \
                    next_sigma_arr, next_len_arr = cal_mu_and_sigma_arr(signal_arr, state_arr, border_arr)
                else:
                    mu_arr, sigma_arr, len_arr, prev_mu_arr, prev_sigma_arr, prev_len_arr, next_mu_arr, \
                    next_sigma_arr, next_len_arr = [-1], [-1], [-1], [-1], [-1], [-1], [-1], [-1], [-1]
            else:
                mu_arr, sigma_arr, len_arr, prev_mu_arr, prev_sigma_arr, prev_len_arr, next_mu_arr, \
                next_sigma_arr, next_len_arr = [-1], [-1], [-1], [-1], [-1], [-1], [-1], [-1], [-1]

#             assert read_info_dict[idx][3] == 1

            write_to_file(os.path.join(save_dir, "curr_mu.csv"), mu_arr)
            write_to_file(os.path.join(save_dir, "curr_sigma.csv"), sigma_arr)
            write_to_file(os.path.join(save_dir, "curr_len.csv"), len_arr)
            write_to_file(os.path.join(save_dir, "prev_mu.csv"), prev_mu_arr)
            write_to_file(os.path.join(save_dir, "prev_sigma.csv"), prev_sigma_arr)
            write_to_file(os.path.join(save_dir, "prev_len.csv"), prev_len_arr)
            write_to_file(os.path.join(save_dir, "next_mu.csv"), next_mu_arr)
            write_to_file(os.path.join(save_dir, "next_sigma.csv"), next_sigma_arr)
            write_to_file(os.path.join(save_dir, "next_len.csv"), next_len_arr)


def check_final(hmm_init_dir, hmm_input_dir, hmm_final_dir, sample):
    read_info_file = os.path.join(hmm_init_dir, sample + "-readname.csv")
    signal_file = os.path.join(hmm_input_dir, sample, "signal.csv")
    curr_len_file = os.path.join(hmm_final_dir, sample, "curr_len.csv")
    curr_mu_file = os.path.join(hmm_final_dir, sample, "curr_mu.csv")
    curr_sigma_file = os.path.join(hmm_final_dir, sample, "curr_sigma.csv")

    try:
        signal_file_lines = wc_count(signal_file)
        curr_len_file_lines = wc_count(curr_len_file)
        read_info_file_lines = wc_count(read_info_file) - 1
        curr_mu_file_lines = wc_count(curr_mu_file)
        curr_sigma_file_lines = wc_count(curr_sigma_file)

        if signal_file_lines == curr_len_file_lines == read_info_file_lines == \
                curr_mu_file_lines == curr_sigma_file_lines:
            print(" * ", sample, " success !")

        else:
            print(" * ", sample, " error !")
    except:
        print(" * ", sample, " error !")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    args = parser.parse_args()

    root_dir = os.path.join(args.root_dir, "4_hhmm")
    hmm_init_dir = os.path.join(root_dir, "hhmm_init")
    hmm_input_dir = os.path.join(root_dir, "hhmm_input")
    hmm_output_dir = os.path.join(root_dir, "hhmm_output")
    hmm_final_dir = os.path.join(root_dir, "hhmm_final")

    signal_file = os.path.join(hmm_input_dir, args.sample, "signal.csv")
    init_border_file = os.path.join(hmm_input_dir, args.sample, "init_border.csv")
    read_info_file = os.path.join(hmm_init_dir, args.sample + "-readname.csv")
    res_border_file = os.path.join(hmm_output_dir, args.sample, "res_border.csv")
    res_state_file = os.path.join(hmm_output_dir, args.sample, "res_state.csv")

    signal_file_lines = wc_count(signal_file)
    init_border_file_lines = wc_count(init_border_file)
    read_info_file_lines = wc_count(read_info_file) - 1
    res_border_file_lines = wc_count(res_border_file)
    res_state_file_lines = wc_count(res_state_file)

    if signal_file_lines != init_border_file_lines and signal_file_lines != read_info_file_lines and \
            signal_file_lines != res_border_file_lines and signal_file_lines != res_state_file_lines:
        print(args.folder, " signal line and init border line are not the same ! ")
        print(signal_file_lines, init_border_file_lines, read_info_file_lines,
              res_border_file_lines, res_state_file_lines)
        sys.exit()

    save_dir = os.path.join(hmm_final_dir, args.sample)
    if os.path.exists(save_dir):
        import shutil
        shutil.rmtree(save_dir)
    os.makedirs(save_dir)

    print(" * starting to generate the final res ...")
    generate_all_mu_sigma(signal_file, res_border_file, res_state_file, read_info_file, save_dir)

    check_final(hmm_init_dir, hmm_input_dir, hmm_final_dir, args.sample)
