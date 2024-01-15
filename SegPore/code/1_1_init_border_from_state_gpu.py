# -*- coding: UTF-8 -*-
"""
@Project : base_calling_workflow 
@File    : 1_1_generate_border_from_state_gpu.py
@Date    : 2022/12/13 20:33 
@Github  : https://github.com/guangzhaocs/SegPore
"""
import os
import csv
import argparse
from csv import reader
import numpy as np


def state2border(state_arr):
    """
    :param state_arr:
    :param idx:
    :return:
    """
    t_final_st_sig_ind_arr = np.where(np.diff(state_arr) == 1)[0] + 1
    t_final_en_sig_ind_arr = np.where(np.diff(state_arr) == -1)[0] + 1

    border_arr = np.array([0, len(state_arr)])
    border_arr = np.sort(np.append(border_arr, np.append(t_final_st_sig_ind_arr, t_final_en_sig_ind_arr)))
    return border_arr


def write_file(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def generate_borders_from_state(_state_file_name, _border_file_name):

    with open(_state_file_name, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for i, row in enumerate(csv_reader):
            state_arr = np.array([int(r) for r in row])
            assert state_arr[0] == i
            state_arr = state_arr[1:]

            if state_arr[0] == 1:
                for en in range(1, len(state_arr)):
                    if state_arr[en - 1] == 1 and state_arr[en] == 0:
                        break
                state_arr[0: en] = 0

            if state_arr[-1] == 1:
                for st in np.arange(len(state_arr) - 1, 1, -1):
                    if state_arr[st - 1] == 0 and state_arr[st] == 1:
                        break
                state_arr[st: len(state_arr)] = 0

            assert state_arr[0] == state_arr[-1] == 0
            border_arr = state2border(state_arr)
            assert len(border_arr) % 2 == 0
            write_file(_border_file_name, border_arr)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--peak-distance', type=int, default=10)
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    args = parser.parse_args()

    state_file_name = os.path.join(args.root_dir, "4_hhmm", "hhmm_init", args.sample + "-state.csv")
    border_file_name = os.path.join(args.root_dir, "4_hhmm", "hhmm_init", args.sample + "-border.csv")

    if os.path.exists(border_file_name):
        os.remove(border_file_name)
    generate_borders_from_state(state_file_name, border_file_name)
    print('* Generating init border from state (GPU) is over. ')