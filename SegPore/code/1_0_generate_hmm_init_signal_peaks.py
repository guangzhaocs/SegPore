# -*- coding: UTF-8 -*-
"""
@Project : base_calling_workflow 
@File    : 1_0_generate_hmm_input_signal.py
@Date    : 2022/12/13 20:33 
# @Github  : https://github.com/guangzhaocs/SegPore
"""
import os
import h5py
import csv
import argparse
import numpy as np
from csv import reader
from scipy.signal import find_peaks


def write_file(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def cal_slope0(y):
    # y is a numpy array
    x0 = np.full_like(y, 1)
    x1 = np.arange(len(x0))
    x = np.vstack((x1, x0))    # 2xn
    xx = np.linalg.inv(np.matmul(x, x.T))   # 2x2
    xy = np.matmul(x, y)
    beta = np.matmul(xx, xy) # 1x2
    return beta[0]  # k


def run_slope(y_arr, win_size_arr=[1, 2]):
    res = np.full_like(y_arr, 0)
    ny = len(y_arr)
    for w in win_size_arr:
        for i in range(w, ny-w):
            res[i] += cal_slope0(y_arr[i-w:i+w+1])
    return res


def generate_peaks(signal_file_name, save_peaks_file_name, peak_distance):

    with open(signal_file_name, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for i, row in enumerate(csv_reader):
            signal_arr = np.array([float(r) for r in row])
            slope_arr = np.abs(run_slope(signal_arr))
            slope_arr = smooth(slope_arr, 3)
            peaks_arr, _ = find_peaks(slope_arr, distance=peak_distance)
            write_file(save_peaks_file_name, peaks_arr)


def rm_exit_file(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def read_normalized_fast5(fast5_file, save_signal_file_name, save_readname_file_name):
    write_file(save_readname_file_name, ['read_name', 'read_index', 'contig', 'trans_position_start',
                                         'trans_position_end', 'start_idx', 'end_idx', 'seq'])

    with h5py.File(fast5_file, 'r+') as fast5_data:
        for read_key in fast5_data:
            read_name = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
            assert read_name == read_key.split("_")[1]
            if "NanopolishEvent" in fast5_data[read_key] and "Normalized" in fast5_data[read_key]:

                normalized_signal = fast5_data[read_key]["Normalized"]['Signal'][:]
                write_file(save_signal_file_name, normalized_signal)

                NanopolishEvent = fast5_data[read_key]['NanopolishEvent']
                read_index = NanopolishEvent.attrs['read_index']
                contig = NanopolishEvent.attrs['contig']
                # order = NanopolishEvent.attrs['order']
                trans_position_start = NanopolishEvent.attrs['trans_position_start']
                trans_position_end = NanopolishEvent.attrs['trans_position_end']
                start_idx = NanopolishEvent.attrs['start_idx']
                end_idx = NanopolishEvent.attrs['end_idx']
                seq = fast5_data[read_key]["NanopolishEvent"]['Reference'][()].decode("utf-8")
                write_file(save_readname_file_name, [read_name, int(read_index), contig,
                                                     int(trans_position_start), int(trans_position_end),
                                                     int(start_idx), int(end_idx), seq])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm signal.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--peak-distance', type=int, default=10)
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    args = parser.parse_args()

    fast5_file_name = os.path.join(args.root_dir, "1_fast5", args.sample + ".standardized.fast5")
    save_signal_dir = os.path.join(args.root_dir, "4_hhmm", "hhmm_init")
    save_signal_file_name = os.path.join(save_signal_dir, args.sample + "-signal.csv")
    save_readname_file_name = os.path.join(save_signal_dir, args.sample + "-readname.csv")
    save_peaks_file_name = os.path.join(save_signal_dir, args.sample + "-peaks.csv")

    rm_exit_file(save_signal_file_name)
    rm_exit_file(save_readname_file_name)
    rm_exit_file(save_peaks_file_name)

    print('* Starting to process for hhmm init ... ')
    read_normalized_fast5(fast5_file_name, save_signal_file_name, save_readname_file_name)
    print('* Starting to generate peaks ... ')
    generate_peaks(save_signal_file_name, save_peaks_file_name, args.peak_distance)
    print('* Generating signal, read_name, and peaks are over. ')