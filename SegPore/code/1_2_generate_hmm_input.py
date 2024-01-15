# -*- coding: UTF-8 -*-
"""
@Project : base_calling 
@File    : generate_hmm_input.py
@Date    : 2022/11/2 12:43 
@Github  : https://github.com/guangzhaocs/SegPore
"""
import os
import argparse


def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])


def check_exist_border_file(hmm_init_dir):

    all_files_list = os.listdir(hmm_init_dir)
    files_dict = dict()

    for file_name in all_files_list:
        file_name = file_name.split('-')[0]
        if file_name not in files_dict:
            files_dict[file_name] = 1
        else:
            files_dict[file_name] += 1

    files_with_border_list = list()
    for file_name in files_dict:
        if files_dict[file_name] == 3:
            files_with_border_list.append(file_name)

    print(hmm_init_dir)
    print("files_dict len : ", len(files_dict))
    print("files_with_border_list len : ", len(files_with_border_list))

    return files_with_border_list


def check_lines(hmm_init_dir, files_with_border_list):

    files_with_same_lines_list = list()

    for i, file_name in enumerate(files_with_border_list):
        signal_file_name = os.path.join(hmm_init_dir, file_name + "-signal.csv")
        border_file_name = os.path.join(hmm_init_dir, file_name + "-border.csv")
        readname_file_name = os.path.join(hmm_init_dir, file_name + "-readname.csv")

        signal_lines = wc_count(signal_file_name)
        border_lines = wc_count(border_file_name)
        reads_lines = wc_count(readname_file_name)

        if signal_lines == border_lines and signal_lines == reads_lines:
            print(i, " ", file_name)
            files_with_same_lines_list.append(file_name)

    print("files_with_same_lines_list len : ", len(files_with_same_lines_list))
    return files_with_same_lines_list


def save_files_to_hmm_input_dir(hmm_init_dir, hmm_input_dir, hmm_output_dir):

    idx = 0
    for file_name in os.listdir(hmm_init_dir):
        if "readname" in file_name:
            file_name = file_name.split("-")[0]
            print(idx, file_name)
            idx += 1
            signal_file_name = os.path.join(hmm_init_dir, file_name + "-signal.csv")
            border_file_name = os.path.join(hmm_init_dir, file_name + "-border.csv")

            save_dir = os.path.join(hmm_input_dir, file_name)
            if not os.path.exists(save_dir):
                os.mkdir(save_dir)

            out_dir = os.path.join(hmm_output_dir, file_name)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            save_signal_file_name = os.path.join(save_dir, "signal.csv")
            save_border_file_name = os.path.join(save_dir, "init_border.csv")

            if os.path.exists(save_signal_file_name):
                os.remove(save_signal_file_name)

            if os.path.exists(save_border_file_name):
                os.remove(save_border_file_name)

            os.system('cp ' + signal_file_name + " " + save_signal_file_name)
            os.system('cp ' + border_file_name + " " + save_border_file_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    args = parser.parse_args()

    root_dir = os.path.join(args.root_dir, "4_hhmm")

    hmm_init_dir = os.path.join(root_dir, "hhmm_init")
    hmm_input_dir = os.path.join(root_dir, "hhmm_input")
    hmm_output_dir = os.path.join(root_dir, "hhmm_output")
    hmm_final_dir = os.path.join(root_dir, "hhmm_final")

    save_files_to_hmm_input_dir(hmm_init_dir, hmm_input_dir, hmm_output_dir)