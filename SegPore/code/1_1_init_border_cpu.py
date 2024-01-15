# -*- coding: UTF-8 -*-
"""
@Project : base_calling 
@File    : 01_read_orignal_and_init_border.py
@Date    : 2022/11/2 12:43 
@Github  : https://github.com/guangzhaocs/SegPore
"""
# -------------------------------------------------------------------------------------------------
#  original file : init_border_parallel.py
# -------------------------------------------------------------------------------------------------
import os, csv
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import norm
import itertools
import argparse
from csv import reader
from tqdm import tqdm


def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])


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


def ker_smooth(y_arr, b=1.0):
    res = np.full_like(y_arr, 0)
    ny = len(y_arr)
    w_arr = -np.arange(-3*b, 3*b+1)**2/(2*b*b)
    w_arr = np.exp(w_arr)  # Gaussian smoothing
    w_arr_sum = np.sum(w_arr)

    win_size = int(3*b)
    for i in range(win_size):
        st, en = i-win_size, i+win_size
        tmpinds = np.arange(st, en+1) >= 0
        res[i] = np.sum(w_arr[tmpinds] * y_arr[0:en+1]) / np.sum(w_arr[tmpinds])
    for i in range(win_size, ny-win_size):
        st, en = i - win_size, i + win_size
        res[i] = np.sum(w_arr * y_arr[st:en+1]) / w_arr_sum
    for i in range(ny-win_size, ny):
        st, en = i-win_size, i+win_size
        tmpinds = np.arange(st, en+1) < ny
        res[i] = np.sum(w_arr[tmpinds] * y_arr[st:ny]) / np.sum(w_arr[tmpinds])
    return res


# calculate the loglik for each point using a uniform prior for the slope
def cal_lr_loglik(y_arr, min_slope=5, n_sig=2.0,
                  slope_flag=False, sig_flag=False,
                  horizontal_line=False, mean_flag=True):
    def cal_loglik(slope, n_sig, mean_flag=mean_flag):
        b = np.mean(y_arr - slope*x_arr)
        res = y_arr - slope*x_arr - b
        if sig_flag:
            n_sig = np.std(res)
        if mean_flag:
            return np.mean(norm.logpdf(res, loc=0, scale=n_sig))
        else:
            return np.sum(norm.logpdf(res, loc=0, scale=n_sig))
    ny = len(y_arr)
    x_arr = np.arange(ny)

    if ny == 0:
        return 0

    if horizontal_line:
        return cal_loglik(0, n_sig)

    exp_slope = cal_slope0(y_arr)
    if not slope_flag:
        if 0.0 <= exp_slope < min_slope:
            exp_slope = min_slope
        elif -min_slope < exp_slope < 0.0:
            exp_slope = -min_slope

    return cal_loglik(exp_slope, n_sig)


def get_peak_win(y_arr, ny, p, max_w=8):
    p_arr = np.arange(p-max_w, p+max_w+1)
    p_arr = p_arr[p_arr >= 0]
    p_arr = p_arr[p_arr < ny]
    tmpdic = {}
    for st, en in itertools.combinations(p_arr, 2):
        if en-st <= 2 or p < st+1 or p > en-1:
            continue
        tmpdic[(st, en)] = cal_lr_loglik(y_arr[st:en + 1]) + \
                           cal_lr_loglik(y_arr[p_arr[0]:st], horizontal_line=True) + \
                           cal_lr_loglik(y_arr[en+1:p_arr[-1]], horizontal_line=True)
    max_win, max_val = None, -np.inf
    for k in tmpdic:
        if tmpdic[k] > max_val:
            max_win, max_val = k, tmpdic[k]
    return max_win


def get_peak_win_(para):
    return get_peak_win(para[0], para[1], para[2], max_w=8)


def get_init_t_wins(y_arr=None, debug=False, peak_distance=10):
    ny = len(y_arr)
    res = np.full_like(y_arr, 0)
    s_arr = np.abs(run_slope(y_arr))
    s_arr = smooth(s_arr, 3)
    # s_arr = ker_smooth(s_arr)
    peak_inds, _ = find_peaks(s_arr, distance=peak_distance)

    # border_res = list()
    # for p in peak_inds:
    #     if debug: plt.axvline(p)
    #     border_res.append(get_peak_win(y_arr, ny, p))

    peak_data = []
    for p in peak_inds[1:-1]:
        peak_data.append((y_arr, ny, p))

    from multiprocessing import Pool
    with Pool() as p:
        border_res = p.map(get_peak_win_, peak_data)

    for b in border_res:
        st = b[0]
        en = b[1]
        res[st:en+1] = 1

    if debug:
        plt.plot(y_arr, color="k")
        plt.plot(s_arr, color="r")
        plt.show()
    return res


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


def init_border_region(obs_arr: np.ndarray, peak_distance=10):
    """

    :param obs_arr: [110.20637099307125, 110.05821547788733, 109.02112687159988, 106.35432759828929, .... ]
    :param init_plot:  True or False
    :return: state_arr:  [0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0 ... ]
                         0 represents Base Region
                         1 represents Trans Region
    """

    state_arr = get_init_t_wins(y_arr=obs_arr, peak_distance=peak_distance)

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

    return state_arr, border_arr


def write_to_file(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--peak-distance', type=int, default=10)
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    args = parser.parse_args()

    signal_file = os.path.join(args.root_dir, "4_hhmm", "hhmm_init", args.sample + "-signal.csv")
    save_border_file = os.path.join(args.root_dir, "4_hhmm", "hhmm_init", args.sample + "-border.csv")

    if os.path.exists(save_border_file):
        os.remove(save_border_file)

    print('* Start to processing ... ')
    num_lines = wc_count(signal_file)

    with open(signal_file, 'r') as read_obj:
        csv_reader = reader(read_obj)
        for i, row in enumerate(tqdm(csv_reader, total=num_lines)):
            row = np.array([float(r) for r in row])
            _, border_arr = init_border_region(obs_arr=row)
            write_to_file(save_border_file, border_arr)
    print('* Generating init border is over. ')