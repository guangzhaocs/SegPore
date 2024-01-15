# -*- coding: utf-8 -*-
# @Time    : 2021/12/28 15:34
# @FileName: plot_density.py
# @Github  : https://github.com/guangzhaocs/SegPore
import os
import csv
import numpy as np
import pandas as pd
import argparse
import math
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from GMM_fix_KNN import GaussianMixtureModel
from scipy.stats import norm


def write(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def plot_all_site_density(nano_df, segpore_df, binwidth=0.5):
    def subplot():
        plt.ylim((0, 0.12))
        plt.xlim((103, 137))

    nano_df['label'] = 1
    segpore_df['label'] = 2

    plt.figure(figsize=(12, 4))
    pal = sns.color_palette("pastel")

    ax1 = plt.subplot(1, 2, 1)
    ax1_1 = sns.histplot(
        data=nano_df, x="mean", hue="label", kde=True, stat="density",
        palette=[pal[0]], legend=False,
        alpha=.7, linewidth=0,  binwidth=binwidth
    )
    ax1_1.set_xlabel("")
    ax1_1.set_ylabel("")
    plt.title("Nanopolish")
    subplot()
    # ax1.axis('off')

    ax2 = plt.subplot(1, 2, 2)
    ax2_1 = sns.histplot(
        data=segpore_df, x="mean", hue="label", kde=True, stat="density",
        palette=[pal[0]], legend=False,
        alpha=.7, linewidth=0,  binwidth=binwidth
    )
    ax2_1.set_xlabel("")
    ax2_1.set_ylabel("")
    subplot()
    plt.title("SegPore eventalign (GGACT)")

    # ax2.axis('off')
    # plt.subplots_adjust(left=0.03, right=0.99, top=0.95, bottom=0.01)
    # plt.show()
    plt.savefig("SegPore_GGACT.jpg")


def plot_density_site_contrast(title, my_res, nano_res, save_folder):
    binwidth = 0.5
    pal = sns.color_palette("pastel")

    plt.figure(figsize=(17, 7))

    plt.subplot(1, 2, 1)
    sns.histplot(nano_res, kde=True, stat="density", linewidth=0, binwidth=binwidth)
    subplot("Nanopolish")

    plt.subplot(1, 2, 2)
    sns.histplot(my_res, kde=True, stat="density", linewidth=0, binwidth=binwidth)
    subplot("Ours")

    plt.suptitle(title)
    # plt.show()
    plt.savefig(os.path.join(save_folder, title+".jpg"))
    plt.clf()


def subplot_for_two_com(mean_1, sigma_1, w_1, mean_2, sigma_2, w_2):

    plt.axvline(x=mean_1, ls='--', color="black")
    plt.axvline(x=mean_2, ls='--', color="red")

    plot_x_arr = np.linspace(mean_1 - 3 * sigma_1, mean_1 + 3 * sigma_1)
    plt.plot(plot_x_arr, w_1 * norm.pdf(plot_x_arr, mean_1, sigma_1), color="red")

    plot_x_arr = np.linspace(mean_2 - 3 * sigma_2, mean_2 + 3 * sigma_2)
    plt.plot(plot_x_arr, w_2 * norm.pdf(plot_x_arr, mean_2, sigma_2), color="red")
    

def plot_density_site_only_one_sample(data, mean_1, sigma_1, w_1, mean_2, sigma_2, w_2):
    binwidth = 0.5
    pal = sns.color_palette("pastel")

    def get_label(x):
        if norm.logpdf(x, loc=mean_1, scale=sigma_1) > norm.logpdf(x, loc=mean_2, scale=sigma_2):
            return 1
        else:
            return 2

    # data['label'] = data['mean'].apply(lambda x: get_label(x))
    data['label'] = 0

    plt.figure()

    sns.histplot(
        data=data, x="mean", hue="label", kde=False, stat="density",
        palette=pal, legend=False,
        alpha=.7, linewidth=0, binwidth=binwidth
    )
    subplot_for_two_com(mean_1, sigma_1, w_1, mean_2, sigma_2, w_2)

    plt.suptitle("GMM with one component mean fixed")
    # plt.savefig("2.jpg")
    plt.show()


def subplot(title):
    plt.title(title)
    plt.ylim((0, 0.20))
    plt.xlim((105, 135))


def fit_GMM(data, fix_mean):
    data = np.expand_dims(np.array(data), 1)

    gm2 = GaussianMixtureModel(n_components=2, random_state=1, fix_mu_flag=True, fix_mu=fix_mean, max_sigma=4.0)
    gm2.fit(data)

    pred_mean = gm2.means_
    mean_1 = round(pred_mean[0][0], 2)
    mean_2 = round(pred_mean[1][0], 2)

    pred_sigma = gm2.covariances_
    sigma_1 = round(np.sqrt(pred_sigma[0][0][0]), 2)
    sigma_2 = round(np.sqrt(pred_sigma[1][0][0]), 2)

    w_1 = round(gm2.weights_[0], 2)
    w_2 = round(gm2.weights_[1], 2)

    bic = gm2.bic(data)

    return mean_1, sigma_1, w_1, mean_2, sigma_2, w_2, bic


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate hhmm init border.')
    parser.add_argument('--root-dir', type=str, default='demo')
    parser.add_argument('--sample', type=str, default='HEK293T_WT_rep1_FAK27249_demo_0')
    parser.add_argument('--eventalign', type=str, default="segpore_eventalign_2D_combined_GGACT.txt")
    args = parser.parse_args()

    segpore_df = pd.read_csv(args.eventalign, header=None, sep='\t')
    segpore_df.columns = ['read_index', 'contig', 'pos', 'kmer', 'kmer_idx', 'mean', 'start_idx', 'end_idx', 'len']

    mean_1, sigma_1, w_1, mean_2, sigma_2, w_2, bic = fit_GMM(list(segpore_df['mean']), 123.83)
    print(mean_1, sigma_1, w_1, mean_2, sigma_2, w_2)

    # plot_density_site_only_one_sample(segpore_df, mean_1, sigma_1, w_1, mean_2, sigma_2, w_2)