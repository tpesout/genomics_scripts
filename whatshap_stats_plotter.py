#!/usr/bin/env python3
import argparse
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import pysam
import os
import math

plt.style.use('ggplot')
text_fontsize = 8
# plt.rcParams['ytick.labelsize']=text_fontsize+4
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

FILE_NAME = "file_name"
BLOCK_N50 = "block_n50"
ALL_SWITCH_RATE = "all_switch_rate"
ALL_HAMMING_RATE = "blockwise_hamming_rate"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from tpesout's whatshap wdl workflow")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                        help='TSV file from whatshap wdl workflow with stats and pairwise info')

    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')

    return parser.parse_args() if args is None else parser.parse_args(args)



def get_label(filename):
    if ".whatshap_phased" in filename.lower():
        tool = "Whatshap"
    elif ".phased" in filename.lower() or "margin_phased":
        tool = "Margin"
    else:
        tool = "UNKNOWN"

    if "ccs" in filename.lower():
        data = "HiFi"
    elif "25x" in filename.lower():
        data = "ONT 25x"
    elif "50x" in filename.lower():
        data = "ONT 50x"
    elif "75x" in filename.lower():
        data = "ONT 75x"
    else:
        data = "UNKNOWN"

    return data + " " + tool


def get_color(filename):
    Margin_COLOR = [(238 / 256.0, 39 / 256.0, 8 / 256.0, 0.7), (238 / 256.0, 39 / 256.0, 8 / 256.0, 0.8),
                    (238 / 256.0, 39 / 256.0, 8 / 256.0, 1.0)]
    WhatsHap_COLOR = [(10/256.0, 118/256.0, 148/256.0, 0.7), (10/256.0, 118/256.0, 148/256.0, 0.8), (10/256.0, 118/256.0, 148/256.0, 1.0)]
    if "whatshap" in filename.lower():
        return WhatsHap_COLOR[0]
    elif "margin" in filename.lower():
        return Margin_COLOR[0]
    else:
        return "grey"


def plot_two_stats(stats, key_x_fn, key_y_fn, figsize=(3.5, 3.5), label_x=None, label_y=None, figure_name=None):
    fig, ((ax1)) = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    for filename in stats.keys():
        file_stats = stats[filename]

        # get marker
        if "UL" in filename:
            marker = '*'
        elif "75x" in filename:
            marker = 'p'
        elif "50x" in filename:
            marker = 's'
        elif "25x" in filename:
            marker = '^'
        else:
            marker = "o"

        color=get_color(filename)
        label=get_label(filename)

        plt.scatter(key_x_fn(file_stats), key_y_fn(file_stats), marker=marker, color=color, label=label, alpha=.3)

    plt.xlim(0, max(map(key_x_fn, stats.values())) * 1.1)
    plt.ylim(0, max(map(key_y_fn, stats.values())) * 1.1)
    if label_x is not None: plt.xlabel(label_x)
    if label_y is not None: plt.ylabel(label_y)
    plt.ticklabel_format(style='plain')
    plt.legend()
    plt.tight_layout()
    if figure_name is not None:
        if not figure_name.endswith(".png"):
            figure_name += ".png"
        figure_name += ".pdf"
        plt.savefig(figure_name, format='pdf', dpi=300)

    plt.show()
    plt.close()


def main():
    args = parse_args()

    raw_list_info = list()
    stats_dict = dict()
    header = None
    with open(args.input) as f:
        for line in f:
            if header is None:
                header = line.strip().lstrip("#").split("\t")
                continue
            if len(line.strip()) == 0: continue
            parts = line.strip().split("\t")
            assert(len(parts) == len(header))
            sample_stats = dict()
            for i, key in enumerate(header):
                sample_stats[key] = parts[i]
            stats_dict[sample_stats[FILE_NAME]] = sample_stats

    print("Got {} records from {}".format(len(stats_dict), args.input))

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(os.path.basename(args.input).split(".")[:-1])


    plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50])/1000000.0, lambda x: float(x[ALL_SWITCH_RATE]), label_x="N50 (Mb)", label_y="Switch Rate",
                   figure_name=None if fig_name is None else fig_name+".switch_rate.1")
    plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50])/1000000.0, lambda x: float(x[ALL_SWITCH_RATE]), label_x="N50 (Mb)", label_y="Switch Rate",
                   figure_name=None if fig_name is None else fig_name+".switch_rate.2", figsize=(2.5, 3.5))
    plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50])/1000000.0, lambda x: float(x[ALL_HAMMING_RATE]), label_x="N50 (Mb)", label_y="Hamming Rate",
                   figure_name=None if fig_name is None else fig_name+".hamming_rate.1")

    # plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50]), lambda x: float(x[ALL_SWITCH_RATE]), "N50", "Switch Rate",
    #                figure_name=None if fig_name is None else fig_name+".switch_rate")
    # plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50]), lambda x: float(x[ALL_HAMMING_RATE]), "N50", "Hamming Rate",
    #                figure_name=None if fig_name is None else fig_name+".hamming_rate")


if __name__ == "__main__":
    main()
