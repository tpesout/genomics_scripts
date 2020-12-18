#!/usr/bin/env python3
import argparse
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import pysam
import os
import math

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


def plot_two_stats(stats, key_x_fn, key_y_fn, label_x=None, label_y=None, figure_name=None):

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

        if "whatshap_phased" in filename:
            color="blue"
        elif "margin_phased" in filename:
            color="red"
        else:
            color="grey"

        label = filename.replace("HG001_", "")\
            .replace("_2_GRCh37_pepper_margin_deepvariant.", " ")\
            .replace("_guppy422_GRCh37_pepper_margin_deepvariant.", " ")\
            .replace("_phased.vcf", "")
        # print("{}: {}".format(label_x, list(key_x_fn(file_stats))))
        # print("{}: {}".format(label_y, list(key_y_fn(file_stats))))
        plt.scatter(key_x_fn(file_stats), key_y_fn(file_stats), marker=marker, color=color, label=label, alpha=.3)

    if label_x is not None: plt.xlabel(label_x)
    if label_y is not None: plt.ylabel(label_y)
    plt.ticklabel_format(style='plain')
    plt.legend()
    plt.tight_layout()
    if figure_name is not None:
        if not figure_name.endswith(".png"):
            figure_name += ".png"
        plt.savefig(figure_name)

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

    plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50]), lambda x: float(x[ALL_SWITCH_RATE]), "N50", "Switch Rate",
                   figure_name=None if fig_name is None else fig_name+".switch_rate")
    plot_two_stats(stats_dict, lambda x: int(x[BLOCK_N50]), lambda x: float(x[ALL_HAMMING_RATE]), "N50", "Hamming Rate",
                   figure_name=None if fig_name is None else fig_name+".hamming_rate")


if __name__ == "__main__":
    main()
