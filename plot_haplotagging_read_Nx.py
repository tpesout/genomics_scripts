#!/usr/bin/env python3
import argparse
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import pysam
import math
import pandas
import haplotagging_stats
import os
import collections


def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from haplotagging_stats tsv")
    parser.add_argument('--input_read_lengths', '-i', dest='input_read_lengths', default=None, required=True, type=str, action='append',
                        help='Read length file generated from "haplotagging_stats.py".  Can set multiple times, all values will be plotted')
    parser.add_argument('--figure_name', '-f', dest='figure_name', required=True, type=str,
                       help='Figure name')

    return parser.parse_args() if args is None else parser.parse_args(args)


def get_color(filename):
    if "margin" in filename.lower():
        return "red"
    if "whatshap" in filename.lower():
        return "blue"
    return "grey"

def get_linestyle(filename):
    if "untagged" in filename.lower():
        return "dotted"
    elif "tagged" in filename.lower():
        return "solid"
    return "dashdot"

def get_label(filename):
    if "margin" in filename.lower():
        if "untagged" in filename.lower():
            return "Margin Untagged"
        elif "tagged" in filename.lower():
            return "Margin Tagged"
    if "whatshap" in filename.lower():
        if "untagged" in filename.lower():
            return "Whatshap Untagged"
        elif "tagged" in filename.lower():
            return "Whatshap Tagged"
    return "Unknown"


def main():
    args = parse_args()
    read_lengths = dict()
    for len_file in args.input_read_lengths:
        lengths = np.genfromtxt(len_file, dtype=np.int32)
        lengths.sort(kind='stable')
        lengths = np.flip(lengths)
        assert(lengths[0] >= lengths[-1])
        read_lengths[os.path.basename(len_file)] = lengths

    length_bin_size = 1000
    cumulative_cov_bin_size = 1000000

    fig, ((ax1)) = plt.subplots(nrows=1,ncols=1)

    for i, filename in enumerate(read_lengths.keys()):
        lengths = read_lengths[filename]
        n50 = None
        half_total = sum(lengths) / 2
        total = 0
        curr_length = None
        curr_length_start_x = None
        first = True
        for length in lengths:
            if curr_length is None:
                curr_length = length // length_bin_size
                curr_length_start_x = 0

            # we need to plot something
            if length // length_bin_size != curr_length:
                # plot
                curr_length_end_x = total // cumulative_cov_bin_size
                if (first):
                    ax1.hlines(curr_length, curr_length_start_x, curr_length_end_x, alpha=.5, linestyles=get_linestyle(filename),
                               color=get_color(filename), label=get_label(filename))
                    ax1.vlines(curr_length_end_x, curr_length, length // length_bin_size,alpha=.5,
                               linestyles=get_linestyle(filename),  color=get_color(filename))
                    first = False
                else:
                    ax1.hlines(curr_length, curr_length_start_x, curr_length_end_x, linestyles=get_linestyle(filename),
                               color=get_color(filename), alpha=.5)
                    ax1.vlines(curr_length_end_x, curr_length, length // length_bin_size, alpha=.5,
                               linestyles=get_linestyle(filename),  color=get_color(filename))

                # iterate
                curr_length_start_x = curr_length_end_x
                curr_length = length // length_bin_size

            total += length
            if n50 is None and total >= half_total:
                n50 = length
        # ax1.vlines(half_total // length_bin_size, 0, n50 // length_bin_size, linestyles='dashed', color=get_color(filename))
        # ax1.annotate("{}".format(n50), (half_total // length_bin_size, n50 // length_bin_size), color=get_color(filename))

    ax1.ticklabel_format(axis='both', style='plain')
    ax1.set_ylim(0,100)
    ax1.set_ylabel("Read Length (kb)")
    ax1.set_xlabel("Cumulative Coverage (Mb)")

    # ax1.set_yscale('log')
    # ax2.set_yscale('log')
    plt.legend()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_stats_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show()


if __name__ == "__main__":
    main()
