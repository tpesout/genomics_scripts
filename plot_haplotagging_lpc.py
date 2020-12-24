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

DECAY = "decay"
LENGTH_SCALE = "length_scale_num_vars"
LENGTH_SCALE_BP = "length_scale_bps"
WEIGHTED_MEAN = "weighted_mean"
APPROX = "approx_"
TOTAL_EFFECTIVE_SIZE = "total_eff_size"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from margin's calcLocalPhasingCorrectness ")
    parser.add_argument('--input', '-i', dest='input_tsvs', default=None, required=True, type=str, action='append',
                        help='TSV file generated from calcLocalPhasingCorrectness.  Can set multiple times, all values will be plotted')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')
    parser.add_argument('--length_scale_bps', '-b', dest='length_scale_bps', default=False, required=False, action='store_true',
                       help='Lenght scale should use {} instead of {}'.format(LENGTH_SCALE_BP, LENGTH_SCALE))
    parser.add_argument('--second_plot_max_length', '-l', dest='second_plot_max_length', default=None, required=False, action='store', type=int,
                       help='Produce second plot with different max value')
    parser.add_argument('--third_plot_max_length', '-L', dest='third_plot_max_length', default=None, required=False, action='store', type=int,
                       help='Produce third plot with different max value')


    return parser.parse_args() if args is None else parser.parse_args(args)

def get_label(filename):
    if ".whatshap_phased" in filename.lower():
        tool = "Whatshap"
    elif ".phased" in filename.lower():
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
    if "whatshap" in filename.lower():
        return "blue"
    elif "margin" in filename.lower():
        return "red"
    else:
        return "grey"


def get_linestyle(filename):
    if "ccs" in filename.lower():
        return "solid"
    elif "25x" in filename.lower():
        return "dotted"
    elif "50x" in filename.lower():
        return "dashdot"
    elif "75x" in filename.lower():
        return "dashed"
    else:
        return "solid"


def main():
    args = parse_args()
    dataframes = dict()
    for tsv in args.input_tsvs:
        dataframes[os.path.basename(tsv)] = pandas.read_csv(tsv, sep='\t')


    if args.second_plot_max_length is None:
        fig, ((ax1)) = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    elif args.third_plot_max_length is None:
        fig, ((ax1, ax2)) = plt.subplots(nrows=1,ncols=2, figsize=(12, 6))
    else:
        # fig = plt.figure(figsize=(12, 6))
        fig = plt.figure(constrained_layout=True, figsize=(12, 6))
        gs = fig.add_gridspec(2, 2)
        ax1 = fig.add_subplot(gs[0:2, 0:1])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 1])

    y_lim_min_1 = 1.0
    y_lim_min_2 = 1.0
    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        length_scale = list()
        decay = list()
        weighted_mean = list()
        total_effective_size = list()
        tes_by_length = list()
        for x in dataframe.iterrows():
            length_scale.append(x[1][LENGTH_SCALE_BP if args.length_scale_bps else LENGTH_SCALE])
            decay.append(x[1][DECAY])
            weighted_mean.append(x[1][WEIGHTED_MEAN])
            # total_effective_size.append(x[1][TOTAL_EFFECTIVE_SIZE])

        ax1.plot(length_scale, weighted_mean, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)
        if args.second_plot_max_length is not None:
            ax2.plot(length_scale, weighted_mean, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)
            for idx in range(1, len(weighted_mean)):
                if length_scale[idx] > args.second_plot_max_length:
                    y_lim_min_1 = min(y_lim_min_1, weighted_mean[idx-1])
                    break
            if args.third_plot_max_length is not None:
                ax3.plot(length_scale, weighted_mean, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)
                for idx in range(1, len(weighted_mean)):
                    if length_scale[idx] > args.third_plot_max_length:
                        y_lim_min_2 = min(y_lim_min_2, weighted_mean[idx-1])
                        break

        # ax2.plot(length_scale, total_effective_size, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)
        # ax3.plot(total_effective_size, weighted_mean, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)

    ax1.set_xscale('log')
    ax1.set_ylabel('Score')
    ax1.set_xlabel('Length (bp)' if args.length_scale_bps else 'Length (variant count)')

    if args.second_plot_max_length is not None:
        ax2.set_xscale('log')
        if args.third_plot_max_length is None:
            ax2.set_xlabel('Length (bp)' if args.length_scale_bps else 'Length (variant count)')
        ax2.set_xlim(0, args.second_plot_max_length)
        ax2.set_ylim(y_lim_min_1 * .999, 1.001)
        if args.third_plot_max_length is not None:
            ax3.set_xscale('log')
            ax3.set_xlabel('Length (bp)' if args.length_scale_bps else 'Length (variant count)')
            ax3.set_xlim(0, args.third_plot_max_length)
            ax3.set_ylim(y_lim_min_2 * .9999, 1.0001)

    plt.legend()
    fig.tight_layout()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show()


if __name__ == "__main__":
    main()
