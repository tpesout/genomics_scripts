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

plt.style.use('ggplot')
text_fontsize = 6
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

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
    elif ".margin_phased" in filename.lower():
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


def plot_with_parts(dataframes, fig_name, figsize, args, max_length=None, legend=False, xlabel=False, ylabel=False):
    fig, ((ax1)) = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    y_lim_min = 1.0
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
        if max_length is not None:
            for idx in range(1, len(weighted_mean)):
                if length_scale[idx] > max_length:
                    y_lim_min = min(y_lim_min, weighted_mean[idx-1])
                    break

        # ax2.plot(length_scale, total_effective_size, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)
        # ax3.plot(total_effective_size, weighted_mean, label=get_label(filename), color=get_color(filename), linestyle=get_linestyle(filename), alpha=.5)

    ax1.set_xscale('log')
    if ylabel:
        ax1.set_ylabel('Score')
    if xlabel:
        ax1.set_xlabel('Length (bp)' if args.length_scale_bps else 'Length (variant count)')

    if max_length is not None:
        ax1.set_xlim(0, max_length)
        if max_length <= 10:
            ax1.set_ylim(y_lim_min * .9999, 1.0001)
        elif max_length <= 100000:
            ax1.set_ylim(y_lim_min * .999, 1.001)

    fig.tight_layout()
    if legend:
        ax1.legend()

    if fig_name is not None:
        fig_name += ".pdf"
        plt.savefig(fig_name, format='pdf', dpi=300)
    plt.show()


def main():
    args = parse_args()
    dataframes = dict()
    for tsv in args.input_tsvs:
        dataframes[os.path.basename(tsv)] = pandas.read_csv(tsv, sep='\t')

    if (args.second_plot_max_length is not None):
        plot_with_parts(dataframes, args.figure_name + ".1", (2.5, 3.5), args,
                        max_length=None, legend=True, xlabel=True, ylabel=True)
        if args.third_plot_max_length is not None:
            plot_with_parts(dataframes, args.figure_name + ".2", (2.5, 1.75), args,
                            max_length=args.second_plot_max_length, legend=False, xlabel=True, ylabel=True)
            plot_with_parts(dataframes, args.figure_name + ".3", (2.5, 1.75), args,
                            max_length=args.third_plot_max_length, legend=False, xlabel=False, ylabel=True)
        else:
            plot_with_parts(dataframes, args.figure_name + ".2", (2.5, 3.5), args,
                            max_length=args.second_plot_max_length, legend=False, xlabel=True, ylabel=False)
    else:
        plot_with_parts(dataframes, args.figure_name, (3.5, 3.5), args,
                        max_length=None, legend=True, xlabel=True, ylabel=True)
    sys.exit()

    if args.second_plot_max_length is None:
        fig, ((ax1)) = plt.subplots(nrows=1, ncols=1, figsize=(3.75, 3.5))
    elif args.third_plot_max_length is None:
        fig, ((ax1, ax2)) = plt.subplots(nrows=1,ncols=2, figsize=(12, 6))
    else:
        # fig = plt.figure(figsize=(12, 6))
        fig = plt.figure(constrained_layout=True, figsize=(4.5, 3.5))
        gs = fig.add_gridspec(2, 2)
        ax1 = fig.add_subplot(gs[0:2, 1:2])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[0, 0])

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
            ax2.set_xlabel('Length (bp)' if args.length_scale_bps else 'Length (variant count)')
            ax3.set_xlim(0, args.third_plot_max_length)
            ax3.set_ylim(y_lim_min_2 * .9999, 1.0001)

    ax1.legend()
    fig.tight_layout()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        fig_name += ".pdf"
        ax1.savefig(fig_name, format='pdf', dpi=300)
        ax2.savefig(fig_name, format='pdf', dpi=300)
        ax3.savefig(fig_name, format='pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
