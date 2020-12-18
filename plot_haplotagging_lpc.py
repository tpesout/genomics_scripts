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

    return parser.parse_args() if args is None else parser.parse_args(args)


def main():
    args = parse_args()
    dataframes = dict()
    for tsv in args.input_tsvs:
        dataframes[os.path.basename(tsv)] = pandas.read_csv(tsv, sep='\t')

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2,ncols=2)

    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        length_scale = list()
        decay = list()
        weighted_mean = list()
        for x in dataframe.iterrows():
            length_scale.append(x[1][LENGTH_SCALE_BP if args.length_scale_bps else LENGTH_SCALE])
            decay.append(x[1][DECAY])
            weighted_mean.append(x[1][WEIGHTED_MEAN])

        ax1.plot(length_scale, weighted_mean, label=filename)
        ax2.plot(length_scale, weighted_mean, label=filename)
        ax3.plot(decay, weighted_mean, label=filename)
        ax4.plot(decay, weighted_mean, label=filename)

    # ax1.set_yscale('log')
    # ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax4.set_xscale('log')

    ax1.set_ylabel('Length Scale')
    ax3.set_ylabel('Decay')
    ax3.set_xlabel('Linear')
    ax4.set_xlabel('Logarithmic')

    plt.legend()
    # fig.tight_layout()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show()


if __name__ == "__main__":
    main()
