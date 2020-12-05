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

FILE_NAME = "file_name"
BLOCK_N50 = "block_n50"
ALL_SWITCH_RATE = "all_switch_rate"
ALL_HAMMING_RATE = "blockwise_hamming_rate"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from tpesout's whatshap wdl workflow")
    parser.add_argument('--input', '-i', dest='input_stats_tsvs', default=None, required=True, type=str, action='append',
                        help='TSV file generated from "haplotagging_stats.py".  Can set multiple times, all values will be plotted')


    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')

    return parser.parse_args() if args is None else parser.parse_args(args)


def main():
    args = parse_args()
    dataframes = dict()
    for tsv in args.input_stats_tsvs:
        dataframes[os.path.basename(tsv)] = pandas.read_csv(tsv, sep='\t')

    fig, ((ax1, ax2)) = plt.subplots(nrows=1,ncols=2)

    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        accuracies = collections.defaultdict(lambda: 0)
        unclassified = collections.defaultdict(lambda: 0)
        for x in dataframe.iterrows():
            if x[1][haplotagging_stats.TOTAL_DEPTH] < 10: continue
            correct_ratio = int(100 * x[1][haplotagging_stats.CORRECT_RATIO])
            accuracies[correct_ratio] += 1
            unclassified_ratio = int(100 * x[1][haplotagging_stats.UNCLASSIFIED_RATIO])
            unclassified[unclassified_ratio] += 1

        ax1.bar(list(range(50,101)), list(accuracies[x] for x in range(50, 101)), label=filename, alpha=.5, color='red' if i == 0 else 'blue')
        ax2.bar(list(range(0,101)), list(unclassified[x] for x in range(0, 101)), label=filename, alpha=.5, color='red' if i == 0 else 'blue')

    ax1.set_yscale('log')
    ax2.set_yscale('log')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
