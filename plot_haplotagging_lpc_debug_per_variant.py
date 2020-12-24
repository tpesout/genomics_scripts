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
import seaborn as sns

DECAY = "decay"
LENGTH_SCALE = "length_scale_num_vars"
LENGTH_SCALE_BP = "length_scale_bps"
WEIGHTED_MEAN = "weighted_mean"
APPROX = "approx_"
TOTAL_EFFECTIVE_SIZE = "total_eff_size"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from margin's calcLocalPhasingCorrectness ")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str, action='store',
                        help='TSV file generated from calcLocalPhasingCorrectness.  Can set multiple times, all values will be plotted')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')
    parser.add_argument('--length_scale_bps', '-b', dest='length_scale_bps', default=False, required=False, action='store_true',
                       help='Lenght scale should use {} instead of {}'.format(LENGTH_SCALE_BP, LENGTH_SCALE))

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
    dataframes[os.path.basename(args.input)] = pandas.read_csv(args.input, sep='\t')

    df_records = list()
    length_scale = "Length Scale (bp)" if args.length_scale_bps else "Length Scale (Var)"
    df_columns = ["Decay", length_scale, "Score", "Count", "Log Count"]


    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        for p, x in enumerate(dataframe.iterrows()):
            decay = x[1].values[0]
            if decay == 0: continue
            if math.isinf(x[1].values[1]): continue
            # length_bp = float("{:.3f}".format(x[1].values[1]))
            # length_var = float("{:.3f}".format(x[1].values[2]))
            # length_var = int(x[1].values[1])
            # length_bp = int(x[1].values[2])
            # length_var = x[1].values[1]
            # length_bp = x[1].values[2]
            # length = length_bp if args.length_scale_bps else length_var
            # if length <= 10:
            length_var = float("{:.3f}".format(x[1].values[1]))
            length_bp = float("{:.3f}".format(x[1].values[2]))
            length = length_bp if args.length_scale_bps else length_var

            scores = collections.defaultdict(lambda: 0)
            for i in range(3, x[1].size):
                score = x[1].values[i]
                if math.isnan(score): continue
                score = int(score * 100)
                scores[score] += 1

            for i in range(0, 101):
                df_records.append([decay, length, i, scores[i], 0 if scores[i] == 0 else math.log(scores[i])])

    df = pandas.DataFrame(df_records, columns=df_columns)
    df_pivoted = df.pivot(length_scale, "Score", "Log Count")


    plt.figure(figsize=(18, 12))
    sns.heatmap(df_pivoted)

    # plt.legend()
    plt.tight_layout()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show()


if __name__ == "__main__":
    main()
