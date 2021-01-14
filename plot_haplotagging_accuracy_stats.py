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
import seaborn as sns
import collections

plt.style.use('ggplot')
text_fontsize = 8
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

FILE_NAME = "file_name"
BLOCK_N50 = "block_n50"
ALL_SWITCH_RATE = "all_switch_rate"
ALL_HAMMING_RATE = "blockwise_hamming_rate"
UNCLASSIFIED_RATIO = 'unclassified_ratio'
CHROM="#chrom"
POS="position"

CLASSIFICATION_KEY = "classifications"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from haplotagging_stats tsv")
    parser.add_argument('--input', '-i', dest='input_stats_tsvs', default=None, required=True, type=str, action='append',
                        help='TSV file generated from "haplotagging_stats.py".  Can set multiple times, all values will be plotted')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str, action='store',
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')
    parser.add_argument('--plot_classified_not_accuracy', dest='plot_classified_not_accuracy', default=False, required=False, action='store_true',
                       help='Will plot the classified ratio instead of accuracy')
    parser.add_argument('--classification', '-c', dest='classification', default=None, required=True, type=str, action='append',
                        help='Specification of the form "Classification Name:/path/to/classification.bed".  If set, will not plot based on depth profiles')

    return parser.parse_args() if args is None else parser.parse_args(args)


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


def get_label(filename):
    if "margin" in filename.lower():
        return "Margin"
    if "whatshap" in filename.lower():
        return "Whatshap"
    return "Unknown"


def get_classifications(classification_list):
    classification_position_map=collections.defaultdict(lambda : collections.defaultdict(lambda : list()))
    classification_position_map[CLASSIFICATION_KEY] = list()

    for classification_specification in classification_list:
        parts = classification_specification.split(":")
        assert(len(parts) == 2)
        classification_name = parts[0].replace("\\n", "\n")
        classification_bed = parts[1]
        classification_position_map[CLASSIFICATION_KEY].append(classification_name)
        with open(classification_bed) as cin:
            for line in cin:
                if line.startswith("#"): continue
                bed_parts = line.split("\t")
                assert(len(bed_parts) >= 3)
                chr = bed_parts[0]
                start = int(bed_parts[1])
                end = int(bed_parts[2])
                if (end - start) < 1000: continue
                for i in range(start//1000, end//1000 + 1):
                    classification_position_map[chr][i*1000].append(classification_name)


    return classification_position_map


def plot_hist(dataframes, args):

    fig, ((ax1)) = plt.subplots(nrows=1,ncols=1)

    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        accuracies = list()
        for x in dataframe.iterrows():
            if x[1][haplotagging_stats.TOTAL_DEPTH] < 10: continue
            if x[1][haplotagging_stats.CORRECT_RATIO] <= 0: continue
            # accuracies.append(min(.99, x[1][haplotagging_stats.CORRECT_RATIO]))
            accuracies.append(x[1][haplotagging_stats.CORRECT_RATIO])

        n, bins, patches = ax1.hist(accuracies, 17, density=True, histtype='stepfilled', cumulative=False, label=get_label(filename), alpha=.5,
                 color=get_color(filename))

    # ax1.set_yscale('log')
    # ax2.set_yscale('log')
    plt.title("Haplotagging Accuracy")
    ax1.set_xlabel("Haplotagging Accuracy (1kb granularity)")
    ax1.set_ylabel("Relative Occurrence")
    plt.legend()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_stats_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show()


def plot_violin_with_depths(dataframes, args):

    merged_dataframe = []
    columns = ["Depth", "Accuracy", "Tool", "Classified Ratio"]
    total_depths = collections.defaultdict(lambda : 0)

    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        for x in dataframe.iterrows():
            if x[1][haplotagging_stats.TOTAL_DEPTH] == 0: continue

            # get depths
            td = x[1][haplotagging_stats.TOTAL_DEPTH]
            if td < 20:
                depth = "<20"
            elif 20 <= td < 40:
                depth = "20-40"
            elif 40 <= td < 60:
                depth = "40-60"
            elif 60 <= td < 80:
                depth = "60-80"
            elif 80 <= td < 100:
                depth = "80-100"
            elif td >= 100:
                depth = ">100"
            else:
                assert(False)

            # quit for non-classified ratio
            cr = x[1][haplotagging_stats.CORRECT_RATIO]
            if not args.plot_classified_not_accuracy and cr <= 0:
                continue

            # save rows
            row = [depth, cr, get_label(filename), 1.0 - x[1][UNCLASSIFIED_RATIO]]
            merged_dataframe.append(row)
            row = ["All", cr, get_label(filename), 1.0 - x[1][UNCLASSIFIED_RATIO]]
            merged_dataframe.append(row)

            # track total depths
            total_depths[depth] += 1

    total_rows = sum(total_depths.values())

    col_order = ["<20", "20-40", "40-60", "60-80", "80-100", ">100", "All"]
    df = pandas.DataFrame(merged_dataframe, columns=columns)
    ax = sns.violinplot(x="Depth", y="Classified Ratio" if args.plot_classified_not_accuracy else "Accuracy",
                        hue="Tool", data=df, split=True, scale="count", inner="quartile", order=col_order, cut=0)

    ax.set_xticklabels(list(map(lambda x: "{}\n{:.1f}%".format(x, 100 if x == "All" else 100.0 * total_depths[x] / total_rows), col_order)))

    plt.legend()
    plt.tight_layout()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_stats_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show()



def plot_violin_with_classifications(dataframes, args, classifications):

    merged_dataframe = []
    columns = ["Stratification", "Accuracy", "Tool", "Classified Ratio"]
    total_depths = collections.defaultdict(lambda : 0)

    plt.figure(figsize=(3.5, 3.5))
    colors = {}
    for i, filename in enumerate(dataframes.keys()):
        dataframe = dataframes[filename]
        colors[get_label(filename)] = get_color(filename)
        for x in dataframe.iterrows():
            if x[1][haplotagging_stats.TOTAL_DEPTH] == 0: continue

            # get position
            chr = x[1][CHROM]
            pos = x[1][POS]

            # quit for non-classified ratio
            cr = x[1][haplotagging_stats.CORRECT_RATIO]
            if not args.plot_classified_not_accuracy and cr <= 0:
                continue

            # save rows
            row = ["All", cr, get_label(filename), 1.0 - x[1][UNCLASSIFIED_RATIO]]
            merged_dataframe.append(row)

            for classification in classifications[chr][pos]:
                row = [classification, cr, get_label(filename), 1.0 - x[1][UNCLASSIFIED_RATIO]]
                merged_dataframe.append(row)

    classifications[CLASSIFICATION_KEY].append("All")
    df = pandas.DataFrame(merged_dataframe, columns=columns)
    ax = sns.violinplot(x="Stratification", y="Classified Ratio" if args.plot_classified_not_accuracy else "Accuracy",
                        hue="Tool", data=df, split=True, scale="count", palette=colors, linewidth=.75, bw=.1,
                        order=classifications[CLASSIFICATION_KEY], cut=0)
    # ax = sns.violinplot(x="Stratification", y="Classified Ratio" if args.plot_classified_not_accuracy else "Accuracy",
    #                     hue="Tool", data=df, split=True, scale="count", inner="quartile",
    #                     order=classifications[CLASSIFICATION_KEY], cut=0)
    ax.set_ylim(.79, 1.01)

    plt.legend()
    plt.tight_layout()

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(list(map(os.path.basename, args.input_stats_tsvs)))
        fig_name += ".png"
    if fig_name is not None:
        plt.savefig(fig_name + ".pdf", format='pdf', dpi=300 )
    plt.show()


def main():
    args = parse_args()
    dataframes = dict()
    for tsv in args.input_stats_tsvs:
        dataframes[os.path.basename(tsv)] = pandas.read_csv(tsv, sep='\t')

    # plot_hist(dataframes, args)
    if args.classification is None or len(args.classification) == 0:
        plot_violin_with_depths(dataframes, args)
    else:
        classification_data = get_classifications(args.classification)
        plot_violin_with_classifications(dataframes, args, classification_data)

if __name__ == "__main__":
    main()
