#!/usr/bin/env python3
import argparse
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pysam
import pandas
import seaborn as sns
import os
from collections import defaultdict as defaultdict
import math

text_fontsize = 6
plt.rcParams.update({'font.size': text_fontsize})
plt.style.use('ggplot')
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

TYPE="Type"
GENE_REGION="Subset"
F1_SCORE="METRIC.F1_Score"
PRECISION="METRIC.Precision"
RECALL="METRIC.Recall"
REGION_SIZE="Subset.Size"
REGION_HICONF_SIZE="Subset.IS_CONF.Size"
TP_T="TRUTH.TP"
TP_Q="QUERY.TP"
FN="TRUTH.FN"
FP="QUERY.FP"

COLUMNS=[TYPE, GENE_REGION, RECALL, PRECISION, F1_SCORE, REGION_SIZE, REGION_HICONF_SIZE, TP_T, TP_Q, FN, FP]


def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots custom information")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                        help='TSV file with desired hap.py results')
    parser.add_argument('--high_conf_threshold', '-t', dest='high_conf_threshold', default=.8, required=False, type=float,
                        help='High confidence coverage threshold below which results will be discarded')

    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')

    return parser.parse_args() if args is None else parser.parse_args(args)


def plot_manual_data(df, figsize=(3.5, 3.5), figure_name=None):
    fig, ((ax1, ax2)) = plt.subplots(nrows=2, ncols=1, figsize=figsize)

    plt.legend()
    plt.tight_layout()

    if figure_name is not None:
        if not figure_name.endswith(".png"):
            figure_name += ".png"
        figure_name += ".pdf"
        fig.savefig(figure_name, format='pdf', dpi=300)
    plt.show()


def main():
    args = parse_args()

    data = defaultdict(lambda: defaultdict(lambda: 0.0))
    column_map = None
    with open(args.input) as tsv_in:
        for line in tsv_in:
            line = line.strip().split("\t")
            if column_map is None:
                column_map = {h:i for (i,h) in filter(lambda x: x[1] in COLUMNS, enumerate(line))}
                continue
            line_data = {c:line[column_map[c]] for c in COLUMNS}
            size = int(float(line_data[REGION_SIZE]))
            hc_size = int(float(line_data[REGION_HICONF_SIZE]))
            tp_t = int(float(line_data[TP_T]))
            tp_q = int(float(line_data[TP_Q]))
            fp = int(float(line_data[FP]))
            fn = int(float(line_data[FN]))
            id = line_data[GENE_REGION]
            data[id][REGION_SIZE] = size
            data[id][REGION_HICONF_SIZE] = hc_size
            data[id][TP_T] += tp_t
            data[id][TP_Q] += tp_q
            data[id][FP] += fp
            data[id][FN] += fn
            # if line_data[TP] != line_data[TP2] and line_data[F1_SCORE] not in ("1.0", ""):
            #     pass

    total_genes = 0
    total_hiconf_sequence = 0
    below_hiconf_threshold = 0
    above_hiconf_threshold = 0
    had_no_variants = 0
    had_variants = 0
    had_variants_but_no_errors = 0
    no_errors = 0
    for (id, gene_data) in data.items():
        total_genes += 1
        if gene_data[REGION_HICONF_SIZE]/max(gene_data[REGION_SIZE], 1) < args.high_conf_threshold:
            below_hiconf_threshold += 1
            continue
        above_hiconf_threshold += 1
        total_hiconf_sequence += gene_data[REGION_HICONF_SIZE]
        if gene_data[TP_T] + gene_data[TP_Q] == 0:
            had_no_variants += 1
        else:
            had_variants += 1
            if gene_data[FP] + gene_data[FN] == 0:
                had_variants_but_no_errors += 1
        if gene_data[FP] + gene_data[FN] == 0:
            no_errors += 1

    print("Stats for {}:".format(args.input))
    print("\tGene bodies:              {} ".format(total_genes))
    print("\t  Above {:.2f} HighConf:    {} ({:.5f})".format(args.high_conf_threshold, above_hiconf_threshold,
                                                                 above_hiconf_threshold/total_genes))
    print("\t    Total HighConf Seq:   {}".format(total_hiconf_sequence))
    print("\t    Without variants:     {} ({:.5f})".format(had_no_variants, had_no_variants/above_hiconf_threshold))
    print("\t    With variants:        {} ({:.5f})".format(had_variants, had_variants/above_hiconf_threshold))
    print("\t      With no errors:     {} ({:.5f})".format(had_variants_but_no_errors, had_variants_but_no_errors/had_variants))
    print("\t    With no errors:       {} ({:.5f})".format(no_errors, no_errors/above_hiconf_threshold))





    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(os.path.basename(args.input).split(".")[:-1])


if __name__ == "__main__":
    main()
