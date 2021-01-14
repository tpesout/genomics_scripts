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
import math

text_fontsize = 6
plt.rcParams.update({'font.size': text_fontsize})
plt.style.use('ggplot')
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

DATA="Data"
TOOL="Tool"
MODULE="Module"
RUNTIME="Runtime (min)"
COST="Cost"
TOOL_MODULE="Tool Module"

Margin_COLOR = [(238 / 256.0, 39 / 256.0, 8 / 256.0, 0.27), (238 / 256.0, 39 / 256.0, 8 / 256.0, 0.8),
                (238 / 256.0, 39 / 256.0, 8 / 256.0, 1.0)]
WhatsHap_COLOR = [(10 / 256.0, 118 / 256.0, 148 / 256.0, 0.27), (10 / 256.0, 118 / 256.0, 148 / 256.0, 0.8),
                  (10 / 256.0, 118 / 256.0, 148 / 256.0, 1.0)]
COLOR_PALLATE={
    "MP Haplotag": Margin_COLOR[0],
    "MP Phase": [x*.7 for x in Margin_COLOR[0]],
    "WH Haplotag": WhatsHap_COLOR[0],
    "WH Phase": [x*.7 for x in WhatsHap_COLOR[0]],
}


def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots custom information")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                        help='TSV file describing runtime info')

    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_input', '-F', dest='figure_name_input', default=False, required=False, action='store_true',
                       help='Figure name should be based off of input file (-f overrides this)')

    return parser.parse_args() if args is None else parser.parse_args(args)


def plot_runtime(df, figsize=(3.5, 3.5), figure_name=None):
    fig, ((ax1, ax2)) = plt.subplots(nrows=2, ncols=1, figsize=figsize)

    g = sns.barplot(x=RUNTIME, y=DATA, hue=TOOL_MODULE, ax=ax1, data=df, palette=COLOR_PALLATE)
    g.legend_.remove()
    g = sns.barplot(x=COST, y=DATA, hue=TOOL_MODULE, ax=ax2, data=df, palette=COLOR_PALLATE)
    g.legend_.remove()

    ax1.set_ylabel("")
    ax2.set_ylabel("")
    ax2.xaxis.set_major_formatter(mtick.StrMethodFormatter('${x:.2f}'))

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

    df = pandas.read_csv(args.input, sep='\t')
    df[TOOL_MODULE] = df.apply(lambda x: "{} {}".format("MP" if x[TOOL] == "Margin" else "WH", x[MODULE]), axis=1)

    fig_name = args.figure_name
    if fig_name is None and args.figure_name_input:
        fig_name = ".".join(os.path.basename(args.input).split(".")[:-1])

    plot_runtime(df, figure_name=fig_name)

if __name__ == "__main__":
    main()
