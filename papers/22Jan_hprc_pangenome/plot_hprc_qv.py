#!/usr/bin/env python3
import argparse
from glob import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import pysam
import math
import pandas as pd
import haplotagging_stats
import os
import collections
import seaborn as sns

MAT_QV="mat_qv"
PAT_QV="pat_qv"


# cat Y1_assemblies_v2_genbank_QC.csv | sed 's/,/\t/g' | awk '{print $1,$16,$17,"Maternal","\n",$1,$19,$20,"Paternal"}' | sed 's/^ //' | sed 's/ $//' | sed 's/ /,/g' | sed 's/mat_//g' | sed 's/pat_//g' >Y1_assemblies_v2_genbank_QC.switch_hamming.csv
def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from haplotagging_stats tsv")
    parser.add_argument('--input_csv', '-i', dest='input_csv', default=None, required=True, type=str,
                        help='CSV file holding data')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default="HPRC_qv", required=False, type=str,
                       help='Figure name')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg):
    print(msg, file=sys.stderr)


def get_color(filename):
    if "maternal" in filename.lower():
        return "darkred"
    if "paternal" in filename.lower():
        return "darkblue"
    return "black"


def main():
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    print(df.head())

    sns.scatterplot(data=df, x=MAT_QV, y=PAT_QV)

    plt.ylabel("Paternal QV")
    plt.xlabel("Maternal QV")
    plt.tight_layout()
    # plt.set_size_inches(12, 12)
    #
    if args.figure_name is not None:
        plt.savefig(args.figure_name+".scatter.png", format='png', dpi=200)
        plt.savefig(args.figure_name+".scatter.pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()

    ax = sns.jointplot(data=df, x=MAT_QV, y=PAT_QV, hue="x")

    ax.ax_joint.legend([],[], frameon=False)
    ax.ax_joint.set_ylabel("Paternal QV")
    ax.ax_joint.set_xlabel("Maternal QV")
    plt.tight_layout()
    # plt.set_size_inches(12, 12)
    #
    if args.figure_name is not None:
        plt.savefig(args.figure_name+".joint.png", format='png', dpi=200)
        plt.savefig(args.figure_name+".joint.pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()


if __name__ == "__main__":
    main()
