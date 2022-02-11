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

NUM_CONTIGS="num_contigs"
TOTAL_LEN="total_len"
HAPLOTYPE="haplotype"
HAPLO_SEX="haplotype-sex"
SEX="sex"

# cat Y1_assemblies_v2_genbank_QC.csv | sed 's/,/\t/g' | awk '{print $1,$2,$3,"Maternal","\n",$1,$6,$7,"Paternal"}' | sed 's/^ //' | sed 's/ $//' | sed 's/ /,/g' | sed 's/mat_//g' | sed 's/pat_//g' >Y1_assemblies_v2_genbank_QC.contig_stats.csv
# cat Y1_assemblies_v2_genbank_QC.full.csv | sed 's/,/\t/g' | awk '{print $1,$2,$3,"Maternal","Maternal-",$23,$23,"\n",$1,$6,$7,"Paternal","Paternal-",$23,$23}' | sed 's/- /-/g' | sed 's/^ //' | sed 's/ $//' | sed 's/ /,/g' | sed 's/mat_//g' | sed 's/pat_//g' >Y1_assemblies_v2_genbank_QC.contig_stats.csv
def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from haplotagging_stats tsv")
    parser.add_argument('--input_csv', '-i', dest='input_csv', default=None, required=True, type=str,
                        help='CSV file holding data')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default="HPRC_contig_stats", required=False, type=str,
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
    # sns.set_palette(sns.color_palette(["darkred", "darkblue"]))

    # sns.boxplot(x=HAPLOTYPE, y=NUM_CONTIGS, data=df)#, palette={"Maternal":"darkred","Paternal":"darkblue"})
    # spax = sns.swarmplot(x=HAPLOTYPE, y=NUM_CONTIGS, hue=SEX, data=df, palette={"Female":"fuchsia","Male":"cyan"}) #color="fuchsia")

    sns.boxplot(x=HAPLO_SEX, y=NUM_CONTIGS, data=df,  order=["Maternal-Female", "Maternal-Male", "Paternal-Female", "Paternal-Male"],
                palette={"Maternal-Male":"darkred","Maternal-Female":"darkred","Paternal-Male":"darkblue","Paternal-Female":"darkblue"})
    spax = sns.swarmplot(x=HAPLO_SEX, y=NUM_CONTIGS, data=df,  order=["Maternal-Female", "Maternal-Male", "Paternal-Female", "Paternal-Male"],
                         palette={"Maternal-Male":"royalblue","Maternal-Female":"crimson","Paternal-Male":"royalblue","Paternal-Female":"crimson"})


    plt.title("")
    plt.ylabel("Contig Count")
    plt.xlabel("Haplotype")
    plt.tight_layout()
    # plt.set_size_inches(12, 12)
    #
    if args.figure_name is not None:
        plt.savefig(args.figure_name+".contig_count.png", format='png', dpi=200)
        plt.savefig(args.figure_name+".contig_count.pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()

    # sns.boxplot(x=HAPLOTYPE, y=TOTAL_LEN, data=df)#, palette={"Maternal":"darkred","Paternal":"darkblue"})
    # spax = sns.swarmplot(x=HAPLOTYPE, y=TOTAL_LEN, hue=SEX, data=df,  palette={"Female":"fuchsia","Male":"cyan"}) #color="fuchsia")

    sns.boxplot(x=HAPLO_SEX, y=TOTAL_LEN, data=df, order=["Maternal-Female", "Maternal-Male", "Paternal-Female", "Paternal-Male"],
                palette={"Maternal-Male":"darkred","Maternal-Female":"darkred","Paternal-Male":"darkblue","Paternal-Female":"darkblue"})
    spax = sns.swarmplot(x=HAPLO_SEX, y=TOTAL_LEN, data=df, order=["Maternal-Female", "Maternal-Male", "Paternal-Female", "Paternal-Male"],
                palette={"Maternal-Male":"royalblue","Maternal-Female":"crimson","Paternal-Male":"royalblue","Paternal-Female":"crimson"})



    plt.title("")
    plt.ylabel("Total Length")
    plt.xlabel("Haplotype")
    plt.tight_layout()
    # plt.set_size_inches(12, 12)
    #
    if args.figure_name is not None:
        plt.savefig(args.figure_name+".total_len.png", format='png', dpi=200)
        plt.savefig(args.figure_name+".total_len.pdf", format='pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
