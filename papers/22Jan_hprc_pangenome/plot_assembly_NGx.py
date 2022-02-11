#!/usr/bin/env python3
import argparse
from glob import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import pysam
import math
import pandas
import haplotagging_stats
import os
import collections

REF_LINEWIDTH=2
LINEWIDTH=1
REF_ALPHA=1.0
ALPHA=.33

MALE_SAMPLES=[
    "HG002",
    "HG005",
    "HG00621",
    "HG00673",
    "HG01106",
    "HG01109",
    "HG01243",
    "HG01258",
    "HG01358",
    "HG01928",
    "HG01952",
    "HG02055",
    "HG02145",
    "HG02486",
    "HG02572",
    "HG02717",
    "HG03098",
    "HG03492",
    "HG03579"]

#NG50, black line smaller, include


def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from haplotagging_stats tsv")
    parser.add_argument('--input_fai_glob', '-i', dest='input_fai_glob', default=None, required=True, type=str,
                        help='Fasta indices')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name')
    parser.add_argument('--reference_fai', '-r', dest='reference_fai', required=False, default=list(), type=str,
                       action="append", help='Reference faidx')
    parser.add_argument('--reference_name', '-R', dest='reference_name', required=False, default=list(), type=str,
                        action="append", help='Reference name, if specified will label reference lines')
    parser.add_argument('--genome_size', '-g', dest='genome_size', required=False, default=3100000000, type=int,
                       help='Genome size in bp')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg):
    print(msg, file=sys.stderr)


def get_length_mod(max_length):
    if max_length < 1000:
        return 1.0, "bp"
    if max_length < 1000000:
        return 1000.0, "kb"
    if max_length < 1000000000:
        return 1000000, "Mb"
    return 1000000000.0, "Gb"


def is_male(filename):
    for m in MALE_SAMPLES:
        if m in filename: return True
    return False


def get_color(filename):
    if "paternal" in filename:
        if is_male(filename):
            return "mediumblue"
        return "darkorchid"
    if "maternal" in filename:
        return "darkorchid"
    return "black"


# def get_color(filename):
#     if "paternal" in filename:
#         return "darkblue"
#     if "maternal" in filename:
#         return "darkred"
#     return "black"


def main():
    args = parse_args()
    files = glob(args.input_fai_glob)
    if len(files) == 0:
        log("Found no files matching {}".format(args.input_fai_glob))
        sys.exit()
    else:
        log("Found {} files matching {}".format(len(files), args.input_fai_glob))

    # get reference faidx if appropriate
    for ref in args.reference_fai:
        tmp = [ref]
        tmp.extend(files)
        files = tmp

    if len(args.reference_name) != len(args.reference_fai):
        log("Got {} reference fais and {} reference names".format(len(args.reference_fai), len(args.reference_name)))

    # for tracking all lengths
    all_file_lengths = list()
    max_contig_length = 0
    index_to_filename = dict()


    # get all read data
    for i, file in enumerate(files):
        file_lengths = list()
        all_file_lengths.append(file_lengths)
        index_to_filename[i] = file
        with open(file) as file_in:
            for line in file_in:
                parts = line.split()
                length = int(parts[1])
                file_lengths.append(length)
                if length > max_contig_length:
                    max_contig_length = length

    # setup
    fig, ((ax1)) = plt.subplots(nrows=1,ncols=1)
    n50_size = args.genome_size / 2
    n50s = list()
    read_length_mod, read_length_mod_id = get_length_mod(max_contig_length)
    cumulative_length_mod, cumulative_length_mod_id = get_length_mod(args.genome_size)

    log("\nPlotting Contig Lengths")
    for i, file_lengths in enumerate(all_file_lengths):
        is_reference = i < len(args.reference_fai)
        file_lengths.sort(reverse=True)
        total = 0
        n50 = None
        prev_length = None
        for length in file_lengths:

            # plot
            new_total = total+length
            if prev_length is not None and prev_length != length:
                ax1.vlines(total / cumulative_length_mod, prev_length / read_length_mod, length / read_length_mod,
                           alpha=REF_ALPHA if is_reference else ALPHA, color=get_color(index_to_filename[i]), linewidth=REF_LINEWIDTH if is_reference else LINEWIDTH)
            ax1.hlines(length / read_length_mod, total / cumulative_length_mod, new_total / cumulative_length_mod,
                       alpha=REF_ALPHA if is_reference else ALPHA, color=get_color(index_to_filename[i]), linewidth=REF_LINEWIDTH if is_reference else LINEWIDTH)

            # iterate
            prev_length = length
            total = new_total
            if n50 is None and total >= n50_size:
                n50 = length
                if not is_reference:
                    n50s.append(n50)
                    log("\tFile {} got NG50 of {}".format(index_to_filename[i], n50))

        # finish
        if prev_length is not None:
            ax1.vlines(total / cumulative_length_mod, prev_length / read_length_mod, 0,
                       alpha=REF_ALPHA if is_reference else ALPHA, color=get_color(index_to_filename[i]),
                       linewidth=REF_LINEWIDTH if is_reference else LINEWIDTH)
        if total >= args.genome_size:
            log("\tFile {} has total length {} >= genome size {}".format(index_to_filename[i], total, args.genome_size))

        #n50 edge case
        if n50 is None and not is_reference:
            n50s.append(0)
            log("\tFile {} got N50 of {}".format(index_to_filename[i], 0))
        if is_reference and len(args.reference_name) > i:
            ax1.annotate(list(reversed(args.reference_name))[i], (0, file_lengths[0] / read_length_mod + 1), fontfamily='monospace', fontsize=12,weight="bold")

    log("")
    avg_n50 = np.mean(n50s)
    log("Average N50: {}".format(avg_n50))
    log("Genome Size: {}".format(args.genome_size))
    ax1.vlines(n50_size / cumulative_length_mod, max_contig_length / read_length_mod + 1, 0, color="black", linestyles="dashed",
               linewidth=LINEWIDTH)
    n50_mod, n50_mod_id = get_length_mod(avg_n50)
    ax1.annotate("Avg NG50: {:d} {}".format(int(avg_n50 / n50_mod), n50_mod_id),
                 (1.1 * n50_size / cumulative_length_mod, max_contig_length / read_length_mod), fontfamily='monospace', fontsize=12,weight="bold")

    ax1.ticklabel_format(axis='both', style='plain')
    # ax1.set_ylim(0,100)
    ax1.set_ylabel("Contig Length ({})".format(read_length_mod_id))
    ax1.set_xlabel("Cumulative Coverage ({})".format(cumulative_length_mod_id))

    fig.tight_layout()
    fig.set_size_inches(12, 12)

    if args.figure_name is not None:
        plt.savefig(args.figure_name+".png", format='png', dpi=200)
        plt.savefig(args.figure_name+".pdf", format='pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
