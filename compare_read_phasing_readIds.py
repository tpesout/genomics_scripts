#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import matplotlib.pyplot as plt

MARGIN_HAP1="hap1"
MARGIN_HAP2="hap2"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Compares phasing for reads haplotyped by margin")
    parser.add_argument('--margin_input_glob', '-i', dest='margin_input_glob', default=None, required=True, type=str,
                        help='Glob matching input files.  For each chunk, the only difference between fileNames '
                             'must be ' + MARGIN_HAP1 + ' and ' + MARGIN_HAP2)
    parser.add_argument('--truth_hap1', '-1', dest='truth_hap1', default=None, required=True, type=str,
                       help='Truth Hap1 readset')
    parser.add_argument('--truth_hap2', '-2', dest='truth_hap2', default=None, required=True, type=str,
                       help='Truth Hap2 readset')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')

    # parser.add_argument('--generic_stats', '-g', dest='generic_stats', action='store_true', default=False,
    #                     help='Print generic stats for all files')
    # parser.add_argument('--read_length', '-l', dest='read_length', action='store_true', default=False,
    #                     help='Print statistics on read length for all files')
    # parser.add_argument('--read_depth', '-d', dest='read_depth', action='store_true', default=False,
    #                     help='Print statistics on read depth for all files')
    # parser.add_argument('--genome_only', dest='genome_only', action='store_true', default=False,
    #                     help='Print only statistics for the whole genome (do not print stats for individual chroms)')
    # parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
    #                     help='Print histograms for length and depth')
    # parser.add_argument('--silent', '-V', dest='silent', action='store_true', default=False,
    #                     help='Print nothing')
    # parser.add_argument('--depth_spacing', '-s', dest='depth_spacing', action='store', default=1000, type=int,
    #                     help='How far to sample read data')
    # parser.add_argument('--depth_range', '-r', dest='depth_range', action='store', default=None,
    #                     help='Whether to only calculate depth within a range, ie: \'100000-200000\'')
    # parser.add_argument('--filter_secondary', dest='filter_secondary', action='store_true', default=False,
    #                     help='Filter secondary alignments out')
    # parser.add_argument('--filter_supplementary', dest='filter_supplementary', action='store_true', default=False,
    #                     help='Filter supplemenary alignments out')
    # parser.add_argument('--filter_read_length_min', dest='read_length_min', action='store', default=None, type=int,
    #                     help='Removes reads with length below this')
    # parser.add_argument('--filter_read_length_max', dest='read_length_max', action='store', default=None, type=int,
    #                     help='Removes reads with length above this')
    # parser.add_argument('--filter_alignment_threshold_min', dest='min_alignment_threshold', action='store',
    #                     default=None, type=int, help='Minimum alignment quality threshold')
    #
    # parser.add_argument('--produce_read_length_tsv', dest='read_length_tsv', action='store',
    #                     default=None, type=str, help='Produce a TSV with read lengths, named as this parameter')
    # parser.add_argument('--read_length_bucket_size', dest='read_length_bucket_size', action='store',
    #                     default=50000, type=int, help='Bucket size for read length TSV')

    return parser.parse_args() if args is None else parser.parse_args(args)

def log(msg):
    print(msg, file=sys.stderr)


def plottit(x, correct, haplotyped=None, avg_correct=None, avg_haplotyped=None, figName=None):
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Chunk IDX')
    ax1.set_ylabel('Phasing Accuracy', color=color)
    ax1.scatter(x, correct, color=color, s=5)
    ax1.plot(x, correct, color=color, linewidth=1, alpha=.5)
    ax1.tick_params(axis='y', labelcolor=color)
    if avg_correct is not None:
        ax1.plot(x, [avg_correct for _ in x], color=color, alpha=.5)

    if haplotyped is not None:
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('Truth Haplotyped Ratio', color=color)  # we already handled the x-label with ax1
        ax2.plot(x, haplotyped, color=color, linewidth=1, alpha=.3)
        ax2.tick_params(axis='y', labelcolor=color)
        if avg_haplotyped is not None:
            ax2.plot(x, [avg_haplotyped for _ in x], color=color, alpha=.3)
    else:
        ax1.set_ylim(.45,1.05)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.set_size_inches(18, 6)
    if figName is not None:
        plt.savefig(figName)
    plt.show()

def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # get truth reads
    truth_h1 = set()
    truth_h2 = set()
    with open(args.truth_hap1, 'r') as fin:
        for line in fin:
            truth_h1.add(line.split(',')[0].strip())
    with open(args.truth_hap2, 'r') as fin:
        for line in fin:
            truth_h2.add(line.split(',')[0].strip())
    log("Found {} truth H1 reads and {} truth H2 reads".format(len(truth_h1), len(truth_h2)))
    truth_all = truth_h1.union(truth_h2)

    # get read phasing pairs
    margin_hap_files = glob.glob(args.margin_input_glob)
    margin_haps = list()
    if len(margin_hap_files) == 0:
        log("Found no files matching {}".format(args.margin_input_glob))
        sys.exit(1)
    if len(margin_hap_files) % 2 != 0:
        log("Found odd number of files ({}) matching {}".format(len(margin_hap_files), args.margin_input_glob))
        sys.exit(1)
    margin_hap_files.sort()
    margin_hap_files.reverse()
    while len(margin_hap_files) > 0:
        hapA = margin_hap_files.pop()
        hapB = margin_hap_files.pop()
        if hapA.replace(MARGIN_HAP2, '').replace(MARGIN_HAP1, '') != hapB.replace(MARGIN_HAP2, '').replace(MARGIN_HAP1, ''):
            log("Found unmatched haplotype sets:\n\t{}\n\t{}".format(hapA, hapB))
            sys.exit(1)
        margin_haps.append([hapA, hapB])

    chunk_idxs = []
    correct_raitos = []
    haplotyped_ratios = []

    total_right = 0
    total_wrong = 0
    total_unknown = 0
    idx = 0
    for hap in margin_haps:
        file_hapA = hap[0]
        file_hapB = hap[1]
        reads_hapA = set()
        reads_hapB = set()

        # read files
        with open(file_hapA) as fin:
            next(fin)
            for line in fin:
                reads_hapA.add(line.split(',')[0].strip())
        with open(file_hapB) as fin:
            next(fin)
            for line in fin:
                reads_hapB.add(line.split(',')[0].strip())

        inter_H1_HA = truth_h1.intersection(reads_hapA)
        inter_H2_HA = truth_h2.intersection(reads_hapA)
        inter_H1_HB = truth_h1.intersection(reads_hapB)
        inter_H2_HB = truth_h2.intersection(reads_hapB)
        untagged_HA = reads_hapA.difference(truth_all)
        untagged_HB = reads_hapB.difference(truth_all)

        H1_HA = len(inter_H1_HA)
        H2_HA = len(inter_H2_HA)
        HU_HA = len(untagged_HA)
        H1_HB = len(inter_H1_HB)
        H2_HB = len(inter_H2_HB)
        HU_HB = len(untagged_HB)

        cis_support = H1_HA + H2_HB
        trans_support = H1_HB + H2_HA

        if (cis_support > trans_support):
            right_count = cis_support
            wrong_count = trans_support
        else:
            right_count = trans_support
            wrong_count = cis_support
        unknown_count = HU_HA + HU_HB

        total_right += right_count
        total_wrong += wrong_count
        total_unknown += unknown_count

        right = 1.0 * right_count / max(right_count + wrong_count, 1)
        haplotyped_ratio= 1.0 * (right_count + wrong_count) / max(right_count + wrong_count + unknown_count, 1)


        print("C{:05d} {:s}\tcorrect_ratio: {:1.6f} \tright: {:5d} \twrong: {:5d} \tunknown: {:5d} \thaplotyped_ratio: {:1.6f} \ttotal_hapA: {:5d} \ttotal_hapB: {:5d}\t\t{:s}".format(
            idx, "*" if right < .9 else " ", right, right_count, wrong_count, unknown_count, haplotyped_ratio, len(reads_hapA), len(reads_hapB), file_hapA.replace(MARGIN_HAP1, '').replace(MARGIN_HAP2, '')))

        chunk_idxs.append(idx)
        correct_raitos.append(right)
        haplotyped_ratios.append(haplotyped_ratio)
        idx += 1

    right = 1.0 * total_right / max(total_right + total_wrong, 1)
    haplotyped_ratio = 1.0 * (total_right + total_wrong) / max(total_right + total_wrong + total_unknown, 1)

    print()
    print(
        "TOTAL \tcorrect_ratio: {:1.6f} \tright: {:5d} \twrong: {:5d} \tunknown: {:5d} \thaplotyped_ratio: {:1.6f}".format(
            right, total_right, total_wrong, total_unknown, haplotyped_ratio))

    # plottit(chunk_idxs, correct_raitos, None, right, None, args.figure_name)
    plottit(chunk_idxs, correct_raitos, haplotyped_ratios, right, haplotyped_ratio, args.figure_name)



if __name__ == "__main__":
    main()











