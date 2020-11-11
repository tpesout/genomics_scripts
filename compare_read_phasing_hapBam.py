#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import matplotlib.pyplot as plt
import pysam
import collections
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr

HP_TAG = "HP"
UNCLASSIFIED = 'u'
CORRECT = 'c'
INCORRECT = 'i'
HETS = 'h'
FP = 'p'
FN = 'n'
SPACING = 1000

def parse_args(args = None):
    parser = argparse.ArgumentParser("Compares phasing for reads haplotyped by margin")
    parser.add_argument('--margin_input_bam', '-i', dest='margin_input_bam', default=None, required=True, type=str,
                        help='BAM file haplotagged by margin. ')
    parser.add_argument('--truth_hap1', '-1', dest='truth_hap1', default=None, required=True, type=str,
                       help='Truth Hap1 readset')
    parser.add_argument('--truth_hap2', '-2', dest='truth_hap2', default=None, required=True, type=str,
                       help='Truth Hap2 readset')
    parser.add_argument('--het_vcf', '-v', dest='het_vcf', default=None, required=False, type=str,
                       help='VCF containing only heterozygous sites')
    parser.add_argument('--result_vcf', '-r', dest='result_vcf', default=None, required=False, type=str,
                       help='VCF containing TP/FP/FN classified sites')
    parser.add_argument('--chunks', '-c', dest='chunks', default=None, required=False, type=str,
                       help='File describing chunk positions')
    parser.add_argument('--phaseset_bed', '-p', dest='phaseset_bed', default=None, required=False, type=str,
                       help='File describing phase sets')
    parser.add_argument('--title', '-t', dest='title', default=None, required=False, type=str,
                       help='Figure title')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_bam', '-F', dest='figure_name_bam', default=False, required=False, action='store_true',
                       help='Figure name should be based off of BAM file (-f overrides this)')
    parser.add_argument('--untagged_only', '-u', dest='untagged_only', default=False, required=False, action='store_true',
                        help='Only plot untagged reads')

    return parser.parse_args() if args is None else parser.parse_args(args)

def log(msg):
    print(msg, file=sys.stderr)


def smooth_values(values,size=5,should_sum=False):
    new_values = []
    for i in range(0,len(values)):
        s=max(0,i-size)
        e=min(len(values),i+size)
        new_values.append(sum(values[s:e]) if should_sum else np.mean(values[s:e]))
    return new_values


def plottit(classification_data, figName=None, phasesets=None, has_het_vcf=False, has_result_vcf=False, chunk_boundaries=None,
            title=None):
    start_idx = min(classification_data.keys())
    end_idx = max(classification_data.keys())
    x = []
    right = []
    rong = []
    total = []
    unclassified = []
    correct_ratio = []
    raw_hets = []
    fp = []
    fn = []
    for i in range(start_idx, end_idx + 1):
        x.append(i)
        ri = classification_data[i][CORRECT]
        ro = classification_data[i][INCORRECT]
        un = classification_data[i][UNCLASSIFIED]
        ht = classification_data[i][HETS]

        right.append(ri)
        rong.append(-1 * ro)
        unclassified.append(un)
        total.append(ri + ro + un)
        correct_ratio.append(None if ri + ro == 0 else 100*abs(max(ri,ro)/(ri + ro)))
        raw_hets.append(ht)
        fp.append(classification_data[i][FP])
        fn.append(classification_data[i][FN])
    unclassified = smooth_values(unclassified, size=20)
    total = smooth_values(total, size=20)
    hets = smooth_values(raw_hets, size=5, should_sum=True)
    log_hets = [0 if h == 0 else 1 if h == 1 else np.log2(h) for h in hets]
    raw_fpfn = list(map(sum, zip(fp, fn)))
    fpfn = smooth_values(raw_fpfn, size=5, should_sum=True)
    log_fpfn = [0 if h == 0 else 1 if h == 1 else np.log2(h) for h in fpfn]
    avg_correct = np.mean(list(filter(lambda x: x is not None, correct_ratio)))
    std_correct = np.std(list(filter(lambda x: x is not None, correct_ratio)))
    avg_depth = np.mean(total)
    top_y_coord = avg_depth * 2

    fig, ((ax1, ax2, ax3)) = plt.subplots(nrows=3,ncols=1,sharex='all',gridspec_kw={'height_ratios': [2, 1, 1]})
    lw = .5

    ax1.set_ylabel('Phasing Partitions')
    ax1.fill_between(x, right, color='tab:red', linewidth=lw)
    ax1.fill_between(x, rong, color='tab:blue', linewidth=lw)
    ax1.set_ylim(-1 * top_y_coord, top_y_coord)
    if has_het_vcf:
        sizes = [(1 if ht == 0 else min(256, lw/16*(4**ht))) for ht in log_hets]
        ax1.scatter(x, [0 for _ in x], color='black', s=sizes, alpha=.05)
    if chunk_boundaries is not None:
        for cb in chunk_boundaries:
            ax1.axvline(x=int(cb/SPACING), linewidth=lw, alpha=.5, color="black", linestyle="dotted")

    ax2.set_ylabel('Correct Ratio')
    if has_het_vcf or has_result_vcf:
        ax2point5 = ax2.twinx()
        means = []
        if has_het_vcf:
            ax2point5.plot(x, log_hets, color='black', alpha=.2, linewidth=lw)
            means.append(np.mean(list(filter(lambda x: x!=0, log_hets))))
        if has_result_vcf:
            ax2point5.plot(x, log_fpfn, color='purple', alpha=.2, linewidth=lw)
            means.append(np.mean(list(filter(lambda x: x!=0, log_fpfn))))
        ax2point5.set_ylim(0, 3 * max(means))
    c = ax2.scatter(x, correct_ratio, c=correct_ratio, s=lw, cmap=plt.cm.get_cmap('RdYlGn'), alpha=.5)
    ax2.plot(x, [avg_correct for _ in x], color='black', alpha=.5, linewidth=lw)
    log("Correct ratio:\n\tAvg: {}\n\tStd: {}".format(avg_correct, std_correct))
    ax2.annotate("{:5.2f}".format(avg_correct), (x[0], avg_correct+1), fontfamily='monospace', fontsize=12,weight="bold")
    if chunk_boundaries is not None:
        for cb in chunk_boundaries:
            ax2.axvline(x=int(cb/SPACING), linewidth=lw, alpha=.5, color="black", linestyle="dotted")
    ax2.set_ylim(45, 105)

    if phasesets is not None:
        top = True
        for ps in phasesets:
            ax2.plot(range(ps[0] // SPACING, ps[1] // SPACING),
                     [48 if top else 46 for _ in range(ps[0] // SPACING, ps[1] // SPACING)],
                     color='black', alpha=.65, linewidth=2)
            ax1.plot(range(ps[0] // SPACING, ps[1] // SPACING),
                     [2 if top else -2 for _ in range(ps[0] // SPACING, ps[1] // SPACING)],
                     color='black', alpha=.65, linewidth=2)
            if ps[0] // SPACING != ps[1] // SPACING:
                top = not top


    ax3.set_ylabel('Classified Depth')
    ax3.plot(x, unclassified, color='lightgrey', linewidth=lw)
    ax3.plot(x, total, color='grey', linewidth=lw)

    avg_unclassified = np.mean(list(filter(lambda x: x != 0, unclassified)))
    ax3.plot(x, [avg_unclassified for _ in x], color='black', alpha=.5, linewidth=lw)
    log("Total Unknown:\n\tAvg: {}\n\tStd: {}".format(avg_unclassified, np.std(unclassified)))
    ax3.annotate("{:3d}x".format(int(avg_unclassified)), (x[0], avg_unclassified+1), fontfamily='monospace', fontsize=12,weight="bold")

    avg_classified = np.mean(list(filter(lambda x: x != 0, total)))
    ax3.plot(x, [avg_classified for _ in x], color='black', alpha=.5, linewidth=lw)
    log("Total Classified:\n\tAvg: {}\n\tStd: {}".format(avg_classified, np.std(total)))
    ax3.annotate("{:3d}x".format(int(avg_classified)), (x[0], avg_classified+1), fontfamily='monospace', fontsize=12,weight="bold")

    ax3.set_ylim(-.05 * top_y_coord, top_y_coord)

    if title is not None:
        ax1.set_title(title)
    fig.tight_layout()
    fig.set_size_inches(18, 12)
    if figName is not None:
        plt.savefig(figName)
    plt.show()
    plt.close()


def save_het_counts(vcf_file, position_classifications):
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"): continue
            pos = int(line.split("\t")[1])
            position_classifications[int(pos/SPACING)][HETS] += 1


def save_fp_fn_counts(vcf_file, position_classifications):
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"): continue
            pos = int(line.split("\t")[1])
            if "FP" in line:
                position_classifications[int(pos/SPACING)][FP] += 1
            # if "FN" in line:
            #     position_classifications[int(pos/SPACING)][FN] += 1


def read_phaseset_bed(bed_file):
    phasesets = list()
    with open(bed_file) as bed:
        for line in bed:
            parts = line.split("\t")
            assert(len(parts) >= 3)
            phasesets.append((int(parts[1]), int(parts[2])))
    return phasesets


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

    # get chunks
    chunk_boundaries = list()
    if args.chunks is not None:
        with open(args.chunks) as cin:
            chunk_idx = 0
            for line in cin:
                parts = line.split(",")
                if chunk_idx != 0:
                    chunk_boundaries.append(int(parts[3]))
                chunk_idx += 1


    # get read phasing pairs
    samfile = None
    read_count = 0
    missing_hp_count = 0
    position_classifications = collections.defaultdict(
        lambda : collections.defaultdict(lambda : 0)
    )
    analyzed_lengths = []
    try:
        samfile = pysam.AlignmentFile(args.margin_input_bam, 'rb' if args.margin_input_bam.endswith("bam") else 'r')
        for read in samfile.fetch():
            read_count += 1
            if not read.has_tag(HP_TAG):
                missing_hp_count += 1
                continue

            hp = read.get_tag(HP_TAG)
            id = read.query_name
            spos = read.reference_start
            epos = read.reference_end
            classifier = UNCLASSIFIED
            if args.untagged_only:
                if hp == 0:
                    classifier = CORRECT
            elif id in truth_all and hp != 0:
                if hp == 1 and id in truth_h1:
                    classifier = CORRECT
                elif hp == 2 and id in truth_h2:
                    classifier = CORRECT
                else:
                    classifier = INCORRECT

            if classifier != UNCLASSIFIED:
                analyzed_lengths.append(epos - spos)

            while spos <= epos:
                pos = int(spos / SPACING)
                position_classifications[pos][classifier] += 1
                spos += SPACING
    finally:
        if samfile is not None: samfile.close()


    if args.het_vcf is not None:
        save_het_counts(args.het_vcf, position_classifications)

    if args.result_vcf is not None:
        save_fp_fn_counts(args.result_vcf, position_classifications)

    log("Classified Read Lengths:")
    log("\tmean:   {}".format(np.mean(analyzed_lengths)))
    log("\tmedain: {}".format(np.median(analyzed_lengths)))
    analyzed_lengths.sort()
    len_total = sum(analyzed_lengths)
    len_curr = 0
    for l in analyzed_lengths:
        len_curr += l
        if len_curr > len_total/2:
            log("\tN50:    {}".format(l))
            break

    figName = args.figure_name
    if figName is None and args.figure_name_bam:
        figName = args.margin_input_bam + ".png"

    phasesets = None if args.phaseset_bed is None else read_phaseset_bed(args.phaseset_bed)
    plottit(position_classifications, chunk_boundaries=chunk_boundaries, phasesets=phasesets, figName=figName,
            has_het_vcf=args.het_vcf is not None, has_result_vcf=args.result_vcf is not None,
            title=args.margin_input_bam if args.title is None and args.figure_name_bam else args.title)





if __name__ == "__main__":
    main()











