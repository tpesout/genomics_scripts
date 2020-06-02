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
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')

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


def plottit(classification_data, figName=None, has_het_vcf=False, has_result_vcf=False, chunk_boundaries=None):
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
    # divider = make_axes_locatable(ax2)
    # cax = divider.append_axes("right", size="3%", pad=0.5)
    # fig.colorbar(c, cax=cax, orientation='vertical')
    ax2.plot(x, [avg_correct for _ in x], color='black', alpha=.5, linewidth=lw)
    if chunk_boundaries is not None:
        for cb in chunk_boundaries:
            ax2.axvline(x=int(cb/SPACING), linewidth=lw, alpha=.5, color="black", linestyle="dotted")
    ax2.set_ylim(45, 105)

    ax3.set_ylabel('Depth, Unclassified')
    ax3.plot(x, total, color='black', linewidth=lw)
    ax3.plot(x, unclassified, color='grey', linewidth=lw)
    ax3.set_ylim(-.05 * top_y_coord, top_y_coord)

    fig.tight_layout()
    fig.set_size_inches(18, 12)
    if figName is not None:
        plt.savefig(figName)
    plt.show()
    plt.close()

    # second one
    # df = pd.DataFrame.from_dict({'het_sites':raw_hets,'nearby_hets':hets,'log_nearby_hets':log_hets,'correctness':correct_ratio})
    # df = df[df['nearby_hets'] > 0]
    # sns.jointplot(x='nearby_hets',y='correctness',data=df,stat_func=pearsonr,alpha=.1)
    # plt.show()
    # plt.close()
    #
    # sns.jointplot(x='log_nearby_hets',y='correctness',data=df,stat_func=pearsonr,alpha=.1)
    # plt.show()

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
            if id in truth_all:
                if hp == 1 and id in truth_h1:
                    classifier = CORRECT
                elif hp == 2 and id in truth_h2:
                    classifier = CORRECT
                else:
                    classifier = INCORRECT

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


    plottit(position_classifications, chunk_boundaries=chunk_boundaries, figName=args.figure_name,
            has_het_vcf=args.het_vcf is not None, has_result_vcf=args.result_vcf is not None)





if __name__ == "__main__":
    main()











