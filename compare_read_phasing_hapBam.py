#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import matplotlib.pyplot as plt
import pysam
import collections
import numpy as np
import os

plt.style.use('ggplot')
text_fontsize = 8
# plt.rcParams['ytick.labelsize']=text_fontsize+4
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')


HP_TAG = "HP"
UNCLASSIFIED = 'u'
CORRECT = 'c'
INCORRECT = 'i'
UNKNOWN = 'k'
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
    parser.add_argument('--high_conf_bed', '-C', dest='high_conf_bed', default=None, required=False, type=str,
                       help='If set, will only count regions within this BED file')
    parser.add_argument('--title', '-t', dest='title', default=None, required=False, type=str,
                       help='Figure title')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                       help='Figure name (will save if set)')
    parser.add_argument('--figure_name_bam', '-F', dest='figure_name_bam', default=False, required=False, action='store_true',
                       help='Figure name should be based off of BAM file (-f overrides this)')
    parser.add_argument('--max_depth', '-d', dest='max_depth', default=72, required=False, type=int,
                       help='What depth should be used on y axes')
    parser.add_argument('--only_natural_switch', '-N', dest='only_natural_switch', default=False, required=False, action='store_true',
                       help='Only plot the natural switch')


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


def plotOnlyNaturalSwitch(classification_data, args, phasesets=None, figName=None):

    start_idx = min(classification_data.keys())
    end_idx = max(classification_data.keys())
    x = []
    right = []
    rong = []
    for i in range(start_idx, end_idx + 1):
        x.append(i / (1000000.0 / SPACING))
        ri = classification_data[i][CORRECT]
        ro = classification_data[i][INCORRECT]

        right.append(min(args.max_depth, ri))
        rong.append(-1 * min(args.max_depth, ro))

    # get plots
    # fig, (ax1) = plt.subplots(nrows=1,ncols=1)
    fig, (ax1) = plt.subplots(nrows=1,ncols=1, figsize=(8, 1.75))
    lw = .5

    ax1.set_ylabel('Concordant/Discordant\nDepth')
    ax1.set_xlabel('Chr1 Positions (Mb)')
    ax1.fill_between(x, right, color='tab:red', linewidth=lw)
    ax1.fill_between(x, rong, color='tab:blue', linewidth=lw)
    ax1.set_ylim(-1 * args.max_depth, args.max_depth)

    if phasesets is not None:
        top = True
        for ps in phasesets:
            modifier = 1000000.0
            start = ps[0] / modifier
            end = ps[1] / modifier
            # ax1.plot(range(start, end), [2 if top else -2 for _ in range(start, end)],
            #          color='black', alpha=.65, linewidth=2)
            # start = ps[0] // SPACING
            # end = ps[1] // SPACING
            ps_range = list(map(lambda x: x/modifier, range(int(start*modifier), int(end*modifier), 10000)))
            ax1.plot(ps_range, [2 if top else -2 for _ in ps_range], color='black', alpha=1, linewidth=.75)
            top = not top

    fig.tight_layout()
    # fig.set_size_inches(12, 3)
    if figName is not None:
        plt.savefig(figName+".svg", format='svg', dpi=50)
        plt.savefig(figName+".pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()



def plottit(classification_data, args, figName=None, phasesets=None, has_het_vcf=False, has_result_vcf=False, chunk_boundaries=None,
            title=None, highconf_positions=None):
    start_idx = min(classification_data.keys())
    end_idx = max(classification_data.keys())
    x = []
    right = []
    rong = []
    total_classified = []
    total_reads = []
    unknown = []
    unclassified = []
    correct_ratio = []
    raw_hets = []
    fp = []
    fn = []
    for i in range(start_idx, end_idx + 1):
        x.append(i)
        ri = classification_data[i][CORRECT]
        ro = classification_data[i][INCORRECT]
        unc = classification_data[i][UNCLASSIFIED]
        unk = classification_data[i][UNKNOWN]
        ht = classification_data[i][HETS]
        if highconf_positions is not None and i not in highconf_positions:
            ri = 0
            ro = 0
            unc = 0
            unk = 0
        total = ri+ro+unc+unk


        right.append(ri)
        rong.append(-1 * ro)
        unknown.append(unk)
        unclassified.append(unk)
        total_classified.append(ri + ro + unc)
        total_reads.append(total)
        correct_ratio.append(None if ri + ro == 0 else 100*abs(max(ri,ro)/(ri + ro)))
        raw_hets.append(ht)
        fp.append(classification_data[i][FP])
        fn.append(classification_data[i][FN])

    # get averages
    avg_unknown = np.mean(list(map(lambda y: y[1], filter(lambda x: x[0] != 0, zip(total_reads, unknown)))))
    avg_classified = np.mean(list(map(lambda y: y[1], filter(lambda x: x[0] != 0, zip(total_reads, total_classified)))))
    avg_total_reads = np.mean(list(filter(lambda x: x != 0, total_reads)))
    avg_correct = np.mean(list(filter(lambda x: x is not None, correct_ratio)))
    std_correct = np.std(list(filter(lambda x: x is not None, correct_ratio)))

    # smooth values
    smoothed_unknown = smooth_values(unknown, size=20)
    smoothed_total_classified = smooth_values(total_classified, size=20)
    smoothed_total_reads = smooth_values(total_reads, size=20)
    smoothed_hets = smooth_values(raw_hets, size=5, should_sum=True)
    smoothed_fpfn = smooth_values(list(map(sum, zip(fp, fn))), size=5, should_sum=True)

    # log scale stuff (for variants)
    log_hets = [0 if h == 0 else 1 if h == 1 else np.log2(h) for h in smoothed_hets]
    log_fpfn = [0 if h == 0 else 1 if h == 1 else np.log2(h) for h in smoothed_fpfn]

    # get plots
    fig, ((ax1, ax2, ax3)) = plt.subplots(nrows=3,ncols=1,sharex='all',gridspec_kw={'height_ratios': [2, 1, 1]})
    lw = .5

    ax1.set_ylabel('Phasing Partitions')
    ax1.fill_between(x, right, color='tab:red', linewidth=lw)
    ax1.fill_between(x, rong, color='tab:blue', linewidth=lw)
    ax1.set_ylim(-1 * args.max_depth, args.max_depth)
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
        concordancy_in_phasesets = list()
        top = True
        for ps in phasesets:
            start = ps[0] // SPACING
            end = ps[1] // SPACING
            if not any([highconf_positions is None or x in highconf_positions for x in range(start, end) ]):
                continue
            ax2.plot(range(start, end), [48 if top else 46 for _ in range(start, end)],
                     color='black', alpha=.65, linewidth=2)
            ax1.plot(range(start, end), [2 if top else -2 for _ in range(start, end)],
                     color='black', alpha=.65, linewidth=2)
            if start != end:
                top = not top
                concordancy_in_phasesets.extend(list(filter(lambda x: x is not None, correct_ratio[start:end])))
        avg_ps_correct = np.mean(concordancy_in_phasesets)
        log("Correct ratio in phasesets:\n\tAvg: {}\n\tStd: {}".format(avg_ps_correct,
                                                                       np.std(concordancy_in_phasesets)))
        # ax2.plot(x, [avg_ps_correct for _ in x], color='black', alpha=.5, linewidth=lw)
        # ax2.annotate("{:5.2f} (Phaseset)".format(avg_ps_correct), (x[0], avg_ps_correct+1), fontfamily='monospace', fontsize=12,weight="bold")

    ax3.set_ylabel('Classified Depth')
    ax3.plot(x, smoothed_unknown, color='red', alpha=.25, linewidth=lw)
    ax3.plot(x, smoothed_total_classified, color='blue', alpha=.25, linewidth=lw)
    ax3.plot(x, smoothed_total_reads, color='black', alpha=.25, linewidth=lw)

    ax3.plot(x, [avg_unknown for _ in x], color='red', alpha=.5, linewidth=lw)
    log("Total Unknown:\n\tAvg: {}\n\tStd: {}".format(avg_unknown, np.std(unknown)))
    ax3.annotate("{:3d}x Unknown".format(int(avg_unknown)), (x[0], avg_unknown+1), fontfamily='monospace', color='darkred', fontsize=12, weight="bold")

    ax3.plot(x, [avg_classified for _ in x], color='blue', alpha=.5, linewidth=lw)
    log("Total Classified:\n\tAvg: {}\n\tStd: {}".format(avg_classified, np.std(total_classified)))
    ax3.annotate("{:3d}x Classified".format(int(avg_classified)), (x[0], avg_classified+1), fontfamily='monospace', color='darkblue', fontsize=12, weight="bold")

    ax3.plot(x, [avg_total_reads for _ in x], color='black', alpha=.5, linewidth=lw)
    log("Total Reads:\n\tAvg: {}\n\tStd: {}".format(avg_total_reads, np.std(total_reads)))
    ax3.annotate("{:3d}x Total Reads".format(int(avg_total_reads)), (x[0], max(avg_total_reads, avg_classified+8)+1), fontfamily='monospace', color='black', fontsize=12,weight="bold")

    ax3.set_ylim(-.05 * args.max_depth, args.max_depth)

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


def get_highconf_positions(bed_file):
    highconf_positions = set()
    with open(bed_file) as bed:
        for line in bed:
            parts = line.split("\t")
            assert(len(parts) >= 3)
            for i in range(int(int(parts[1]) / SPACING), int(int(parts[2]) / SPACING) + 1):
                highconf_positions.add(i)
    return highconf_positions


def get_position_classifications(bam_location, truth_h1_ids, truth_h2_ids, region=None, verbose=True):
    # get read phasing pairs
    samfile = None
    read_count = 0
    missing_hp_count = 0
    position_classifications = collections.defaultdict(
        lambda : collections.defaultdict(lambda : 0)
    )
    analyzed_lengths = []
    try:
        samfile = pysam.AlignmentFile(bam_location, 'rb' if bam_location.endswith("bam") else 'r')
        for read in samfile.fetch(region=region):
            read_count += 1
            id = read.query_name
            spos = read.reference_start
            epos = read.reference_end

            if not read.has_tag(HP_TAG):
                missing_hp_count += 1
                classifier = UNCLASSIFIED
            else:
                hp = read.get_tag(HP_TAG)
                if hp == 0:
                    classifier = UNCLASSIFIED
                elif id not in truth_h1_ids and id not in truth_h2_ids:
                    classifier = UNKNOWN
                elif hp == 1 and id in truth_h1_ids:
                    classifier = CORRECT
                elif hp == 2 and id in truth_h2_ids:
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

    if verbose:
        log("Classified Read Lengths{}:".format("" if region is None else " for {}".format(region)))
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

    return position_classifications


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

    # classifiy positions for reads
    position_classifications = get_position_classifications(args.margin_input_bam, truth_h1, truth_h2)

    highconf_positions = None
    if args.high_conf_bed is not None:
        highconf_positions = get_highconf_positions(args.high_conf_bed)

    if args.het_vcf is not None:
        save_het_counts(args.het_vcf, position_classifications)

    if args.result_vcf is not None:
        save_fp_fn_counts(args.result_vcf, position_classifications)

    figName = args.figure_name
    if figName is None and args.figure_name_bam:
        figName = os.path.basename(args.margin_input_bam) + (".natural_switch" if args.only_natural_switch else ("" if args.high_conf_bed is None else ".highconf")) + ".png"

    phasesets = None if args.phaseset_bed is None else read_phaseset_bed(args.phaseset_bed)
    if (args.only_natural_switch):
        plotOnlyNaturalSwitch(position_classifications, args, phasesets, figName)
    else:
        plottit(position_classifications, args=args, chunk_boundaries=chunk_boundaries, phasesets=phasesets, figName=figName,
                has_het_vcf=args.het_vcf is not None, has_result_vcf=args.result_vcf is not None,
                title=args.margin_input_bam if args.title is None and args.figure_name_bam else args.title,
                highconf_positions=highconf_positions)





if __name__ == "__main__":
    main()











