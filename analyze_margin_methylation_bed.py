#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pysam
import collections
import numpy as np
import os
import seaborn as sns
import pandas
import math
from scipy.stats import pearsonr

plt.style.use('ggplot')
text_fontsize = 8
# plt.rcParams['ytick.labelsize']=text_fontsize+4
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

CHR_IDX = 0
START_IDX = 1
END_IDX = 2
STRAND_IDX = 3
METHYL_COUNT_IDX = 4
COVERAGE_IDX = 5
FWD = "+"
REV = "-"

CPG_START_IDX = 0
CPG_END_IDX = 1
CPG_METHYLS_IDX = 2
CPG_METHYL_RATIO_IDX = 3

class MethylRecord:
    def __init__(self, chr, pos, strand, coverage, methyl_count):
        self.chr = chr
        self.pos = int(pos)
        self.strand = strand
        self.coverage = int(coverage)
        self.methyl_count = int(methyl_count)
        self.methyl_ratio = self.methyl_count / self.coverage

    def __str__(self):
        return "MethylRecord({}:{} {} {}/{})".format(self.chr, self.pos, self.strand, self.methyl_count, self.coverage)


class MethylLocus:
    def __init__(self, chr, cpg, fwd_methyl: MethylRecord, rev_methyl: MethylRecord):
        self.fwd_methyl = fwd_methyl
        self.rev_methyl = rev_methyl
        self.chr = chr
        self.start_pos = cpg[CPG_START_IDX]
        self.end_pos = cpg[CPG_END_IDX]
        self.explicit_coverage = (0 if fwd_methyl is None else fwd_methyl.coverage) + (0 if rev_methyl is None else rev_methyl.coverage)
        self.estimated_coverage = 0 if fwd_methyl is None and rev_methyl is None else (self.explicit_coverage if None not in [fwd_methyl, rev_methyl] else 2 * self.explicit_coverage)
        self.total_methyl = (0 if fwd_methyl is None else fwd_methyl.methyl_count) + (0 if rev_methyl is None else rev_methyl.methyl_count)
        self.avg_methyl_ratio = self.total_methyl / max(1, self.estimated_coverage)
        self.cpg = cpg


def parse_args():
    parser = argparse.ArgumentParser("Produce stats and refactoring for methylation bed file based from bisulfite data")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                       help='methylation bed input file')
    parser.add_argument('--cpg_bed', '-c', dest='cpg_bed', required=True, type=str,
                       help='bed file describing CpG sites')

    parser.add_argument('--min_coverage', '-d', dest='min_coverage', required=False, default=8, type=int,
                       help='Minimum depth/coverage for each directional record to be used in analysis')
    parser.add_argument('--boolean_truth_methyl_threshold', '-T', dest='boolean_truth_methyl_threshold', required=False, default=.95, type=float,
                       help='Threshold used to quantify boolean "methylated" or "not" in truth')
    parser.add_argument('--boolean_query_methyl_threshold', '-t', dest='boolean_query_methyl_threshold', required=False, default=.8, type=float,
                       help='Threshold used to quantify boolean "methylated" or "not" in query')

    parser.add_argument('--output_base', '-o', dest='output_base', required=False, default=None, type=str,
                       help='Write output files, otherwise will use input bed name as prefix')
    parser.add_argument('--write_full_bed', '-1', dest='write_full_bed', default=False, required=False, action='store_true',
                       help='Write filtered and merged methylation loci')
    parser.add_argument('--write_quartile_beds', '-4', dest='write_quartile_beds', default=False, required=False, action='store_true',
                       help='Stratify calls into quartiles')
    parser.add_argument('--write_quintile_beds', '-5', dest='write_quintile_beds', default=False, required=False, action='store_true',
                       help='Stratify calls into quintiles')
    parser.add_argument('--plot', '-p', dest='plot', default=False, required=False, action='store_true',
                       help='Produce plots describing results')

    return parser.parse_args()


def log(msg):
    print(msg, file=sys.stderr)


def print_histogram(data, idx_str_function, size = 32):
    max_size = max(1, max(data))
    total = max(1, sum(data))
    for i, x in enumerate(data):
        count = int(size * x / max_size)
        log("  {:>10s}: {}{} {:>8d} {:.3f}".format(idx_str_function(i), "#" * count, " " * (size - count), x, x/total))


def print_methyl_loci_differences(loci_map):
    # prep
    size_factor = 10
    size = size_factor * 2 + 1
    coverage_value_differences = [0 for x in range(size)]
    coverage_ratio_differences = [0 for x in range(size)]
    methyl_value_differences = [0 for x in range(size)]
    methyl_ratio_differences = [0 for x in range(size)]

    for chr_loci in loci_map.values():
        for locus in chr_loci:
            fwd = locus.fwd_methyl
            rev = locus.rev_methyl
            if None in [fwd, rev]: continue

            # values + idx
            coverage_value = fwd.coverage - rev.coverage
            coverage_value_idx = min(size_factor * 2 , max(0, size_factor + coverage_value))
            coverage_ratio = coverage_value / ((fwd.coverage + rev.coverage)/2)
            coverage_ratio_idx = min(size_factor * 2 , max(0, size_factor + int(size_factor * coverage_ratio)))
            methyl_value = fwd.methyl_count - rev.methyl_count
            methyl_value_idx = min(size_factor * 2 , max(0, size_factor + methyl_value))
            methyl_ratio = methyl_value / max(0.0001, ((fwd.methyl_count + rev.methyl_count)/2))
            methyl_ratio_idx = min(size_factor * 2 , max(0, size_factor + int(size_factor * methyl_ratio)))

            coverage_value_differences[coverage_value_idx] += 1
            coverage_ratio_differences[coverage_ratio_idx] += 1
            methyl_value_differences[methyl_value_idx] += 1
            methyl_ratio_differences[methyl_ratio_idx] += 1

    # print it
    log("\nCoverage Value:")
    print_histogram(coverage_value_differences, lambda x: "<= -{}".format(size_factor) if x == 0 else \
        (">= +{}".format(size_factor) if x == size_factor * 2 else "{:+d}".format(x - size_factor)))
    log("\nCoverage Ratio:")
    print_histogram(coverage_ratio_differences, lambda x: "<= {:+d}%".format((0 - size_factor) * size_factor) if x == 0 else \
        (">= {:+d}%".format(size_factor * size_factor) if x == size_factor * 2 else "{:+d}%".format(int((x - size_factor) * size_factor))))
    log("\nMethyl Value:")
    print_histogram(methyl_value_differences, lambda x: "<= -{}".format(size_factor) if x == 0 else \
        (">= +{}".format(size_factor) if x == size_factor * 2 else "{:+d}".format(x - size_factor)))
    log("\nMethyl Ratio:")
    print_histogram(methyl_ratio_differences, lambda x: "<= {:+d}%".format((0 - size_factor) * size_factor) if x == 0 else \
        (">= {:+d}%".format(size_factor * size_factor) if x == size_factor * 2 else "{:+d}%".format(int((x - size_factor) * size_factor))))


def print_methyl_accuracy_differences(differences):
    # prep
    size_factor = 10
    size = size_factor * 2 + 1
    buckets = [0 for _ in range(size)]
    bucket_tostring = lambda x: "{:+d}%".format(int((x - size_factor) * size_factor))
    for diff in differences:
        idx = int(np.round(diff * size_factor)) + size_factor
        # print("{:+.5f} -> {:2d} -> {}".format(diff, idx, bucket_tostring(idx)))
        buckets[idx] += 1
    log("Methyl accuracy differences (T-Q):")
    print_histogram(buckets, bucket_tostring)


def plot_truth_query_heatmap(truth_query_pairs, output_base=None, size=10, fig_size=(8, 8)):

    # get data
    get_idx=lambda x: min(size-1, int(x*size))
    pair_maps = [[0 for y in range(size)] for x in range(size)]
    for tqp in truth_query_pairs:
        pair_maps[size - 1 - get_idx(tqp[1])][get_idx(tqp[0])] += 1
    log_pair_maps = [[0 if x == 0 else math.log10(x) for x in y] for y in pair_maps]

    fig, ax = plt.subplots(figsize=fig_size)
    im = ax.imshow(log_pair_maps)
    cbar = plt.colorbar(im,fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels([int(10**x) for x in cbar.ax.get_yticks()])

    # label ticks
    ax.set_xticks([x for x in range(size)])
    ax.set_yticks([x for x in range(size)])
    ax.set_xticklabels(["{}-{}%".format(int(100.0 * x / size), int(100.0 * (x+1) / size)) for x in range(size)])
    ax.set_yticklabels(reversed(["{}-{}%".format(int(100.0 * x / size), int(100.0 * (x+1) / size)) for x in range(size)]))

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # add labels
    for i in range(size):
        for j in range(size):
            text = ax.text(j, i, pair_maps[i][j],
                           ha="center", va="center", color="w")

    plt.suptitle("BS vs Guppy Methyl Ratio", y=1)
    ax.set_ylabel("Guppy Methyl Ratio")
    ax.set_xlabel("Bisulfite Methyl Ratio")
    fig.tight_layout()
    if output_base is not None:
        filename="{}.bs_vs_guppy_heatmap.png".format(output_base)
        log("\tSaving plot to {}".format(filename))
        plt.savefig(filename)
    plt.show()
    plt.close()


def plot_jointplot_hist(truth_query_pairs, output_base=None, bins=5, fig_size=8):
    # df_prep = []
    # for tqp in truth_query_pairs:
    #     df_prep.append([tqp[0], tqp[1]])

    columns = ["Bisulfite Methyl Ratio", "Guppy Methyl Ratio"]
    df = pandas.DataFrame(truth_query_pairs, columns=columns)

    sns.jointplot(data=df, x=columns[0], y=columns[1], bins=bins, height=fig_size, kind="hist", marginal_kws=dict(bins=bins))

    plt.suptitle("BS vs Guppy Methyl Ratio", y=1)
    plt.tight_layout()
    if output_base is not None:
        filename="{}.bs_vs_guppy_jointplot_hist_{}.png".format(output_base, bins)
        log("\tSaving plot to {}".format(filename))
        plt.savefig(filename)
    plt.show()
    plt.close()


def plot_jointplot_reg(truth_query_pairs, output_base=None, fig_size=8, bins=20):
    # df_prep = []
    # for tqp in truth_query_pairs:
    #     df_prep.append([tqp[0], tqp[1]])

    columns = ["Bisulfite Methyl Ratio", "Guppy Methyl Ratio"]
    num_of_datapoint = 25000
    if len(truth_query_pairs) > num_of_datapoint:
        ratio = num_of_datapoint / len(truth_query_pairs)
        truth_query_pairs = list(filter(lambda x: np.random.random() < ratio, truth_query_pairs))
    df = pandas.DataFrame(truth_query_pairs, columns=columns)

    avg_points_per_bin = len(truth_query_pairs) / (bins * bins)
    # alpha = min(1.0, 10.0 / avg_points_per_bin)
    # log("\tappb: {}, alpha: {}".format(avg_points_per_bin, alpha))
    p = sns.jointplot(data=df, x=columns[0], y=columns[1], marker='o', height=fig_size, kind="reg",
                  marginal_kws={'bins': 20}, joint_kws = { 'line_kws':{'color':'red'}, 'scatter_kws':{'alpha': .05, 'linewidth':0}})
    p.ax_joint.text(0.05, 0.95, "Pearson's r: {:.3f}".format(pearsonr(df[columns[0]], df[columns[1]])[0]), horizontalalignment='left', size='medium', color='black', weight='semibold')

    # plt.suptitle("BS vs Guppy Methyl Ratio", y=1)
    plt.suptitle("Bisulfite vs Guppy 4.5.4 Methylation Ratio (chr20)", y=.975)
    # plt.tight_layout()
    if output_base is not None:
        filename="{}.bs_vs_guppy_jointplot_reg.png".format(output_base)
        log("\tSaving plot to {}".format(filename))
        plt.savefig(filename.replace(".png", ".pdf"), format='pdf', dpi=300)
        # plt.savefig(filename)
    plt.show()
    plt.close()



def main():
    args = parse_args()

    # data we want
    all_methyl_records = collections.defaultdict(lambda: list())
    all_methyl_loci = collections.defaultdict(lambda: list())
    all_cpg_sites = collections.defaultdict(lambda: list())

    # read input file
    log("Reading margin records from {}".format(args.input))
    with open(args.input) as fin:
        for linenr, line in enumerate(fin):
            if line.startswith("#") or len(line.strip()) == 0: continue
            parts = line.strip().split("\t")
            assert len(parts) >= COVERAGE_IDX, "Malformed line {}: '{}'".format(linenr, line)
            chr, pos, strand, coverage, methyl_count = parts[CHR_IDX], parts[START_IDX], parts[STRAND_IDX], \
                                                       parts[COVERAGE_IDX], parts[METHYL_COUNT_IDX]
            record = MethylRecord(chr, pos, strand, coverage, methyl_count)
            all_methyl_records[chr].append(record)
    for methyl_list in all_methyl_records.values():
        methyl_list.sort(key=lambda x: x.pos)
    total_methyl_records = sum(list(map(len, all_methyl_records.values())))
    log("Got {} total methyl records".format(total_methyl_records))

    # get cpg sites
    log("Reading CpG sites from {}".format(args.cpg_bed))
    with open(args.cpg_bed) as fin:
        for linenr, line in enumerate(fin):
            if line.startswith("#") or len(line.strip()) == 0: continue
            parts = line.split("\t")
            assert len(parts) >= 3, "Malformed line {}: '{}'".format(linenr, line)
            chr, start, end = parts[0], int(parts[1]), int(parts[2])
            methyl_ratio = None if len(parts) < 5 else float(parts[4])
            assert (end - start) % 2 == 0, "CpG site with odd size at line {}: '{}'".format(linenr, line)
            while start < end:
                all_cpg_sites[chr].append((start, start+2, [], methyl_ratio))
                start += 2
    for cpg_list in all_cpg_sites.values():
        cpg_list.sort(key=lambda x: x[START_IDX])
    log("Got {} CpG sites".format(sum(list(map(len, all_cpg_sites.values())))))

    # attach methyls to cpg sites
    skipped_for_not_cpg = 0
    for chr in all_cpg_sites.keys():
        cpg_iter = iter(all_cpg_sites[chr])
        methyl_iter = iter(all_methyl_records[chr])
        curr_cpg = next(cpg_iter, None)
        curr_methyl = next(methyl_iter, None)
        step = 0
        while curr_cpg is not None and curr_methyl is not None:
            if curr_cpg[CPG_START_IDX] <= curr_methyl.pos < curr_cpg[CPG_END_IDX]:
                curr_cpg[CPG_METHYLS_IDX].append(curr_methyl)
                curr_methyl = next(methyl_iter, None)
            elif curr_methyl.pos < curr_cpg[CPG_START_IDX]:
                skipped_for_not_cpg += 1
                curr_methyl = next(methyl_iter, None)
            elif curr_cpg[CPG_END_IDX] <= curr_methyl.pos:
                curr_cpg = next(cpg_iter, None)
            step += 1
    log("Skipped {} ({:.5f}) methyl records for not falling within CpG sites".format(skipped_for_not_cpg, 1.0 * skipped_for_not_cpg / total_methyl_records))

    # accuracy stats between bed and margin methyls (if available)
    threshold_tp = 0
    threshold_fp = 0
    threshold_tn = 0
    threshold_fn = 0
    total_percentages_off = []
    truth_query_pairs = []

    # get methyl cpg sites
    skipped_for_wrong_cpg = 0
    skipped_for_missing_depth = 0
    for chr in all_cpg_sites.keys():
        for cpg in all_cpg_sites[chr]:
            methyls = cpg[CPG_METHYLS_IDX]
            fwd_records = list(filter(lambda x: x.strand == FWD and x.pos == cpg[CPG_START_IDX], methyls))
            rev_records = list(filter(lambda x: x.strand == REV and x.pos == cpg[CPG_START_IDX] + 1, methyls))
            assert len(fwd_records) <= 1 and len(rev_records) <= 1, "Got multiple FWD or REV records for CpG {}".format(cpg)
            fwd = None if len(fwd_records) == 0 else fwd_records[0]
            rev = None if len(rev_records) == 0 else rev_records[0]
            skipped_for_wrong_cpg += len(methyls) - (1 if fwd is not None else 0) - (1 if rev is not None else 0)

            # build locus
            methyl_locus = MethylLocus(chr, cpg, fwd, rev)

            # should save locus?
            if None in [fwd, rev] and methyl_locus.estimated_coverage < 2 * args.min_coverage:
                skipped_for_missing_depth += 1
                continue
            if None not in [fwd, rev] and (fwd.coverage < args.min_coverage or rev.coverage < args.min_coverage):
                skipped_for_missing_depth += 1
                continue
            all_methyl_loci[chr].append(methyl_locus)

            # quantify accuracy
            truth_methyl_ratio = round(cpg[CPG_METHYL_RATIO_IDX], 5)
            query_methyl_ratio = round(methyl_locus.avg_methyl_ratio, 5)
            if truth_methyl_ratio is None:
                continue
            if truth_methyl_ratio >= args.boolean_truth_methyl_threshold:
                if query_methyl_ratio >= args.boolean_query_methyl_threshold:
                    threshold_tp += 1
                else:
                    threshold_fn += 1
            else:
                if query_methyl_ratio >= args.boolean_query_methyl_threshold:
                    threshold_fp += 1
                else:
                    threshold_tn += 1
            total_percentages_off.append(truth_methyl_ratio - query_methyl_ratio)
            truth_query_pairs.append((truth_methyl_ratio, query_methyl_ratio))

    # loggit
    log("Got {} methyl loci, skipped {} ({:.5f}) methyl sites because of wrong CpG setup".format(
        sum(list(map(len, all_methyl_loci.values()))), skipped_for_wrong_cpg, 1.0 * skipped_for_wrong_cpg / total_methyl_records))
    print_methyl_loci_differences(all_methyl_loci)

    # log accuracy
    log("Boolean accuracies based on threshold of T{}/Q{}:".format(args.boolean_truth_methyl_threshold, args.boolean_query_methyl_threshold))
    log("\tTP: {}".format(threshold_tp))
    log("\tFP: {}".format(threshold_fp))
    log("\tFN: {}".format(threshold_fn))
    log("\tTN: {}".format(threshold_tn))
    log("\tPrecision:   {}".format(threshold_tp / max(1, threshold_tp + threshold_fp)))
    log("\tRecall:      {}".format(threshold_tp / max(1, threshold_tp + threshold_fn)))
    log("\tSpecificity: {}".format(threshold_tn / max(1, threshold_tn + threshold_fp)))
    if len(total_percentages_off) > 0:
        log("Accuracies based off distance (T-Q):")
        log("\tMean distance: {}".format(np.mean(total_percentages_off)))
        log("\tMean positive distance: {}".format(np.mean(list(filter(lambda x: x > 0, total_percentages_off)))))
        log("\tMean negative distance: {}".format(np.mean(list(filter(lambda x: x < 0, total_percentages_off)))))
    print_methyl_accuracy_differences(total_percentages_off)

    # for outputs
    output_base = args.output_base
    if output_base is None:
        output_base = ".".join(os.path.basename(args.input).split(".")[:-1])

    # plot
    if args.plot:
        log("Plotting data")
        plot_jointplot_reg(truth_query_pairs, output_base=output_base)
        sys.exit()
        plot_jointplot_hist(truth_query_pairs, output_base=output_base)
        plot_jointplot_hist(truth_query_pairs, output_base=output_base, bins=10)
        plot_jointplot_hist(truth_query_pairs, output_base=output_base, bins=20)
        plot_truth_query_heatmap(truth_query_pairs, output_base=output_base)

    # write
    should_write = True in [args.write_quartile_beds, args.write_quintile_beds, args.write_full_bed]
    if should_write:
        log("\nWriting output to files with base {}:".format(output_base))

        full_file_name = "{}.full.bed".format(output_base)
        quartile_file_names = list(map(lambda x: "{}.{}.bed".format(output_base, x), ["0-25", "25-50", "50-75", "75-100"]))
        quintile_file_names = list(map(lambda x: "{}.{}.bed".format(output_base, x), ["0-20", "20-40", "40-60", "60-80", "80-100"]))
        full_out = None
        quartile_outs = None
        quintile_outs = None

        try:
            header = "#chr\tstart_pos\tend_pos\testimated_coverage\tquery_methyl_ratio\ttruth_methyl_ratio\n"
            # open files
            output_files = []
            if args.write_full_bed:
                full_out = open(full_file_name, 'w')
                full_out.write(header)
                output_files.append(full_file_name)
            if args.write_quartile_beds:
                quartile_outs = list(map(lambda x: open(x, 'w'), quartile_file_names))
                for o in quartile_outs:
                    o.write(header)
                output_files.extend(quartile_file_names)
            if args.write_quintile_beds:
                quintile_outs = list(map(lambda x: open(x, 'w'), quintile_file_names))
                for o in quintile_outs:
                    o.write(header)
                output_files.extend(quintile_file_names)
            # finish logging
            for fn in output_files:
                log("\t{}".format(fn))

            # look at all loci and write out
            for chr in sorted(list(all_methyl_loci.keys())):
                for locus in all_methyl_loci[chr]:
                    bed_data = [locus.chr, locus.start_pos, locus.end_pos, locus.estimated_coverage, round(locus.avg_methyl_ratio, 5), locus.cpg[CPG_METHYL_RATIO_IDX]]
                    bed_line = "\t".join(map(str, bed_data)) + "\n"

                    r = locus.avg_methyl_ratio
                    if args.write_full_bed:
                        full_out.write(bed_line)
                    if args.write_quartile_beds:
                        out_idx = 0 if r < .25 else 1 if r < .5 else 2 if r < .75 else 3
                        quartile_outs[out_idx].write(bed_line)
                    if args.write_quintile_beds:
                        out_idx = 0 if r < .2 else 1 if r < .4 else 2 if r < .6 else 3 if r < .8 else 4
                        quintile_outs[out_idx].write(bed_line)

        finally:
            # close file
            if full_out is not None: full_out.close()
            if quartile_outs is not None: [x.close() for x in quartile_outs]
            if quintile_outs is not None: [x.close() for x in quintile_outs]


if __name__ == "__main__":
    main()
