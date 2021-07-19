#!/usr/bin/env python
from __future__ import print_function
import argparse
import gzip
import sys
import os
import numpy as np
import random as rand
import math



def parse_args():
    parser = argparse.ArgumentParser("Downsample FASTQ based on length and total coverage")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                       help='FASTQ input file')
    parser.add_argument('--output', '-o', dest='output', required=False, default=None, type=str,
                       help='Write output to file (otherwise stdout)')
    parser.add_argument('--stats_only', '-s', dest='stats_only', action='store_true', default=False,
                       help='Only print stats on input FQ')

    # filtering
    parser.add_argument('--min_read_length', '-l', dest='min_read_length', default=None, type=int,
                       help='Only consider reads with this length or greater')
    parser.add_argument('--max_read_length', '-L', dest='max_read_length', default=None, type=int,
                       help='Only consider reads with this length or lesser')
    parser.add_argument('--min_quality', '-q', dest='min_quality', default=None, type=int,
                       help='Only include reads with this average quality or greater')

    # whether to keep reads
    parser.add_argument('--read_ratio', '-a', dest='read_ratio', default=None, type=float,
                       help='Keep approximately this ratio of the reads (must be between 0 and 1)')
    parser.add_argument('--total_bases', '-t', dest='total_bases', default=None, type=str,
                       help='Keep reads totaling around this amount of total nucleotides')
    parser.add_argument('--coverage_depth', '-c', dest='coverage_depth', default=None, type=float,
                       help='Keep reads such that roughly this depth is achieved (--reference_length parameter is required)')
    parser.add_argument('--reference_length', '-r', dest='reference_length', default=None, type=str,
                       help='Length of reference')

    return parser.parse_args()


def log(msg):
    print(msg, file=sys.stderr)


def human2comp(size_str):
    if size_str is None: return None
    size_str = size_str.upper()
    if not size_str.endswith("B"): return int(size_str)
    if size_str.endswith("KB"): return 1000 * int(size_str.rstrip("KB"))
    if size_str.endswith("MB"): return 1000 * 1000 * int(size_str.rstrip("MB"))
    if size_str.endswith("GB"): return 1000 * 1000 * 1000 * int(size_str.rstrip("GB"))
    if size_str.endswith("TB"): return 1000 * 1000 * 1000 * 1000 * int(size_str.rstrip("TB"))
    raise Exception("Unsupported size suffix: {}".format(size_str))


def bin_and_print_data(data, non_log_bin_size=None, log_base=None, indent_count=3):
    if (non_log_bin_size is None and log_base is None) or (non_log_bin_size is not None and log_base is not None):
        log("{}ERROR: bin_and_print_data invoked incorrectly".format("  " * indent_count))
        return
    # fill bins
    bins = dict()
    get_bin = lambda x: int(math.log(x, log_base) if log_base is not None else (x / non_log_bin_size))
    for d in data:
        b = get_bin(d)
        if b not in bins: bins[b] = 0
        bins[b] += 1
    # print buckets
    min_bin = min(bins.keys())
    max_bin = max(bins.keys())
    max_value = max(bins.values())
    last_print_skipped = False
    for b in range(min_bin, max_bin + 1):
        id = ("%d^%3d:\t" % (log_base, b)) if log_base is not None else ("%d* %4d:\t" % (non_log_bin_size, b))
        count = int(bins[b]) if b in bins else 0
        if count * 500 < max_value and log_base is None:
            if not last_print_skipped:
                log("{} ...".format('  '*indent_count))
            last_print_skipped = True
            continue
        pound_count = int(32.0 * count / max_value)
        log("{} {} {}{} {:6d}".format('  '*(indent_count), id, "#"*pound_count, " "*(32 - pound_count), count))
        last_print_skipped = False


def main():
    args = parse_args()

    # determine depth mechanism
    if not args.stats_only:
        if args.output is None:
            raise Exception("--output parameter is required if --stats_only is not set")
        inclusion_mechanisms = len(list(filter(lambda x: x is not None, [args.read_ratio, args.total_bases, args.coverage_depth])))
        if inclusion_mechanisms == 0:
            raise Exception("Must specify mechanism to keep reads: read_ratio, total_bases, or coverage_depth")
        elif inclusion_mechanisms > 1:
            raise Exception("Must specify only one mechanism to keep reads: read_ratio, total_bases, or coverage_depth")
        elif args.coverage_depth is not None and args.reference_length is None:
            raise Exception("If coverage_depth mechanism is used, reference_length must also be specified")
        elif args.read_ratio is not None and not (0.0 < args.read_ratio < 1.0):
            raise Exception("Parameter read_ratio must be between 0.0 and 1.0 (exclusive)")


    # get input parameters
    input = args.input
    if not os.path.isfile(input):
        raise Exception("Input file {} does not exist".format(input))

    # filters
    filters = list()
    if args.min_read_length is not None: filters.append("min length {}".format(args.min_read_length))
    if args.max_read_length is not None: filters.append("max length {}".format(args.max_read_length))
    def filter_read(read_data):
        assert type(read_data) == list
        read = read_data[1]
        read_length = len(read)
        if args.min_read_length is not None and read_length < args.min_read_length: return False
        if args.max_read_length is not None and read_length > args.max_read_length: return False
        return True

    # preprocessing
    log("Preprocessing {}:".format(input))
    if len(filters) != 0: log("  Filters: {}".format(", ".join(filters)))
    all_read_lengths = list()
    filtered_read_lengths = list()
    def handle_read_analysis(read, linenr):
        # format sanity check
        if not (read[0].startswith("@") and read[2].startswith("+") and
                len(read[1]) == len(read[3])):
            log("Error at approx line {} in {}: does not appear to be in FQ format.".format(linenr, input))
            sys.exit(1)
        all_read_lengths.append(len(read[1]))
        if filter_read(read):
            filtered_read_lengths.append(len(read[1]))

    # read file
    current_read = list()
    linenr = 0
    with open(input) as instream:
        for line in instream:
            linenr += 1
            current_read.append(line.strip())
            # have we finished the current read?
            if len(current_read) == 4:
                handle_read_analysis(current_read, linenr)
                current_read = list()

    # stats
    reference_size = human2comp(args.reference_length)
    def print_stats(read_lengths):
        log("    Total reads: {}".format(len(read_lengths)))
        log("    Total bases: {} {}".format(sum(read_lengths),
                                            "" if reference_size is None else
                                            ("({:.2f}x coverage)".format(1.0 * sum(read_lengths) / reference_size))))
        log("    Avg length:  {}".format(int(np.mean(read_lengths))))
        log("    Med length:  {}".format(int(np.median(read_lengths))))
        log("    Min length:  {}".format(min(read_lengths)))
        log("    Max length:  {}".format(max(read_lengths)))
        half_total_read_length = sum(read_lengths) / 2
        read_lengths.sort()
        for l in read_lengths:
            half_total_read_length -= l
            if half_total_read_length <= 0:
                log("    Read N50:    {}".format(l))
                break
        if (args.stats_only):
            log("    Lengths in log space:")
            bin_and_print_data(read_lengths, log_base=2)
            log("    Lengths in real space:")
            bin_and_print_data(read_lengths, non_log_bin_size=1000)


    log("  Read {} lines".format(linenr))
    if len(all_read_lengths) == 0:
        print("No reads found!")
        sys.exit(0)
    if reference_size is not None: log("  Assuming reference size of {}".format(reference_size))
    log("  All reads:")
    print_stats(all_read_lengths)
    if len(all_read_lengths) == len(filtered_read_lengths):
        log("  No reads filtered")
    else:
        log("  Reads passing filter:")
        print_stats(filtered_read_lengths)

    # quit early if only printing stats
    if args.stats_only:
        log("Only printing stats.")
        sys.exit(0)

    # determine metric for accepting a read
    likelihood_per_read = None
    if args.read_ratio is not None:
        likelihood_per_read = args.read_ratio
        log("Using read ratio of {} for read selection".format(likelihood_per_read))
    elif args.total_bases is not None:
        likelihood_per_read = 1.0 * human2comp(args.total_bases) / np.mean(filtered_read_lengths) / len(filtered_read_lengths)
        log("Selecting {} total bases, with a per-read likelihood of {}".format(args.total_bases, likelihood_per_read))
    elif args.coverage_depth is not None:
        likelihood_per_read = 1.0 * args.coverage_depth * reference_size / np.mean(filtered_read_lengths) / len(filtered_read_lengths)
        log("Selecting for a coverage depth of {} for a reference of size {} with a per-read likelihood of {}".format(args.coverage_depth, reference_size, likelihood_per_read))
    should_keep_read = lambda x: rand.random() <= likelihood_per_read

    # downsample the fq
    log("Downsampling input FASTQ")
    saved_read_lengths = list()
    with open(input) as instream, open(args.output, 'w') as outstream:
        for line in instream:
            current_read.append(line)
            # have we finished the current read?
            if len(current_read) == 4:
                if filter_read(current_read) and should_keep_read(current_read):
                    for r in current_read:
                        outstream.write(r)
                    saved_read_lengths.append(len(current_read[1]))
                current_read = list()

    # print stats
    log("  Saved reads:")
    print_stats(saved_read_lengths)
    log("Fin.")


if __name__ == "__main__":
    main()