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


CHR_IDX = 0
START_IDX = 1
END_IDX = 2
SCORE_IDX = 4
STRAND_IDX = 5
COVERAGE_IDX = 9
METHYL_PERCENTAGE_IDX = 10
FWD = "+"
REV = "-"

class MethylRecord:
    def __init__(self, chr, pos, score, strand, coverage, methyl_percent):
        self.chr = chr
        self.pos = int(pos)
        self.strand = strand
        self.coverage = int(coverage)
        self.methyl_ratio = int(methyl_percent) / 100.0
        self.methyl_count = round(self.coverage * self.methyl_ratio)


class MethylLocus:
    def __init__(self, fwd_methyl: MethylRecord, rev_methyl: MethylRecord, sanity_check = True):
        if sanity_check:
            assert fwd_methyl.chr == rev_methyl.chr, "MethylLocus chr mismatch"
            assert fwd_methyl.pos == rev_methyl.pos - 1, "MethylLocus pos mismatch"
            assert fwd_methyl.strand == FWD and rev_methyl.strand == REV, "MethyLocus strand mismatch"
        self.fwd_methyl = fwd_methyl
        self.rev_methyl = rev_methyl
        self.chr = fwd_methyl.chr
        self.start_pos = fwd_methyl.pos
        self.end_pos = rev_methyl.pos
        self.coverage = fwd_methyl.coverage + rev_methyl.coverage
        self.avg_methyl_ratio = (fwd_methyl.methyl_ratio + rev_methyl.methyl_ratio) / 2


def parse_args():
    parser = argparse.ArgumentParser("Produce stats and refactoring for methylation bed file based from bisulfite data")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                       help='methylation bed input file')

    parser.add_argument('--min_coverage', '-d', dest='min_coverage', required=False, default=8, type=int,
                       help='Minimum depth/coverage for record to be used in analysis')

    parser.add_argument('--output_base', '-o', dest='output_base', required=False, default=None, type=str,
                       help='Write output files, otherwise will use input bed name as prefix')
    parser.add_argument('--write_full_bed', '-1', dest='write_full_bed', default=False, required=False, action='store_true',
                       help='Write filtered and merged methylation loci')
    parser.add_argument('--write_quartile_beds', '-4', dest='write_quartile_beds', default=False, required=False, action='store_true',
                       help='Stratify calls into quartiles')
    parser.add_argument('--write_quintile_beds', '-5', dest='write_quintile_beds', default=False, required=False, action='store_true',
                       help='Stratify calls into quintiles')

    return parser.parse_args()


def log(msg):
    print(msg, file=sys.stderr)


def print_histogram(data, idx_str_function, size = 32):
    max_size = max(data)
    total = sum(data)
    for i, x in enumerate(data):
        count = int(size * x / max_size)
        log("{:>10s}: {}{} {:>8d} {:.3f}".format(idx_str_function(i), "#" * count, " " * (size - count), x, x/total))


def print_methyl_loci_differences(loci):
    # prep
    size_factor = 10
    size = size_factor * 2 + 1
    coverage_value_differences = [0 for x in range(size)]
    coverage_ratio_differences = [0 for x in range(size)]
    methyl_value_differences = [0 for x in range(size)]
    methyl_ratio_differences = [0 for x in range(size)]

    for locus in loci:
        fwd = locus.fwd_methyl
        rev = locus.rev_methyl

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


def main():
    args = parse_args()

    # data we want
    all_methyl_records = list()
    all_methyl_loci = list()

    # read input file
    log("Reading records from {}".format(args.input))
    with open(args.input) as fin:
        for linenr, line in enumerate(fin):
            if line.startswith("#") or len(line.strip()) == 0: continue
            parts = line.split("\t")
            assert len(parts) > METHYL_PERCENTAGE_IDX, "Malformed line {}: '{}'".format(linenr, line)
            chr, pos, score, strand, coverage, methyl_percentage = parts[CHR_IDX], parts[START_IDX], parts[SCORE_IDX], \
                                                                   parts[STRAND_IDX], parts[COVERAGE_IDX], \
                                                                   parts[METHYL_PERCENTAGE_IDX]
            record = MethylRecord(chr, pos, score, strand, coverage, methyl_percentage)
            all_methyl_records.append(record)

    log("Got {} total methyl records".format(len(all_methyl_records)))
    assert(len(all_methyl_records) % 2 == 0)

    # construct loci
    skipped_for_depth = 0
    fwd_record = None
    rev_record = None
    for idx, record in enumerate(all_methyl_records):
        # iteration
        if rev_record is not None:
            fwd_record = None
            rev_record = None

        # identify fwd/rev
        if fwd_record is None:
            assert record.strand == FWD, "Expected {} record at idx {}: {}".format(FWD, idx, record)
            fwd_record = record
            continue
        else:
            assert record.strand == REV, "Expected {} record at idx {}: {}".format(REV, idx, record)
            rev_record = record

        # filter and build methyl locus
        if fwd_record.coverage < args.min_coverage or rev_record.coverage < args.min_coverage:
            skipped_for_depth += 1
            continue

        # build locus
        methyl_locus = MethylLocus(fwd_record, rev_record)
        all_methyl_loci.append(methyl_locus)

    # loggit
    log("Got {} methyl loci, skipped {} because of missing depths".format(len(all_methyl_loci), skipped_for_depth))
    print_methyl_loci_differences(all_methyl_loci)

    # write
    should_write = True in [args.write_quartile_beds, args.write_quintile_beds, args.write_full_bed]
    if should_write:
        output_base = args.output_base
        if output_base is None:
            output_base = ".".join(os.path.basename(args.input).split(".")[:-1])
        log("\nWriting output to files with base {}:".format(output_base))

        full_file_name = "{}.full.bed".format(output_base)
        quartile_file_names = list(map(lambda x: "{}.{}.bed".format(output_base, x), ["0-25", "25-50", "50-75", "75-100"]))
        quintile_file_names = list(map(lambda x: "{}.{}.bed".format(output_base, x), ["0-20", "20-40", "40-60", "60-80", "80-100"]))
        full_out = None
        quartile_outs = None
        quintile_outs = None

        try:
            # open files
            output_files = []
            if args.write_full_bed:
                full_out = open(full_file_name, 'w')
                output_files.append(full_file_name)
            if args.write_quartile_beds:
                quartile_outs = list(map(lambda x: open(x, 'w'), quartile_file_names))
                output_files.extend(quartile_file_names)
            if args.write_quintile_beds:
                quintile_outs = list(map(lambda x: open(x, 'w'), quintile_file_names))
                output_files.extend(quintile_file_names)
            # finish logging
            for fn in output_files:
                log("\t{}".format(fn))

            # look at all loci and write out
            for locus in all_methyl_loci:
                bed_data = [locus.chr, locus.start_pos, locus.end_pos + 1, locus.coverage, round(locus.avg_methyl_ratio, 5)]
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
