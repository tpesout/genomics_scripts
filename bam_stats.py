#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import numpy as np
import pysam
import math

# read keys
R_ID = "id"
R_START_POS = "start_pos"
R_END_POS = "end_pos"
R_LENGTH = "length"
R_SECONDARY = "secondary_alignment"
R_SUPPLEMENTARY = "supplementary_alignment"
R_MAPPING_QUALITY = "mapping_quality"
R_CHROMOSOME = "chromosome"

#bam summary
B_READ_COUNT = "read_count"
B_SECONDARY_COUNT = "secondary_count"
B_SUPPLEMENTARY_COUNT = "supplementary_count"
B_FILTERED_READ_COUNT = "filtered_read_count"
B_CHROMOSOME = "chromosome"
B_MEDIAN_QUAL = "median_quality"

#length summary
L_LOG_LENGTH_BUCKETS = "log_length_buckets"
L_LOG_LENGTH_BUCKETS_ALT = "log_length_buckets_alt"
L_MIN = "min_length"
L_MAX = "max_length"
L_AVG = "avg_length"
L_MED = 'median_len'
L_STD = "std_lenght"
L_N50 = "N50"
L_LOG_BASE = "log_base"
L_LOG_BASE_ALT = "log_base_alt"
L_LOG_MAX = "log_max"
L_ALL_LENGTHS = "all_lengths"

# depth summary
D_MAX = "max_depth"
D_MIN = "min_depth"
D_MED = "median_depth"
D_AVG = "avg_depth"
D_STD = "std_depth"
D_ALL_DEPTHS = "all_depths"
D_ALL_DEPTH_POSITIONS = "all_depth_positions"
D_ALL_DEPTH_MAP = "all_depth_map"
D_LOG_DEPTH_BINS = "all_depth_bins"
D_DEPTH_BINS = "depth_bins"
D_SPACING = "depth_spacing"
D_START_IDX = "depth_start_idx"
D_RANGE = "region"

# misc
GENOME_KEY = "genome"
LENGTH_LOG_BASE_DEFAULT=2
LENGTH_LOG_BASE_ALT=10
LENGTH_LOG_MAX_DEFAULT=32
percent = lambda small, big: int(100.0 * small / big) if big != 0 else 0.0
log = lambda value, base: math.log(value, base) if value > 0 else 0


def parse_args(args = None):
    parser = argparse.ArgumentParser("Provides statistics on a BAM/SAM file")
    parser.add_argument('--input_glob', '-i', dest='input_glob', default=None, required=True, type=str,
                        help='Glob matching SAM or BAM file(s)')
    parser.add_argument('--output_file', '-o', dest='output_file', default=None, required=False, type=str,
                       help='Write output to file instead of stdout')
    parser.add_argument('--generic_stats', '-g', dest='generic_stats', action='store_true', default=False,
                        help='Print generic stats for all files')
    parser.add_argument('--read_length', '-l', dest='read_length', action='store_true', default=False,
                        help='Print statistics on read length for all files')
    parser.add_argument('--read_depth', '-d', dest='read_depth', action='store_true', default=False,
                        help='Print statistics on read depth for all files')
    parser.add_argument('--genome_only', dest='genome_only', action='store_true', default=False,
                        help='Print only statistics for the whole genome (do not print stats for individual chroms)')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                        help='Print histograms for length and depth')
    parser.add_argument('--silent', '-V', dest='silent', action='store_true', default=False,
                        help='Print nothing')
    parser.add_argument('--depth_spacing', '-s', dest='depth_spacing', action='store', default=1000, type=int,
                        help='How far to sample read data')
    parser.add_argument('--region', '-r', dest='region', action='store', default=None,
                        help='Whether to only calculate depth within a range, ie: \'100000-200000\'')
    parser.add_argument('--filter_secondary', dest='filter_secondary', action='store_true', default=False,
                        help='Filter secondary alignments out')
    parser.add_argument('--filter_supplementary', dest='filter_supplementary', action='store_true', default=False,
                        help='Filter supplemenary alignments out')
    parser.add_argument('--filter_read_length_min', dest='read_length_min', action='store', default=None, type=int,
                        help='Removes reads with length below this')
    parser.add_argument('--filter_read_length_max', dest='read_length_max', action='store', default=None, type=int,
                        help='Removes reads with length above this')
    parser.add_argument('--filter_alignment_threshold_min', dest='min_alignment_threshold', action='store',
                        default=None, type=int, help='Minimum alignment quality threshold')

    parser.add_argument('--produce_read_length_tsv', dest='read_length_tsv', action='store',
                        default=None, type=str, help='Produce a TSV with read lengths, named as this parameter')
    parser.add_argument('--read_length_bucket_size', dest='read_length_bucket_size', action='store',
                        default=50000, type=int, help='Bucket size for read length TSV')

    return parser.parse_args() if args is None else parser.parse_args(args)


def get_read_summary(read):
    summary = {
        R_START_POS: read.reference_start,
        R_END_POS: read.reference_end,
        R_LENGTH: read.reference_length,
        R_ID: read.query_name,
        R_SECONDARY: read.is_secondary,
        R_SUPPLEMENTARY: read.is_supplementary,
        R_MAPPING_QUALITY: read.mapping_quality,
        R_CHROMOSOME: read.reference_name
    }
    return summary


def caluclate_n50(lengths):
    half_total_length = sum(lengths) / 2
    lengths.sort(reverse=True)
    for l in lengths:
        half_total_length -= l
        if half_total_length <= 0:
            return l


def get_read_length_summary(read_summaries, length_log_base=LENGTH_LOG_BASE_DEFAULT,
                            length_log_max=LENGTH_LOG_MAX_DEFAULT, length_log_base_alt=LENGTH_LOG_BASE_ALT):
    # what we look for
    log_length_bins = [0 for _ in range(length_log_max)]
    log_lenght_alt_bins = [0 for _ in range(length_log_base_alt)] if length_log_base_alt is not None else None
    all_lengths = []

    for read_summary in read_summaries:
        if read_summary[R_LENGTH] is None:
            #todo
            continue
        log_length_bins[int(log(read_summary[R_LENGTH], length_log_base))] += 1.0
        if length_log_base_alt is not None: log_lenght_alt_bins[int(log(read_summary[R_LENGTH], length_log_base_alt))] += 1.0
        all_lengths.append(read_summary[R_LENGTH])

    summary = {
        L_LOG_LENGTH_BUCKETS: log_length_bins,
        L_LOG_LENGTH_BUCKETS_ALT: log_lenght_alt_bins,
        L_MAX: max(all_lengths) if len(all_lengths) != 0 else 0,
        L_MIN: min(all_lengths) if len(all_lengths) != 0 else 0,
        L_AVG: np.mean(all_lengths) if len(all_lengths) != 0 else 0,
        L_MED: np.median(all_lengths) if len(all_lengths) != 0 else 0,
        L_STD: np.std(all_lengths) if len(all_lengths) != 0 else 0,
        L_N50: caluclate_n50(all_lengths) if len(all_lengths) != 0 else 0,
        L_LOG_BASE: length_log_base,
        L_LOG_BASE_ALT: length_log_base_alt,
        L_LOG_MAX: length_log_max,
        L_ALL_LENGTHS: all_lengths
    }

    return summary


def graph_read_length_summary(summary, title, save_name=None):
    log_lengths = list(summary.keys())
    log_lengths.sort()

    x = log_lengths
    y = [summary[l] for l in log_lengths]

    import matplotlib.pyplot as plt
    plt.bar(x, y, color='g')
    plt.xlabel("Read Length (Log {})".format(LENGTH_LOG_BASE_DEFAULT))
    plt.ylabel("Count")
    plt.title(title)
    if save_name is not None:
        plt.savefig(save_name)
    else:
        plt.show()


def print_generic_read_stats(summary, output, verbose=False, genome_only=False):
    keys = list(summary.keys())
    keys.sort()
    #print genome last
    if GENOME_KEY in keys:
        keys.remove(GENOME_KEY)
        keys.append(GENOME_KEY)
    #print them all
    for chrom in keys:

        # if we only have one chrom, skip genome reporting
        if len(keys) == 2 and chrom == GENOME_KEY and not genome_only: continue
        # if genome_only, skip chroms
        if genome_only and chrom != GENOME_KEY: continue

        B_READ_COUNT = "read_count"
        B_SECONDARY_COUNT = "secondary_count"
        B_SUPPLEMENTARY_COUNT = "supplementary_count"
        print("\tGENERIC STATS: {}".format(chrom), file=output)
        print("\t\tcount        : {}".format(summary[chrom][B_READ_COUNT]), file=output)
        print("\t\tsecondary    : {} ({}%)".format(summary[chrom][B_SECONDARY_COUNT],
                                                   percent(summary[chrom][B_SECONDARY_COUNT],
                                                           summary[chrom][B_READ_COUNT])), file=output)
        print("\t\tsupplenatary : {} ({}%)".format(summary[chrom][B_SUPPLEMENTARY_COUNT],
                                                   percent(summary[chrom][B_SUPPLEMENTARY_COUNT],
                                                           summary[chrom][B_READ_COUNT])), file=output)
        print("\t\tmedian qual  : {}".format(summary[chrom][B_MEDIAN_QUAL]), file=output)


def print_log_binned_data(log_bins, output, indent_count=3):
    extant_log_bins = list(filter(lambda x: log_bins[x] != 0, [x for x in range(len(log_bins))]))
    if len(extant_log_bins) == 0:
        print("{} [No Data]".format('\t'*indent_count), file=output)
        return
    max_bucket = max(extant_log_bins)
    min_bucket = min(extant_log_bins)
    max_bucket_size = max(log_bins)
    total_bucket_size = sum(log_bins)
    total_bucket_size_left = total_bucket_size
    for bucket in range(min_bucket, max_bucket + 1):
        id = "%3d:" % bucket
        count = int(log_bins[bucket])
        pound_count = int(32.0 * count / max_bucket_size)
        of_total = 1.0 * count / total_bucket_size
        at_least = 1.0 * total_bucket_size_left / total_bucket_size
        total_bucket_size_left -= count
        print("{} {} {}{} {:6d}\t({:.3f} of total)\t({:.3f} at least)".format(
            '\t'*indent_count, id, "#"*pound_count, " "*(32 - pound_count), count, of_total, at_least), file=output)


def print_read_length_summary(summary, output, verbose=False, genome_only=False):
    keys = list(summary.keys())
    keys.sort()
    #print genome last
    if GENOME_KEY in keys:
        keys.remove(GENOME_KEY)
        keys.append(GENOME_KEY)
    #print them all
    for chrom in keys:

        # if we only have one chrom, skip genome reporting
        if len(keys) == 2 and chrom == GENOME_KEY and not genome_only: continue
        # if genome_only, skip chroms
        if genome_only and chrom != GENOME_KEY: continue

        print("\tREAD LENGTHS: {}".format(chrom), file=output)
        print("\t\tmin: {}".format(summary[chrom][L_MIN]), file=output)
        print("\t\tmax: {}".format(summary[chrom][L_MAX]), file=output)
        print("\t\tmed: {}".format(summary[chrom][L_MED]), file=output)
        print("\t\tavg: {}".format(int(summary[chrom][L_AVG])), file=output)
        print("\t\tstd: {}".format(int(summary[chrom][L_STD])), file=output)
        print("\t\tN50: {}".format(int(summary[chrom][L_N50])), file=output)

        if verbose:
            print("\t\tread length log_{}:".format(summary[chrom][L_LOG_BASE]), file=output)
            print_log_binned_data(summary[chrom][L_LOG_LENGTH_BUCKETS], output)
            if L_LOG_LENGTH_BUCKETS_ALT in summary[chrom] and summary[chrom][L_LOG_LENGTH_BUCKETS_ALT] is not None:
                print("\t\tread length log_{}:".format(summary[chrom][L_LOG_BASE_ALT]), file=output)
                print_log_binned_data(summary[chrom][L_LOG_LENGTH_BUCKETS_ALT], output)


def  get_read_depth_summary(read_summaries, spacing, region=None):
    S, E = 's', 'e'

    # get reads which start or end on spacing interval
    positions = []
    for summary in read_summaries:
        if summary[R_LENGTH] is None:
            #todo
            continue
        positions.append((S, int(summary[R_START_POS]/spacing)))
        positions.append((E, int(summary[R_END_POS]/spacing)))

    # sort them: we iterate from the high end by popping off
    positions.sort(key=lambda x: x[1])
    start_idx = positions[0][1]
    end_idx = positions[-1][1]

    # data we care about
    depths = [0 for _ in range(end_idx - start_idx + 1)]
    depth_positions = []

    # iterate over all read starts and ends
    depth = 0
    idx = end_idx
    while idx >= start_idx:
        curr = positions.pop()
        while curr[1] == idx:
            if curr[0] == E: depth += 1
            if curr[0] == S: depth -= 1
            # iterate
            if len(positions) == 0: break
            else: curr = positions.pop()
        positions.append(curr)
        # save and iterate
        pos = idx - start_idx
        depths[pos] = depth
        depth_positions.append(idx)
        idx -= 1

    #todo I don't like that I don't know why I need to do this
    depth_positions.reverse()

    assert depth == 0
    assert len(positions) == 1
    assert len(depths) == len(depth_positions)

    depth_map = {pos: depth for pos, depth in zip(depth_positions, depths)}

    # check range before outputting summary
    included_range = None
    if region is not None:
        # get range
        included_range = list(map(int, region.split(':')[-1].split("-")))
        if len(included_range) != 2:
            raise Exception("Malformed depth range: '{}'".format("-".join(map(str, included_range))))
        range_start = int(included_range[0]/spacing)
        range_end = int(included_range[1]/spacing)
        # sanity check
        if range_start > end_idx or range_end < start_idx or range_start >= range_end:
            raise Exception("Range {} outside of bounds of chunks: {}".format("-".join(map(str, included_range)),
                                                                              "-".join(map(str, [start_idx*spacing, end_idx*spacing]))))

        new_depths = list()
        new_depth_positions = list()
        new_depth_map = dict()
        for i in range(range_start, range_end):
            new_depth_positions.append(i)
            depth = depth_map[i] if i in depth_map else 0
            new_depths.append(depth)
            new_depth_map[i] = depth

        # update values
        depths = new_depths
        depth_positions = new_depth_positions
        depth_map = new_depth_map
        start_idx = max(start_idx, range_start)
        assert len(depths) > 0
        assert len(new_depths) == len(new_depth_positions)

    # get read depth bins
    log_depth_bins = [0 for _ in range(int(log(max(1,max(depths)), 2)) + 1)]
    depth_bins = [0 for _ in range(int(max(1,max(depths))) + 1)]
    for depth in depths:
        depth_bins[depth] += 1
        if depth == 0:
            log_depth_bins[0] += 1
        else:
            log_depth_bins[int(log(depth, 2))] += 1

    # get depth summary
    summary = {
        D_MAX: max(depths),
        D_MIN: min(depths),
        D_MED: np.median(depths),
        D_AVG: np.mean(depths),
        D_STD: np.std(depths),
        D_ALL_DEPTHS: depths,
        D_ALL_DEPTH_POSITIONS: depth_positions,
        D_ALL_DEPTH_MAP: depth_map,
        D_LOG_DEPTH_BINS: log_depth_bins,
        D_DEPTH_BINS: depth_bins,
        D_SPACING: spacing,
        D_START_IDX: start_idx
    }
    if included_range is not None:
        summary[D_RANGE] = included_range

    return summary


def get_genome_depth_summary(summaries):
    depths = list()
    for summary in summaries.values():
        depths.extend(summary[D_ALL_DEPTHS])

    summary = {
        D_MAX: max(depths) if len(depths) > 0 else 0,
        D_MIN: min(depths) if len(depths) > 0 else 0,
        D_MED: np.median(depths) if len(depths) > 0 else 0,
        D_AVG: np.mean(depths) if len(depths) > 0 else 0,
        D_STD: np.std(depths) if len(depths) > 0 else 0,
        D_ALL_DEPTHS: None,
        D_ALL_DEPTH_POSITIONS: None,
        D_ALL_DEPTH_MAP: None,
        D_LOG_DEPTH_BINS: None,
        D_SPACING: None,
        D_START_IDX: None,
        D_RANGE: None
    }
    return summary


def write_read_length_tsv(reads, filename, bucket_size=50000):
    length_to_bucket = lambda x: int(1.0 * x / bucket_size)
    read_lengths = dict()

    for read in reads:
        bucket = length_to_bucket(read[R_LENGTH])
        while len(read_lengths) <= bucket:
            read_lengths[len(read_lengths)] = 0
        read_lengths[bucket] += 1

    with open(filename, 'w') as output:
        output.write("#min_length\tmax_length\tread_count\n")
        started = False
        for i in range(len(read_lengths)):
            if not started and read_lengths[i] == 0:
                continue
            started = True
            output.write("{}\t{}\t{}\n".format(bucket_size * i, bucket_size * i + bucket_size - 1, read_lengths[i]))


def print_read_depth_summary(summary, output, verbose=False, genome_only=False):
    keys = list(summary.keys())
    keys.sort()
    #print genome last
    if GENOME_KEY in keys:
        keys.remove(GENOME_KEY)
        keys.append(GENOME_KEY)

    for chrom in keys:
        # if we only have one chrom, skip genome reporting
        if len(keys) == 2 and chrom == GENOME_KEY and not genome_only: continue
        # if genome_only, skip chroms
        if genome_only and chrom != GENOME_KEY: continue

        print("\tREAD DEPTHS: {}".format(chrom), file=output)
        print("\t\tmax: {}".format(summary[chrom][D_MAX]), file=output)
        print("\t\tmin: {}".format(summary[chrom][D_MIN]), file=output)
        print("\t\tmed: {}".format(summary[chrom][D_MED]), file=output)
        print("\t\tavg: {}".format(summary[chrom][D_AVG]), file=output)
        print("\t\tstd: {}".format(summary[chrom][D_STD]), file=output)

        if chrom != GENOME_KEY and summary[chrom][D_LOG_DEPTH_BINS] is not None:
            log_depth_bins = summary[chrom][D_LOG_DEPTH_BINS]
            total_depths = sum(log_depth_bins)
            log_depth_pairs = [(i, log_depth_bins[i]) for i in range(len(log_depth_bins))]
            log_depth_pairs.sort(key=lambda x: x[1], reverse=True)
            print("\t\tmost frequent read depths [floor(log2(depth))]:", file=output)
            for i in range(0,min(len(list(filter(lambda x: x[1] != 0, log_depth_pairs))), 3)):
                print("\t\t\t#{}: depth:{}({}) count:{} ({}%)".format(i + 1, log_depth_pairs[i][0],
                                                                      2**log_depth_pairs[i][0], log_depth_pairs[i][1],
                                                                      int(100.0 * log_depth_pairs[i][1] / total_depths)),
                      file=output)

        if chrom != GENOME_KEY and D_RANGE in summary[chrom] and summary[chrom][D_RANGE] is not None:
            print("\t\tdepths with spacing {}{}:".format(summary[chrom][D_SPACING],
                                                         "" if summary[chrom][D_RANGE] is None else
                                                         ", and range {}".format(summary[chrom][D_RANGE])), file=output)
            for idx in summary[chrom][D_ALL_DEPTH_POSITIONS]:
                depth = summary[chrom][D_ALL_DEPTH_MAP][idx]
                id = "{:12d}:".format(idx * summary[chrom][D_SPACING])
                pound_count = int(32.0 * depth / summary[chrom][D_MAX]) if summary[chrom][D_MAX] != 0 else 0
                print("\t\t\t{} {} {}".format(id, '#' * pound_count, depth), file=output)

        if verbose:
            if chrom != GENOME_KEY and summary[chrom][D_DEPTH_BINS] is not None:
                print("\t\tread depth at above intervals:", file=output)
                print_log_binned_data(summary[chrom][D_DEPTH_BINS], output)

            if chrom != GENOME_KEY and summary[chrom][D_LOG_DEPTH_BINS] is not None:
                print("\t\tread depth log_2 at above intervals:", file=output)
                print_log_binned_data(log_depth_bins, output)
                # max_bucket = max(list(filter(lambda x: log_depth_bins[x] != 0, [x for x in range(16)])))
                # min_bucket = min(list(filter(lambda x: log_depth_bins[x] != 0, [x for x in range(16)])))
                # max_bucket_size = max(log_depth_bins)
                # for bucket in range(min_bucket, max_bucket + 1):
                #     id = "%3d:" % bucket
                #     count = log_depth_bins[bucket]
                #     pound_count = int(32.0 * count / max_bucket_size)
                #     print("\t\t\t{} {} {}".format(id, "#" * pound_count, count))


def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)
    if True not in [args.generic_stats, args.read_depth, args.read_length]:
        args.generic_stats, args.read_depth, args.read_length = True, True, True

    # get filenames, sanity check
    in_alignments = glob.glob(args.input_glob)
    if len(in_alignments) == 0:
        print("No files matching {}".format(args.input_glob))
        return 1
    else:
        if not args.silent: print("Analyzing {} files".format(len(in_alignments)))

    # data we care about
    bam_summaries = dict()
    length_summaries = dict()
    depth_summaries = dict()
    all_read_summaries = list()

    # iterate over all alignments
    for alignment_filename in in_alignments:
        # sanity check
        if not (alignment_filename.endswith("sam") or alignment_filename.endswith("bam")):
            print("Matched file {} has unexpected filetype".format(alignment_filename))
            continue

        # prep
        bam_summaries[alignment_filename] = {}
        length_summaries[alignment_filename] = {}
        depth_summaries[alignment_filename] = {}

        # data we're gathering
        read_summaries = list()
        chromosomes =  set()

        # get read data we care about
        samfile = None
        read_count = 0
        try:
            if not args.silent: print("Read {}:".format(alignment_filename))
            samfile = pysam.AlignmentFile(alignment_filename, 'rb' if alignment_filename.endswith("bam") else 'r')
            region_chrom=None if args.region is None else args.region.split(":")[0]
            region_start=None if args.region is None else int(args.region.split(":")[1].split("-")[0])
            region_end=None if args.region is None else int(args.region.split(":")[1].split("-")[1])
            for read in (samfile.fetch() if args.region is None else samfile.fetch(region_chrom, region_start, region_end)):
                read_count += 1
                summary = get_read_summary(read)
                read_summaries.append(summary)
                chromosomes.add(read.reference_name)
        finally:
            if samfile is not None: samfile.close()

        bad_read_count = len(list(filter(lambda x: x[R_LENGTH] is None, read_summaries)))
        if bad_read_count > 0 and not args.silent:
            print("\tGot {}/{} ({}%) bad reads in {}. Filtering out."
                  .format(bad_read_count, len(read_summaries), int(100.0 * bad_read_count / len(read_summaries)),
                          alignment_filename), file=sys.stderr)
        read_summaries = list(filter(lambda x: x[R_LENGTH] is not None, read_summaries))

        # filter if appropriate
        did_filter = False
        if args.filter_secondary:
            if not args.silent: print("\tFiltering secondary reads")
            read_summaries = list(filter(lambda x: not x[R_SECONDARY], read_summaries))
            did_filter = True
        if args.filter_supplementary:
            if not args.silent: print("\tFiltering supplementary reads")
            read_summaries = list(filter(lambda x: not x[R_SUPPLEMENTARY], read_summaries))
            did_filter = True
        if args.min_alignment_threshold is not None:
            if not args.silent: print("\tFiltering reads below map quality {}".format(args.min_alignment_threshold))
            read_summaries = list(filter(lambda x: x[R_MAPPING_QUALITY] >= args.min_alignment_threshold, read_summaries))
            did_filter = True
        if args.read_length_min is not None:
            if not args.silent: print("\tFiltering reads below length {}".format(args.read_length_min))
            read_summaries = list(filter(lambda x: x[R_LENGTH] >= args.read_length_min, read_summaries))
            did_filter = True
        if args.read_length_max is not None:
            if not args.silent: print("\tFiltering reads above length {}".format(args.read_length_max))
            read_summaries = list(filter(lambda x: x[R_LENGTH] <= args.read_length_max, read_summaries))
            did_filter = True
        if did_filter:
            filtered_read_count = len(read_summaries)
            if not args.silent:
                print("\tFiltering removed {}/{} reads ({}% remaining) "
                      .format((read_count - filtered_read_count), read_count, 100 * filtered_read_count / read_count))


        # summarize
        for chromosome in chromosomes:

            #prep
            chromosome_reads = list()
            chrom_read_count = 0
            chrom_sec_count = 0
            chrom_sup_count = 0

            # analyze
            for read in read_summaries:
                if read[R_CHROMOSOME] == chromosome:
                    chromosome_reads.append(read)
                    chrom_read_count += 1
                    if read[R_SECONDARY]: chrom_sec_count += 1
                    if read[R_SUPPLEMENTARY]: chrom_sup_count += 1

            # filtered out all reads
            if len(chromosome_reads) == 0: continue

            # summarize
            bam_summaries[alignment_filename][chromosome] = {
                B_READ_COUNT: chrom_read_count,
                B_SECONDARY_COUNT: chrom_sec_count,
                B_SUPPLEMENTARY_COUNT: chrom_sup_count,
                B_MEDIAN_QUAL: np.median(list(map(lambda x: x[R_MAPPING_QUALITY], chromosome_reads))),
                B_FILTERED_READ_COUNT:len(chromosome_reads),
                B_CHROMOSOME: chromosome
            }
            length_summaries[alignment_filename][chromosome] = get_read_length_summary(chromosome_reads)
            depth_summaries[alignment_filename][chromosome] = get_read_depth_summary(chromosome_reads,
                                                                                     spacing=args.depth_spacing,
                                                                                     region=args.region)

        # whole file summaries
        bam_summaries[alignment_filename][GENOME_KEY] = {
            B_READ_COUNT: read_count,
            B_SECONDARY_COUNT: len(list(filter(lambda x: x[R_SECONDARY], read_summaries))),
            B_SUPPLEMENTARY_COUNT: len(list(filter(lambda x: x[R_SUPPLEMENTARY], read_summaries))),
            B_MEDIAN_QUAL: np.median(list(map(lambda x: x[R_MAPPING_QUALITY], read_summaries))),
            B_FILTERED_READ_COUNT: len(read_summaries),
            B_CHROMOSOME: GENOME_KEY
        }
        length_summaries[alignment_filename][GENOME_KEY] = get_read_length_summary(read_summaries)
        depth_summaries[alignment_filename][GENOME_KEY] = get_genome_depth_summary(depth_summaries[alignment_filename])

        # print
        try:
            output = sys.stdout
            if args.output_file is not None:
                output = open(args.output_file, 'w+')
            if args.generic_stats:
                if not args.silent: print_generic_read_stats(bam_summaries[alignment_filename], output,
                                                             verbose=args.verbose, genome_only=args.genome_only)
            if args.read_length:
                if not args.silent: print_read_length_summary(length_summaries[alignment_filename], output,
                                                             verbose=args.verbose, genome_only=args.genome_only)
            if args.read_depth:
                if not args.silent: print_read_depth_summary(depth_summaries[alignment_filename], output,
                                                             verbose=args.verbose, genome_only=args.genome_only)
        finally:
            if args.output_file is not None:
                output.close()

        # save
        all_read_summaries.extend(read_summaries)

    # do whole run analysis
    if args.read_length and not args.silent and len(in_alignments) > 1:
        print_read_length_summary({'ALL_FILES':get_read_length_summary(all_read_summaries)}, verbose=args.verbose)

    # tsv
    if args.read_length_tsv is not None:
        write_read_length_tsv(all_read_summaries, args.read_length_tsv, args.read_length_bucket_size)

    return bam_summaries, length_summaries, depth_summaries




if __name__ == "__main__":
    main()
