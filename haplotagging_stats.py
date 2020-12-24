#!/usr/bin/env python3
import argparse
import sys
import os
import pysam
import collections
import numpy as np
from multithread import *

HP_TAG = "HP"
UNCLASSIFIED = 'unclassified'
CORRECT = 'correct'
INCORRECT = 'incorrect'
UNKNOWN = 'unknown'
CORRECT_RATIO = 'correct_ratio'
UNCLASSIFIED_RATIO = 'unclassified_ratio'
TOTAL_DEPTH = 'total_depth'
CHR = "chrom"
POS = "position"

OUTPUT_ORDER = [CHR, POS, TOTAL_DEPTH, CORRECT, INCORRECT, UNKNOWN, UNCLASSIFIED, CORRECT_RATIO, UNCLASSIFIED_RATIO]
OUTPUT_TYPES = {
    CHR: str,
    POS: int,
    TOTAL_DEPTH: int,
    CORRECT: int,
    INCORRECT: int,
    UNKNOWN: int,
    UNCLASSIFIED: int,
    CORRECT_RATIO: float,
    UNCLASSIFIED_RATIO: float
}

BAM_LOCATION='bam_location'
CHROM='chrom'

UNTAGGED_READS="UTR:"
TAGGED_READS="TR:"

SPACING = 1000

def parse_args(args = None):
    parser = argparse.ArgumentParser("Compares phasing for reads haplotyped by margin")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                        help='Haplotagged BAM file ')
    parser.add_argument('--truth_hap1', '-1', dest='truth_hap1', default=None, required=True, type=str,
                       help='Truth Hap1 readset')
    parser.add_argument('--truth_hap2', '-2', dest='truth_hap2', default=None, required=True, type=str,
                       help='Truth Hap2 readset')
    parser.add_argument('--output_filename', '-o', dest='output_filename', default=None, required=False, type=str,
                       help='Base of output files (will write to filename based off of input bam if not set)')
    parser.add_argument('--threads', '-t', dest='threads', default=1, required=False, type=int,
                       help='Thread count')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg):
    print(msg, file=sys.stderr)


def get_position_classifications_and_lengths(bam_location, truth_h1_ids, truth_h2_ids, region=None, verbose=True):
    # get read phasing pairs
    samfile = None
    read_count = 0
    missing_hp_count = 0
    position_classifications = collections.defaultdict(
        lambda : collections.defaultdict(lambda : 0)
    )
    tagged_lengths = []
    untagged_lengths = []
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

            if classifier == UNCLASSIFIED:
                untagged_lengths.append(epos - spos)
            else:
                tagged_lengths.append(epos - spos)

            while spos <= epos:
                pos = int(spos / SPACING)
                position_classifications[pos][classifier] += 1
                spos += SPACING
    finally:
        if samfile is not None: samfile.close()

    if verbose:
        log("Classified Read Lengths{}:".format("" if region is None else " for {}".format(region)))
        log("\tmean:   {}".format(np.mean(tagged_lengths)))
        log("\tmedian: {}".format(np.median(tagged_lengths)))
        tagged_lengths.sort()
        len_total = sum(tagged_lengths)
        len_curr = 0
        for l in tagged_lengths:
            len_curr += l
            if len_curr > len_total/2:
                log("\tN50:    {}".format(l))
                break
        total_class = sum(tagged_lengths)
        total_unclass = sum(untagged_lengths)
        log("Total Read Length: {} ({} classified, {} unclassified)".format(total_class + total_unclass, total_class,
                                                                            total_unclass))

    return position_classifications, tagged_lengths, untagged_lengths


def classify_chrom_haplotagging_service(work_queue, done_queue, bam_location=None, truth_hap1_loc=None, truth_hap2_loc=None,
                                        service_name="classify_chrom_haplotagging_service"):
    # sanity
    assert(truth_hap1_loc is not None and truth_hap2_loc is not None)

    # prep
    total_handled = 0
    failure_count = 0

    truth_h1 = set()
    truth_h2 = set()
    with open(truth_hap1_loc, 'r') as fin:
        for line in fin:
            truth_h1.add(line.split(',')[0].strip())
    with open(truth_hap2_loc, 'r') as fin:
        for line in fin:
            truth_h2.add(line.split(',')[0].strip())

    #catch overall exceptions
    fout=None
    output_name = "temp_haplotagging.{}.{}.tsv".format(os.path.basename(bam_location), "_".join(current_process().name.split()))
    try:
        fout = open(output_name, 'w')
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                # logging
                log("[{}] '{}' processing {}".format(service_name, current_process().name, f))

                # get read data
                classifications, tagged_lengths, untagged_lengths = get_position_classifications_and_lengths(
                    bam_location, truth_h1, truth_h2, region=f[CHROM])

                # write classifications
                for bucket in sorted(classifications.keys()):
                    classification = classifications[bucket]
                    # data
                    ri = classification[CORRECT]
                    ro = classification[INCORRECT]
                    unk = classification[UNKNOWN]
                    unc = classification[UNCLASSIFIED]
                    # derived data
                    td = ri + ro + unk + unc
                    classification[TOTAL_DEPTH] = td
                    classification[CORRECT_RATIO] = -1.0 if ri + ro == 0 else 1.0 * max(ri, ro) / (ri + ro)
                    classification[UNCLASSIFIED_RATIO] = -1.0 if td == 0 else 1.0 * unc / td
                    # position data
                    classification[CHR] = f[CHROM]
                    classification[POS] = SPACING * bucket

                    # write it
                    print("\t".join(map(str, [classification[x] for x in OUTPUT_ORDER])), file=fout)

                # write length data
                for tl in tagged_lengths:
                    print("{}{}".format(TAGGED_READS, tl), file=fout)
                for utl in untagged_lengths:
                    print("{}{}".format(UNTAGGED_READS, utl), file=fout)

            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                log("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        fout.close()
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))
        done_queue.put(output_name)


def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # get chroms
    samfile = pysam.AlignmentFile(args.input)
    index_stats = samfile.get_index_statistics()
    chroms = {x.contig: x.mapped for x in list(filter(lambda x: x.mapped > 0, index_stats))}
    chroms = [x[0] for x in reversed(sorted(chroms.items(), key=lambda x: x[1]))]

    # classifiy positions for reads
    service_args = {
        'bam_location': args.input,
        'truth_hap1_loc': args.truth_hap1,
        'truth_hap2_loc': args.truth_hap2,
    }

    log("Running service getting position classifications")
    total, failure, messages = run_service(classify_chrom_haplotagging_service, chroms, {}, CHROM, args.threads, service_args, log)
    log("Finished position classification service over {} entries with {} failures".format(total, failure))

    output_base = args.output_filename
    if output_base is None:
        output_base = os.path.basename(args.input).rstrip(".bam")
    haplotagging_filename="{}.haplotagging_stats.tsv".format(output_base)
    tagged_read_lengths_filename="{}.tagged_read_lengths.txt".format(output_base)
    untagged_read_lengths_filename="{}.untagged_read_lengths.txt".format(output_base)

    out_ht = None
    out_tr = None
    out_utr = None
    try:
        out_ht = open(haplotagging_filename, 'w')
        out_tr = open(tagged_read_lengths_filename, 'w')
        out_utr = open(untagged_read_lengths_filename, 'w')

        # header
        print("#" + "\t".join(OUTPUT_ORDER), file=out_ht)
        for filename in messages:
            with open(filename) as fin:
                for line in fin:
                    if line.startswith(UNTAGGED_READS):
                        out_utr.write(line.lstrip(UNTAGGED_READS))
                    elif line.startswith(TAGGED_READS):
                        out_tr.write(line.lstrip(TAGGED_READS))
                    else:
                        out_ht.write(line)
            os.remove(filename)


    finally:
        if out_ht is not None:
            out_ht.close()
        if out_tr is not None:
            out_tr.close()
        if out_utr is not None:
            out_utr.close()






if __name__ == "__main__":
    main()











