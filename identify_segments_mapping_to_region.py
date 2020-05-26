#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import sys
import numpy as np
import pysam
import time
import tempfile
import os
from subprocess import run
from contextlib import closing
from multithread import *
import subprocess

FORWARD="FORWARD"
BACKWARD="BACKWARD"
FILE_LOCATION = 'file_location'

def parse_args(args = None):
    parser = argparse.ArgumentParser("Identifies segments from an assembly which map to a truth region")
    parser.add_argument('--assembly', '-a', dest='assembly', required=True, type=str,
                        help='Input assembly')
    parser.add_argument('--reference', '-r', dest='reference', required=True, type=str,
                        help='Input reference')
    parser.add_argument('--alignment', '-l', dest='alignment', required=False, type=str, default=None,
                        help='Sorted alignment file in BAM format of assembly to reference (if provided, will use this bam instead of aligning)')

    parser.add_argument('--threads', '-t', dest='threads', required=False, type=int, default=1,
                        help='Number of threads to use')
    parser.add_argument('--minimap_preset', '-p', dest='minimap_preset', required=False, type=str, default='asm5',
                        help='Alignment preset to use for minimap')
    parser.add_argument('--include_supplementary', '-s', dest='include_supplementary', required=False, type=str, default='asm5',
                        help='Include supplementary alignments')


    return parser.parse_args() if args is None else parser.parse_args(args)


def direction_service(work_queue, done_queue, service_name="direction_service", reference_location=None, tmp_dir=None):
    # sanity check
    assert reference_location is not None
    assert tmp_dir is not None

    # prep
    total_handled = 0
    failure_count = 0

    #catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                # logging
                print("[{}] '{}' processing {}".format(service_name, current_process().name, f))
                file_location = f[FILE_LOCATION]

                forward = determine_strand(file_location, reference_location, tmp_dir)

                if forward is not None:
                    done_queue.put("{}:{}".format(FORWARD if forward else BACKWARD, file_location))

            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                print("[{}] ".format(service_name) + error)
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
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))


def determine_strand(segment_filename, reference_location, tmp_dir):

    print("Determining strand for {}".format(segment_filename))
    args = ["minimap2", "-a", "-x", "asm20", reference_location, segment_filename]
    out_filename = os.path.join(tmp_dir, "{}.sam".format(os.path.basename(segment_filename)))

    print("\t" + (" ".join(args)) + " > {}".format(out_filename))
    with open(out_filename, "w") as output_file, open(os.devnull, 'w') as devnull:
        run(args, stdout=output_file, stderr=devnull, check=True)

    # get alignment stats
    aligned_count = 0
    read_length = 0
    forward_aligned_length = 0
    backward_aligned_length = 0
    forward_primary = None
    samfile = pysam.AlignmentFile(out_filename, 'r')
    for read in samfile.fetch():
        if read.is_unmapped:
            print("\tunmapped read")
            continue

        # read length
        read_len = read.infer_read_length()
        if read_length != 0 and read_length != read_len:
            print("\tunexpected length! ({} / {})".format(read_length, read_len))
        read_length = max(read_len, read_length)
        aligned_len = read.query_alignment_length

        # soft/hard clipped
        cig_stats = read.get_cigar_stats()
        clipped = 0 #cig_stats[]

        # forward
        fwd = not read.is_reverse
        primary = not (read.is_secondary or read.is_supplementary)

        # document
        aligned_count += 1
        if fwd: forward_aligned_length += aligned_len
        else: backward_aligned_length += aligned_len
        if primary:
            if forward_primary is not None:
                print("\tfound multiple primary mappings!")
            forward_primary = fwd

    # loggit
    print("\talignments:{} len:{} fwd_len:{} bkwd_len:{} fwd_pri:{}".format(
        aligned_count, read_length, forward_aligned_length, backward_aligned_length, forward_primary))

    # recommend
    if aligned_count == 0:
        print("\tstrand could not be determined")
        return None
    is_forward = None
    if forward_aligned_length > backward_aligned_length:
        is_forward = True
    elif backward_aligned_length > forward_aligned_length:
        is_forward = False
    else:
        is_forward = forward_primary
    print("\tis_forward: {}".format(is_forward))

    return is_forward

def get_fasta_name(name):
    parts = os.path.basename(name).split(".")
    while parts[-1] in ["fasta", "fa"]:
        parts = parts[:-1]
    return "_".join(parts)


def main():
    args = parse_args()

    # sanity check
    assembly_file = args.assembly
    if not os.path.isfile(assembly_file):
        raise Exception("Input assembly does not exist: {}".format(assembly_file))
    reference_file = args.reference
    if not os.path.isfile(reference_file):
        raise Exception("Input reference does not exist: {}".format(reference_file))
    alignment_file = args.alignment
    if alignment_file is not None and not os.path.isfile(alignment_file):
        raise Exception("Alignment file was specified but does not exist: {}".format(alignment_file))

    # names
    assembly_id = get_fasta_name(assembly_file)
    reference_id = get_fasta_name(reference_file)

    # configuration
    minimap_preset = args.minimap_preset
    threads = args.threads
    include_supplementary = args.include_supplementary

    # do alignment if necessary
    if alignment_file is None:
        unsorted_file = "REF_{}.READS_{}.unsorted.sam"
        alignment_file = "REF_{}.READS_{}.bam"

        # align
        with open(unsorted_file, 'w') as out:
            subprocess.check_call(['minimap2', '-a', '-x', minimap_preset, '-t', str(threads), reference_file, assembly_file], stdout=out)

        # sort
        subprocess.check_call(['samtools', 'sort', '-@', str(threads), '-o', alignment_file, unsorted_file])
        if not os.path.isfile(alignment_file):
            raise Exception("Alignment file was not created: {}".format(alignment_file))

        os.remove(unsorted_file)

    # index
    if not os.path.isfile("{}.bai"):
        subprocess.check_call(["samtools", "index", "-@", str(threads), alignment_file])

    bamfile = pysam.AlignmentFile(alignment_file, 'r')
    for read in bamfile.fetch():
        if read.is_unmapped:
            print("\tunmapped read")
            continue

        # read length
        read_len = read.infer_read_length()
        if read_length != 0 and read_length != read_len:
            print("\tunexpected length! ({} / {})".format(read_length, read_len))
        read_length = max(read_len, read_length)
        aligned_len = read.query_alignment_length

        # soft/hard clipped
        cig_stats = read.get_cigar_stats()
        clipped = 0  # cig_stats[]

        # forward
        fwd = not read.is_reverse
        primary = not (read.is_secondary or read.is_supplementary)

        # document
        aligned_count += 1
        if fwd:
            forward_aligned_length += aligned_len
        else:
            backward_aligned_length += aligned_len
        if primary:
            if forward_primary is not None:
                print("\tfound multiple primary mappings!")
            forward_primary = fwd







if __name__ == "__main__":
    main()