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

FORWARD="FORWARD"
BACKWARD="BACKWARD"
FILE_LOCATION = 'file_location'

def parse_args(args = None):
    parser = argparse.ArgumentParser("Separates a shasta assembly into individual contigs")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                        help='Input assembly')
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=False, default=None, type=str,
                        help='Prefix for output files.  Directories will be made if "/" character exists.  '
                             'Defaults to input filename')
    parser.add_argument('--strandify', '-s', dest='strandify', required=False, action='store_true', default=False,
                        help='Attempts to separate into forward and reverse strands with minimap (requires true ref)')
    parser.add_argument('--reference', '-r', dest='reference', required=False, type=str, default=None,
                        help='Input assembly')
    parser.add_argument('--threads', '-t', dest='threads', required=False, type=int, default=1,
                        help='Number of threads to use for strandification')

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


def main():
    args = parse_args()

    # mandatory prep
    input_file = args.input
    if not os.path.isfile(input_file):
        raise Exception("Input file does not exist: {}".format(input_file))
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = ".".join(input_file.split(".")[:-1])

    # optional prep
    if args.strandify:
        if args.reference is None or not os.path.isfile(args.reference):
            raise Exception("Need --reference parameter to sort reads into strands")
        args.tmp = tempfile.gettempdir()

    # ensure output directory
    if "/" in output_prefix:
        curr_dir = os.curdir
        for dir in output_prefix.split('/')[:-1]:
            curr_dir = os.path.join(curr_dir, dir)
            if not os.path.isdir(curr_dir): os.mkdir(curr_dir)

    # prep for stranding
    all_segment_filenames = list()

    # split files
    segment_filename = None
    segment_out = None
    try:
        with open(input_file, 'r') as fa_in:
            for line in fa_in:
                if line.startswith(">"):
                    # new segment
                    segment_name = line.lstrip(">").split()[0].rstrip()

                    # file handling
                    if segment_out is not None:
                        segment_out.close()
                        all_segment_filenames.append(segment_filename)
                    segment_filename = "{}.{}.fa".format(output_prefix, segment_name)
                    print("Writing segment {} to {}".format(segment_name, segment_filename))
                    segment_out = open(segment_filename, 'w')

                # santity check
                if segment_out is None:
                    raise Exception("No output stream was opened.  Perhaps a malformed input file?")

                # write
                segment_out.write(line)

    finally:
        if segment_out is not None:
            segment_out.close()
            all_segment_filenames.append(segment_filename)



    if args.strandify:
        print("Starting strand determination service")
        service_args = {'reference_location': args.reference, 'tmp_dir': args.tmp}
        iterable_args = {}
        total, failure, messages = run_service(direction_service, all_segment_filenames, iterable_args, FILE_LOCATION,
                                               args.threads, service_arguments=service_args)

        fwd_filename = "{}.forward.fa".format(output_prefix)
        bkwd_filename = "{}.backward.fa".format(output_prefix)
        print("Writing forward segments to: {}".format(fwd_filename))
        print("Writing backward segments to: {}".format(bkwd_filename))
        fwd_out = None
        bkwd_out = None
        try:
            fwd_out = open(fwd_filename, 'w')
            bkwd_out = open(bkwd_filename, 'w')

            for msg in messages:
                msg = msg.split(":")
                if msg[0] == FORWARD:
                    print("Writing {} as FWD".format(msg[1]))
                    with open(msg[1]) as file_in:
                        for line in file_in:
                            fwd_out.write(line)
                elif msg[0] == BACKWARD:
                    print("Writing {} as BKWD".format(msg[1]))
                    with open(msg[1]) as file_in:
                        for line in file_in:
                            bkwd_out.write(line)

        finally:
            if fwd_out is not None:
                fwd_out.close()
            if bkwd_out is not None:
                bkwd_out.close()


if __name__ == "__main__":
    main()