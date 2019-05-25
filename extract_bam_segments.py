#!/usr/bin/env python3
from __future__ import print_function
import argparse
import os
import subprocess
from multithread import *

INPUT="input"
OUTPUT="output"
SEGMENT="segment"

PIECE="piece"

def parse_args():
    parser = argparse.ArgumentParser("Extracts BAM segments")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                       help='BAM  input file')
    parser.add_argument('--segment_file', '-s', dest='segment_file', required=True, type=str,
                       help='File with segments descriptors')
    parser.add_argument('--output', '-o', dest='output', required=True, type=str,
                       help='Write output to file')

    parser.add_argument('--threads', '-t', dest='threads', required=False, default=1, type=int,
                       help='Write output to file')

    return parser.parse_args()



def get_segments(work_queue, done_queue, service_name="example_service"):
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

                segment_file = "{}.{}.bam".format(f[OUTPUT], f[SEGMENT])
                print("Extracting reads into {}".format(segment_file))
                with open(segment_file, 'w') as out_fh:
                    subprocess.check_call(['samtools', 'view', '-hb', f[INPUT], f[SEGMENT]], stdout=out_fh)

                done_queue.put("{}:{}".format(PIECE, segment_file))

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


def main():
    args = parse_args()

    # get segments
    print("Reading segments from {}".format(args.segment_file))
    segments = list()
    with open(args.segment_file) as seg_f:
        for line in seg_f:
            segments.append(line.strip())
    print("Found {} segments".format(len(segments)))

    # get individual pieces pieces
    pieces = list()
    print("Reading from {}".format(args.input))
    iter_args = {
        INPUT: args.input,
        OUTPUT: args.output
    }

    total, failure, messages = run_service(get_segments, segments, iter_args, SEGMENT, args.threads)
    print("Had {} failures out of {}".format(failure, total))
    for msg in messages:
        msg = msg.split(":")
        if msg[0] == PIECE:
            pieces.append(msg[1])

    # merge
    all_pieces = []
    all_pieces.extend(pieces)
    print("Merging into {}".format(args.output))
    if len(pieces) > 1000:
        large_piece_chunks = []
        curr_idx = 0
        while curr_idx < len(pieces):
            large_piece_chunk = "{}.{}.bam"
            large_piece_chunks.append(large_piece_chunk)
            current_pieces = pieces[curr_idx:min(curr_idx + 1000, len(pieces))]
            merge_cmd = ['samtools', 'merge', '-@', str(args.threads), large_piece_chunk]
            merge_cmd.extend(current_pieces)
            subprocess.check_call(merge_cmd)
            curr_idx += 1000

        all_pieces.extend(large_piece_chunks)
        pieces = large_piece_chunks

    # final merge
    merge_cmd = ['samtools', 'merge', '-@', str(args.threads), args.output]
    merge_cmd.extend(pieces)
    subprocess.check_call(merge_cmd)

    #cleanup
    print("Cleaning up tmp files")
    for file in all_pieces:
        os.remove(file)

    print("Fin.")


if __name__ == "__main__":
    main()