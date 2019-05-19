#!/usr/bin/env python3
from __future__ import print_function
import argparse
import os
import subprocess



def parse_args():
    parser = argparse.ArgumentParser("Extracts BAM segments")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                       help='BAM  input file')
    parser.add_argument('--segment_file', '-s', dest='segment_file', required=True, type=str,
                       help='File with segments descriptors')
    parser.add_argument('--output', '-o', dest='output', required=True, type=str,
                       help='Write output to file')

    return parser.parse_args()


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
    for segment in segments:
        segment_file = "{}.{}.bam".format(args.output, segment)
        pieces.append(segment_file)
        print("Extracting reads into {}".format(segment_file))
        with open(segment_file, 'w') as out_fh:
            subprocess.check_call(['samtools', 'view', '-hb', args.input, segment], stdout=out_fh)

    # merge
    print("Merging into {}".format(args.output))
    merge_cmd = ['samtools', 'merge', args.output]
    merge_cmd.extend(pieces)
    subprocess.check_call(merge_cmd)

    #cleanup
    print("Cleaning up tmp files")
    for file in pieces:
        os.remove(file)

    print("Fin.")


if __name__ == "__main__":
    main()