#!/usr/bin/env python3
from __future__ import print_function
import argparse
import os
import subprocess



def parse_args():
    parser = argparse.ArgumentParser("Extracts FASTA/FASTQ segments from another fasta")
    parser.add_argument('--input', '-i', dest='input', required=True, type=str,
                       help='FASTA/FASTQ input file')
    parser.add_argument('--segment_file', '-s', dest='segment_file', required=True, type=str,
                       help='File with segments descriptors')
    parser.add_argument('--output', '-o', dest='output', required=True, type=str,
                       help='Write output to file')
    parser.add_argument('--invert', '-v', action='store_true', dest='invert', default=False,
                        help="Keep segments NOT in segment list")

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

    if args.invert:
        inverted_filename = "{}.inverted_segments.txt".format(args.output)
        print("Inverting segments.  Writing to {}".format(inverted_filename))
        inverted_segments = []
        with open(inverted_filename, 'w') as invert_out, open(args.input, 'r') as fa_in:
            for line in fa_in:
                if line.startswith(">"):
                    seg = line[1:].split()[0]
                    if seg not in segments:
                        inverted_segments.append(seg)
                        print(seg, file=invert_out)
        print("Found {} segments (inverted)".format(len(inverted_segments)))
        segments = inverted_segments


    # ensure faidx
    if not os.path.isfile("{}.fai".format(args.input)):
        print("Generating index on {}".format(args.input))
        subprocess.check_call(['samtools', 'faidx', args.input])

    # get pieces
    print("Reading from {}".format(args.input))
    with open(args.output, 'w') as out_fh:
        for segment in segments:
            print("Extracting {}".format(segment))
            subprocess.check_call(['samtools', 'faidx', args.input, segment], stdout=out_fh)

    print("Fin.")


if __name__ == "__main__":
    main()