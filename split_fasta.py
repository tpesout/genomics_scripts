#!/usr/bin/env python

import argparse


def parse_args():
    parser = argparse.ArgumentParser("Produces BED file of bam where read depth is at or above a threshold")
    parser.add_argument('--input_fasta', '-i', dest='input_fasta', required=True, type=str,
                        help='Input fasta file')
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', default="output_", type=str,
                       help="Output will be of form ${OUTPUT_PREFIX}${IDX}${OUTPUT_SUFFIX}")
    parser.add_argument('--output_suffix', '-O', dest='output_suffix', default=".fasta", type=str,
                       help="Output will be of form ${OUTPUT_PREFIX}${IDX}${OUTPUT_SUFFIX}")

    parser.add_argument('--minimum_lines', '-l', dest='minimum_lines', required=True, type=int,
                        help='Minimum number of fasta lines to keep per file')
    parser.add_argument('--index_format', '-d', dest='index_format', default="%03d",
                        help='Format string for chunk index')

    args = parser.parse_args()
    return args


def main():
    # prep
    args = parse_args()

    with open(args.input_fasta, 'r') as infile:
        current_line = 0
        current_idx = 0
        current_outfile = None
        try:
            for line in infile:
                if current_outfile is None or (line.startswith(">") and current_line >= args.minimum_lines):
                    if current_outfile is not None:
                        current_outfile.close()
                    outfile_idx = args.index_format % current_idx
                    outfilename = args.output_prefix + outfile_idx + args.output_suffix
                    print("Writing to file {} after writing {} lines".format(outfilename, current_line))
                    current_outfile = open(outfilename, 'w')
                    current_line = 0
                    current_idx += 1
                current_outfile.write(line)
                current_line += 1

        finally:
            if current_outfile is not None: current_outfile.close()


if __name__ == "__main__":
    main()
