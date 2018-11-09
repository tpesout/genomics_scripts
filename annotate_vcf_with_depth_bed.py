#!/usr/bin/env python
from __future__ import print_function
import argparse
import gzip
import sys
import os

START = 0
STOP = 1
DEPTH = 2

def parse_args():
    parser = argparse.ArgumentParser("Produce a BED based on VCF calls")
    parser.add_argument('--input_vcf', '-v', dest='input_vcf', required=True, type=str,
                       help='VCF to apply depths to')
    parser.add_argument('--input_bed', '-b', dest='input_bed', required=True, type=str,
                       help='Depth BED to annotate calls with')
    parser.add_argument('--output_vcf', '-o', dest='output_vcf', required=False, default=None, type=str,
                       help='Write output to file (otherwise stdout)')
    parser.add_argument('--depth_tag', '-t', dest='depth_tag', default="DP", type=str,
                       help='Tag to use when annotating depth')

    return parser.parse_args()


def main():
    args = parse_args()

    # get input parameters
    input_vcf = args.input_vcf
    if not os.path.isfile(input_vcf):
        raise Exception("Input VCF {} does not exist".format(input_vcf))
    input_bed = args.input_bed
    if not os.path.isfile(input_bed):
        raise Exception("Input BED {} does not exist".format(input_bed))
    depth_tag = args.depth_tag

    # assumes everything is sorted
    depths = dict()
    with open(input_bed) as input:
        # validation tracking
        linenr = 0
        last_chr = None
        last_stop = None

        # examine all bed lines
        for line in input:
            parts = line.split()
            if len(parts) < 4: raise Exception("Line {} in {} was malformed: {}".format(linenr, input_bed, line))

            # get data we care about
            chr = parts[0]
            start = int(parts[1])
            stop = int(parts[2])
            depth = int(parts[3])

            # sanity checks
            if chr not in depths:
                depths[chr] = list()
                last_chr = chr
                last_stop = 0
            if chr != last_chr: raise Exception("Encountered previously found contig {} at line {} (current {})".format(
                chr, linenr, list(depths.keys())))
            if stop <= last_stop or start < last_stop: raise Exception("BED appears to be unsorted at line {}".format(
                linenr))

            # save
            depths[chr].append([start, stop, depth])

            # iterate
            last_stop = stop
            linenr += 1
    print("Found {} contigs and {} depths in {}".format(len(depths), sum(map(len, depths.values())), input_vcf), file=sys.stderr)

    # now open input and output
    output = None
    try:
        # get output file
        if args.output_vcf is not None:
            output = open(args.output_vcf, 'w')
        else:
            output = sys.stdout

        # prep
        found_tag_format = False
        def get_depth_tag_index(tag_parts):
            for i, t in enumerate(tag_parts):
                if t == depth_tag: return i
            return None
        missing_contigs = set()

        # prep for iteration
        curr_chr = None
        curr_bed = None
        curr_bed_itor = None

        # read input, write to output
        with open(args.input_vcf) as input:
            linenr = -1
            for line in input:
                linenr += 1
                if line.startswith("#"):
                    if line.startswith("##FORMAT=<ID={}".format(depth_tag)):
                        found_tag_format = True
                    if not line.startswith("##") and not found_tag_format:
                        output.write("##FORMAT=<ID={},Number=1,Type=Integer,Description=\"Read depth (custom)\">\n".format(depth_tag))
                    output.write(line)
                    continue

                # get parts
                parts = line.rstrip().split("\t")
                if len(parts) < 10: raise Exception("Line {} in {} appears to be malformed: {}".format(linenr, input_vcf, line))

                # get data
                chr = parts[0]
                pos = int(parts[1])
                tags = parts[8]

                # maybe init depth iterator
                if chr not in depths:
                    if chr not in missing_contigs:
                        print("Contig {} missing from bed!".format(chr), file=sys.stderr)
                        missing_contigs.add(chr)
                    continue
                if curr_chr != chr:
                    curr_bed_itor = iter(depths[chr])
                    curr_bed = next(curr_bed_itor, None)
                while (curr_bed is not None and curr_bed[STOP] <= pos):
                    curr_bed = next(curr_bed_itor, None)

                # get depth
                if curr_bed is None:
                    depth = '0'
                elif pos >= curr_bed[START] and pos < curr_bed[STOP]:
                    depth = str(curr_bed[DEPTH])
                else:
                    raise Exception("Unexpected depth state:\n\t{}\n\tcurr_bed:{}\n".format("\\t".join(parts), curr_bed))

                # save in appropriate spot
                tags_parts = tags.split(":")
                dti = get_depth_tag_index(tags_parts)
                if dti is None:
                    tags_parts.append(depth_tag)
                    for i in range(9,len(parts)):
                        values_parts = parts[i].split(":")
                        values_parts.append(depth)
                        parts[i] = ":".join(values_parts)
                else:
                    for i in range(9,len(parts)):
                        values_parts = parts[i].split(":")
                        values_parts[dti] = depth
                        parts[i] = ":".join(values_parts)

                # save the values
                parts[8] = ":".join(tags_parts)
                output.write("{}\n".format("\t".join(parts)))

    # make sure output is closed (if appropriate)
    finally:
        if output is not None and args.output_vcf is not None: output.close()

    # fin
    print("Fin.", file=sys.stderr)


if __name__ == "__main__":
    main()