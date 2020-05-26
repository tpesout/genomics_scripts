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
from collections import defaultdict
import shutil

import subprocess

is_secondary = lambda x: int(x) & 256 > 0
is_supplementary = lambda x: int(x) & 2048 > 0

def parse_args(args = None):
    parser = argparse.ArgumentParser("Separates a shasta assembly into individual contigs")
    parser.add_argument('--alignment', '-l', dest='alignment', required=True, type=str,
                        help='Assembly aligned to true reference (SAM format required)')
    parser.add_argument('--reference', '-r', dest='reference', required=True, type=str,
                        help='True reference fa')
    parser.add_argument('--assembly', '-s', dest='assembly', required=True, type=str,
                        help='Assembly fa')
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True, type=str,
                        help='Prefix for output files.')
    parser.add_argument("--ng", dest='ng', default=75, type=int,
                        help='ng parameter for juipiter')
    parser.add_argument("--min_bundle_size", dest='min_bundle_size', default=100000, type=int,
                        help='minBundleSize parameter for jupiter')
    parser.add_argument("--no_secondary", dest='include_secondary', default=True, action='store_false',
                        help="do not include secondary alignments")
    parser.add_argument("--no_supplementary", dest='include_supplementary', default=True, action='store_false',
                        help="do not include supplementary alignments")

    return parser.parse_args() if args is None else parser.parse_args(args)



def main():
    args = parse_args()

    # sanity check
    missing = list(filter(lambda x: not os.path.isfile(x), [args.alignment, args.reference, args.assembly]))
    if len(missing) > 0:
        raise Exception("Missing files: {}".format(missing))

    # ensure output
    if "/" in args.output_prefix:
        parts = args.output_prefix.split("/")
        curr_out_dir = "."
        for part in parts[:-1]:
            curr_out_dir = "/".join([curr_out_dir, part])
            if not os.path.isdir(curr_out_dir): os.mkdir(curr_out_dir)

    # build maps:
    header_lines = list()
    ref_contig_to_assembly_segment = defaultdict(set)
    assembly_segment_to_ref_contig = defaultdict(set)

    # get all contigs and segments
    print("Getting contigs and ref sequences from {}".format(args.alignment))
    with open(args.alignment, 'r') as align_in:
        for line in align_in:
            if line.startswith("@"):
                header_lines.append(line)
                continue
            line = line.split()
            ref_contig = line[2].strip()
            assembly_segments = line[0].strip()
            secondary = is_secondary(line[1])
            supplementary = is_supplementary(line[1])
            # no seconds or sups if configured
            if secondary and not args.include_secondary: continue
            if supplementary and not args.include_supplementary: continue
            # skip reads less than MBS
            if len(line[9]) < args.min_bundle_size: continue
            # otherwise save it
            ref_contig_to_assembly_segment[ref_contig].add(assembly_segments)
            assembly_segment_to_ref_contig[assembly_segments].add(ref_contig)

    # write to files split by ref contig
    contig_handles = dict()
    contig_filenames = dict()
    print("Splitting SAM into files: {}".format(args.output_prefix))
    try:
        # open files and write header
        for ref_contig in ref_contig_to_assembly_segment.keys():
            filename = "{}.{}.sam".format(args.output_prefix, ref_contig)
            contig_filenames[ref_contig] = filename
            # skip files that are already created
            if os.path.isfile(filename):
                contig_handles[ref_contig] = None
                continue
            # write the header
            contig_handles[ref_contig] = open(filename, 'w')
            for header in header_lines:
                contig_handles[ref_contig].write(header)

        # write lines to handles
        with open(args.alignment, 'r') as align_in:
            for line in align_in:
                # skip header lines
                if line.startswith("@"):
                    continue
                # for every contig that this segment aligns to, write this line to that contig's bam
                line_parts = line.split()
                secondary = is_secondary(line_parts[1])
                supplementary = is_supplementary(line_parts[1])
                # no seconds or sups if configured
                if secondary and not args.include_secondary: continue
                if supplementary and not args.include_supplementary: continue
                # skip reads less than MBS
                if len(line_parts[9]) < args.min_bundle_size: continue

                assembly_segment = line_parts[0]
                for rc in assembly_segment_to_ref_contig[assembly_segment]:
                    # skip it if the file already exists
                    if contig_handles[rc] is None:
                        continue
                    contig_handles[rc].write(line)

    finally:
        # close everything
        for ref_contig in ref_contig_to_assembly_segment.keys():
            #..unless it was already present
            if ref_contig in contig_handles and contig_handles[ref_contig] is not None:
                contig_handles[ref_contig].close()

    # run jupiter on them
    print("\nRunning jupiter")
    for ref_contig in contig_filenames.keys():
        contig_filename = contig_filenames[ref_contig]
        jupiter_dir ="{}.jupiter".format(contig_filename)
        if not os.path.isdir(jupiter_dir): os.mkdir(jupiter_dir)
        jupiter_cmd = [
            "jupiter",
            "name={}".format(ref_contig),
            "fa={}".format(os.path.abspath(args.assembly)),
            "ref={}".format(os.path.abspath(args.reference)),
            "sam={}".format(os.path.abspath(contig_filename)),
            "ng={}".format(args.ng),
            "minBundleSize={}".format(args.min_bundle_size)
        ]
        print("\nRUNNING:\n\tin:{}\n\t{}".format(jupiter_dir, " ".join(jupiter_cmd)))
        subprocess.run(jupiter_cmd, stdout=sys.stdout, stderr=sys.stderr, cwd=jupiter_dir)
        image_loc = "{}.jupiter/{}.png".format(contig_filename, ref_contig)
        if os.path.isfile(image_loc):
            shutil.copy(image_loc, "{}.jupiter.{}.png".format(contig_filename, ref_contig))

    print("\nFin.")






if __name__ == "__main__":
    main()