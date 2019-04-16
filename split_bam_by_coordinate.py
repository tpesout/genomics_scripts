#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import sys

import bam_stats



CHR = "c"
START = "s"
END = "e"
DESC = "d"

def parse_args():
    parser = argparse.ArgumentParser("Split BAM by region")
    parser.add_argument('--input_bam_glob', '-i', dest='input_bam_glob', required=True, type=str,
                       help='Glob matching input BAMs (will perform for all bams)')
    parser.add_argument('--coordinate_tsv', '-c', dest='coordinate_tsv', required=True, type=str,
                       help='Coordinates for splitting ($CHROM\\t$START\\t$END)')
    parser.add_argument('--output_location', '-o', dest='output_location', default=".", type=str,
                       help='Location where output files are put')
    parser.add_argument('--description_column', '-d', dest='description_column', default=None, type=int,
                       help='0-based index of description field in TSV (not required)')
    parser.add_argument('--produce_bam_stats', '-s', dest='produce_bam_stats', default=False, action='store_true',
                       help='produce bam_stats output afterwards')

    parser.add_argument('--flag_by_depth', dest='min_depth_threshold', default=None, action='store', type=int,
                       help='Will flag output if mean or median depth exceeds this threshold')

    return parser.parse_args()


def get_output_filename(input_file_location, coordinates):

    input_file_name = os.path.basename(input_file_location)
    input_file_parts = input_file_name.split(".")
    coord_description = coordinates[CHR] if len(coordinates) < 3 else "{}_{}-{}".format(
        coordinates[CHR], coordinates[START], coordinates[END])
    output_file_name = "{}.{}".format(".".join(input_file_parts[0:-1]), coord_description)
    if DESC in coordinates:
        output_file_name += "." + coordinates[DESC]
    output_file_name += ".bam"
    return output_file_name


def main():
    args = parse_args()
    assert False not in [len(args.input_bam_glob) > 0, os.path.isfile(args.coordinate_tsv), os.path.isdir(args.output_location)]

    coords = list()
    with open(args.coordinate_tsv) as tsv_in:
        for line in tsv_in:
            # skip
            if line.startswith("#"): continue
            if len(line.strip()) == 0: continue

            # get line parts
            line = line.split("\t")

            # save relevant data
            coord = {CHR:line[0]}
            if len(line) >= 3:
                coord[START] = int(line[1])
                coord[END] = int(line[2])
            if args.description_column is not None:
                coord[DESC] = "_".join(line[args.description_column].split())
            coords.append(coord)

    for file in glob.glob(args.input_bam_glob):
        for coord in coords:
            out_filename = get_output_filename(file, coord)
            out_location = os.path.join(args.output_location, out_filename)
            region = coord[CHR] if len(coord) < 3 else "{}:{}-{}".format(coord[CHR], coord[START], coord[END])

            # log what's happening
            print("{}:\n\tloc:  {}\n\tdesc: {}\n\tout:  {}".format(
                file, region, coord[DESC] if DESC in coord else '', out_location), file=sys.stderr)

            # get the region and index
            samtools_args = ['samtools', 'view', '-hb', file, region]
            with open(out_location, 'w') as output:
                subprocess.check_call(samtools_args, stdout=output)
            subprocess.check_call(['samtools', 'index', out_location])

            # maybe make stats
            if args.produce_bam_stats:
                stats_filename = "{}.stats.txt".format(out_filename)
                stats_location = os.path.join(args.output_location, stats_filename)
                bam_stats_args = ['-i', out_location, '-g', '-l', '-d', '-v', '-o', stats_location]
                if len(coord) >= 3:
                    bam_stats_args.extend(['-r', '{}-{}'.format(coord[START], coord[END])])

                # get summaries
                generic_summaries, _, depth_summaries = bam_stats.main(bam_stats_args)

                # are we flagging output based on depth?
                if args.min_depth_threshold is not None:

                    # find depth
                    generic_summary = generic_summaries[out_location][bam_stats.GENOME_KEY]
                    depth_summary = depth_summaries[out_location][bam_stats.GENOME_KEY]
                    depth_summary_value = int(max(depth_summary[bam_stats.D_MED], depth_summary[bam_stats.D_AVG]))
                    print("\tsection_depth:{}".format(depth_summary_value), file=sys.stderr)

                    # handle min depth
                    if args.min_depth_threshold <= depth_summary_value:
                        print("\tFLAGGED!", file=sys.stderr)
                        # also try to get stats for only primary reads
                        bam_stats_args_prim = ['-i', out_location, '-d', '-V', '--filter_secondary',
                                               '--filter_supplementary', '-r', '{}-{}'.format(coord[START], coord[END])]
                        _, _, prim_depth_summaries = bam_stats.main(bam_stats_args_prim)
                        prim_depth_summary = prim_depth_summaries[out_location][bam_stats.GENOME_KEY]
                        prim_depth_summary_value = int(max(prim_depth_summary[bam_stats.D_MED],
                                                           prim_depth_summary[bam_stats.D_AVG]))

                        # document depths in flagged filename
                        flag_filename = "FLAGGED.DEPTH_{:04d}_p{:04d}.{}.stats.txt".format(
                            depth_summary_value, prim_depth_summary_value, out_filename)
                        with open(os.path.join(args.output_location, flag_filename), 'w') as flag_out:
                            print("####################################", file=flag_out)
                            print("bam_file:{}".format(out_filename), file=flag_out)
                            print("bam_stats_file:{}".format(stats_filename), file=flag_out)
                            print("depth_threshold:{}".format(args.min_depth_threshold), file=flag_out)
                            print("median_depth:{}".format(depth_summary[bam_stats.D_MED]), file=flag_out)
                            print("mean_depth:{}".format(depth_summary[bam_stats.D_AVG]), file=flag_out)
                            print("total_reads:{}".format(generic_summary[bam_stats.B_READ_COUNT]), file=flag_out)
                            print("filtered_reads:{}".format(
                                generic_summary[bam_stats.B_SECONDARY_COUNT] + generic_summary[bam_stats.B_SUPPLEMENTARY_COUNT]),
                                file=flag_out)
                            print("primary_median_depth:{}".format(prim_depth_summary[bam_stats.D_MED]), file=flag_out)
                            print("primary_mean_depth:{}".format(prim_depth_summary[bam_stats.D_AVG]), file=flag_out)

    print("Fin.", file=sys.stderr)






if __name__ == "__main__":
    main()