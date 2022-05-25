#!/usr/bin/env python

import sys
import collections
import numpy as np


def log(msg):
    print(msg, file=sys.stderr)


def save_phaseset(phasesets, chr, pos, ps_id):

    phasesets[chr][ps_id][0] = min(phasesets[chr][ps_id][0], pos)
    phasesets[chr][ps_id][1] = max(phasesets[chr][ps_id][1], pos)
    phasesets[chr][ps_id][2] += 1


def main():

    phasesets = collections.defaultdict(
        lambda : collections.defaultdict(
            lambda: [sys.maxsize, 0, 0]
        )
    )
    phased_variants = 0

    for i, line in enumerate(sys.stdin):

        try:
            # header
            if line.startswith("#"): continue

            # get elements
            parts = line.split()
            if len(parts) < 10:
                log("Line {} was malformed: {}".format(i, line))
                sys.exit(1)

            # get info we care about
            chr = parts[0]
            pos = int(parts[1]) - 1 # vcf is 1-based (ew)
            # find where (or if) ps is annotated
            format_parts = parts[8].split(":")
            ps_idx = None
            for j, f in enumerate(format_parts):
                if f == "PS":
                    ps_idx=j
                    break
            # it might not be defined
            if ps_idx is None:
                continue
            ps_id = parts[9].split(":")[ps_idx]

            # save data
            phased_variants += 1
            save_phaseset(phasesets, chr, pos, ps_id)

        except Exception as e:
            log("Got error reading line {}: {}".format(i, line))
            raise e

    # loggit
    if len(phasesets) == 0:
        log("Found no phasesets!")
        sys.exit(0)

    # get phasesets
    phaseset_lengths = list()
    print("#chr\tstart_pos\tend_pos\tphaseset_id\tcontained_variants")
    for chr in sorted(list(phasesets.keys())):
        for ps_id in sorted(list(phasesets[chr].keys()), key=lambda x: phasesets[chr][x][0]):
            ps_info = phasesets[chr][ps_id]
            print("{}\t{}\t{}\t{}\t{}".format(chr, ps_info[0], ps_info[1], ps_id, ps_info[2]))
            phaseset_lengths.append(ps_info[1] - ps_info[0])

    # log info
    log("Got {} phasesets from {} phased variants with:".format(len(phaseset_lengths), phased_variants))
    log("\tMean:    {}".format(int(np.mean(phaseset_lengths))))
    log("\tMedian:  {}".format(int(np.median(phaseset_lengths))))
    log("\tMax:     {}".format(max(phaseset_lengths)))
    half_total = sum(phaseset_lengths) / 2
    for x in sorted(phaseset_lengths, reverse=True):
        half_total -= x
        if half_total <= 0:
            log("\tN50:     {}".format(x))
            break


if __name__ == "__main__":
    main()