#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas
from matplotlib.pyplot import cm
import os
import seaborn as sns
from collections import defaultdict

SEQUENCE_IDENTITY_IDX = 13
ALIGNMENT_IDENTITY_IDX = 14

SAMPLE = "Sample"
SEQUENCE_IDENTITY = "Sequence Identity"
ALIGNMENT_IDENTITY = "Alignment Identity"

READ_LENGTH_GRANULARITY = 100

def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from margin's calcLocalPhasingCorrectness ")
    parser.add_argument('--input', '-i', dest='input_csvs', default=None, required=True, type=str, action='append',
                        help='Input read identity CSV files (can list multiple)')
    parser.add_argument('--identifier', '-I', dest='identifiers', default=None, required=False, type=str, action='append',
                        help='Input identifiers (can list multiple)')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default="output", required=False, type=str,
                       help='Figure name')

    return parser.parse_args() if args is None else parser.parse_args(args)


def main():
    args = parse_args()
    # read_lengths = defaultdict(lambda: defaultdict(lambda: 0))
    id_list = []
    all_identities = list()

    # for read length plot
    fig, ax = plt.subplots(figsize=(6, 6))
    colors = iter(cm.rainbow(np.linspace(0, 1, len(args.input_csvs))))

    for i, csv in enumerate(args.input_csvs):
        id = csv
        if args.identifiers is not None and i < len(args.identifiers):
            id = args.identifiers[i]
        id_list.append(id)

        with open(csv) as fin:
            read_lengths = defaultdict(lambda: 0)
            for j, line in enumerate(fin):
                line = line.strip().split(sep=",")
                if j == 0:
                    if line[SEQUENCE_IDENTITY_IDX] != "sequence_identity" or line[ALIGNMENT_IDENTITY_IDX] != "alignment_identity":
                        raise Exception("Unexpected identity headers: {}".format(line))
                    continue
                if len(line) in (0,2): continue
                seq_iden = float(line[SEQUENCE_IDENTITY_IDX])
                aln_iden = float(line[ALIGNMENT_IDENTITY_IDX])

                # save rows
                row = [id, seq_iden, aln_iden]
                all_identities.append(row)
                # row = ["All", seq_iden, aln_iden]
                # all_identities.append(row)

                # read length
                length = abs(int(line[4]) - int(line[3]))
                read_lengths[int(length/READ_LENGTH_GRANULARITY)] += 1

            color = next(colors)
            total_coverage = 0
            curr_len = max(read_lengths.keys())
            first = True
            while curr_len > 0:

                current_len_sequence = curr_len * read_lengths[curr_len]
                new_total_coverage = total_coverage + current_len_sequence
                ax.hlines(y=curr_len*READ_LENGTH_GRANULARITY/1000, xmin=total_coverage/1000000.0, xmax=new_total_coverage/1000000.0, color=color)

                if (first):
                    ax.hlines(curr_len*READ_LENGTH_GRANULARITY/1000, total_coverage, new_total_coverage, alpha=.5,
                               color=color, label=id)
                    ax.vlines(new_total_coverage, curr_len*READ_LENGTH_GRANULARITY/1000, (curr_len-1)*READ_LENGTH_GRANULARITY/1000, alpha=.5, color=color)
                    first = False
                else:
                    ax.hlines(curr_len*READ_LENGTH_GRANULARITY/1000, total_coverage, new_total_coverage, alpha=.5,
                               color=color)
                    ax.vlines(new_total_coverage, curr_len*READ_LENGTH_GRANULARITY/1000, (curr_len-1)*READ_LENGTH_GRANULARITY/1000, alpha=.5, color=color)

                total_coverage = new_total_coverage
                curr_len -= 1

    plt.legend()
    ax.set_xlabel("Cumulative Coverage (Gb)")
    ax.set_ylabel("Read Length (kb)")
    plt.savefig("{}.read_nx.png".format(args.figure_name))
    plt.show()
    plt.close()



    columns = [SAMPLE, SEQUENCE_IDENTITY, ALIGNMENT_IDENTITY]
    median = np.median(list(map(lambda x: x[2], all_identities)))
    mean = np.mean(list(map(lambda x: x[2], all_identities)))

    fig, ax = plt.subplots(figsize=(3*len(id_list), 6))
    df = pandas.DataFrame(all_identities, columns=columns)
    ax = sns.violinplot(x=SAMPLE, y=ALIGNMENT_IDENTITY,
                        data=df, order=id_list, linewidth=0)

    ax.annotate("Median: {:.5f}".format(median), xy=(-0.25, median+.005), color='black', size=12, fontweight='bold', font='Courier New')
    ax.axhline(median, color='black', lw=2, linestyle='dashed')
    ax.annotate("Mean:   {:.5f}".format(mean), xy=(-0.25, mean-.015), color='black', size=12, fontweight='bold', font='Courier New')
    ax.axhline(mean, color='black', lw=2, linestyle='dotted')

    ax.set_ylim(0.6, 1.025)

    plt.legend()

    plt.savefig("{}.identity.png".format(args.figure_name))
    plt.show()


if __name__ == "__main__":
    main()
