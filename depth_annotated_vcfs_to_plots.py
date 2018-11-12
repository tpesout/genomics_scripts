#!/usr/bin/env python
from __future__ import print_function
import argparse
import gzip
import sys
import os
import math

TP='tp'
FP='fp'
FN='fn'

def parse_args():
    parser = argparse.ArgumentParser("Produce a BED based on VCF calls")
    parser.add_argument('--tp_vcf', '-tp', dest='tp_vcf', required=True, type=str,
                       help='true positive VCF')
    parser.add_argument('--fp_vcf', '-fp', dest='fp_vcf', required=True, type=str,
                       help='false positive VCF')
    parser.add_argument('--fn_vcf', '-fn', dest='fn_vcf', required=True, type=str,
                       help='false negative VCF')
    parser.add_argument('--max_depth', '-d', dest='max_depth', required=False, default=64, type=int,
                       help='Maximum VCF (all values here and above are consolidated to this value)')
    parser.add_argument('--depth_tag', '-t', dest='depth_tag', default="DP", type=str,
                       help='Tag to use when annotating depth')
    parser.add_argument('--title', dest='title', default=None, type=str,
                       help='plot title')

    return parser.parse_args()


def save_depths_and_counts(args, input_vcf, depth_values, type_key):
    # prep
    depth_tag = args.depth_tag
    max_depth = args.max_depth
    input_open_fcn = open if not input_vcf.endswith("gz") else gzip.open
    def get_depth_tag_index(tag_parts):
        for i, t in enumerate(tag_parts):
            if t == depth_tag: return i
        return None

    # read it all
    with input_open_fcn(input_vcf) as input:
        linenr = -1
        for line in input:
            linenr += 1
            if line.startswith("#"): continue

            # get parts
            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                raise Exception( "Line {} in {} appears to be malformed: {}".format(linenr, input_vcf, line))

            # get data
            tags_parts = parts[8].split(":")
            values_parts = parts[9].split(":")
            dti = get_depth_tag_index(tags_parts)
            if dti is None:
                print("Call from {} does not have depth tag {}: {}".format(input_vcf, depth_tag, "\\t".join(parts)))
                continue
            depth = int(values_parts[dti])
            depth = min(depth, max_depth)

            # save the values
            depth_values[depth][type_key] += 1


def plot_it(depths, precisions, recalls, fmeasures, total_calls=None, title=None):
    import matplotlib.pyplot as plt

    max_f = max(fmeasures)
    max_f_depths = list(filter(lambda i: fmeasures[i[0]] == max_f, enumerate(depths)))

    plt.plot(depths, precisions, color='blue', label='Precision')
    plt.plot(depths, recalls, color='red', label='Recall')
    plt.plot(depths, fmeasures, color='purple', label='F-Measure')
    for d in max_f_depths:
        plt.axvline(x=d[1], color='purple')
    plt.legend(loc=0)

    if total_calls is not None:
        sizes = list(map(lambda x: 4 ** math.log10(1 if x == 0 else x), total_calls))

        plt.scatter(depths, fmeasures, s=sizes, alpha=.25, color='purple')
        plt.plot(depths, precisions, color='blue')
        plt.plot(depths, recalls, color='red')
        plt.plot(depths, fmeasures, color='purple')

    if title is not None: plt.title("MarginPhase Nanopore: " + title)
    plt.show()


def plot_it_all(depths, d_precisions, d_recalls, d_fmeasures, precisions, recalls, fmeasures,
                d_total_calls=None, title=None):
    import matplotlib.pyplot as plt

    # at depths
    plt.plot(depths, d_precisions, linestyle=":", color='blue', label='Precision (At Depth)')
    plt.plot(depths, d_recalls, linestyle=":", color='red', label='Recall (At Depth)')
    plt.plot(depths, d_fmeasures, linestyle=":", color='purple', label='F-Measure (At Depth)')

    # at or above
    plt.plot(depths, precisions, color='blue', label='Precision')
    plt.plot(depths, recalls, color='red', label='Recall')
    plt.plot(depths, fmeasures, color='purple', label='F-Measure')

    # depths
    if d_total_calls is not None:
        sizes = list(map(lambda x: 4 ** math.log10(1 if x == 0 else x), d_total_calls))
        plt.scatter(depths, d_fmeasures, s=sizes, alpha=.2, color='grey', label="Total Calls at Depth")

    plt.legend(loc=0)

    # best depth (at or above)
    max_f = max(fmeasures)
    max_f_depths = list(filter(lambda i: fmeasures[i[0]] == max_f, enumerate(depths)))
    for d in max_f_depths:
        plt.axvline(x=d[1], color='purple')

    if title is not None: plt.title(title)
    plt.show()



def main():
    args = parse_args()

    # get input parameters
    if not os.path.isfile(args.tp_vcf): raise Exception("TP VCF {} does not exist".format(args.tp_vcf))
    if not os.path.isfile(args.fp_vcf): raise Exception("FP VCF {} does not exist".format(args.fp_vcf))
    if not os.path.isfile(args.fn_vcf): raise Exception("FN VCF {} does not exist".format(args.fn_vcf))

    # init data
    depth_values = {d:{TP:0,FP:0,FN:0} for d in range(args.max_depth + 1)}

    save_depths_and_counts(args, args.tp_vcf, depth_values, TP)
    save_depths_and_counts(args, args.fp_vcf, depth_values, FP)
    save_depths_and_counts(args, args.fn_vcf, depth_values, FN)

    # functions
    precision_fcn = lambda tp, fp, fn: 0.0 if (tp + fp) == 0 else 1.0 * tp / (tp + fp)
    recall_fcn = lambda tp, fp, fn: 0.0 if (tp + fn) == 0 else 1.0 * tp / (tp + fn)
    fmeasure_fcn = lambda p, r: 0.0 if p+r == 0 else (2.0 * p * r / (p + r))

    # plotting
    depths = list()
    d_precisions = list()
    d_recalls = list()
    d_fmeasures = list()
    d_totals = list()
    t_precisions = list()
    t_recalls = list()
    t_fmeasures = list()
    t_totals = list()

    # data we track
    curr_depth = args.max_depth
    tp_total = 0
    fp_total = 0
    fn_total = 0
    total = 0
    while curr_depth >= 0:
        dtp = depth_values[curr_depth][TP]
        dfp = depth_values[curr_depth][FP]
        dfn = depth_values[curr_depth][FN]
        dtotal = dtp + dfp + dfn

        dp = precision_fcn(dtp, dfp, dfn)
        dr = recall_fcn(dtp, dfp, dfn)
        df = fmeasure_fcn(dp, dr)

        tp_total += dtp
        fp_total += dfp
        fn_total += dfn
        total += dtotal

        total_p = precision_fcn(tp_total, fp_total, fn_total)
        total_r = recall_fcn(tp_total, fp_total, fn_total)
        total_f = fmeasure_fcn(total_p, total_r)

        print(("{:3d}:\tprecision: {:0.4f}\trecall: {:0.4f}\tf-measure: {:0.4f}\tcalls_at_depth:{:8d}\t\t"
              "Total Above Depth:\tprecision: {:0.4f}\trecall: {:0.4f}\tf-measure: {:0.4f}\ttotal_calls:{:8d}").format(
            curr_depth, dp, dr, df, dtotal, total_p, total_r, total_f, total))

        depths.append(curr_depth)
        d_precisions.append(dp)
        d_recalls.append(dr)
        d_fmeasures.append(df)
        d_totals.append(dtotal)
        t_precisions.append(total_p)
        t_recalls.append(total_r)
        t_fmeasures.append(total_f)
        t_totals.append(total)

        curr_depth -= 1

    # plots
    plot_it(depths, d_precisions, d_recalls, d_fmeasures, total_calls=d_totals,
            title=None if not args.title else (args.title + ": At Depth"))
    plot_it(depths, t_precisions, t_recalls, t_fmeasures,
            title=None if not args.title else (args.title + ": At or Above Depth"))
    plot_it_all(depths, d_precisions, d_recalls, d_fmeasures, t_precisions, t_recalls, t_fmeasures,
                d_total_calls=d_totals, title=args.title)


    # fin
    print("Fin.", file=sys.stderr)


if __name__ == "__main__":
    main()