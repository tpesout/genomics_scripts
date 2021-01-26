#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import matplotlib.pyplot as plt
import pysam
import collections
import numpy as np
import os

GENE_GTF_CHROM_IDX=0
GENE_GTF_TYPE_IDX=2
GENE_GTF_START_POS_IDX=3
GENE_GTF_END_POS_IDX=4
GENE_GTF_INFO_IDX=8

TYPE_GENE="g"
TYPE_HICONF="h"
TYPE_SWITCH="s"
TYPE_PHASESET="p"

WHOLLY="wholly_covered"
PARTIALLY="partially_covered"
NOT_COVERED="not_covered"

SNP="snp"
INDEL="indel"
TP="tp"
FP="fp"
FN="fn"

class BedRecord:
    def __init__(self, chrom, start, end, type):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.type = type
        self.overlaps = list()


class GeneInfo(BedRecord):
    def __init__(self, chrom, start, end, info):
        super().__init__(chrom, start, end, TYPE_GENE)
        self.info = dict()
        for inf in info.split(";"):
            key_val = inf.strip().split()
            if len(key_val) == 0: continue
            self.info[key_val[0]] = key_val[1].strip("\"")
        self.precision = None
        self.recall = None
        self.snp_precision = None
        self.snp_recall = None
        self.indel_precision = None
        self.indel_recall = None
        self.snp_tp = 0
        self.snp_fp = 0
        self.snp_fn = 0
        self.indel_tp = 0
        self.indel_fp = 0
        self.indel_fn = 0

    def calculate_precision_recall(self, snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn):
        self.snp_tp = snp_tp
        self.snp_fp = snp_fp
        self.snp_fn = snp_fn
        self.indel_tp = indel_tp
        self.indel_fp = indel_fp
        self.indel_fn = indel_fn
        if self.get_total_variants() == 0:
            return False
        self.precision = 1.0 * (snp_tp + indel_tp) / max(1, snp_tp+snp_fp+indel_tp+indel_fp)
        self.recall = 1.0 * (snp_tp + indel_tp) / max(1, snp_tp+snp_fn+indel_tp+indel_fn)
        self.snp_precision = 1.0 * snp_tp / max(1, snp_tp+snp_fp)
        self.snp_recall = 1.0 * snp_tp / max(1, snp_tp+snp_fn)
        self.indel_precision = 1.0 * indel_tp / max(1, indel_tp+indel_fp)
        self.indel_recall = 1.0 * indel_tp / max(1, indel_tp+indel_fn)
        return True

    def get_total_variants(self):
        return sum([self.snp_tp, self.snp_fp, self.snp_fn, self.indel_tp, self.indel_fp, self.indel_fn])

    def get_gene_name(self):
        if "gene_name" in self.info:
            return "{} ({}:{}-{})".format(self.info["gene_name"], self.chrom, self.start, self.end)
        return "{}:{}-{}".format(self.chrom, self.start, self.end)


def parse_args(args = None):
    parser = argparse.ArgumentParser("Analyzes phasing over genes")
    parser.add_argument('--phased_vcf', '-v', dest='phased_vcf', default=None, required=True, type=str,
                       help='Phased VCF being analyzed')
    parser.add_argument('--gene_gtf', '-g', dest='gene_gtf', default=None, required=True, type=str,
                       help='GTF holding gene annotations')
    parser.add_argument('--whatshap_compare_switch_bed', '-w', dest='whatshap_compare_switch_bed', default=None, required=False, type=str,
                       help='switch_error.bed file produced from whatshap ')
    parser.add_argument('--high_confidence_bed', '-c', dest='high_confidence_bed', default=None, required=False, type=str,
                       help='high confidence bed from annotated truth')
    parser.add_argument('--truth_annotated_vcf', '-a', dest='truth_annotated_vcf', default=None, required=False, type=str,
                       help='Annotated truth from hap.py')
    parser.add_argument('--repair_chr_names', '-C', dest='repair_chr_names', default=False, required=False, action='store_true',
                       help='Add "chr" to contig names without (there are sometimes "1" and "chr1" mismatches for grch37)')
    parser.add_argument('--output_base', '-o', dest='output_base', default=None, required=False, type=str,
                       help='Output base file name')
    parser.add_argument('--output_base_from_input', '-O', dest='output_base_from_input', default=False, required=False, action='store_true',
                       help='Output file name will be determined from --phased_vcf parameter')
    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg, end="\n"):
    print(msg, file=sys.stderr, end=end)


def fix_chrom(chrom, args):
    return chrom if not args.repair_chr_names or chrom.startswith("chr") else "chr" + chrom

def unfix_chrom(chrom, args):
    if not args.repair_chr_names:
        return chrom
    if chrom == "chrM":
        return chrom
    return chrom if not chrom.startswith("chr") else chrom.lstrip("chr")


def get_bed_regions(args):
    genes = collections.defaultdict(lambda: list())
    switches = collections.defaultdict(lambda: list())
    high_conf_regions = collections.defaultdict(lambda: list())

    # genes
    with open(args.gene_gtf) as fin:
        for line in fin:
            # errors and invalid records
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) <= GENE_GTF_INFO_IDX:
                log("Unexpected line reading gene GTF: '{}'".format(line))
                continue
            # skip non-gene parts
            type = parts[GENE_GTF_TYPE_IDX]
            if type != 'gene': continue
            # get data
            chrom = fix_chrom(parts[GENE_GTF_CHROM_IDX], args)
            start_pos = int(parts[GENE_GTF_START_POS_IDX])
            end_pos = int(parts[GENE_GTF_END_POS_IDX])
            info = parts[GENE_GTF_INFO_IDX]
            gene = GeneInfo(chrom, start_pos, end_pos, info)
            # save it
            genes[chrom].append(gene)
    log("Found {} genes in {}".format(sum(map(len, genes.values())), args.gene_gtf))
    for gv in genes.values():
        gv.sort(key=lambda x: x.start)

    # switches
    if args.whatshap_compare_switch_bed is not None:
        with open(args.whatshap_compare_switch_bed) as fin:
            for line in fin:
                # errors and invalid records
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    log("Unexpected line reading switches: '{}'".format(line))
                    continue
                chrom = fix_chrom(parts[0], args)
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                switches[chrom].append(BedRecord(chrom, start_pos, end_pos, TYPE_SWITCH))
        log("Found {} switches in {}".format(sum(map(len, switches.values())), args.whatshap_compare_switch_bed))
        for sw in switches.values():
            sw.sort(key=lambda x: x.start)
    else:
        switches = None

    # high conf
    if args.high_confidence_bed is not None:
        with open(args.high_confidence_bed) as fin:
            for line in fin:
                # errors and invalid records
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    log("Unexpected line reading high conf bed: '{}'".format(line))
                    continue
                chrom = fix_chrom(parts[0], args)
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                high_conf_regions[chrom].append(BedRecord(chrom, start_pos, end_pos, TYPE_HICONF))
        log("Found {} high confidence regions in {}".format(sum(map(len, high_conf_regions.values())),
                                                            args.high_confidence_bed))
        for hc in high_conf_regions.values():
            hc.sort(key=lambda x: x.start)
    else:
        high_conf_regions = None

    return genes, switches, high_conf_regions


def get_phasesets_from_vcf(vcf_loc, args):
    phasesets = collections.defaultdict(lambda : list())
    vcf = pysam.VariantFile(vcf_loc)
    if len(vcf.header.samples) != 1:
        log("Got {} samples in phased VCF {}, expect only 1".format(len(vcf.header.samples), vcf_loc))
    sample = vcf.header.samples[0]
    current_phaseset = None
    current_phaseset_start = None
    current_phaseset_end = None
    for record in vcf.fetch():
        chrom = record.contig
        pos = record.pos
        phaseset = record.format.get("PS")
        # phaseset = list(record.format.iteritems())

        # skip if not phased
        if phaseset is None:
            continue
        phaseset = record.samples[0]["PS"]
        # phaseset init
        if current_phaseset is None:
            current_phaseset = phaseset
            current_phaseset_start = pos
        # phaseset record
        if phaseset != current_phaseset:
            phasesets[fix_chrom(chrom, args)].append(BedRecord(chrom, current_phaseset_start, current_phaseset_end, TYPE_PHASESET))
            current_phaseset = phaseset
            current_phaseset_start = pos
        # always update end of the current ps
        current_phaseset_end = pos

    # record final phaseset
    if current_phaseset is not None:
        phasesets[fix_chrom(chrom, args)].append(BedRecord(chrom, current_phaseset_start, current_phaseset_end, TYPE_PHASESET))

    #loggit
    log("Found {} phasesets in {}".format(sum(map(len, phasesets.values())), vcf_loc))
    return phasesets


def get_precision_recall_for_genes(genes, args):
    if (args.truth_annotated_vcf) is None:
        log("No truth annotated vcf is present")
        return

    vcf = pysam.VariantFile(args.truth_annotated_vcf)
    if len(vcf.header.samples) != 2:
        log("Got {} samples in phased VCF {}, expect 2".format(len(vcf.header.samples), args.truth_annotated_vcf))
        sys.exit(1)

    # get samples
    query_sample_idx = -1
    truth_sample_idx = -1
    for i,sample in enumerate(vcf.header.samples):
        if sample.upper() == "TRUTH":
            truth_sample_idx = i
        elif sample.upper() == "QUERY":
            query_sample_idx = i
    if query_sample_idx == -1 or truth_sample_idx == -1:
        log("Expected TRUTH and QUERY samples in {}, got {}".format(args.truth_annotated_vcf, vcf.header.samples))
        sys.exit(1)

    # for sanity checks
    update_zero_count = 0
    update_one_count = 0
    update_multiple_count = 0
    high_conf_variants = 0
    not_high_conf_variants = 0
    missing_tag_variants = 0
    genes_with_variants = 0
    genes_without_variants = 0

    # get variants for all genes
    for chrom in sorted(list(genes.keys())):
        for gene in genes[chrom]:
            snp_tp = 0
            snp_fp = 0
            snp_fn = 0
            indel_tp = 0
            indel_fp = 0
            indel_fn = 0


            # we only track in hiconf
            higconf_beds = list(filter(lambda x: x.type == TYPE_HICONF, gene.overlaps))
            try:
                vcf_records = vcf.fetch(unfix_chrom(chrom, args), gene.start, gene.end)
            except:
                log("\nFAILED TO GET RECORDS FOR {}:{}-{}\n".format(unfix_chrom(chrom, args), gene.start, gene.end))
                continue
            for record in vcf_records:

                # only get hiconf
                pos = record.pos
                hiconf = False
                for hcb in higconf_beds:
                    if hcb.start <= pos < hcb.end:
                        hiconf = True
                        break
                if not hiconf:
                    not_high_conf_variants += 1
                    continue
                high_conf_variants += 1
                # must have BD and BVT tags
                if record.format.get("BD") is None or record.format.get("BVT") is None:
                    missing_tag_variants += 1
                    continue

                truth_bd = record.samples[truth_sample_idx]["BD"]
                truth_bvt = record.samples[truth_sample_idx]["BVT"]
                query_bd = record.samples[query_sample_idx]["BD"]
                query_bvt = record.samples[query_sample_idx]["BVT"]

                pre_total = sum([snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn])
                if truth_bd == "FN":
                    if truth_bvt == "INDEL":
                        indel_fn += 1
                    elif truth_bvt == "SNP":
                        snp_fn += 1
                    else:
                        pass
                if query_bd == "FP":
                    if query_bvt == "INDEL":
                        indel_fp += 1
                    elif query_bvt == "SNP":
                        snp_fp += 1
                    else:
                        pass
                if query_bd == "TP" or truth_bd == "TP":
                    if query_bvt == "INDEL":
                        indel_tp += 1
                    elif query_bvt == "SNP":
                        snp_tp += 1
                    elif truth_bvt == "INDEL":
                        indel_tp += 1
                    elif truth_bvt == "SNP":
                        snp_tp += 1

                # qc
                post_total = sum([snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn])
                diff = post_total - pre_total
                if diff == 0:
                    update_zero_count += 1
                elif diff == 1:
                    update_one_count += 1
                else:
                    update_multiple_count += 1

            gene.calculate_precision_recall(snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn)
            if gene.get_total_variants() == 0:
                genes_without_variants += 1
            else:
                genes_with_variants += 1

    log("Analyzed Genes:")
    log("\tVariant classification counts:   0:{} 1:{} 2+:{}".format(update_zero_count, update_one_count, update_multiple_count))
    log("\tVariants in HiConf:              {}".format(high_conf_variants))
    log("\tVariants not in HiConf:          {}".format(not_high_conf_variants))
    log("\tVariants missing tags:           {}".format(missing_tag_variants))
    log("\tGenes with variants:             {}".format(genes_with_variants))
    log("\tGenes without variants:          {}".format(genes_without_variants))


def search_bed_list_for_overlapping_beds(start, end, bed_list):
    overlaps = list()

    # binary search to find starting bed pos
    start_idx = 0
    end_idx = len(bed_list)
    curr_idx = None
    while start_idx < end_idx:
        curr_idx = int(start_idx + end_idx) // 2
        curr_bed = bed_list[curr_idx]
        if curr_bed.start <= start < curr_bed.end:
            break
        if curr_bed.start > start:
            end_idx = curr_idx
        elif curr_bed.start < start:
            start_idx = curr_idx + 1
        else:
            assert(False)
        curr_idx = None

    # return all overlapping beds
    if curr_idx is None:
        return overlaps
    while curr_idx < len(bed_list) and bed_list[curr_idx].start < end:
        overlaps.append(bed_list[curr_idx])
        curr_idx += 1
    return overlaps


def get_genes_stratified_by_bed(genes, bed_regions, bed_desc):

    genes_wholly_covered = collections.defaultdict(lambda: list())
    genes_partially_covered = collections.defaultdict(lambda: list())
    genes_not_covered = collections.defaultdict(lambda: list())

    if bed_regions is None:
        log("Missing bed for {}".format(bed_desc))
        return genes_wholly_covered, genes_partially_covered, genes
    log("Calculating gene overlaps for {}".format(bed_desc))

    for chr in set(genes.keys()).intersection(set(bed_regions.keys())):
        chr_genes = genes[chr]
        chr_region = bed_regions[chr]
        for gene in chr_genes:
            hiconf_overlaps = search_bed_list_for_overlapping_beds(gene.start, gene.end, chr_region)
            gene.overlaps.extend(hiconf_overlaps)
            if len(hiconf_overlaps) == 0:
                genes_not_covered[chr].append(gene)
            elif len(hiconf_overlaps) == 1:
                if hiconf_overlaps[0].start <= gene.start and hiconf_overlaps[0].end >= gene.end:
                    genes_wholly_covered[chr].append(gene)
                else:
                    genes_partially_covered[chr].append(gene)
            else:
                genes_partially_covered[chr].append(gene)
        log("\tFor {}, of {} total genes got {} wholly covered, {} partially covered, and {} not covered".format(
            chr, len(chr_genes), len(genes_wholly_covered[chr]), len(genes_partially_covered[chr]), len(genes_not_covered[chr])))

    return genes_wholly_covered, genes_partially_covered, genes_not_covered


def get_coverage_state(wholly, partially, nott, chrom, gene):
    if gene in wholly[chrom]:
        return WHOLLY
    if gene in partially[chrom]:
        return PARTIALLY
    if gene in nott[chrom]:
        return NOT_COVERED
    return None


def plot_coverage_donut(coverage_breakdown, output_base):

    # data we plot
    wholeH = 0
    wholeH_wholeP = 0
    wholeH_wholeP_wholepartS = 0
    wholeH_wholeP_notS = 0
    wholeH_partP = 0
    wholeH_partP_wholepartS = 0
    wholeH_partP_notS = 0
    wholeH_notP = 0
    partH = 0
    partH_wholeP = 0
    partH_wholeP_wholepartS = 0
    partH_wholeP_notS = 0
    partH_partP = 0
    partH_partP_wholepartS = 0
    partH_partP_notS = 0
    partH_notP = 0
    notH = 0
    notH_wholeP = 0
    notH_partP = 0
    notH_notP = 0

    # fill in the data
    for highconf in [WHOLLY, PARTIALLY, NOT_COVERED]:
        for phased in [WHOLLY, PARTIALLY, NOT_COVERED]:
            for swichd in [WHOLLY, PARTIALLY, NOT_COVERED]:
                value = coverage_breakdown[highconf][phased][swichd]
                if value == 0: continue
                if highconf == WHOLLY:
                    wholeH += value
                    if phased == WHOLLY:
                        wholeH_wholeP += value
                        if swichd == NOT_COVERED:
                            wholeH_wholeP_notS += value
                        else:
                            wholeH_wholeP_wholepartS += value
                    elif phased == PARTIALLY:
                        wholeH_partP += value
                        if swichd == NOT_COVERED:
                            wholeH_partP_notS += value
                        else:
                            wholeH_partP_wholepartS += value
                    elif phased == NOT_COVERED:
                        wholeH_notP += value
                elif highconf == PARTIALLY:
                    partH += value
                    if phased == WHOLLY:
                        partH_wholeP += value
                        if swichd == NOT_COVERED:
                            partH_wholeP_notS += value
                        else:
                            partH_wholeP_wholepartS += value
                    elif phased == PARTIALLY:
                        partH_partP += value
                        if swichd == NOT_COVERED:
                            partH_partP_notS += value
                        else:
                            partH_partP_wholepartS += value
                    elif phased == NOT_COVERED:
                        partH_notP += value
                elif highconf == NOT_COVERED:
                    notH += value
                    if phased == WHOLLY:
                        notH_wholeP += value
                    elif phased == PARTIALLY:
                        notH_partP += value
                    elif phased == NOT_COVERED:
                        notH_notP += value
    # sanity
    assert(wholeH_wholeP_wholepartS + wholeH_wholeP_notS == wholeH_wholeP)
    assert(wholeH_partP_wholepartS + wholeH_partP_notS == wholeH_partP)
    assert(wholeH_wholeP + wholeH_partP + wholeH_notP == wholeH)
    assert(partH_wholeP_wholepartS + partH_wholeP_notS == partH_wholeP)
    assert(partH_partP_wholepartS + partH_partP_notS == partH_partP)
    assert(partH_wholeP + partH_partP + partH_notP == partH)
    assert(notH_wholeP + notH_partP + notH_notP == notH)
    total = wholeH + partH + notH

    # Colors
    a, b, c = [plt.cm.Blues, plt.cm.Greens, plt.cm.Reds]
    w = 'w'

    # hiconf
    group_names = ['Fully HighConf\n{}\n{:5.2f}%'.format(wholeH, 100.0*wholeH/total),
                   'Partially HighConf\n{}\n{:5.2f}%'.format(partH, 100.0*partH/total),
                   'Not HighConf\n{}\n{:5.2f}%'.format(notH, 100.0*notH/total)]
    group_size = [wholeH, partH, notH]
    group_colors = [a(0.9), b(0.9), c(0.9)]

    # phased
    subgroup_names = ['Fully Phased\n{}\n{:5.2f}%'.format(wholeH_wholeP, 100.0*wholeH_wholeP/max(1,wholeH)), '', ''] + \
                     ['Fully Phased\n{}\n{:5.2f}%'.format(partH_wholeP, 100.0*partH_wholeP/max(1,partH)), '', ''] + \
                     ['Fully Phased\n{}\n{:5.2f}%'.format(notH_wholeP, 100.0*notH_wholeP/max(1,notH)), '', '']
    subgroup_size = [wholeH_wholeP, wholeH_partP, wholeH_notP] + \
                    [partH_wholeP, partH_partP, partH_notP] + \
                    [notH_wholeP, notH_partP, notH_notP]
    # subgroup_colors = [a(.7), a(.4), a(.1)] + \
    #                   [b(.7), b(.4), b(.1)] + \
    #                   [c(.7), c(.4), c(.1)]
    subgroup_colors = [a(.6), a(.3), w] + \
                      [b(.6), b(.3), w] + \
                      [c(.6), c(.3), w]

    # switch
    switch_names = ['No Switch Error\n{}\n{:5.2f}%'.format(wholeH_wholeP_notS, 100.0*wholeH_wholeP_notS/max(1,wholeH_wholeP)), ''] + ['', ''] + ['']+ \
                   ['No Switch Error\n{}\n{:5.2f}%'.format(partH_wholeP_notS, 100.0*partH_wholeP_notS/max(1,partH_wholeP)), ''] + ['', ''] + [''] + \
                   ['', '', '']
    switch_size = [wholeH_wholeP_notS, wholeH_wholeP_wholepartS] + [wholeH_partP_notS, wholeH_partP_wholepartS] + [wholeH_notP] + \
                    [partH_wholeP_notS, partH_wholeP_wholepartS] + [partH_partP_notS, partH_partP_wholepartS] + [partH_notP] + \
                    [notH_wholeP, notH_partP, notH_notP]
    # switch_colors = [a(.6), a(.5)] + [a(.3), a(.2)] + [a(.1)] + \
    #                 [b(.6), b(.5)] + [b(.3), b(.2)] + [b(.1)] + \
    #                 [c(.7), c(.4), c(.1)]
    switch_colors = [a(.5), a(.4)] + [w, w] + [w] + \
                    [b(.5), b(.4)] + [w, w] + [w] + \
                    [w, w, w]



    # First Ring (outside)
    fig, ax = plt.subplots()
    ax.axis('equal')
    mypie, _ = ax.pie(group_size, radius=1.3, labels=group_names, colors=group_colors)
    plt.setp(mypie, width=0.3, edgecolor='white')

    # Second Ring (Middle)
    mypie2, _ = ax.pie(subgroup_size, radius=1.3 - 0.3, labels=subgroup_names, labeldistance=0.7,
                       colors=subgroup_colors)
    plt.setp(mypie2, width=0.4, edgecolor='white')

    # Third Ring (Inside)
    mypie3, _ = ax.pie(switch_size, radius=1.3 - 0.3 - 0.4, labels=switch_names, labeldistance=0.3,
                       colors=switch_colors)
    plt.setp(mypie3, width=0.6, edgecolor='white')


    plt.margins(0, 0)

    # show it
    if output_base is not None:
        plt.savefig("{}.coverage_donut.png".format(output_base))
    plt.show()


def plot_coverage_grouped_bar(coverage_breakdown, output_base):

    # data we plot
    wholeH = 0
    wholeH_wholeP = 0
    wholeH_wholeP_wholepartS = 0
    wholeH_wholeP_notS = 0
    wholeH_partP = 0
    wholeH_partP_wholepartS = 0
    wholeH_partP_notS = 0
    wholeH_notP = 0
    partH = 0
    partH_wholeP = 0
    partH_wholeP_wholepartS = 0
    partH_wholeP_notS = 0
    partH_partP = 0
    partH_partP_wholepartS = 0
    partH_partP_notS = 0
    partH_notP = 0
    notH = 0
    notH_wholeP = 0
    notH_partP = 0
    notH_notP = 0

    # fill in the data
    for highconf in [WHOLLY, PARTIALLY, NOT_COVERED]:
        for phased in [WHOLLY, PARTIALLY, NOT_COVERED]:
            for swichd in [WHOLLY, PARTIALLY, NOT_COVERED]:
                value = coverage_breakdown[highconf][phased][swichd]
                if value == 0: continue
                if highconf == WHOLLY:
                    wholeH += value
                    if phased == WHOLLY:
                        wholeH_wholeP += value
                        if swichd == NOT_COVERED:
                            wholeH_wholeP_notS += value
                        else:
                            wholeH_wholeP_wholepartS += value
                    elif phased == PARTIALLY:
                        wholeH_partP += value
                        if swichd == NOT_COVERED:
                            wholeH_partP_notS += value
                        else:
                            wholeH_partP_wholepartS += value
                    elif phased == NOT_COVERED:
                        wholeH_notP += value
                elif highconf == PARTIALLY:
                    partH += value
                    if phased == WHOLLY:
                        partH_wholeP += value
                        if swichd == NOT_COVERED:
                            partH_wholeP_notS += value
                        else:
                            partH_wholeP_wholepartS += value
                    elif phased == PARTIALLY:
                        partH_partP += value
                        if swichd == NOT_COVERED:
                            partH_partP_notS += value
                        else:
                            partH_partP_wholepartS += value
                    elif phased == NOT_COVERED:
                        partH_notP += value
                elif highconf == NOT_COVERED:
                    notH += value
                    if phased == WHOLLY:
                        notH_wholeP += value
                    elif phased == PARTIALLY:
                        notH_partP += value
                    elif phased == NOT_COVERED:
                        notH_notP += value
    # sanity
    assert(wholeH_wholeP_wholepartS + wholeH_wholeP_notS == wholeH_wholeP)
    assert(wholeH_partP_wholepartS + wholeH_partP_notS == wholeH_partP)
    assert(wholeH_wholeP + wholeH_partP + wholeH_notP == wholeH)
    assert(partH_wholeP_wholepartS + partH_wholeP_notS == partH_wholeP)
    assert(partH_partP_wholepartS + partH_partP_notS == partH_partP)
    assert(partH_wholeP + partH_partP + partH_notP == partH)
    assert(notH_wholeP + notH_partP + notH_notP == notH)
    total = wholeH + partH + notH

    # set width of bar
    barWidth = 0.25

    # set height of bar
    # bars1 = [wholeH, wholeH_wholeP, wholeH_wholeP_notS]
    # bars2 = [partH, partH_wholeP, partH_wholeP_notS]
    # bars3 = [notH, notH_wholeP, 0]
    bars1 = [wholeH, partH, notH]
    bars2 = [wholeH_wholeP, partH_wholeP, notH_wholeP]
    bars3 = [wholeH_wholeP_notS, partH_wholeP_notS, 0]

    a, b, c = [plt.cm.Blues, plt.cm.Greens, plt.cm.Reds]

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1,ncols=1)
    rect1 = ax1.bar(r1, bars1, color=[a(0.6), b(0.6), c(0.6)], hatch='', width=barWidth, edgecolor='white', label='Total Genes')
    rect2 = ax1.bar(r2, bars2, color=[a(0.45), b(0.45), c(0.45)], hatch='-',  width=barWidth, edgecolor='white', label='Wholly Phased')
    rect3 = ax1.bar(r3, bars3, color=[a(0.3), b(0.3), c(0.3)], hatch='x', width=barWidth, edgecolor='white', label='No Switch Errors')

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], ['Wholly\nin High\nConfidence', 'Partly\nin High\nConfidence', 'Not\nin High\nConfidence'])

    def autolabel(rects, percentage_denominator):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{} ({:4.1f}%)'.format(height, 100.0*height/max(1,percentage_denominator[i])),
                        xy=(rect.get_x() + rect.get_width() / 2, height - .35 * max(bars1)),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90,
                        ha='center', va='bottom')#, weight='bold')

    autolabel(rect1, [total, total, total])
    autolabel(rect2, bars1)
    autolabel(rect3, bars2)

    # Create legend & Show graphic
    plt.legend()

    # show it
    if output_base is not None:
        plt.savefig("{}.coverage_grouped_bar.png".format(output_base))
    plt.show()



def plot_coverage_grouped_bar_manual_for_files(output_base):

    # data we plot
    total_genes = 60656
    samples = ["Total Genes\n(GRCh38)", "HG003", "HG004"]
    coverage = [60656, 53817, 55234]


    def autolabel(rects, percentage_denominator):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{} ({:4.1f}%)'.format(height, 100.0*height/max(1,percentage_denominator[i])),
                        xy=(rect.get_x() + rect.get_width() / 2, height - .35 * total_genes),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90,
                        ha='center', va='bottom')#, weight='bold')

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1,ncols=1)
    for i in range(len(samples)):
        rect = ax1.bar([i], [coverage[i]], color=plt.cm.Reds(0.6), hatch='' if i == 0 else '-', edgecolor='white')
        autolabel(rect, [total_genes])

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r for r in range(len(samples))], samples)

    # show it
    if output_base is not None:
        plt.savefig("{}.coverage_grouped_bar_manual.png".format(output_base))
    plt.show()

def main():
    args = parse_args()

    files = [args.phased_vcf, args.gene_gtf, args.whatshap_compare_switch_bed, args.truth_annotated_vcf, args.high_confidence_bed]
    missing_files = list(filter(lambda x: x is not None and not os.path.isfile, files))
    if len(missing_files) != 0:
        log("Could not find specified files: {}".format(missing_files))
        sys.exit(1)

    # prep for writing output
    output_base = args.output_base
    if output_base is None and args.output_base_from_input:
        output_base = os.path.basename(args.phased_vcf).rstrip(".gz").rstrip(".vcf")

    # one-off plot for the paper
    # plot_coverage_grouped_bar_manual_for_files(output_base)
    # sys.exit()

    # get data
    phasesets = get_phasesets_from_vcf(args.phased_vcf, args)
    genes, switches, high_conf_regions = get_bed_regions(args)
    genes_wholly_hiconf, genes_partially_hiconf, genes_not_hiconf = get_genes_stratified_by_bed(genes, high_conf_regions, "HIGHCONF")
    genes_wholly_phased, genes_partially_phased, genes_not_phased = get_genes_stratified_by_bed(genes, phasesets, "PHASESET")
    genes_wholly_swichd, genes_partially_swichd, genes_not_swichd = get_genes_stratified_by_bed(genes, switches, "SWITCHES")
    get_precision_recall_for_genes(genes, args)

    # output data
    output_data_tsv = None
    if output_base is not None:
        output_data_tsv = "{}.gene_phasing_stats.tsv".format(output_base)
    outputstream = None

    # info over all genes
    coverage_breakdown = collections.defaultdict(lambda : collections.defaultdict(lambda : collections.defaultdict(lambda : 0)))
    highconf_variants = collections.defaultdict(lambda : collections.defaultdict(lambda : 0))
    total_genes = 0

    # write data out
    try:
        if output_data_tsv is not None:
            outputstream = open(output_data_tsv, 'w')

        # write out determined data
        for chrom in sorted(list(genes.keys())):
            for gene in genes[chrom]:
                # get coverages
                hiconf_cov = get_coverage_state(genes_wholly_hiconf, genes_partially_hiconf, genes_not_hiconf, chrom, gene)
                phased_cov = get_coverage_state(genes_wholly_phased, genes_partially_phased, genes_not_phased, chrom, gene)
                swichd_cov = get_coverage_state(genes_wholly_swichd, genes_partially_swichd, genes_not_swichd, chrom, gene)

                # describe overlaps
                hiconf_overlaps = []
                phased_overlaps = []
                swichd_overlaps = []
                for bed in gene.overlaps:
                    dest = (hiconf_overlaps if bed.type == TYPE_HICONF else
                        (phased_overlaps if bed.type == TYPE_PHASESET else
                        (swichd_overlaps if bed.type == TYPE_SWITCH else None)))
                    dest.append("{}:{}-{}".format(bed.chrom, bed.end, bed.start))
                hiconf_overlaps = "[\"" + "\",\"".join(hiconf_overlaps) + "\"]"
                phased_overlaps = "[\"" + "\",\"".join(phased_overlaps) + "\"]"
                swichd_overlaps = "[\"" + "\",\"".join(swichd_overlaps) + "\"]"

                # get precision/recall
                total_variants = gene.get_total_variants()
                snp_precision = "-1" if total_variants == 0 else "{:.5f}".format(gene.snp_precision)
                snp_recall = "-1" if total_variants == 0 else "{:.5f}".format(gene.snp_recall)
                snp_tp, snp_fp, snp_fn = gene.snp_tp, gene.snp_fp, gene.snp_fn
                indel_precision = "-1" if total_variants == 0 else "{:.5f}".format(gene.indel_precision)
                indel_recall = "-1" if total_variants == 0 else "{:.5f}".format(gene.indel_recall)
                indel_tp, indel_fp, indel_fn = gene.indel_tp, gene.indel_fp, gene.indel_fn

                # for whole stats
                total_genes += 1
                coverage_breakdown[hiconf_cov][phased_cov][swichd_cov] += 1
                if hiconf_cov == WHOLLY:
                    highconf_variants[SNP][TP] += snp_tp
                    highconf_variants[SNP][FP] += snp_fp
                    highconf_variants[SNP][FN] += snp_fn
                    highconf_variants[INDEL][TP] += indel_tp
                    highconf_variants[INDEL][FP] += indel_fp
                    highconf_variants[INDEL][FN] += indel_fn


                # todo
                # write to outputstream

    finally:
        if output_data_tsv is not None and outputstream is not None:
            outputstream.close()

    log("")
    log("Coverage Breakdown:")
    for highconf in [WHOLLY, PARTIALLY, NOT_COVERED]:
        hicf_printed = False
        hicf_total = 0
        for phased in [WHOLLY, PARTIALLY, NOT_COVERED]:
            phsd_printed = False
            phsd_total = 0
            for swichd in [WHOLLY, PARTIALLY, NOT_COVERED]:
                value = coverage_breakdown[highconf][phased][swichd]
                if value == 0: continue
                hicf_str = "HICONF {}  ".format(highconf) if not hicf_printed else "" #26
                phsd_str = "PHASED {}  ".format(phased) if not phsd_printed else ""
                swch_str = "SWITCH {}  ".format(swichd)
                log("\t{:26s} {:26s} {:26s} {:5d}  {:5.1f}%  {:5.1f}%".format( hicf_str, phsd_str, swch_str, value,
                    100.0 * value / sum(coverage_breakdown[highconf][phased].values()), 100.0 * value / total_genes))
                hicf_printed = True
                phsd_printed = True
                hicf_total += value
                phsd_total += value
            if phsd_total != 0 and len(list(filter(lambda x: x > 0, coverage_breakdown[highconf][phased].values()))) > 1:
                log("\t{:26s} {:26s} {:26s} {:5d}  {:6s}  {:5.1f}%".format("", "", "       TOTAL", phsd_total, "", 100.0 * phsd_total / total_genes))
        if hicf_total != 0:
            log("\t{:26s} {:26s} {:26s} {:5d}  {:6s}  {:5.1f}%".format("", "       TOTAL", "", hicf_total, "", 100.0 * hicf_total / total_genes))
    plot_coverage_donut(coverage_breakdown, output_base)
    plot_coverage_grouped_bar(coverage_breakdown, output_base)



    snp_tp = highconf_variants[SNP][TP]
    snp_fp = highconf_variants[SNP][FP]
    snp_fn = highconf_variants[SNP][FN]
    indel_tp = highconf_variants[INDEL][TP]
    indel_fp = highconf_variants[INDEL][FP]
    indel_fn = highconf_variants[INDEL][FN]
    log("")
    log("Wholly Covered HiConf Gene Variant Stats:")
    log("\tSNP Precision:   {:.5f} {:8}/{}+{}".format(snp_tp/max(1, snp_tp+snp_fp), snp_tp, snp_tp, snp_fp))
    log("\tSNP Recall:      {:.5f} {:8}/{}+{}".format(snp_tp/max(1, snp_tp+snp_fn), snp_tp, snp_tp, snp_fn))
    log("\tINDEL Precision: {:.5f} {:8}/{}+{}".format(indel_tp/max(1, indel_tp+indel_fp), indel_tp, indel_tp, indel_fp))
    log("\tINDEL Recall:    {:.5f} {:8}/{}+{}".format(indel_tp/max(1, indel_tp+indel_fn), indel_tp, indel_tp, indel_fn))
    log("\tALL Precision:   {:.5f} {:8}/{}+{}".format((snp_tp+indel_tp)/max(1, indel_tp+indel_fp+snp_tp+snp_fp), snp_tp+indel_tp, indel_tp+snp_tp, indel_fp+snp_fp))
    log("\tALL Recall:      {:.5f} {:8}/{}+{}".format((snp_tp+indel_tp)/max(1, indel_tp+indel_fn+snp_tp+snp_fn), snp_tp+indel_tp, indel_tp+snp_tp, indel_fn+snp_fn))


    pass









if __name__ == "__main__":
    main()

