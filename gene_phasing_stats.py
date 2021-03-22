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


plt.style.use('ggplot')
text_fontsize = 5
# plt.rcParams['ytick.labelsize']=text_fontsize+4
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')


GENE_GTF_CHROM_IDX=0
GENE_GTF_TYPE_IDX=2
GENE_GTF_START_POS_IDX=3
GENE_GTF_END_POS_IDX=4
GENE_GTF_INFO_IDX=8

GENE_TYPE_GENE="gene"
GENE_TYPE_EXON="exon"
GENE_TYPE_CDS="CDS"
GENE_TYPE_START_CODON="start_codon"
GENE_TYPE_END_CODON="end_codon"

TYPE_GENE="g"
TYPE_HICONF="h"
TYPE_SWITCH="s"
TYPE_PHASESET="p"
TYPE_EXON="x"
TYPE_CODING="c"

WHOLLY="wholly"
PARTIALLY="partially"
NOT_COVERED="not"

ERROR="error"
NO_ERROR="no error"


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
        self.exon_frameshift_tp = 0
        self.exon_frameshift_fp = 0
        self.exon_frameshift_fn = 0
        self.coding_frameshift_tp = 0
        self.coding_frameshift_fp = 0
        self.coding_frameshift_fn = 0
        self.highconf_overlap = 0.0

    def calculate_precision_recall(self, snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn):
        self.snp_tp = snp_tp
        self.snp_fp = snp_fp
        self.snp_fn = snp_fn
        self.indel_tp = indel_tp
        self.indel_fp = indel_fp
        self.indel_fn = indel_fn
        if self.get_total_variants() == 0:
            return False
        self.snp_precision = 1.0 * snp_tp / max(1, snp_tp+snp_fp)
        self.snp_recall = 1.0 * snp_tp / max(1, snp_tp+snp_fn)
        self.indel_precision = 1.0 * indel_tp / max(1, indel_tp+indel_fp)
        self.indel_recall = 1.0 * indel_tp / max(1, indel_tp+indel_fn)
        return True

    def get_total_variants(self):
        return sum([self.snp_tp, self.snp_fp, self.snp_fn, self.indel_tp, self.indel_fp, self.indel_fn])

    def get_gene_custom_identifier(self):
        if "gene_name" in self.info:
            return "{} ({}:{}-{})".format(self.info["gene_name"], self.chrom, self.start, self.end)
        return "{}:{}-{}".format(self.chrom, self.start, self.end)

    def get_gene_id(self):
        return self.info["gene_id"] if "gene_id" in self.info else None



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
    parser.add_argument('--stratification_tsv', '-s', dest='stratification_tsv', default=None, required=False, type=str,
                       help='Custom-created TSV from hap.py with one subtype and no duplicate records other than ' +\
                       'INDEL and SNP. Tool will analyze genes without errors if set')
    parser.add_argument('--gene_ids_of_interest', '-n', dest='gene_ids_of_interest', default=None, required=False, type=str,
                       help='Only analyze these genes (if set)')
    parser.add_argument('--repair_chr_names', '-C', dest='repair_chr_names', default=False, required=False, action='store_true',
                       help='Add "chr" to contig names without (there are sometimes "1" and "chr1" mismatches for grch37)')

    parser.add_argument('--only_protien_coding', '-p', dest='only_protien_coding', default=False, required=False, action='store_true',
                       help='Only run for genes with type "protein_coding"')
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


def get_bed_regions(args, genes_of_interest=None):
    genes = collections.defaultdict(lambda: list())
    switches = collections.defaultdict(lambda: list())
    high_conf_regions = collections.defaultdict(lambda: list())

    # genes
    with open(args.gene_gtf) as fin:
        current_gene = None
        for line in fin:
            # errors and invalid records
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) <= GENE_GTF_INFO_IDX:
                log("Unexpected line reading gene GTF: '{}'".format(line))
                continue
            # skip non-gene parts
            type = parts[GENE_GTF_TYPE_IDX]
            chrom = fix_chrom(parts[GENE_GTF_CHROM_IDX], args)
            start_pos = int(parts[GENE_GTF_START_POS_IDX])
            end_pos = int(parts[GENE_GTF_END_POS_IDX])
            if type == GENE_TYPE_GENE:
                # get data
                info = parts[GENE_GTF_INFO_IDX]
                current_gene = GeneInfo(chrom, start_pos, end_pos, info)
                # always construct, but only save it if it's protein_coding (or we want all genes)
                if not args.only_protien_coding or current_gene.info["gene_type"] == 'protein_coding':
                    # further filtering based on restricted gene set (if appropriate)
                    if genes_of_interest is None or current_gene.get_gene_id() in genes_of_interest:
                        genes[chrom].append(current_gene)
            elif type == GENE_TYPE_EXON:
                current_gene.overlaps.append(BedRecord(chrom, start_pos, end_pos, TYPE_EXON))
            elif type in (GENE_TYPE_CDS, GENE_TYPE_START_CODON, GENE_TYPE_END_CODON):
                current_gene.overlaps.append(BedRecord(chrom, start_pos, end_pos, TYPE_CODING))

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


def has_frameshift_allele(record: pysam.VariantRecord):
    ref_size = len(record.ref)
    for allele in record.alts:
        if abs(len(allele) - ref_size) % 3 != 0:
            return True
    return False


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
    genes_with_exon_frameshift = 0
    genes_with_exon_frameshift_err = 0
    genes_with_coding_frameshift = 0
    genes_with_coding_frameshift_err = 0
    total_gene_sequence = 0
    total_exon_sequence = 0
    total_coding_sequence = 0

    # get variants for all genes
    for chrom in sorted(list(genes.keys())):
        for gene in genes[chrom]:
            snp_tp = 0
            snp_fp = 0
            snp_fn = 0
            indel_tp = 0
            indel_fp = 0
            indel_fn = 0

            frameshift_exon_tp = 0
            frameshift_exon_fp = 0
            frameshift_exon_fn = 0
            frameshift_coding_tp = 0
            frameshift_coding_fp = 0
            frameshift_coding_fn = 0

            # we only track in hiconf
            higconf_beds = list(filter(lambda x: x.type == TYPE_HICONF, gene.overlaps))
            exon_beds = list(filter(lambda x: x.type == TYPE_EXON, gene.overlaps))
            coding_beds = list(filter(lambda x: x.type == TYPE_CODING, gene.overlaps))

            total_gene_sequence += gene.end - gene.start
            total_exon_sequence += sum(map(lambda x: x.end - x.start, exon_beds))
            total_coding_sequence += sum(map(lambda x: x.end - x.start, coding_beds))
            try:
                vcf_records = vcf.fetch(unfix_chrom(chrom, args), gene.start, gene.end)
            except:
                if unfix_chrom(chrom, args) != "chrM":
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

                # exon or coding
                exon_region = False
                coding_region = False
                for e in exon_beds:
                    if e.start <= pos < e.end:
                        exon_region = True
                        break
                for c in coding_beds:
                    if c.start <= pos < c.end:
                        coding_region = True
                        break

                # variant info
                truth_bd = record.samples[truth_sample_idx]["BD"]
                truth_bvt = record.samples[truth_sample_idx]["BVT"]
                query_bd = record.samples[query_sample_idx]["BD"]
                query_bvt = record.samples[query_sample_idx]["BVT"]
                is_frameshift = has_frameshift_allele(record)

                pre_total = sum([snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn])
                if truth_bd == "FN":
                    if truth_bvt == "INDEL":
                        indel_fn += 1
                        if is_frameshift:
                            if exon_region: frameshift_exon_fn += 1
                            if coding_region: frameshift_coding_fn += 1
                    elif truth_bvt == "SNP":
                        snp_fn += 1
                if query_bd == "FP":
                    if query_bvt == "INDEL":
                        indel_fp += 1
                        if is_frameshift:
                            if exon_region: frameshift_exon_fp += 1
                            if coding_region: frameshift_coding_fp += 1
                    elif query_bvt == "SNP":
                        snp_fp += 1
                if query_bd == "TP" or truth_bd == "TP":
                    if query_bvt == "INDEL" or truth_bvt == "INDEL":
                        indel_tp += 1
                        if is_frameshift:
                            if exon_region: frameshift_exon_tp += 1
                            if coding_region: frameshift_coding_tp += 1
                    if query_bvt == "SNP" or truth_bvt == "SNP":
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

            # precision and recal
            gene.calculate_precision_recall(snp_tp, snp_fp, snp_fn, indel_tp, indel_fp, indel_fn)
            if gene.get_total_variants() == 0:
                genes_without_variants += 1
            else:
                genes_with_variants += 1

            # frameshift tracking
            gene.exon_frameshift_tp = frameshift_exon_tp
            gene.exon_frameshift_fp = frameshift_exon_fp
            gene.exon_frameshift_fn = frameshift_exon_fn
            if frameshift_exon_fn + frameshift_exon_fp + frameshift_exon_tp > 0:
                genes_with_exon_frameshift += 1
                if frameshift_exon_fn + frameshift_exon_fp > 0:
                    genes_with_exon_frameshift_err += 1
            gene.coding_frameshift_tp = frameshift_coding_tp
            gene.coding_frameshift_fp = frameshift_coding_fp
            gene.coding_frameshift_fn = frameshift_coding_fn
            if frameshift_coding_fp + frameshift_coding_fn + frameshift_coding_tp > 0:
                genes_with_coding_frameshift += 1
                if frameshift_coding_fp + frameshift_coding_fn > 0:
                    genes_with_coding_frameshift_err += 1

    log("")
    log("Analyzed Genes:")
    log("\tTotal gene sequence               {}".format(total_gene_sequence))
    log("\tTotal exon sequence               {} ({:6.4f})".format(total_exon_sequence, total_exon_sequence/total_gene_sequence))
    log("\tTotal coding sequence             {} ({:6.4f})".format(total_coding_sequence, total_coding_sequence/total_gene_sequence))
    log("\tVariant classification counts:    0:{} 1:{} 2+:{}".format(update_zero_count, update_one_count, update_multiple_count))
    log("\tVariants in HiConf:               {}".format(high_conf_variants))
    log("\tVariants not in HiConf:           {}".format(not_high_conf_variants))
    log("\tVariants missing tags:            {}".format(missing_tag_variants))
    log("\tGenes without variants:           {} ({:.3f})".format(
        genes_without_variants, genes_without_variants / (genes_with_variants + genes_without_variants)))
    log("\tGenes with variants:              {} ({:.3f})".format(
        genes_with_variants, genes_with_variants / (genes_with_variants + genes_without_variants)))
    log("\tGenes with exon frameshift:       {} ({:.3f})".format(
        genes_with_exon_frameshift, genes_with_exon_frameshift / (genes_with_variants)))
    log("\tGenes with exon frameshift err:   {} ({:.3f})".format(
        genes_with_exon_frameshift_err, genes_with_exon_frameshift_err / genes_with_exon_frameshift))
    log("\tGenes with coding frameshift:     {} ({:.3f})".format(
        genes_with_coding_frameshift, genes_with_coding_frameshift / (genes_with_variants)))
    log("\tGenes with coding frameshift err: {} ({:.3f})".format(
        genes_with_coding_frameshift_err, genes_with_coding_frameshift_err / genes_with_coding_frameshift))


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


def get_genes_stratified_by_bed(genes, bed_regions, bed_desc, update_gene_hiconf_coverage=False):

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
            if update_gene_hiconf_coverage:
                total_coverage = 0.0
                for overlap in hiconf_overlaps:
                    total_coverage += min(gene.end, overlap.end) - max(gene.start, overlap.start)
                gene.highconf_overlap = total_coverage / (gene.end - gene.start)
                if gene.highconf_overlap > 1.0:
                    log("Gene {} got highconf coverage {}".format(gene.get_gene_custom_identifier(), gene.highconf_overlap))
        log("\tFor {}, of {} total genes got {} wholly covered, {} partially covered, and {} not covered".format(
            chr, len(chr_genes), len(genes_wholly_covered[chr]), len(genes_partially_covered[chr]), len(genes_not_covered[chr])))

    return genes_wholly_covered, genes_partially_covered, genes_not_covered


STRAT_TYPE= "Type"
STRAT_GENE_REGION= "Subset"
STRAT_F1_SCORE= "METRIC.F1_Score"
STRAT_PRECISION= "METRIC.Precision"
STRAT_RECALL= "METRIC.Recall"
STRAT_REGION_SIZE= "Subset.Size"
STRAT_REGION_HICONF_SIZE= "Subset.IS_CONF.Size"
STRAT_TP_T= "TRUTH.TP"
STRAT_TP_Q= "QUERY.TP"
STRAT_FN="TRUTH.FN"
STRAT_FP="QUERY.FP"
COLUMNS=[STRAT_TYPE, STRAT_GENE_REGION, STRAT_RECALL, STRAT_PRECISION, STRAT_F1_SCORE, STRAT_REGION_SIZE, STRAT_REGION_HICONF_SIZE, STRAT_TP_T, STRAT_TP_Q, STRAT_FN, STRAT_FP]


def get_gene_stratification_analysis(strat_tsv, high_conf_threshold=.0):

    data = collections.defaultdict(lambda: collections.defaultdict(lambda: 0.0))
    column_map = None
    with open(strat_tsv) as tsv_in:
        for line in tsv_in:
            line = line.strip().split("\t")
            if column_map is None:
                column_map = {h:i for (i,h) in filter(lambda x: x[1] in COLUMNS, enumerate(line))}
                continue
            line_data = {c:line[column_map[c]] for c in COLUMNS}
            size = int(float(line_data[STRAT_REGION_SIZE]))
            hc_size = int(float(line_data[STRAT_REGION_HICONF_SIZE]))
            tp_t = int(float(line_data[STRAT_TP_T]))
            tp_q = int(float(line_data[STRAT_TP_Q]))
            fp = int(float(line_data[STRAT_FP]))
            fn = int(float(line_data[STRAT_FN]))
            id = line_data[STRAT_GENE_REGION]
            data[id][STRAT_REGION_SIZE] = size
            data[id][STRAT_REGION_HICONF_SIZE] = hc_size
            data[id][STRAT_TP_T] += tp_t
            data[id][STRAT_TP_Q] += tp_q
            data[id][STRAT_FP] += fp
            data[id][STRAT_FN] += fn
            data[id][STRAT_GENE_REGION] = id
            # if line_data[TP] != line_data[TP2] and line_data[F1_SCORE] not in ("1.0", ""):
            #     pass

    gene_ids_above_threshold = set()
    no_error_gene_ids_above_threshold = set()
    def get_gene_id(gene_data):
        return gene_data[STRAT_GENE_REGION].split(";")[1]

    total_genes = 0
    below_hiconf_threshold = 0
    above_hiconf_threshold = 0
    had_no_variants = 0
    had_variants = 0
    had_variants_but_no_errors = 0
    no_errors = 0
    for (id, gene_data) in data.items():
        total_genes += 1
        if gene_data[STRAT_REGION_HICONF_SIZE]/max(gene_data[STRAT_REGION_SIZE], 1) < high_conf_threshold:
            below_hiconf_threshold += 1
            continue
        above_hiconf_threshold += 1
        gene_ids_above_threshold.add(get_gene_id(gene_data))
        if gene_data[STRAT_TP_T] + gene_data[STRAT_TP_Q] == 0:
            had_no_variants += 1
        else:
            had_variants += 1
            if gene_data[STRAT_FP] + gene_data[STRAT_FN] == 0:
                had_variants_but_no_errors += 1
        if gene_data[FP] + gene_data[STRAT_FN] == 0:
            no_errors += 1
            no_error_gene_ids_above_threshold.add(get_gene_id(gene_data))

    log("Stats for {}:".format(strat_tsv))
    log("\tGene bodies:              {} ".format(total_genes))
    log("\t  Above {:.2f} HighConf:    {} ({:.5f})".format(high_conf_threshold, above_hiconf_threshold,
                                                                 above_hiconf_threshold/total_genes))
    log("\t    With variants:        {} ({:.5f})".format(had_variants, had_variants/above_hiconf_threshold))
    log("\t      With no errors:     {} ({:.5f})".format(had_variants_but_no_errors, had_variants_but_no_errors/had_variants))
    log("\t    With no errors:       {} ({:.5f})".format(no_errors, no_errors/above_hiconf_threshold))

    return data, gene_ids_above_threshold, no_error_gene_ids_above_threshold


def get_gene_id_to_gene_map(genes):
    gene_map = dict()
    for chrom_genes in genes.values():
        for gene in chrom_genes:
            id = gene.get_gene_id()
            assert(id is not None)
            gene_map[id] = gene
    return gene_map


def print_highconf_coverage_histogram(genes, bucket_count = 20, coverage_threshold = .8):
    buckets = [0 for _ in range(bucket_count)]

    with_exon_frameshift = [0 for _ in range(bucket_count)]
    with_exon_frameshift_tp = [0 for _ in range(bucket_count)]
    with_exon_frameshift_fp = [0 for _ in range(bucket_count)]
    with_exon_frameshift_fn = [0 for _ in range(bucket_count)]

    with_coding_frameshift = [0 for _ in range(bucket_count)]
    with_coding_frameshift_tp = [0 for _ in range(bucket_count)]
    with_coding_frameshift_fp = [0 for _ in range(bucket_count)]
    with_coding_frameshift_fn = [0 for _ in range(bucket_count)]

    for chrom in genes.keys():
        for gene in genes[chrom]:
            gene_coverage = gene.highconf_overlap
            bucket_idx = int(gene_coverage * bucket_count-1)
            buckets[bucket_idx] += 1
            # exon frameshifts
            if gene.exon_frameshift_fn + gene.exon_frameshift_fp > 0:
                with_exon_frameshift[bucket_idx] += 1
                if gene.exon_frameshift_fn > 0:
                    with_exon_frameshift_fn[bucket_idx] += 1
                if gene.exon_frameshift_fp > 0:
                    with_exon_frameshift_fp[bucket_idx] += 1
            if gene.exon_frameshift_tp > 0:
                with_exon_frameshift_tp[bucket_idx] += 1
            # coding frameshift
            if gene.coding_frameshift_fn + gene.coding_frameshift_fp > 0:
                with_coding_frameshift[bucket_idx] += 1
                if gene.coding_frameshift_fn > 0:
                    with_coding_frameshift_fn[bucket_idx] += 1
                if gene.coding_frameshift_fp > 0:
                    with_coding_frameshift_fp[bucket_idx] += 1
            if gene.coding_frameshift_tp > 0:
                with_coding_frameshift_tp[bucket_idx] += 1
    max_size = max(buckets)
    log("High Conf Coverage Ratios over Genes:")
    log("\tCOV   HIST                             COUNT -- EXON     CODING   EXON_TP/FP/FN   CODING_TP/FP/FN")
    for i in reversed(range(bucket_count)):
        log("\t{:4.2f}+ {:32s} {:5d} -- {:6.4f}   {:6.4f}   {:13s}   {:13s}".format(
            i / bucket_count, "#" * int(32 * buckets[i] / max_size), buckets[i],
            with_exon_frameshift[i] / buckets[i], with_coding_frameshift[i] / buckets[i],
            "{}/{}/{}".format(with_exon_frameshift_tp[i], with_exon_frameshift_fp[i], with_exon_frameshift_fn[i]),
            "{}/{}/{}".format(with_coding_frameshift_tp[i], with_coding_frameshift_fp[i], with_coding_frameshift_fn[i])))


def print_frameshift_rates_above_highconf_coverage_threhsold(genes, coverage_threshold=.8):
    total_genes = 0

    exon_total_frameshifts = 0
    exon_frameshift_tp = 0
    exon_frameshift_fp = 0
    exon_frameshift_fn = 0
    exon_genes_with_frameshift = 0
    exon_genes_with_frameshift_fp = 0
    exon_genes_with_frameshift_fn = 0
    exon_genes_with_frameshift_fp_or_fn = 0

    coding_total_frameshifts = 0
    coding_frameshift_tp = 0
    coding_frameshift_fp = 0
    coding_frameshift_fn = 0
    coding_genes_with_frameshift = 0
    coding_genes_with_frameshift_fp = 0
    coding_genes_with_frameshift_fn = 0
    coding_genes_with_frameshift_fp_or_fn = 0

    for chrom in genes.keys():
        for gene in genes[chrom]:
            gene_coverage = gene.highconf_overlap
            if gene_coverage < coverage_threshold:
                continue
            total_genes += 1

            # exon stuff
            exon_total_frameshifts += gene.exon_frameshift_fn + gene.exon_frameshift_tp
            exon_frameshift_tp += gene.exon_frameshift_tp
            exon_frameshift_fp += gene.exon_frameshift_fp
            exon_frameshift_fn += gene.exon_frameshift_fn
            if gene.exon_frameshift_fn + gene.exon_frameshift_tp > 0:
                exon_genes_with_frameshift += 1
            if gene.exon_frameshift_fn > 0:
                exon_genes_with_frameshift_fn += 1
            if gene.exon_frameshift_fp > 0:
                exon_genes_with_frameshift_fp += 1
            if gene.exon_frameshift_fn + gene.exon_frameshift_fp > 0:
                exon_genes_with_frameshift_fp_or_fn += 1

            # coding stuff
            coding_total_frameshifts += gene.coding_frameshift_fn + gene.coding_frameshift_tp
            coding_frameshift_tp += gene.coding_frameshift_tp
            coding_frameshift_fp += gene.coding_frameshift_fp
            coding_frameshift_fn += gene.coding_frameshift_fn
            if gene.coding_frameshift_fn + gene.coding_frameshift_tp > 0:
                coding_genes_with_frameshift += 1
            if gene.coding_frameshift_fn > 0:
                coding_genes_with_frameshift_fn += 1
            if gene.coding_frameshift_fp > 0:
                coding_genes_with_frameshift_fp += 1
            if gene.coding_frameshift_fn + gene.coding_frameshift_fp > 0:
                coding_genes_with_frameshift_fp_or_fn += 1

    log("Frameshift Rates in genes covered spanned by {:.1f}%+ high confidence:".format(coverage_threshold*100.0))
    log("\tExon frameshift TP/FP/FN:         {}/{}/{}".format(exon_frameshift_tp, exon_frameshift_fp, exon_frameshift_fn))
    log("\tExon frameshift precision:        {:.5f}".format(exon_frameshift_tp/(exon_frameshift_fp + exon_frameshift_tp)))
    log("\tExon frameshift recall:           {:.5f}".format(exon_frameshift_tp/(exon_frameshift_fn + exon_frameshift_tp)))
    log("\tCoding frameshift TP/FP/FN:       {}/{}/{}".format(coding_frameshift_tp, coding_frameshift_fp, coding_frameshift_fn))
    log("\tCoding frameshift precision:      {:.5f}".format(coding_frameshift_tp/(coding_frameshift_fp + coding_frameshift_tp)))
    log("\tCoding frameshift recall:         {:.5f}".format(coding_frameshift_tp/(coding_frameshift_fn + coding_frameshift_tp)))
    log("\tGenes with exon frameshift:       {} ({:.5f})".format(exon_genes_with_frameshift, exon_genes_with_frameshift / total_genes))
    log("\tGenes with exon frameshift err:   {} ({:.5f})".format(exon_genes_with_frameshift_fp_or_fn, exon_genes_with_frameshift_fp_or_fn / exon_genes_with_frameshift))
    log("\tGenes with coding frameshift:     {} ({:.5f})".format(coding_genes_with_frameshift, coding_genes_with_frameshift / total_genes))
    log("\tGenes with coding frameshift err: {} ({:.5f})".format(coding_genes_with_frameshift_fp_or_fn, coding_genes_with_frameshift_fp_or_fn / coding_genes_with_frameshift))



def get_coverage_state(wholly, partially, nott, chrom, gene):
    if gene in wholly[chrom]:
        return WHOLLY
    if gene in partially[chrom]:
        return PARTIALLY
    if gene in nott[chrom]:
        return NOT_COVERED
    return None


def print_coverage_information(coverage_breakdown, gene_lengths=None, outfile=None, first="HICONF", second="PHASED", third="SWITCH"):
    total_genes = 0
    for x in coverage_breakdown.values():
        for y in x.values():
            for z in y.values():
                total_genes += z

    for highconf in [WHOLLY, PARTIALLY, NOT_COVERED]:
        hicf_printed = False
        hicf_total = 0
        for phased in [WHOLLY, PARTIALLY, NOT_COVERED]:
            phsd_printed = False
            phsd_total = 0
            for swichd in [WHOLLY, PARTIALLY, NOT_COVERED]:
                value = coverage_breakdown[highconf][phased][swichd]
                if value == 0: continue
                hicf_str = "{} {}  ".format(first, highconf) if not hicf_printed else "" #26
                phsd_str = "{} {}  ".format(second, phased) if not phsd_printed else ""
                swch_str = "{} {}  ".format(third, swichd)
                lengths_str = ""
                if gene_lengths is not None  and len(gene_lengths[highconf][phased][swichd]) > 0:
                    lengths = gene_lengths[highconf][phased][swichd]
                    lengths_str = "{:6d} | {:6d} | {:6d}".format(int(np.quantile(lengths, .25)), int(np.quantile(lengths, .5)), int(np.quantile(lengths, .75)))
                log("\t{:26s} {:26s} {:26s} {:5d}  {:5.1f}%  {:5.1f}%  {}".format( hicf_str, phsd_str, swch_str, value,
                    100.0 * value / sum(coverage_breakdown[highconf][phased].values()), 100.0 * value / total_genes, lengths_str))
                hicf_printed = True
                phsd_printed = True
                hicf_total += value
                phsd_total += value
            if phsd_total != 0 and len(list(filter(lambda x: x > 0, coverage_breakdown[highconf][phased].values()))) > 1:
                lengths_str = ""
                if gene_lengths is not None:
                    # lengths = [item for item in gene_lengths[highconf][phased]]
                    lengths = []
                    for x in gene_lengths[highconf][phased].values():
                        for y in x:
                            lengths.append(y)
                    lengths_str = "{:6d} | {:6d} | {:6d}".format(int(np.quantile(lengths, .25)), int(np.quantile(lengths, .5)), int(np.quantile(lengths, .75)))
                log("\t{:26s} {:26s} {:26s} {:5d}  {:6s}  {:5.1f}%  {}".format("", "", "       TOTAL", phsd_total, "", 100.0 * phsd_total / total_genes, lengths_str))
        if hicf_total != 0:
            lengths_str = ""
            if gene_lengths is not None:
                # lengths = [item for sublist in gene_lengths[highconf] for item in sublist]
                lengths = []
                for x in gene_lengths[highconf].values():
                    for y in x.values():
                        for z in y:
                            lengths.append(z)
                lengths_str = "{:6d} | {:6d} | {:6d}".format(int(np.quantile(lengths, .25)), int(np.quantile(lengths, .5)), int(np.quantile(lengths, .75)))
            log("\t{:26s} {:26s} {:26s} {:5d}  {:6s}  {:5.1f}%  {}".format("", "       TOTAL", "", hicf_total, "", 100.0 * hicf_total / total_genes, lengths_str))


def write_coverage_information(coverage_breakdown, gene_lengths, outfile, file_id):

    modified_coverage_breakdown = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: 0)))
    modified_coverage_breakdown_gene_lengths = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: list())))

    for highconf in [WHOLLY, PARTIALLY, NOT_COVERED]:
        for phased in [WHOLLY, PARTIALLY, NOT_COVERED]:
            for swichd in [WHOLLY, PARTIALLY, NOT_COVERED]:
                value = coverage_breakdown[highconf][phased][swichd]
                lengths = gene_lengths[highconf][phased][swichd]
                ps_m = NO_ERROR if phased == NOT_COVERED or highconf == NOT_COVERED else ERROR
                sw_m = NO_ERROR if swichd == NOT_COVERED else ERROR
                modified_coverage_breakdown[highconf][ps_m][sw_m] += value
                for l in lengths:
                    modified_coverage_breakdown_gene_lengths[highconf][ps_m][sw_m].append(l)

    for hiconf in [WHOLLY, PARTIALLY, NOT_COVERED]:
        for phased in [ERROR, NO_ERROR]:
            for swichd in [ERROR, NO_ERROR]:
                lengths = modified_coverage_breakdown_gene_lengths[hiconf][phased][swichd]
                data = [
                    hiconf,
                    phased,
                    swichd,
                    file_id,
                    str(modified_coverage_breakdown[hiconf][phased][swichd]),
                    "--" if len(lengths) == 0 else str(int(np.quantile(lengths, .25))),
                    "--" if len(lengths) == 0 else str(int(np.quantile(lengths, .5))),
                    "--" if len(lengths) == 0 else str(int(np.quantile(lengths, .75)))
                ]
                print("\t".join(data), file=outfile)
                # print("\t".join(data))





def plot_coverage_donut(coverage_breakdown, output_base, figsize=(3.5, 3.5)):

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
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=figsize)
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
        plt.savefig("{}.coverage_donut.pdf".format(output_base), format='pdf', dpi=300)
    plt.show()


def plot_coverage_grouped_bar(coverage_breakdown, output_base, figsize=(3.5, 3.5)):

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

    a, b, c = [plt.cm.Purples, plt.cm.Blues, plt.cm.Greens]

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=figsize)
    rect1 = ax1.bar(r1, bars1, color=[a(0.6), b(0.6), c(0.6)], hatch='', width=barWidth, edgecolor='white', label='Total Genes')
    rect2 = ax1.bar(r2, bars2, color=[a(0.45), b(0.45), c(0.45)], hatch='',  width=barWidth, edgecolor='white', label='Wholly Phased')
    rect3 = ax1.bar(r3, bars3, color=[a(0.3), b(0.3), c(0.3)], hatch='', width=barWidth, edgecolor='white', label='No Switch Errors')

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], ['Wholly\nin High\nConfidence', 'Partly\nin High\nConfidence', 'Not\nin High\nConfidence'])

    def autolabel(rects, percentage_denominator):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{} ({:4.1f}%)'.format(height, 100.0*height/max(1,percentage_denominator[i])),
                        xy=(rect.get_x() + rect.get_width() / 2, height - .25 * max(bars1)),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90,
                        ha='center', va='bottom')#, weight='bold')

    autolabel(rect1, [total, total, total])
    autolabel(rect2, bars1)
    autolabel(rect3, bars2)

    # Create legend & Show graphic
    plt.legend()
    plt.tight_layout()
    ax1.set_ylabel("Gene Count")
    ax1.set_xlabel("Coverage Classification")

    # show it
    if output_base is not None:
        plt.savefig("{}.coverage_grouped_bar.pdf".format(output_base), format='pdf', dpi=300)
    plt.show()


def plot_stratification_grouped_bar(stratification_breakdown, output_base, figsize=(3.5, 3.5)):

    # data we plot
    wholeP = 0
    wholeP_notS = 0
    wholeP_notS_noE = 0
    partP = 0
    partP_notS = 0
    partP_notS_noE = 0
    notP = 0
    notP_notS = 0
    notP_notS_noE = 0

    # fill in the data
    for phased in [WHOLLY, PARTIALLY, NOT_COVERED]:
        for switched in [WHOLLY, PARTIALLY, NOT_COVERED]:
            for errored in [WHOLLY, NOT_COVERED]:
                value = stratification_breakdown[phased][switched][errored]
                if value == 0: continue
                if phased == WHOLLY:
                    wholeP += value
                    if switched == NOT_COVERED:
                        wholeP_notS += value
                        if errored == NOT_COVERED:
                            wholeP_notS_noE += value
                elif phased == PARTIALLY:
                    partP += value
                    if switched == NOT_COVERED:
                        partP_notS += value
                        if errored == NOT_COVERED:
                            partP_notS_noE += value
                elif phased == NOT_COVERED:
                    notP += value
                    if switched == NOT_COVERED:
                        notP_notS += value
                        if errored == NOT_COVERED:
                            notP_notS_noE += value
    total = wholeP + partP + notP

    # set width of bar
    barWidth = 0.25

    # set height of bar
    bars1 = [wholeP, partP, notP]
    bars2 = [wholeP_notS, partP_notS, notP_notS]
    bars3 = [wholeP_notS_noE, partP_notS_noE, notP_notS_noE]

    a, b, c = [plt.cm.Purples, plt.cm.Blues, plt.cm.Greens]

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=figsize)
    rect1 = ax1.bar(r1, bars1, color=[a(0.6), b(0.6), c(0.6)], hatch='', width=barWidth, edgecolor='white', label='Total Genes')
    rect2 = ax1.bar(r2, bars2, color=[a(0.45), b(0.45), c(0.45)], hatch='',  width=barWidth, edgecolor='white', label='No Switch Errors')
    rect3 = ax1.bar(r3, bars3, color=[a(0.3), b(0.3), c(0.3)], hatch='', width=barWidth, edgecolor='white', label='No CDS Variant Errors')

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], ['Wholly\nPhased', 'Partly\nPhased', 'Not\nPhased'])

    def autolabel(rects, percentage_denominator):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{} ({:4.1f}%)'.format(height, 100.0*height/max(1,percentage_denominator[i])),
                        xy=(rect.get_x() + rect.get_width() / 2, height - .25 * max(bars1)),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90,
                        ha='center', va='bottom')#, weight='bold')

    autolabel(rect1, [total, total, total])
    autolabel(rect2, bars1)
    autolabel(rect3, bars2)

    # Create legend & Show graphic
    plt.legend()
    plt.tight_layout()
    ax1.set_ylabel("Gene Count")
    ax1.set_xlabel("Coverage Classification")

    # show it
    if output_base is not None:
        plt.savefig("{}.stratification_grouped_bar.pdf".format(output_base), format='pdf', dpi=300)
    plt.show()


def plot_manual_unerrored_grouped_bar(output_base, figsize=(3.5, 3.5), ont=True):
    if ont:
        total = 3793
        total_phased = 3540
        total_phased_no_switch = 3502
        total_cds_no_error = 7 + 30 + 3481 + 37 + 2 + 213
        total_exon_no_error = 6 + 27 + 3121 + 29 + 2 + 199
        total_gene_no_error = 5 + 10 + 1738 + 6 + 1 + 124
        total_phased_no_switch_cds_no_error = 3481
        total_phased_no_switch_exon_no_error = 3121
        total_phased_no_switch_gene_no_error = 1738
    else:
        total = 2500 + 365 + 928
        total_phased = 2500
        total_phased_no_switch = 2476
        total_cds_no_error = 928 + 365 + 2474 + 19 + 5
        total_exon_no_error = 916 + 1 + 354 + 5 + 2446 + 18 + 5
        total_gene_no_error = 740 + 1 + 190 + 2 + 2086 + 13 + 5
        total_phased_no_switch_cds_no_error = 2474
        total_phased_no_switch_exon_no_error = 2446
        total_phased_no_switch_gene_no_error = 2086

    # set width of bar
    barWidth = 0.25

    # set height of bar
    bars1 = [total, total_phased, total_phased_no_switch]
    bars2 = [total_cds_no_error, total_exon_no_error, total_gene_no_error]
    bars3 = [total_phased_no_switch_cds_no_error, total_phased_no_switch_exon_no_error, total_phased_no_switch_gene_no_error]

    r1_col = plt.cm.Reds if not ont else plt.cm.Oranges
    r2_col = plt.cm.YlOrBr if ont else plt.cm.Greens
    r3_col = plt.cm.Blues if ont else plt.cm.Purples

    r1_col = plt.cm.Reds if ont else plt.cm.YlOrBr
    r2_col = plt.cm.PuRd if ont else plt.cm.Greens
    r3_col = plt.cm.Purples if ont else plt.cm.Blues

    r1_colors = [r1_col(0.55), r1_col(0.45), r1_col(0.35)] #if ont else [r1_col(0.6), r1_col(0.5), r1_col(0.4)]
    r2_colors = [r2_col(0.6), r2_col(0.5), r2_col(0.4)] if ont else [r2_col(0.6), r2_col(0.45), r2_col(0.3)]
    r3_colors = [r3_col(0.6), r3_col(0.45), r3_col(0.3)] #if not ont else [r3_col(0.45), r3_col(0.55), r3_col(0.65)]

    # Set position of bar on X axis
    r1 = [1, 1+barWidth, 1+2*barWidth]
    r2 = [2, 2+barWidth, 2+2*barWidth]
    r3 = [3, 3+barWidth, 3+2*barWidth]

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=figsize)
    rect1 = ax1.bar(r1, bars1, color=r1_colors, hatch='', width=barWidth, edgecolor='white')
    rect2 = ax1.bar(r2, bars2, color=r2_colors, hatch='', width=barWidth, edgecolor='white')
    rect3 = ax1.bar(r3, bars3, color=r3_colors, hatch='', width=barWidth, edgecolor='white')

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r + barWidth for r in range(1,4)], ['Gene\nPhasing', 'Genes\nWithout\nSNP or\nINDEL\nErrors', 'Genes\nWithout\nSNP, INDEL,\nor Switch\nErrors'])
    ylim_min = 0
    ylim_max = total
    def autolabel(rects, percentage_denominator, labels):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{} ({:4.1f}%)'.format(height, 100.0 * height / max(1, percentage_denominator[i])),
                         xy=(rect.get_x() + rect.get_width() / 2, height - .25 * max(bars1)),
                         xytext=(0, 3),  # 3 points vertical offset
                         textcoords="offset points", rotation=90,
                         ha='center', va='bottom')  # , weight='bold')
            ax1.annotate('{}'.format(labels[i]),
                        xy=(rect.get_x() + rect.get_width() / 2, ylim_min + .025 * (ylim_max - ylim_min)),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90, color="black",
                        ha='center', va='bottom')#, weight='bold')

    autolabel(rect1, [total, total, total], ["All Genes", "Wholly Phased", "No Switch"])
    autolabel(rect2, [total, total, total], ["CDS", "Exon", "Gene"])
    autolabel(rect3, [total, total, total], ["CDS", "Exon", "Gene"])

    # Create legend & Show graphic
    ax1.set_ylabel("Gene Count")
    ax1.set_xlabel("Classification")
    plt.tight_layout()

    # show it
    if output_base is not None:
        plt.savefig("{}.{}.errored_grouped_bar.pdf".format(output_base, "ont" if ont else "ccs"), format='pdf', dpi=300)
    plt.show()



def plot_manual_coverage_bar(output_base, figsize=(3.5, 3.5)):
    # data we plot
    total_genes_38 = 60656
    total_genes_37 = 62438

    samples = ["Total\nGenes\nGRCh38", "HG003\nWholly\nPhased\n(GRCh38)", "HG004\nWholly\nPhased\n(GRCh38)",
               "Total\nGenes\nGRCh37", "HG001\nWholly\nPhased\n(GRCh37)", "HG005\nWholly\nPhased\n(GRCh37)",
               "HG006\nWholly\nPhased\n(GRCh37)", "HG007\nWholly\nPhased\n(GRCh37)"]
    samples = ["Total\nGenes\nGRCh38", "HG003\n(GRCh38)", "HG004\n(GRCh38)", "Total\nGenes\nGRCh37", "HG001\n(GRCh37)",
               "HG005\n(GRCh37)", "HG006\n(GRCh37)", "HG007\n(GRCh37)"]
    coverage = [60656, 53817, 55234, 62438, 57175, 53150, 53112, 54116]

    def autolabel(rects, percentage_denominator):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{} ({:4.1f}%)'.format(height, 100.0 * height / max(1, percentage_denominator[i])),
                         xy=(rect.get_x() + rect.get_width() / 2, height - .25 * max(coverage)),
                         xytext=(0, 3),  # 3 points vertical offset
                         textcoords="offset points", rotation=90,
                         ha='center', va='bottom')  # , weight='bold')

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=figsize)
    for i in range(len(samples)):
        rect = ax1.bar([i], [coverage[i]],
                       color=(plt.cm.Oranges if i < 3 else plt.cm.Reds)(0.6 if i in (0, 3) else .45),
                       hatch='' if i == 0 else '', edgecolor='white')
        autolabel(rect, [total_genes_38 if i < 3 else total_genes_37])

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r for r in range(len(samples))], samples)
    ax1.set_ylabel("Gene Count")
    ax1.set_xlabel("Reference / Sample")

    # show it
    plt.tight_layout()
    if output_base is not None:
        plt.savefig("{}.coverage_grouped_bar_manual.pdf".format(output_base), format='pdf', dpi=300)
    plt.show()


def plot_manual_accuracy_bar(output_base, plot_snp, figsize=(3.5, 3.5)):

    # data we plot
    GENOME = "Genome"
    ALL_GENES = "All Genes"
    PC_GENES = "Protein\nCoding\nGenes"
    PC_SEQUENCE = "Protein\nCoding\nSequence\n(CDS)"
    PC_FRAMESHIFT = "CDS\nFrameshift"
    CCS = "PacBio-HiFi"
    ONT = "Nanopore"
    data_snp = {
        GENOME: {
            ONT : [0.996314, 0.998169],
            CCS : [0.998062, 0.999391]},
        ALL_GENES: {
            ONT : [0.996481, 0.998097],
            CCS : [0.998197, 0.999384]},
        PC_GENES: {
            ONT : [0.996751, 0.99799],
            CCS : [0.998441, 0.999382]},
        PC_SEQUENCE: {
            ONT : [0.997887, 0.998641],
            CCS : [0.998994, 0.999446]},
    }
    data_indel = {
        GENOME: {
            ONT : [0.878512, 0.60077],
            CCS : [0.92602, 0.948736]},
        ALL_GENES: {
            ONT : [0.877943, 0.595042],
            CCS : [0.922887, 0.947847]},
        PC_GENES: {
            ONT : [0.876731, 0.584914],
            CCS : [0.918991, 0.946128]},
        PC_SEQUENCE: {
            ONT : [0.926893, 0.799544],
            CCS : [0.909278, 0.984055]},
        PC_FRAMESHIFT: {
            ONT : [0.85161, 0.80488],
            CCS : [0.82292, 0.98137]},
    }


    # set width of bar
    barWidth = 0.15

    # set height of bar
    # bars1 = [wholeH, wholeH_wholeP, wholeH_wholeP_notS]
    # bars2 = [partH, partH_wholeP, partH_wholeP_notS]
    # bars3 = [notH, notH_wholeP, 0]
    data = data_snp if plot_snp else data_indel
    bar_ordering = [GENOME, ALL_GENES, PC_GENES, PC_SEQUENCE] #+ ([] if plot_snp else [PC_FRAMESHIFT])
    bar_ont_precision = [data[x][ONT][0] for x in bar_ordering]
    bar_ont_recall = [data[x][ONT][1] for x in bar_ordering]
    bar_ont_f1 = [2 * data[x][ONT][0] * data[x][ONT][1] / (data[x][ONT][0] + data[x][ONT][1]) for x in bar_ordering]
    bar_ccs_precision = [data[x][CCS][0] for x in bar_ordering]
    bar_ccs_recall = [data[x][CCS][1] for x in bar_ordering]
    bar_ccs_f1 = [2 * data[x][CCS][0] * data[x][CCS][1] / (data[x][CCS][0] + data[x][CCS][1]) for x in bar_ordering]

    C_ONT = plt.cm.Reds if plot_snp else plt.cm.Oranges
    C_CCS = plt.cm.Blues if plot_snp else plt.cm.Purples
    S_P = .6
    S_R = .45
    S_F = .3

    # Set position of bar on X axis
    r_ont_precision = np.arange(len(bar_ordering))
    r_ont_recall = [x + barWidth for x in r_ont_precision]
    r_ont_f1 = [x + barWidth for x in r_ont_recall]
    r_ccs_precision = [x + barWidth for x in r_ont_f1]
    r_ccs_recall = [x + barWidth for x in r_ccs_precision]
    r_ccs_f1 = [x + barWidth for x in r_ccs_recall]

    # Make the plot
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=figsize)
    rect_op = ax1.bar(r_ont_precision, bar_ont_precision, color=C_ONT(S_P), hatch='', width=barWidth, edgecolor='white', label=ONT+' Precision')
    rect_or = ax1.bar(r_ont_recall, bar_ont_recall, color=C_ONT(S_R), hatch='',  width=barWidth, edgecolor='white', label=ONT+' Recall')
    rect_of = ax1.bar(r_ont_f1, bar_ont_f1, color=C_ONT(S_F), hatch='', width=barWidth, edgecolor='white', label=ONT+' F1-Score')
    rect_cp = ax1.bar(r_ccs_precision, bar_ccs_precision, color=C_CCS(S_P), hatch='', width=barWidth, edgecolor='white', label=CCS + ' Precision')
    rect_cr = ax1.bar(r_ccs_recall, bar_ccs_recall, color=C_CCS(S_R), hatch='',  width=barWidth, edgecolor='white', label=CCS + ' Recall')
    rect_cf = ax1.bar(r_ccs_f1, bar_ccs_f1, color=C_CCS(S_F), hatch='', width=barWidth, edgecolor='white', label=CCS + ' F1-Score')

    # Add xticks on the middle of the group bars
    # plt.xlabel('', fontweight='bold')
    plt.xticks([r + 2.5*barWidth for r in range(len(bar_ordering))], bar_ordering)
    ylim_min = .99 if plot_snp else .0
    ylim_max = 1.0

    def autolabel(rects, msg):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for i, rect in enumerate(rects):
            height = rect.get_height()
            if height == 0: continue
            ax1.annotate('{:5.4f}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height - .125 * (ylim_max - ylim_min)),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90, color="black",
                        ha='center', va='bottom')#, weight='bold')
            ax1.annotate('{}'.format(msg),
                        xy=(rect.get_x() + rect.get_width() / 2, ylim_min + .025 * (ylim_max - ylim_min)),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", rotation=90, color="black",
                        ha='center', va='bottom')#, weight='bold')
    autolabel(rect_op, ONT + " Precision")
    autolabel(rect_or, ONT + " Recall")
    autolabel(rect_of, ONT + " F1-Score")
    autolabel(rect_cp, CCS + " Precision")
    autolabel(rect_cr, CCS + " Recall")
    autolabel(rect_cf, CCS + " F1-Score")

    # Create legend & Show graphic
    # plt.legend()
    # plt.tight_layout()
    ax1.set_ylabel("Accuracy")
    ax1.set_xlabel("Gene Stratification")
    ax1.set_ylim(ylim_min, ylim_max)

    # show it
    if output_base is not None:
        plt.savefig("{}.gene_accuracies_{}.pdf".format(output_base, "snp" if plot_snp else "indel"), format='pdf', dpi=300)
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
    # plot_accuracy_grouped_bar_manual_for_files(output_base, True)
    # plot_accuracy_grouped_bar_manual_for_files(output_base, False)
    # plot_coverage_grouped_bar_manual_for_files(output_base)
    # plot_manual_unerrored_grouped_bar(output_base)
    # plot_manual_unerrored_grouped_bar(output_base, ont=False)
    # sys.exit()

    # get data
    restricted_gene_ids = None
    if args.gene_ids_of_interest is not None:
        restricted_gene_ids = set()
        with open(args.gene_ids_of_interest) as gioi:
            for line in gioi:
                restricted_gene_ids.add(line.strip())

    phasesets = get_phasesets_from_vcf(args.phased_vcf, args)
    genes, switches, high_conf_regions = get_bed_regions(args, restricted_gene_ids)
    genes_wholly_hiconf, genes_partially_hiconf, genes_not_hiconf = get_genes_stratified_by_bed(genes, high_conf_regions, "HIGHCONF", True)
    genes_wholly_phased, genes_partially_phased, genes_not_phased = get_genes_stratified_by_bed(genes, phasesets, "PHASESET")
    genes_wholly_swichd, genes_partially_swichd, genes_not_swichd = get_genes_stratified_by_bed(genes, switches, "SWITCHES")
    # get_precision_recall_for_genes(genes, args)

    # info over all genes
    coverage_breakdown = collections.defaultdict(lambda : collections.defaultdict(lambda : collections.defaultdict(lambda : 0)))
    coverage_breakdown_gene_lengths = collections.defaultdict(lambda : collections.defaultdict(lambda : collections.defaultdict(lambda : list())))
    highconf_variants = collections.defaultdict(lambda : collections.defaultdict(lambda : 0))
    total_genes = 0

    strat_coverage_breakdown, gene_ids_above_threshold, no_error_gene_ids_above_threshold = None, None, None
    strat_coverage_breakdown_gene_lengths = None
    if args.stratification_tsv is not None:
        gene_strat_data, gene_ids_above_threshold, no_error_gene_ids_above_threshold = get_gene_stratification_analysis(
            args.stratification_tsv)
        # (only for genes with strat coverage)
        # VCF_PHASED -> SWITCHES -> HAS_ERROR (wholly/not)
        strat_coverage_breakdown = collections.defaultdict(
            lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: 0)))
        strat_coverage_breakdown_gene_lengths = collections.defaultdict(
            lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: list())))

    # get coverage data
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
                if dest is None: continue
                dest.append("{}:{}-{}".format(bed.chrom, bed.end, bed.start))

            # get precision/recall
            total_variants = gene.get_total_variants()
            snp_tp, snp_fp, snp_fn = gene.snp_tp, gene.snp_fp, gene.snp_fn
            indel_tp, indel_fp, indel_fn = gene.indel_tp, gene.indel_fp, gene.indel_fn

            # for whole stats
            total_genes += 1
            coverage_breakdown[hiconf_cov][phased_cov][swichd_cov] += 1
            coverage_breakdown_gene_lengths[hiconf_cov][phased_cov][swichd_cov].append(gene.end - gene.start)
            if hiconf_cov == WHOLLY:
                highconf_variants[SNP][TP] += snp_tp
                highconf_variants[SNP][FP] += snp_fp
                highconf_variants[SNP][FN] += snp_fn
                highconf_variants[INDEL][TP] += indel_tp
                highconf_variants[INDEL][FP] += indel_fp
                highconf_variants[INDEL][FN] += indel_fn

            # stratification info (if available)
            if strat_coverage_breakdown is not None:
                gene_id = gene.get_gene_id()
                if gene_id in gene_ids_above_threshold:
                    error_coverage = NOT_COVERED if gene_id in no_error_gene_ids_above_threshold else WHOLLY
                    strat_coverage_breakdown[phased_cov][swichd_cov][error_coverage] += 1
                    strat_coverage_breakdown_gene_lengths[phased_cov][swichd_cov][error_coverage].append(gene.end - gene.start)


    # log("")
    # print_highconf_coverage_histogram(genes)
    # log("")
    # print_frameshift_rates_above_highconf_coverage_threhsold(genes)
    # log("")
    # print_frameshift_rates_above_highconf_coverage_threhsold(genes, 0.0)
    log("")
    log("Coverage Breakdown:")
    print_coverage_information(coverage_breakdown, coverage_breakdown_gene_lengths, first="HICONF", second="PHASED", third="SWITCH")
    if strat_coverage_breakdown is not None:
        log("")
        log("Stratification Breakdown:")
        print_coverage_information(strat_coverage_breakdown, strat_coverage_breakdown_gene_lengths, first="PHASED", second="SWITCH", third="ERRORD")
        strat_basename = os.path.basename(args.stratification_tsv)
        strat_id = strat_basename
        if "gene" in strat_id:
            strat_id = "gene"
        elif "exon" in strat_id:
            strat_id = "exon"
        elif "cds" in strat_id:
            strat_id = "cds"
        with open("{}.{}.strat_summary.tsv".format(output_base, strat_id), 'w') as outfile:
            write_coverage_information(strat_coverage_breakdown, strat_coverage_breakdown_gene_lengths, outfile, strat_id)

    # plot_coverage_donut(coverage_breakdown, output_base)
    plot_stratification_grouped_bar(strat_coverage_breakdown, output_base)
    plot_coverage_grouped_bar(coverage_breakdown, output_base)

    snp_tp = highconf_variants[SNP][TP]
    snp_fp = highconf_variants[SNP][FP]
    snp_fn = highconf_variants[SNP][FN]
    indel_tp = highconf_variants[INDEL][TP]
    indel_fp = highconf_variants[INDEL][FP]
    indel_fn = highconf_variants[INDEL][FN]
    if sum([snp_tp, snp_fp, snp_fn, indel_tp, indel_fn, indel_fp]) != 0:
        log("")
        log("Wholly Covered HiConf Gene Variant Stats:")
        log("\tSNP Precision:   {:.5f} {:8}/{}+{}".format(snp_tp/max(1, snp_tp+snp_fp), snp_tp, snp_tp, snp_fp))
        log("\tSNP Recall:      {:.5f} {:8}/{}+{}".format(snp_tp/max(1, snp_tp+snp_fn), snp_tp, snp_tp, snp_fn))
        log("\tSNP F1:          {:.5f} {:8}/({}+({}+{})/2)".format(snp_tp/max(1, snp_tp+(snp_fn+snp_fp)/2), snp_tp, snp_tp, snp_fn,snp_fp))
        log("\tINDEL Precision: {:.5f} {:8}/{}+{}".format(indel_tp/max(1, indel_tp+indel_fp), indel_tp, indel_tp, indel_fp))
        log("\tINDEL Recall:    {:.5f} {:8}/{}+{}".format(indel_tp/max(1, indel_tp+indel_fn), indel_tp, indel_tp, indel_fn))
        log("\tINDEL F1:        {:.5f} {:8}/({}+({}+{})/2)".format(indel_tp/max(1, indel_tp+(indel_fn+indel_fp)/2), indel_tp, indel_tp, indel_fn,indel_fp))


    pass









if __name__ == "__main__":
    main()

