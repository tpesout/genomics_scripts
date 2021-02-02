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

TYPE_GENE="gene"
TYPE_EXON="exon"
TYPE_CDS="CDS"
TYPE_START_CODON="start_codon"
TYPE_STOP_CODON="stop_codon"
TYPE__PROTEIN_OCDING="protein_coding"

INFO_GENE_TYPE="gene_type"
INFO_GENE_ID="gene_id"
INFO_GENE_NAME="gene_name"


def parse_args(args = None):
    parser = argparse.ArgumentParser("Produces stratification beds from gencode")
    parser.add_argument('--gencode_gtf', '-g', dest='gencode_gtf', default=None, required=True, type=str,
                       help='GTF holding GENCODE annotations')
    parser.add_argument('--output_dir', '-d', dest='output_dir', default=None, required=False, type=str,
                       help='Output directory')
    parser.add_argument('--output_base', '-o', dest='output_base', default=None, required=False, type=str,
                       help='Output base file name')
    parser.add_argument('--output_base_from_input', '-O', dest='output_base_from_input', default=False, required=False, action='store_true',
                       help='Output directory and file name will be determined from --gencode_gtf parameter')
    parser.add_argument('--protein_coding_only', '-p', dest='protein_coding_only', default=False, required=False, action='store_true',
                       help='Only')
    parser.add_argument('--strip_chr_names', '-C', dest='strip_chr_names', default=False, required=False, action='store_true',
                       help='Remove "chr" from contig names (for grch37)')
    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg, end="\n"):
    print(msg, file=sys.stderr, end=end)


def make_dir(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    if not os.path.isdir(dir):
        log("Could not make directory: {}".format(dir))
        sys.exit(1)


def get_gene_info(gene_data):
    if type(gene_data) == str:
        gene_data = gene_data.split("\t")
    assert(gene_data[GENE_GTF_TYPE_IDX] == TYPE_GENE)
    return {stat.split()[0].strip(): stat.split()[1].strip().strip("\"") for stat in gene_data[GENE_GTF_INFO_IDX].strip().strip(";").split(";")}


def is_gene_line(gene_data):
    if type(gene_data) == str:
        gene_data = gene_data.split("\t")
    return gene_data[GENE_GTF_TYPE_IDX] == TYPE_GENE


def get_gene_filename(gene_info, type):
    return "genes/{}--{}--{}.bed".format(gene_info[INFO_GENE_NAME], gene_info[INFO_GENE_ID], type)


def get_stratification_entry(filename):
    return "{}\t{}\n".format(os.path.basename(filename[:-4] if filename.endswith(".bed") else filename).strip(), filename.strip())


def get_bed_record(record_parts, gene_info, args):
    return "\t".join(
        [(record_parts[i] if i != GENE_GTF_CHROM_IDX or not args.strip_chr_names else record_parts[i].lstrip("chr")).strip()
            for i in [GENE_GTF_CHROM_IDX, GENE_GTF_START_POS_IDX, GENE_GTF_END_POS_IDX]] +
        ["{};{};{}".format(gene_info[INFO_GENE_NAME], gene_info[INFO_GENE_ID], gene_info[INFO_GENE_TYPE])]
    ) + "\n"


def main():
    args = parse_args()

    files = [args.gencode_gtf]
    missing_files = list(filter(lambda x: x is not None and not os.path.isfile, files))
    if len(missing_files) != 0:
        log("Could not find specified files: {}".format(missing_files))
        sys.exit(1)

    # prep for writing output
    output_base = args.output_base
    if output_base is None and args.output_base_from_input:
        output_base = os.path.basename(args.gencode_gtf).rstrip(".gz").rstrip(".gtf")
    if output_base is None:
        log("Need output filename configuration")
        sys.exit(1)

    # output directories
    output_dir = args.output_dir
    if output_dir is None and args.output_base_from_input:
        output_dir = output_base + ".gene_stratification"
    if output_dir is None:
        log("Need output directory configuration")
        sys.exit(1)
    output_genes_dir = os.path.join(output_dir, "genes")
    output_genome_dir = os.path.join(output_dir, "genome")
    make_dir(output_dir)
    make_dir(output_genes_dir)
    make_dir(output_genome_dir)

    # prep for writing
    stratification_filename = os.path.join(output_dir, "{}.stratification.tsv".format(output_base))
    all_gene_filename = os.path.join("genome", "all_genes.bed")
    all_exon_filename = os.path.join("genome", "all_exons.bed")
    all_cds_filename = os.path.join("genome", "all_cds.bed")
    stratification_file = None
    all_gene_file = None
    all_exon_file = None
    all_cds_file = None
    log("Writing stratifications to {}".format(stratification_filename))
    
    # tracking
    total_gencode_records = 0
    total_genes_written = 0
    total_entries_written = 0
    total_stratifications = 0

    # actually write the file
    try:
        stratification_file = open(stratification_filename, 'w')
        all_gene_file = open(os.path.join(output_dir, all_gene_filename), 'w')
        all_exon_file = open(os.path.join(output_dir, all_exon_filename), 'w')
        all_cds_file = open(os.path.join(output_dir, all_cds_filename), 'w')

        all_gene_file.write("# contains all gene records from {}\n".format(args.gencode_gtf))
        all_exon_file.write("# contains all exon records from {}\n".format(args.gencode_gtf))
        all_cds_file.write("# contains all CDS records from {}\n".format(args.gencode_gtf))
        
        # update strat file
        stratification_file.write(get_stratification_entry(all_gene_filename))
        stratification_file.write(get_stratification_entry(all_exon_filename))
        stratification_file.write(get_stratification_entry(all_cds_filename))

        # class for handling gene info
        class GeneInfo():
            def __init__(self, gene_parts):
                if type(gene_parts) == str:
                    gene_parts = gene_parts.split("\t")
                assert (gene_parts[GENE_GTF_TYPE_IDX] == TYPE_GENE)
                self.gene_info = get_gene_info(gene_parts)
                self.gene_bed_record = get_bed_record(gene_parts, self.gene_info, args)
                self.exon_entries = []
                self.cds_entries = []

            def save_entry(self, entry_parts):
                assert (entry_parts[GENE_GTF_TYPE_IDX] != TYPE_GENE)
                if entry_parts[GENE_GTF_TYPE_IDX] == TYPE_EXON:
                    self.exon_entries.append(get_bed_record(entry_parts, self.gene_info, args))
                elif entry_parts[GENE_GTF_TYPE_IDX] in (TYPE_CDS, TYPE_START_CODON, TYPE_STOP_CODON):
                    self.cds_entries.append(get_bed_record(entry_parts, self.gene_info, args))

            def write(self):
                # early exit for protein coding option
                if args.protein_coding_only and self.gene_info[INFO_GENE_TYPE] != TYPE__PROTEIN_OCDING:
                    return 0, 0, 0
                # genes
                gene_filename = get_gene_filename(self.gene_info, "gene")
                with open(os.path.join(output_dir, gene_filename), 'w') as fout:
                    fout.write("# contains a gene record from {}\n".format(args.gencode_gtf))
                    fout.write(self.gene_bed_record)
                all_gene_file.write(self.gene_bed_record)
                stratification_file.write(get_stratification_entry(gene_filename))
                gene_strats_written = 1
                gene_records_written = 1
                # exons
                if len(self.exon_entries) > 0:
                    exon_filename = get_gene_filename(self.gene_info, "exon")
                    with open(os.path.join(output_dir, exon_filename), 'w') as fout:
                        fout.write("# contains exon records from {}\n".format(args.gencode_gtf))
                        for exon_record in self.exon_entries:
                            fout.write(exon_record)
                            all_exon_file.write(exon_record)
                            gene_records_written += 1
                    stratification_file.write(get_stratification_entry(exon_filename))
                    gene_strats_written += 1
                # cds
                if len(self.cds_entries) > 0:
                    cds_filename = get_gene_filename(self.gene_info, "cds")
                    with open(os.path.join(output_dir, cds_filename), 'w') as fout:
                        fout.write("# contains CDS records from {}\n".format(args.gencode_gtf))
                        for cds_record in self.exon_entries:
                            fout.write(cds_record)
                            all_cds_file.write(cds_record)
                            gene_records_written += 1
                    stratification_file.write(get_stratification_entry(cds_filename))
                    gene_strats_written += 1
                return 1, gene_strats_written, gene_records_written

        # do it all
        current_gene = None
        log("Reading GENCODE annotations from {}".format(args.gencode_gtf))
        with open(args.gencode_gtf) as gencode_in:
            for line in gencode_in:
                if line.startswith("#"): continue
                total_gencode_records += 1
                line_parts = line.split("\t")
                
                if is_gene_line(line_parts):
                    if current_gene is not None:
                        genes_written, stratifications_written, entries_written = current_gene.write()
                        total_genes_written += genes_written
                        total_stratifications += stratifications_written
                        total_entries_written += entries_written
                    current_gene = GeneInfo(line_parts)
                else:
                    current_gene.save_entry(line_parts)
                    
            if current_gene is not None:
                genes_written, stratifications_written, entries_written = current_gene.write()
                total_genes_written += genes_written
                total_stratifications += stratifications_written
                total_entries_written += entries_written

    finally:
        if stratification_file is not None:
            stratification_file.close()
        if all_gene_file is not None:
            all_gene_file.close()
        if all_exon_file is not None:
            all_exon_file.close()
        if all_cds_file is not None:
            all_cds_file.close()

    log("From {} gencode entries wrote {} genes, {} stratifications, and {} entries".format(
        total_gencode_records, total_genes_written, total_stratifications, total_entries_written))








if __name__ == "__main__":
    main()

