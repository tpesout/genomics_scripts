
from __future__ import print_function
import sys
import subprocess

# read.0
# read.1
# read.2
# chr
# subs
# ptype

existed = set()

def assert_exists(loc):
    if not loc.startswith("s3://"): return True
    if loc in existed: return True
    if len(subprocess.check_output(['s3cmd', 'ls', loc])) == 0:
        print("DNE: {}".format(loc), file=sys.stderr)
        return False
    existed.add(loc)
    return True

UUID = "NA12878.grch38.{read0}.chr{chr}"
BAM = ""
CHR = "chr{chr}"
FA = "s3://margin-phase/fasta/hg38.chr{chr}.fa"
# PARAMS = "s3://margin-phase/params/pecan/params.{read1}.{sub}.{ptype}.pec_hmm_2.json"
PARAMS = "s3://margin-phase/params/final/params.pacbio.181009.json"

reads = [
    ['pb', 'pb', 's3://margin-phase/bam/grch38/NA12878/pb/NA12878.grch38.pb.mm.q30pri.chr{chr}.bam'],
    # ['np2-q0', 'np', 's3://margin-phase/bam/nanopore2/NA12878.hg38.np2.{chr}.bam'],
    # ['pb-q0', 'pb', 's3://margin-phase/bam/realigned/NA12878.hg38.pb.mm.{chr}.bam'],
    # ['np2-q30', 'np', 's3://margin-phase/bam/nanopore2.q30/NA12878.hg38.np2.q30.{chr}.bam'],
    # ['pb-q30', 'pb', 's3://margin-phase/bam/pacbio.q30/NA12878.hg38.pb.mm.q30.{chr}.bam'],
#     ['np2-q30nsu', 'np', 's3://margin-phase/bam/nanopore2.q30nsu/NA12878.hg38.np2.q30nsu.{chr}.bam'],
#     ['pb-q30nsu', 'pb', 's3://margin-phase/bam/pacbio.q30nsu/NA12878.hg38.pb.mm.q30nsu.{chr}.bam'],
]

chrs = ["{}".format(x) for x in range(1,23)]
chrs.append("X")
# chrs.append("chrY")
# chrs = ['chr6']

# chrs = ["chr{}".format(x) for x in range(1,23)]
# chrs.append("chrX")
# chrs.append("chrY")
# chrs = ['chr6']

# subs = ['998', '9995']
# subs = ['9995']
subs = ['hs9993']
# subs = ['hs9993', '9995']

# ptypes = ['rp', 'plain']
# ptypes = ['rp-train1', 'nofilt']
# ptypes = ['rp', 'rp-train1', 'plain']
# ptypes = ['rp-nofilt']
# ptypes = ['rp-nofilt-train1.v1', 'rp-nofilt-train1.v2', 'rp-nofilt-train2.v1']
ptypes = ['rp-nofilt.v3']
# ptypes = ['plain']

have_everything = True
for read in reads:
    for chr in chrs:
        for sub in subs:
            for ptype in ptypes:
                uuid = UUID.format(read0=read[0], chr=chr, sub=sub, ptype=ptype)
                bam = read[2].format(read0=read[0], chr=chr)
                fa = FA.format(chr=chr)
                contig = CHR.format(chr=chr)
                params = PARAMS.format(read1=read[1], sub=sub, ptype=ptype)

                sample = [uuid, bam, contig, fa, params]
                print("\t".join(sample))
                for s in sample:
                    have_everything = have_everything and assert_exists(s)
                if not have_everything: sys.exit(1)



# read_errs = ['re98', 're998']
#
# have_everything = True
# unique_uuids = set()
# for read in reads:
#     for chr in chrs:
#         for sub in subs:
#             for read_err in read_errs:
#                 for ptype in ptypes:
#                     uuid = UUID.format(read0=read[0], chr=chr, sub=sub, ptype=ptype, read_err=read_err)
#                     bam = read[2].format(read0=read[0], chr=chr)
#                     chr = CHR.format(chr=chr)
#                     fa = FA.format(chr=chr)
#                     params = PARAMS.format(read1=read[1], sub=sub, read_err=read_err, ptype=ptype)
#                     vcf = VCF.format(chr=chr)
#
#                     sample = [uuid, bam, chr, fa, params, vcf]
#                     print("\t".join(sample))
#                     if uuid in unique_uuids:
#                         print("Duplicated UUID: {}".format(uuid), file=sys.stderr)
#                     unique_uuids.add(uuid)
#                     for s in sample:
#                         have_everything = have_everything and assert_exists(s)
#                     if not have_everything: sys.exit(1)