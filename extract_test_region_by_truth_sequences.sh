#!/bin/bash

ASM=$1
ALN=$2
TRUTH_H1=$3
TRUTH_H2=$4
ID=$5

TC=16

if [[ -z "$5" ]] ; then
	echo "usage: $0 ASM ALN TRUTH_H1 TRUTH_H2 IDENTIFIER"
	exit 1
fi

echo -e "ASM: $ASM\nALN: $ALN\nT_1: $TRUTH_H1\nT_2: $TRUTH_H2"
sleep 1

set -o xtrace

minimap2 -ax asm20 -t $TC $ASM $TRUTH_H1 $TRUTH_H2 | samtools view -F 0x904 -q 60 >tmp.truth_to_asm.sam

TC=$(wc -l tmp.truth_to_asm.sam | awk '{print $1}')
if [[ "$TC" != "2" ]] ; then
	echo "Expected 2 truth alignments, got $TC"
	exit 1
fi

H1S=$(head -n1 tmp.truth_to_asm.sam | awk '{print $4}')
H2S=$(tail -n1 tmp.truth_to_asm.sam | awk '{print $4}')
H1L=$(head -n1 tmp.truth_to_asm.sam | awk '{print $10}' | wc -c)
H2L=$(tail -n1 tmp.truth_to_asm.sam | awk '{print $10}' | wc -c)
CHR=$(head -n1 tmp.truth_to_asm.sam | awk '{print $3}')

S_MIN=$(( $H1S > $H2S ? $H2S : $H1S ))
S_MAX=$(( $H1S < $H2S ? $H2S : $H1S ))
L_MIN=$(( $H1L > $H2L ? $H2L : $H1L ))
L_MAX=$(( $H1L < $H2L ? $H2L : $H1L ))

echo "Min start: $S_MIN (< $S_MAX)"
echo "Max length: $L_MAX (> $L_MIN)"

START=$(( $S_MIN - 1000 ))
END=$(( $S_MAX + $L_MAX + 1000 ))

echo "Chrom: $CHR"
echo "Start: $START"
echo "End:   $END"

REGION="$CHR:$START-$END"
RID="$CHR-$START-$END"

# extract reads
samtools view -hb -@ $TC $ALN $REGION >tmp.region.bam
samtools fastq -@ $TC tmp.region.bam >tmp.reads.fastq

# extract region
samtools faidx $ASM
samtools faidx $ASM $REGION | sed 's/:/_/' > $ID.reference.$RID.fasta

# align to region
minimap2 -ax map-ont -t $TC $ID.reference.$RID.fasta tmp.reads.fastq | samtools view -hb -F 0x904 -q 60 | samtools sort >$ID.align.$RID.bam
samtools index -@ $TC $ID.align.$RID.bam

# cleanup
rm tmp.*

echo "Fin."
