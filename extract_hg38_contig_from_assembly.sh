#!/bin/bash

export ASM=$1
export TC=8
export CHR=chrX
export HG38=/home/tpesout/data/reference/hg38.fa
#export HG38=/home/ubuntu/data/human/reference/hg38/hg38.fa

minimap2 -ax asm20 -t $TC $HG38 $ASM | samtools view -hb >$ASM.unsorted.bam
samtools sort -@ $TC $ASM.unsorted.bam | samtools view -hb >$ASM.bam
samtools index -@ $TC $ASM.bam
rm $ASM.unsorted.bam

samtools view -F 0x904 $ASM.bam $CHR | awk '{print $1}' | uniq | sort | uniq >$ASM.hg38_${CHR}_0x904_segments.txt
samtools view -F 0x904 -q 60 $ASM.bam $CHR | awk '{print $1}' | uniq | sort | uniq >$ASM.hg38_${CHR}_0x904q60_segments.txt

extract_fasta_segments.py -i $ASM -s $ASM.hg38_${CHR}_0x904_segments.txt -o $ASM.hg38_${CHR}_0x904.fa
extract_fasta_segments.py -i $ASM -s $ASM.hg38_${CHR}_0x904q60_segments.txt -o $ASM.hg38_${CHR}_0x904q60.fa
