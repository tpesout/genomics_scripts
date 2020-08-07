#!/bin/bash

if [[ -z "$1" ]] ; then
	echo "usage: $0 VARIANT_POS"
	exit 1
fi

BAM=HG002.pFDA_r941g360.grch38_chr1.mvp-pdv-filtfix2-rle.haplotagged.bam
TRUTH_H1=/home/tpesout/work/marginPolish/phasingData/haplotagged/200527_HG002_trio-binned-miten/HG002_haplotype-HG003_ids.fixed.txt
TRUTH_H2=/home/tpesout/work/marginPolish/phasingData/haplotagged/200527_HG002_trio-binned-miten/HG002_haplotype-HG004_ids.fixed.txt
CHUNKS=HG002.pFDA_r941g360.grch38_chr1.mvp-pdv-filtfix2-rle.chunks.csv

CHUNK_ID=$(cat $CHUNKS | sed 's/,/\t/g' | awk '{printf "%d\t%s\n", NR - 1, $0}' | awk -v POS=$1 '{if ($5 <= POS && $6 >= POS ) {print $1 ".chr1-" $3 "-" $4} }')
FILE=$(ls phasingInfo/*$CHUNK_ID*)
REGION=$(echo $CHUNK_ID | sed 's/.*\.//' | sed 's/-/:/')

echo
echo "Pos:    $1"
echo "Region: $REGION"
echo "Chunk:  $CHUNK_ID"
echo "File:   $FILE"
echo

IMG=img/variant-$1-$2-C${CHUNK_ID}.png

if [[ -f $IMG ]] ; then
	echo "Image already exists:"
	echo
	echo $IMG
	sleep 5
fi

python3 ~/dev/marginPolish/scripts/plot_mvp_bubble_data.py -i $FILE -f $IMG -t "$1 $2 $CHUNK_ID" -1 $TRUTH_H1 -2 $TRUTH_H2 &

sleep 7

bam_stats.py -i $BAM -r $REGION -gld -s 5000 2>&1 | tee $(echo $IMG | sed 's/.png/.stats.txt/')

echo "Fin."
