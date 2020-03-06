#!/bin/bash

# usage
if [ -z "$2" ]; then
	echo "$0: read_filename data_descriptor [shasta_descriptor]"
	exit 1
fi

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
#set -u
# echo each line of the script to stdout so we can see what is happening
# to turn off echo do 'set +o xtrace'
set -o xtrace

##################
##### CONFIG #####
##################

# argument config
export READS=$1
export DESC=$2
if [ -z "$3" ]; then
	export SHASTA_DESC=""
else
	export SHASTA_DESC="$3_"
fi
export READS_TO_TRUE=$(basename $READS).reference.bam

# machine config
export TC=36
export TRUE_REF="../Microbial.reference.fasta"
export SPECIES="../Microbial.species.txt"
export SEGMENT_DIR="../segments"
export REFERENCE_DIR="../reference"
export SHASTA_ASM_DIR="shasta"
export SHASTA_CONFIG="/home/ubuntu/software/shasta/build/shasta-install/conf/Nanopore-Dec2019.conf"
export SHASTA_CONFIG="$SHASTA_CONFIG --Assembly.consensusCaller Modal"

# misc parameters
export SEP_FILTER="0x904"
export ALN_FILTER="0x904"
export TRUTH_FILTER="0x904"
export VALIDATE_SPECIES="Staphylococcus_aureus"
export TEST_SPECIES="Escherichia_coli"

# steps to run
export ALIGN_FULL=true
export EXTRACT_FQ=true
export ASM_SHASTA=true
export ALN_SHASTA=true
export TRAIN_TEST=true

##################
##### SCRIPT #####
##################

# align all to reference
if $ALIGN_FULL ; then
  minimap2 -ax map-ont -t $TC $TRUE_REF $READS | samtools view -hb -F $SEP_FILTER >unsorted.bam 
  samtools sort -@ $TC -o $READS_TO_TRUE unsorted.bam 
  samtools index -@ $TC $READS_TO_TRUE
  rm unsorted.bam
fi

# extract into per-species reads
if $EXTRACT_FQ ; then
  mkdir align_by_species
  mkdir fastq_by_species
  cat $SPECIES | xargz bash -c "extract_bam_segments.py -i $READS_TO_TRUE -s $SEGMENT_DIR/{}.reference_segments.txt -o align_by_species/{}.$DESC.bam -t $TC ; samtools fastq align_by_species/{}.$DESC.bam >fastq_by_species/{}.$DESC.fastq"
fi

# shasta
if $ASM_SHASTA ; then
  mkdir tmp $SHASTA_ASM_DIR
  cat $SPECIES | xargz bash -c "shasta --config $SHASTA_CONFIG --input fastq_by_species/{}.$DESC.fastq --assemblyDirectory tmp/{} ; cat tmp/{}/Assembly.fasta | sed 's/^>/>{}_/' >$SHASTA_ASM_DIR/{}.${SHASTA_DESC}shasta.fasta"
fi
 
# align to shasta
if $ALN_SHASTA ; then
  mkdir align_to_shasta
  cat $SPECIES | xargz bash -c "minimap2 -ax map-ont -t $TC $SHASTA_ASM_DIR/{}.${SHASTA_DESC}shasta.fasta fastq_by_species/{}.$DESC.fastq | samtools view -hb -F $ALN_FILTER >align_to_shasta/{}.unsorted.bam ; samtools sort -@ $TC -o align_to_shasta/{}.${DESC}_to_${SHASTA_DESC}shasta.$ALN_FILTER.bam align_to_shasta/{}.unsorted.bam ; samtools index -@ $TC align_to_shasta/{}.${DESC}_to_${SHASTA_DESC}shasta.$ALN_FILTER.bam ; rm align_to_shasta/{}.unsorted.bam"
fi

# make train test validate data
if $TRAIN_TEST ; then
  # prep 
  mkdir train_test
  mkdir train_test/tmp
  cp align_to_shasta/* train_test/tmp

  # validate
  mv train_test/tmp/$VALIDATE_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.$ALN_FILTER.bam train_test/Microbial.validate_$VALIDATE_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.bam
  mv train_test/tmp/$VALIDATE_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.$ALN_FILTER.bam.bai train_test/Microbial.validate_$VALIDATE_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.bam.bai
  cp $REFERENCE_DIR/$VALIDATE_SPECIES.reference.fasta train_test/Microbial.validate_$VALIDATE_SPECIES.truth.fasta

  # test
  mv train_test/tmp/$TEST_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.$ALN_FILTER.bam train_test/Microbial.test_$TEST_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.bam
  mv train_test/tmp/$TEST_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.$ALN_FILTER.bam.bai train_test/Microbial.test_$TEST_SPECIES.${DESC}_to_${SHASTA_DESC}shasta.bam.bai
  
  # train
  samtools merge -@ $TC train_test/Microbial.train.${DESC}_to_${SHASTA_DESC}shasta.bam train_test/tmp/*bam
  samtools index -@ $TC train_test/Microbial.train.${DESC}_to_${SHASTA_DESC}shasta.bam

  # truth align
  cat $SHASTA_ASM_DIR/* >train_test/Microbial.${SHASTA_DESC}shasta.fasta
  minimap2 -ax asm20 -t $TC train_test/Microbial.${SHASTA_DESC}shasta.fasta $TRUE_REF | samtools view -hb -F $TRUTH_FILTER >train_test/Microbial.truth_to_shasta.unsorted.bam
  samtools sort -@ $TC -o train_test/Microbial.truth_to_${SHASTA_DESC}shasta.bam train_test/Microbial.truth_to_shasta.unsorted.bam && rm train_test/Microbial.truth_to_shasta.unsorted.bam

  # cleanup
  rm -r train_test/tmp
fi

# fin
echo "Fin."
