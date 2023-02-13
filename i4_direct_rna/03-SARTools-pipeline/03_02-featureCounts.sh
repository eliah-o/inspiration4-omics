#!/bin/bash
set -euo pipefail

# input:
INDEX_FILE_PATH="$1"
ANNOTATION_FILE_PATH="$2"
GENOMIC_BAM="$3"
# output:
FEATURE_COUNTS="$4"

# parameters:
NUM_THREADS=32

eval featureCounts -T $NUM_OF_THREADS \
    -s 1 -O --fraction -J -L \
    -G "$INDEX_FILE_PATH" \
    -a "$ANNOTATION_FILE_PATH" -t exon -g gene_id \
    -o "$FEATURE_COUNTS" \
    "$GENOMIC_BAM.bam"
