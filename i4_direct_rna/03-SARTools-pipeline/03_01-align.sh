#!/bin/bash
set -euo pipefail

# input:
GENOMIC_REF_FASTA="$1"
GENOMIC_REF_MMI="$2"
READS_FASTQ="$3"
# output:
GENOMIC_BAM="$4"

# parameters:
NUM_THREADS=32
MM2_OPTIONS="-ax splice -uf -k14 --secondary=no"

if [ ! -f "$GENOMIC_REF_MMI" ]; then
    minimap2 -t$NUM_THREADS $MM2_OPTIONS \
        -d "$GENOMIC_REF_MMI" "$GENOMIC_REF_FASTA"
fi

minimap2 -t$NUM_THREADS $MM2_OPTIONS "$GENOMIC_REF_MMI" "$READS_FASTQ" \
    | samtools view -bh \
    > "$GENOMIC_BAM.unsorted"

samtools sort -@$NUM_THREADS -O BAM -o "$GENOMIC_BAM" "$GENOMIC_BAM.unsorted"
rm "$GENOMIC_BAM.unsorted"

samtools index -@$NUM_THREADS "$GENOMIC_BAM"
