#!/bin/bash
set -euo pipefail

# input:
TRANSCRIPTOMIC_REF_FASTA="$1"
TRANSCRIPTOMIC_REF_MMI="$2"
READS_FASTQ="$3"
# output:
TRANSCRIPTOMIC_BAM="$4"

# parameters:
NUM_THREADS=32
MM2_OPTIONS="-ax splice -uf -k14 --secondary=no"

if [ ! -f "$TRANSCRIPTOMIC_REF_MMI" ]; then
    minimap2 -t$NUM_THREADS $MM2_OPTIONS \
        -d "$TRANSCRIPTOMIC_REF_MMI" "$TRANSCRIPTOMIC_REF_FASTA"
fi

minimap2 -t$NUM_THREADS $MM2_OPTIONS "$TRANSCRIPTOMIC_REF_MMI" "$READS_FASTQ" \
    | samtools view -bh \
    > "$TRANSCRIPTOMIC_BAM.unsorted"

samtools sort -@$NUM_THREADS -O BAM \
    -o "$TRANSCRIPTOMIC_BAM" "$TRANSCRIPTOMIC_BAM.unsorted"
rm "$TRANSCRIPTOMIC_BAM.unsorted"

samtools index -@$NUM_THREADS "$TRANSCRIPTOMIC_BAM"
