#!/bin/bash
set -euo pipefail

# input:
TRANSCRIPTOMIC_BAM="$1"
TRANSCRIPTOMIC_REF_FASTA="$2"
GUPPY_SUMMARY_TXT="$3"
READS_FASTQ="$4"
READS_BLOW5="$5"

# output:
EVENTALIGN_TSV="$6"

# parameters (depend on GPU capacity):
F5C_OPTIONS="-K 32768 -B 16M"

f5c_x86_64_linux_cuda eventalign --rna \
        -b "$TRANSCRIPTOMIC_BAM" \
        -g "$TRANSCRIPTOMIC_REF_FASTA" \
        -r "$READS_FASTQ" \
        --slow5 "$READS_FASTQ" \
        --summary "$GUPPY_SUMMARY_TXT" \
        --signal-index --scale-events \
        $F5C_OPTIONS \
    > "$EVENTALIGN_TSV"
