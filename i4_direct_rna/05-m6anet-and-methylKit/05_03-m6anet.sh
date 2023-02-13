#!/bin/bash
set -euo pipefail

# input:
EVENTALIGN_TSV="$1"
# intermediate:
M6ANET_DIR="$2"
# output:
M6ANET_TSV_GZ="$3"

# parameters:
NUM_THREADS=8

mkdir -p "$M6ANET_DIR/dataprep"
mkdir -p "$M6ANET_DIR/inference"

m6anet-dataprep --n_processes $NUM_THREADS \
    --eventalign "$EVENTALIGN_TSV" --out_dir "$M6ANET_DIR/dataprep"
m6anet-run_inference --n_processes $NUM_THREADS \
    --input_dir "$M6ANET_DIR/dataprep" --out_dir "$M6ANET_DIR/inference"
python python/csv2tsv.py \
    "$M6ANET_DIR/inference/data.result.csv" "$M6ANET_TSV_GZ"
