#!/bin/bash
set -euo pipefail

# input:
FAST5_PASS_DIR="$1"
# temporary files:
BLOW5_TEMP_DIR="$2"
FASTQ_TEMP_DIR="$3"
# output:
BLOW5="$4"
FASTQ_GZ="$5"

mkdir -p "$BLOW5_TEMP_DIR"
mkdir -p "$FASTQ_TEMP_DIR"

guppy_basecaller --disable_pings \
    --input_path "$FAST5_PASS_DIR" --recursive \
    --save_path "$FASTQ_TEMP_DIR" --compress_fastq \
    -c rna_r9.4.1_70bps_hac_prom.cfg --calib_detect \
    -x 'cuda:0' --num_callers 8 --gpu_runners_per_device 8 \
    --chunks_per_runner 1024 --chunk_size 1000 \
    --trim_barcodes --trim_adapters --trim_primers \
    --do_read_splitting --detect_mid_strand_adapter

cat "$FASTQ_TEMP_DIR/pass/*.fastq.gz" > "$FASTQ_GZ"
rm "$FASTQ_TEMP_DIR/pass/*.fastq.gz"

slow5tools f2s --to blow5 --iop 16 -d "$BLOW5_TEMP_DIR" "$FAST5_PASS_DIR"
slow5tools merge "$BLOW5_TEMP_DIR" -o "$BLOW5" -t16
rm -r "$BLOW5_TEMP_DIR"
f5c_x86_64_linux_cuda index -t 12 --slow5 "$BLOW5" "$FASTQ_GZ"

