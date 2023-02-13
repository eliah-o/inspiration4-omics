#!/bin/bash
set -euo pipefail

# input:
SUMMARY_TXT="$1"
# output:
PYCOQC_JSON="$2"

pycoQC -f "SUMMARY_TXT" -j "$PYCOQC_JSON"
